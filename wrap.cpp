/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// c++_driver = simple example of how an umbrella program
//              can invoke LAMMPS as a library on some subset of procs
// Syntax: c++_driver P in.lammps
//         P = # of procs to run LAMMPS on
//             must be <= # of procs the driver code itself runs on
//         in.lammps = LAMMPS input script
// See README for compilation instructions

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "iostream"
#include "fstream"


#include "lammps.h"         // these are LAMMPS include files
#include "input.h"
#include "atom.h"
#include "library.h"
#include "group.h"
#include "domain.h"

#include "xrand.h"

//#include "many2one.h"
//#include "one2many.h"
//#include "files.h"
//#include "memory.h"
//#include "error.h"

#define QUOTE_(x) #x
#define QUOTE(x) QUOTE_(x)

#include "lmppath.h"
#include QUOTE(LMPPATH/src/lammps.h)
#include QUOTE(LMPPATH/src/library.h)
#include QUOTE(LMPPATH/src/input.h)
#include QUOTE(LMPPATH/src/modify.h)
#include QUOTE(LMPPATH/src/fix.h)
#include QUOTE(LMPPATH/src/fix_external.h)
#include QUOTE(LMPPATH/src/fix_rigid_nve_small.h)


using namespace LAMMPS_NS;
using namespace std;

struct Info {
  int me;
  LAMMPS *lmp;
};

FixRigidNVESmall *fixrigid;
FixExternal *fixexternal;
//FixRigidNHSmall *fixrigid;

void force_callback(void *, bigint, int, int *, double **, double **);
void simulate(LAMMPS* lmp);
double get_single_variable(LAMMPS* lmp, const char* pippo);
void setup(LAMMPS* lmp, int me);
void create_template_molecule(LAMMPS* lmp, int me, int ngrow);
void update(LAMMPS* lmp, int me, int lammps, int nsteps);
int const Nmax_bacteria = 10000;// when this value set very large, the run speed will slow down
int const N_mesh = 1000;
int create_bacterium(LAMMPS* lmp, double x0, double y0, double z0, double angle, int type, int state, int tolerance, double vx = 0, double vy = 0, double vz = 0, bool check = true);
void grow_molecules(LAMMPS* lmp);
int touch(LAMMPS* lmp, int ns, int ne, double radius);
void trail_record(LAMMPS* lmp);
void mesh(LAMMPS* lmp);
void trail_count(LAMMPS* lmp, FILE *out, double tempo);
void compute_velocity(LAMMPS* lmp);
void visit_frequency(int Ncharge,int Nmesh);
int num_bacteria = 0;
double total_time = 0;
// all times expressed in seconds
// t_run, t_tumble, together with viscosity, velocity and torque, fix the MSD
double t_run = 0.3;
double t_tumble = 0.5;
//double t_run = 0.99;
//double t_tumble = 0.01;
// reproduction time
double t_reproduction = 60*60; // = 1h

double t_stationary = 3;
double t_immobile = 0.001; // typical time before a mobile bacterium becomes immobile
double t_mobile = 0.002; // typical time before an immobile bacterium becomes mobile

bool psl = true;
bool change_motility = false;
bool reproduction = false;
bool only_running = false;// if you want only running, you also need to turn down swith_run_tumble
bool switch_run_tumble = true;
int number_initial = 1;
int Ncharge=1000;//maximum time of visiting
int count2=0;// for count how many times the visit_frequency funtion is used
// probabilites
double p_sibiling_bounded = 0.6;

double Lbox,Frun,Torque,Kn,DiameterPush,viscosity,Velocity,Mesh_size;
double Length_Molecule,Diameter,dt_integration;
int NBeads_Molecule, NumGrowMolecules,Nmesh,Grid_cover;

void write_positions(LAMMPS* lmp, FILE* out, double tempo);
void open_dump(LAMMPS* lmp, int seconds);
void close_dump(LAMMPS* lmp);
void delete_bacterium(LAMMPS* lmp, int m);
//void change_state(LAMMPS* lmp, int n);
void run_and_tumble();
void fix_times(int m);
void compute_msd(LAMMPS* lmp, int lammps, int me);



class Bact{
  public:
  int exist; // -1 don't exist;
  int nstart;
  int nend;
  int state; // 0 run, 1 tumble; -1 = never used
  int type; // growth state
  int mol_id;
  int bodytag;
  double X[2];
  double t_transition;// initial = 0
  double t_mobility;// initial = 0
  double Torque;
}Bacterium[Nmax_bacteria];

class Trail{
public:
  int q;
  int trail_id;//backup
}Trail[N_mesh][N_mesh];

int main(int narg, char **arg){
  system("rm dumps/*");
  init_random(1,1);
  // setup MPI and various communicators
  // driver runs on all procs in MPI_COMM_WORLD
  // comm_lammps only has 1st P procs (could be all or any subset)
  MPI_Init(&narg,&arg);
  MPI_Comm comm = MPI_COMM_WORLD;
  int me,nprocs;
  MPI_Comm_rank(comm,&me);
  MPI_Comm_size(comm,&nprocs);

  int nprocs_lammps = nprocs;
  if (nprocs_lammps > nprocs) {
    if (me == 0)
      printf("ERROR: LAMMPS cannot use more procs than available\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  int lammps;
  if (me < nprocs_lammps) lammps = 1;
  else lammps = MPI_UNDEFINED;
  MPI_Comm comm_lammps;
  MPI_Comm_split(MPI_COMM_WORLD,lammps,0,&comm_lammps);

  //FILE* out;

  // open LAMMPS input script

  // run the input script thru LAMMPS one line at a time until end-of-file
  // driver proc 0 reads a line, Bcasts it to all procs
  // (could just send it to proc 0 of comm_lammps and let it Bcast)
  // all LAMMPS procs call input->one() on the line

  LAMMPS *lmp;
  if (lammps == 1) lmp = new LAMMPS(0,NULL,comm_lammps);

  setup(lmp, me);
  //Mesh_size = Lbox/Nmesh;
  mesh(lmp);

//  total_time = 0;

//  compute_msd(lmp, lammps, me); return 0; // compute the MSD of crawling bacteria

  int angle, type, state;
  create_template_molecule(lmp,me,NumGrowMolecules);
  for(int q = 0; q < number_initial; q++){
    angle = -10; // angle < 0; angle will be random;
    type = 0; // 0,1,2,.... ;
    state = 0; // used to fix run, sticke, tumble,...;
    if(only_running == true) state = 0;
    else{ if(Xrandom() < (t_run/(t_run+t_tumble))) state = 0; else state = 1;}
    create_bacterium(lmp, Xrandom()*Lbox, Xrandom()*Lbox, 0, angle, type, state,0);
//    cout << q+1 << "	" << num_bacteria << endl;
  }
  run_and_tumble(); // fic the time of the transitions

  int time_dumps = 5; // dump ever X second
  int final_time = 1*60*60;  // we simulate 7h
  int total_step = final_time/dt_integration;
  double tau= Mesh_size/Velocity;// this is actually quite nice
  int timestep_record = int(tau/dt_integration);
  int timestep_grow = int(t_reproduction/(dt_integration*NumGrowMolecules));
  int next_grow = timestep_grow;
  int next_record = timestep_record;
//  cout << "Grow " <<  timestep_grow << "	" << timestep_grow*dt_integration << endl; getchar();

  int step_done = 0;
  open_dump(lmp,time_dumps);
  int timestep_todo;

  //ofstream outrecord("recortimes.dat");
  //ns=Bacterium[0].nstart;

  while(step_done < total_step){
    if(next_grow <= next_record){
//      cout << "AAAA " << timestep_todo << endl;
      timestep_todo = next_grow-step_done;
      update(lmp,me,lammps,timestep_todo);
      if(reproduction){
         grow_molecules(lmp);
//         cout << "Done reproduction after " << step_done*dt_integration << " sec " << endl; getchar();
      }
      if(next_grow == next_record){
        trail_record(lmp);
       //x=lmp->atom->x;
        //is = lmp->atom->map(ns);
        //outrecord <<  step_done<<' '<<x[is][0]<<' '<<x[is][1]<<endl;
        compute_velocity(lmp);
        next_record += timestep_record;
      }
      step_done += timestep_todo;
      next_grow += timestep_grow;
    }else{
      timestep_todo = next_record-step_done;
      update(lmp,me,lammps,timestep_todo);
      trail_record(lmp);
      compute_velocity(lmp);
      next_record += timestep_record;
      step_done += timestep_todo;
    }
  }
  close_dump(lmp);
  // for record the trail count
  MPI_Finalize();

  if (lammps == 1) delete lmp;
  return 0;
}

// dump every x seconds
void open_dump(LAMMPS* lmp, int seconds){
  int every = int(seconds/dt_integration);
  int not_every= 100*every;
  char run[1024];
  sprintf(run,"group TOTAL union mobile stationary"); lmp-> input->one(run);

  sprintf(run,"reset_timestep 0");
  lmp->input->one(run);
  sprintf(run,"dump WRITE mobile custom %d dumps/dumpfile.*.txt x y i_index mol",every);
  lmp->input->one(run);
  sprintf(run,"dump_modify WRITE pad 10");
  lmp->input->one(run);
  sprintf(run,"dump_modify WRITE sort id");
  lmp->input->one(run);

  sprintf(run,"dump WRITE2 stationary custom %d dumps/s_dumpfile.*.txt x y i_index mol", every);
  lmp->input->one(run);
  sprintf(run,"dump_modify WRITE2 pad 10");
  lmp->input->one(run);
  sprintf(run,"dump_modify WRITE2 sort id");
  lmp->input->one(run);

  sprintf(run,"dump WRITE3 MESH custom %d mesh/dumpfile.*.txt x y id q",not_every);
  lmp->input->one(run);
  sprintf(run,"dump_modify WRITE3 pad 10");
  lmp->input->one(run);
  sprintf(run,"dump_modify WRITE3 sort id");
  lmp->input->one(run);
}

void close_dump(LAMMPS* lmp){
  char run[1024];
  sprintf(run,"undump  WRITE"); lmp->input->one(run);
  sprintf(run,"undump  WRITE2"); lmp->input->one(run);
}



void update(LAMMPS* lmp, int me, int lammps, int nsteps){

  char run[1024];
  int ifix;
//  sprintf(run,"echo none"); lmp->input->one(run);


  sprintf(run,"fix NVE mobile rigid/nve/small molecule"); lmp->input->one(run);
  sprintf(run,"fix VIS     mobile     viscous %g", viscosity); lmp->input->one(run);
  sprintf(run,"neigh_modify exclude molecule all"); lmp->input->one(run);
  // make info avaiable to callback function
  // the purpose here is to addforce on previous forces per step
  Info info;
  info.me = me;
  info.lmp = lmp;
  // set callback to function inside fix external
  sprintf(run,"fix CALLBACK mobile external pf/callback 1 1"); // change all to mobile
  lmp->input->one(run);
  //don't know how it works here
  ifix = lmp->modify->find_fix("CALLBACK");
  fixexternal = (FixExternal *) lmp->modify->fix[ifix];
  fixexternal->set_callback(force_callback,&info);

  sprintf(run,"run %d", nsteps);
  lmp->input->one(run);
  sprintf(run,"unfix NVE");
  lmp->input->one(run);
  sprintf(run,"unfix CALLBACK");
  lmp->input->one(run);
  sprintf(run,"unfix VIS");
  lmp->input->one(run);
  fixrigid = NULL;


//  is = lmp->atom->map(Bacterium[0].nstart); ie =lmp->atom->map(Bacterium[0].nend);
//  cout << "EEEEE " << is << "	" << ie << "	" << lmp->atom->x[is][0] << "	" << lmp->atom->x[is][1] << "	" <<  "    " << lmp->atom->x[ie][0] << "   " << lmp->atom->x[ie][1] << endl;
//  getchar();
}

void delete_bacterium(LAMMPS* lmp, int m){
  char line[1024];
  int ns = Bacterium[m].nstart;
  int ne = Bacterium[m].nend;
  int k = 0; for(int q = 0; q < Nmax_bacteria; q++) if(Bacterium[q].exist != -1) k++;
  Bacterium[m].exist = -1;  //
  int k1 = 0; for(int q = 0; q < Nmax_bacteria; q++) if(Bacterium[q].exist != -1) k1++;
  sprintf(line,"group deleted id %d:%d", ns, ne); lmp->input->one(line);
  sprintf(line,"neigh_modify exclude group all deleted"); lmp->input->one(line); // no need to compute forces between deleted molecules
  sprintf(line,"group molecule_del id %d:%d", ns, ne); lmp->input->one(line);
  if( (Bacterium[m].state == 0) || (Bacterium[m].state == 1)){
    sprintf(line,"group mobile subtract mobile molecule_del"); lmp->input->one(line);
  }else if(Bacterium[m].state == 2){
    sprintf(line,"group stationary subtract stationary molecule_del"); lmp->input->one(line);
  }
//  sprintf(line,"delete_atoms group molecule_del mol yes compress no"); lmp->input->one(line);
//  we don't delete the molecule, but we shift it so that it does not interact with the others
  sprintf(line,"displace_atoms molecule_del move 0.0 0.0 10.0");lmp->input->one(line);
  sprintf(line,"group molecule_del delete"); lmp->input->one(line);
}

int create_bacterium(LAMMPS* lmp, double x0, double y0, double z0, double angle, int type, int state, int tolerance, double vx, double vy, double vz, bool check){
  int natoms = static_cast<int> (lmp->atom->natoms);
  int m = 0;
  int in_cylinder;
  char line[1024];
  int is, ie, ns, ne, max_id;
  double ** x;
  imageint *image;
  double angle_horizontal;
  double unwrap_ns[3];
  double unwrap_ne[3];
  double dx,dy,dz;

  m = 0;
  bool keep_loop = true;
  bool never_used, correct_type;
  while(keep_loop){
    if( Bacterium[m].exist < 0 ){ // free bacterium
      if( Bacterium[m].type == -1) never_used = true; else never_used = false;
      if( Bacterium[m].type == type) correct_type = true; else correct_type = false;
      if( (never_used) || (!never_used && correct_type) ) keep_loop = false;
    }
    if(keep_loop) m++;
    if(m == Nmax_bacteria){ cout << "Too many bacteria; increase Nmax_bacteria; I am not inserting;" << endl; return 0; };
  }

  // if the particle was not used before; fix the id and create it
  if(never_used){
    max_id = 0; for(int i = 0; i < natoms; i++) max_id = MAX(max_id,lmp->atom->tag[i]);
    ns = max_id+1; // global id
    ne = max_id+NBeads_Molecule; // global id
    Bacterium[m].nstart = ns;
    Bacterium[m].nend = ne;
    // create the particle in the center of the box
    sprintf(line,"create_atoms 0 single %g %g %g mol TemplateBacteria_%d 1", 0.5*Lbox, 0.5*Lbox, 0.0, type); lmp->input->one(line);
  }else{
    ns = Bacterium[m].nstart;
    ne = Bacterium[m].nend;
  }
  // group and image and local id
  sprintf(line,"group molecule id %d:%d", ns, ne); lmp->input->one(line);
  x = lmp->atom->x; // x seems the position for all the bacteria
  image = lmp->atom->image;
  ie = lmp->atom->map(ne); // local id
  is = lmp->atom->map(ns);  // local id
  lmp->domain->unmap(x[is],image[is],unwrap_ns);
  lmp->domain->unmap(x[ie],image[ie],unwrap_ne);
  // move the center of the particle to the center of the box
  dx = x[is][0]+(unwrap_ne[0]-unwrap_ns[0])/2; // center of the molecule, x
  dy = x[is][1]+(unwrap_ne[1]-unwrap_ns[1])/2; // center of the molecule, y
  dz = x[is][2]+(unwrap_ne[2]-unwrap_ns[2])/2; // center of the molecule, z
  dx = -dx+0.5*Lbox;// what is the box here? move every one of the bacterium to the center?
  dy = -dy+0.5*Lbox;// maybe just relevent to the center, to change the distrubution of bacteria to concertrate on the middle
  dz = -dz;
  sprintf(line,"displace_atoms molecule move %g %g %g",dx,dy,dz);lmp->input->one(line);

  if(angle < 0){ // rotate the particle by a random amount
     sprintf(line,"displace_atoms molecule rotate %g %g %g %g %g %g %g",0.5*Lbox,0.5*Lbox,0.0,0.0,0.0,1.0,Xrandom()*360);// Xrandom is a random variable between 0 to 1
     lmp->input->one(line);
     sprintf(line,"velocity molecule set %g %g %g",0.0,0.0,0.0); lmp->input->one(line); // we are growing a molecule; the velocity is the previous one
  }else{ // rotate the particle by a prescribed amount; this is used when a bacterium grows, as we replace a molecule with a larger one in the same position
     angle_horizontal = atan2(unwrap_ne[1]-unwrap_ns[1],unwrap_ne[0]-unwrap_ns[0])*180/M_PI;
     sprintf(line,"displace_atoms molecule rotate %g %g %g %g %g %g %g",0.5*Lbox,0.5*Lbox,0.0,0.0,0.0,1.0,angle-angle_horizontal); lmp->input->one(line);
     sprintf(line,"velocity molecule set %g %g %g",vx,vy,vz); lmp->input->one(line); // we are growing a molecule; the velocity is the previous one
  }
  // place the particle in the correct position: baricenter of the particle will be x0,y0
  sprintf(line,"displace_atoms molecule move %g %g %g",x0-0.5*Lbox,y0-0.5*Lbox,0.0);lmp->input->one(line);
  // set the index of the molecule; the index is an index we use to understand that a bacterium is the same when it grows
  sprintf(line,"set group molecule i_index %d", num_bacteria+1); lmp->input->one(line);
  Bacterium[m].mol_id = num_bacteria+1;
  sprintf(line,"set group molecule charge 1"); lmp->input->one(line);

  if(check) in_cylinder = touch(lmp, ns, ne, Diameter); else in_cylinder = 0; // there is no need to check for overlaps when a molecule is divided in two

  if( in_cylinder <= (NBeads_Molecule+tolerance)  ){
    Bacterium[m].exist = 1;
    Bacterium[m].state = state;
    Bacterium[m].type = type;
    Bacterium[m].X[0] = x[is][0];
    Bacterium[m].X[1] = x[is][1];
    if( (state == 0) || (state == 1) ){
      sprintf(line,"group mobile id %d:%d", ns, ne); lmp->input->one(line);
    }else if(state == 2){
      sprintf(line,"group stationary id %d:%d", ns, ne); lmp->input->one(line);
    }
    if( never_used == false ){
      sprintf(line,"group deleted subtract deleted molecule"); lmp->input->one(line); // eliminate the molecule from the group of the deleted particles
    }
    num_bacteria++;
    // remove the auxiliary group
    sprintf(line,"group molecule delete"); lmp->input->one(line);
    return ns; // one added molecule
  }
  // otherwise, delete the molecule
  delete_bacterium(lmp, m);
  sprintf(line,"group molecule delete"); lmp->input->one(line);
  return 0; // no added molecule
}

// number of sphers that are in contact with a region made of several sphers of radius radius, superiposing particles ns,ne
// we need to check the periodic boundary condition
int touch(LAMMPS* lmp, int ns, int ne, double radius){
  double** x = lmp->atom->x;
  int nspheres = 0;
  int n;
  double x0,y0;
  char line[1000];
  char line2[1000];
  sprintf(line2," ");
  int sx,ex, sy,ey;

  for(int i = ns; i <= ne; i++){
    n = lmp->atom->map(i);
    x0 = x[n][0];
    y0 = x[n][1];
    sx = ex = sy = ey = 0;
    if( x0 <= radius ) ex = 1;
    if( x0 >= (Lbox - radius)) sx = -1;
    if( y0 <= radius ) ey = 1;
    if( y0 >= (Lbox - radius)) sy = -1;
    for(int a = sx; a <= ex; a++){
      for(int b = sy; b <= ey; b++){
        sprintf(line,"region Sphere_%d sphere %g %g %g %g", nspheres, x0 + a*Lbox, y0+b*Lbox, x[n][2], radius);  //  sphere with radius = diameter of each bead
        lmp->input->one(line);
        sprintf(line," Sphere_%d", nspheres);
        strcat(line2,line);
        nspheres++;
      }
    }
  }
  // union of all spheres
  sprintf(line,"region Cylinder union %d %s", nspheres, line2); lmp->input->one(line);
  // group of all spheres in the cylinder

  sprintf(line,"group touch region Cylinder"); lmp->input->one(line);
  sprintf(line,"group touch subtract touch MESH");lmp->input->one(line); //do not account the mesh atoms
  // index of the group
  int igroup = lmp->group->find("touch");
  // number of particles within the group
  int in_cylinder = lmp->group->count(igroup);
  // remove the auxiliary groups and regions
  sprintf(line,"region Cylinder delete"); lmp->input->one(line);
  for(int i = 0; i < nspheres; i++){
    sprintf(line,"region Sphere_%d delete", i); lmp->input->one(line);
  }
  sprintf(line,"group touch delete"); lmp->input->one(line);
  return in_cylinder;
}

void grow_molecules(LAMMPS* lmp){
  int grow[Nmax_bacteria];// arrary of int
  double angle;
  double** x = lmp->atom->x;
  double** v = lmp->atom->v;
  char line[1024];
  sprintf(line,"index");
  int pp;
  int idx = lmp->atom->find_custom(line,pp);
  int* ind = lmp->atom->ivector[idx]; // vector containing the indeces of the atoms
  imageint *image = lmp->atom->image;
  double unwrap_ns[3], unwrap_ne[3], norm, norm2, vel[3], dr[3];
  int new_ns, new_ns2, is, ie, ns, ne;
  double x1,y1,x2,y2, index;
  bool stays;
  for(int m = 0; m < Nmax_bacteria; m++) grow[m] = Bacterium[m].exist;
//  FILE* out;
  int p;
//  bool fix_sibiling, ok_splitting;
  for(int m = 0; m < Nmax_bacteria; m++){
    if((grow[m] != -1)&(Bacterium[m].state !=2)){
      ns = Bacterium[m].nstart; // gloabal id
      ne = Bacterium[m].nend; // global id
//      state = Bacterium[m].state;
      is = lmp->atom->map(ns); // local id
      ie = lmp->atom->map(ne);   // local id
      lmp->domain->unmap(x[is],image[is],unwrap_ns);
      lmp->domain->unmap(x[ie],image[ie],unwrap_ne);
      index = ind[ie]; //  index of the bacterium
      norm2 = 0;
      for(int k = 0; k < 3; k++){
        dr[k] = unwrap_ne[k]-unwrap_ns[k];// size of the molecule
        norm2 += dr[k]*dr[k]; //mode
        vel[k] = v[ie][k];     // record of the velocity, the end of the molecule
      }
      vel[2] = 0;  //only on 2D for the moment
      norm = sqrt(norm2);
      angle = atan2(dr[1]/norm,dr[0]/norm);
      if(angle < 0) angle += 2*M_PI;
      angle *= 180/M_PI;//?
      if(Bacterium[m].type < NumGrowMolecules-1){
        x1 = (unwrap_ne[0]+unwrap_ns[0])/2.0;
        y1 = (unwrap_ne[1]+unwrap_ns[1])/2.0;
        // try adding a longer molecue; if attempt is successfull, then the shorter molecule is removed
        new_ns = create_bacterium(lmp, x1, y1, 0, angle, Bacterium[m].type+1, Bacterium[m].state, NBeads_Molecule, vel[0],vel[1],vel[2]);
        if(new_ns > 0){
           delete_bacterium(lmp, m);
           sprintf(line,"group molecule_new id %d:%d", new_ns, new_ns+NBeads_Molecule-1); lmp->input->one(line);
           sprintf(line,"set group molecule_new i_index %g", index); lmp->input->one(line);
           sprintf(line,"group molecule_new delete"); lmp->input->one(line);
           p = 0; while(Bacterium[p].nstart != new_ns) p++;
           Bacterium[p].t_mobility = Bacterium[m].t_mobility;
           Bacterium[p].state = Bacterium[m].state;
           Bacterium[p].Torque = Bacterium[m].Torque;
        }
      }else{ // split a molecule
      //  return;
//  x & y of the two short molecules
        // randomly select who is the father, who is the child; father will have coordinated x1,y1; child x2,y2
        if(Xrandom() < 0.5){
          x1 = x[is][0]+0.5*(Length_Molecule-Diameter)*(dr[0]/norm);
          y1 = x[is][1]+0.5*(Length_Molecule-Diameter)*(dr[1]/norm);
          x2 = x[ie][0]-0.5*(Length_Molecule-Diameter)*(dr[0]/norm);
          y2 = x[ie][1]-0.5*(Length_Molecule-Diameter)*(dr[1]/norm);
        }else{
          x2 = x[is][0]+0.5*(Length_Molecule-Diameter)*(dr[0]/norm);
          y2 = x[is][1]+0.5*(Length_Molecule-Diameter)*(dr[1]/norm);
          x1 = x[ie][0]-0.5*(Length_Molecule-Diameter)*(dr[0]/norm);
          y1 = x[ie][1]-0.5*(Length_Molecule-Diameter)*(dr[1]/norm);
        }

//        out = fopen("prima.txt","w"); write_positions(lmp,out,0); fclose(out);
        delete_bacterium(lmp,m);

        if(Xrandom() < p_sibiling_bounded) stays = true; else stays = false; // the silbiling stays with some probability


        //father
        new_ns = create_bacterium(lmp, x1, y1, 0, angle, 0, Bacterium[m].state, 0, vel[0],vel[1],vel[2], false);
        p = 0; while(Bacterium[p].nstart != new_ns) p++;
        Bacterium[p].t_transition = Bacterium[m].t_transition;
        Bacterium[p].Torque = Bacterium[m].Torque;
        Bacterium[p].t_mobility = Bacterium[m].t_mobility;
        Bacterium[p].state = Bacterium[m].state;
        sprintf(line,"group molecule_new id %d:%d", new_ns, new_ns+NBeads_Molecule-1); lmp->input->one(line);
        sprintf(line,"set group molecule_new i_index %g", index); lmp->input->one(line);
        sprintf(line,"group molecule_new delete"); lmp->input->one(line);

        if(stays){
          if(Xrandom() < 0.5) angle += 180;
          new_ns2 = create_bacterium(lmp, x2, y2, 0, angle, 0, Bacterium[m].state, 0, 0.0, 0.0,vel[2], false);
          p = 0; while(Bacterium[p].nstart != new_ns2) p++;
          Bacterium[p].Torque = Bacterium[m].Torque/(1+Bacterium[m].type*1.0/NumGrowMolecules);
          if(Bacterium[p].state == 0) Bacterium[p].t_transition = total_time -1.0/t_run*log(Xrandom());
          else if(Bacterium[p].state == 1) Bacterium[p].t_transition = total_time -1.0/t_tumble*log(Xrandom());
          Bacterium[m].t_mobility = total_time - 1.0/t_mobile*log(Xrandom());
        };
      }
    }
  }
}

void write_positions(LAMMPS* lmp, FILE *out, double tempo){
  double **x = lmp->atom->x;
  int is;
  for(int m = 0; m < Nmax_bacteria; m++){
    if(Bacterium[m].exist != -1){
       for(int n = 0; n < NBeads_Molecule; n++){
//       for(int n = 0; n < NBeads_Molecule; n += NBeads_Molecule-1){
         is = lmp->atom->map(Bacterium[m].nstart+n);
         fprintf(out,"%g %g %g\n", x[is][0], x[is][1], x[is][2]);
//         fprintf(out,"\n\n");
       }
       fprintf(out,"\n\n");
    }
  }
}

void mesh(LAMMPS* lmp){
  //this function is to produce the mesh and initial the atoms in the center of the grid
  char line[1064];
//  double x,y,z;
//  int new_id=0;
  //int natoms =static_cast<int>(lmp->atom->natoms);
  //how to really link with atom id?
  for(int i=0;i<Nmesh;i++){
    for(int j=0; j<Nmesh; j++){
//      x=Mesh_size*(i+0.5);
//      y=Mesh_size*(j+0.5);
//      z=0.0;
//      sprintf(line,"create_atoms 2 single %g %g %g", x,y,z);lmp->input->one(line);
//      new_id=new_id+1;
//      cout<<"new_id="<<new_id<<endl;
//      Trail[i][j].trail_id = new_id;
      Trail[i][j].q=0;
    }
  }
  // set all the mesh atom into one groups
  sprintf(line,"group MESH type 2");lmp->input->one(line);
//  sprintf(line,"neigh_modify exclude type 2 2"); lmp->input->one(line); // no need to compute forces between grid atoms
}

void trail_record(LAMMPS* lmp){
  double **x = lmp->atom->x;
  int is,im;
  int i,j;
  char line[1064];
  double xm,ym;





  if(psl){
  for(int m=0; m<Nmax_bacteria; m++){
    if(Bacterium[m].exist!= -1){
      is=lmp->atom->map(Bacterium[m].nstart);
      //im=is+int(NBeads_Molecule/2);
      im=is;
      xm = x[im][0]; if(xm < 0) xm += Lbox; else if(xm > Lbox) xm -= Lbox;
      ym = x[im][1]; if(ym < 0) ym += Lbox; else if(ym > Lbox) ym -= Lbox;
      i=floor(xm/Mesh_size);
      j=floor(ym/Mesh_size);

      //ofstream outij("ij.dat",ios::app);
      //outij<<i<<' '<<j<<endl;
      if(i == Nmesh) i-=1;
      if(i < 0) i = 0;
      if(j == Nmesh) j-=1;
      if(j < 0) j = 0;

      if(Trail[i][j].q == 0){ // site has never been visited, we create an atom
        sprintf(line,"create_atoms 2 single %g %g %g", Mesh_size*(i+0.5),Mesh_size*(j+0.5),0.0);lmp->input->one(line);
        Trail[i][j].trail_id = lmp->atom->natoms;
        sprintf(line,"group MESH type 2");lmp->input->one(line);
        //sprintf(line,"neigh_modify exclude type 2 2"); lmp->input->one(line); // no need to compute forces between grid atoms
      }else if(Trail[i][j].q<Ncharge){
      Trail[i][j].q= Trail[i][j].q+1;
      sprintf(line,"variable nm equal %d",Trail[i][j].trail_id);lmp->input->one(line);
      sprintf(line,"set atom ${nm} charge %d",Trail[i][j].q); lmp->input->one(line);
      }
      //if(Trail[i][j].q == 0){ // site has never been visited, we create an atom
        //sprintf(line,"create_atoms 2 single %g %g %g", Mesh_size*(i+0.5),Mesh_size*(j+0.5),0.0);lmp->input->one(line);
      ///  Trail[i][j].trail_id = lmp->atom->natoms;
      //  sprintf(line,"group MESH type 2");lmp->input->one(line);
      //  sprintf(line,"neigh_modify exclude type 2 2"); lmp->input->one(line); // no need to compute forces between grid atoms
      //}
      //Trail[i][j].q= Trail[i][j].q+1;
      //sprintf(line,"variable nm equal %d",Trail[i][j].trail_id);lmp->input->one(line);
      //sprintf(line,"set atom ${nm} charge %d",Trail[i][j].q); lmp->input->one(line);
    }
  }
  }

}



double get_single_variable(LAMMPS* lmp, const char* pippo){
  char variable[100];
  strcpy(variable, pippo);
  char S_NULL[100]; sprintf(S_NULL,"NULL");
  double* pippone = (double *) lammps_extract_variable(lmp, variable, S_NULL);
  return pippone[0];
}


// run the setupfile, and get some info about the particles
void setup(LAMMPS* lmp, int me){
  lmp->input->file("in.setup");
  Length_Molecule = get_single_variable(lmp,"LL");
  NBeads_Molecule = int(get_single_variable(lmp,"NumAtomMolecule"));
  Diameter = get_single_variable(lmp,"Diameter");
  Lbox = get_single_variable(lmp,"Lbox");
  Frun = get_single_variable(lmp,"Frun");
  Torque = get_single_variable(lmp,"Torque");
  dt_integration = get_single_variable(lmp,"dt_integration");
  Kn = get_single_variable(lmp,"Kn");
  DiameterPush = get_single_variable(lmp,"DiameterPush");
  viscosity = get_single_variable(lmp,"viscosity");
  NumGrowMolecules = get_single_variable(lmp,"NumGrowMolecules");
  Velocity= get_single_variable(lmp,"velocity");
  Nmesh=get_single_variable(lmp,"Nmesh");
  Mesh_size=get_single_variable(lmp,"mesh_size");
  for(int n = 0; n < Nmax_bacteria; n++){
    Bacterium[n].exist = -1;
    Bacterium[n].type = -1;
  }
}




void create_template_molecule(LAMMPS* lmp, int me, int ngrow){
  double N = NBeads_Molecule;
  double LL = Length_Molecule; // initial lenght; the final one will be the double
  double length;
  char line[1024];
  system("rm template.bacteria*"); // eliminate existing templates
  for(int k = 0; k < ngrow; k++){
    sprintf(line,"template.bacteria_%d",k);
    ofstream outfile(line);
    length = LL*(1+k*1.0/(ngrow-1)*1.01);
    outfile << "#template for a bacteria with " << N  << "atoms, of lenght " << length << endl << endl;
    outfile << N << " atoms" << endl;
    outfile << 0 << " bonds" << endl;
    outfile << endl;
    outfile << "Coords" << endl << endl;
    for(int n = 0; n < N; n++){
      outfile << n+1 << "	" << 0.5*Diameter+(length-Diameter)*n*1.0/(N-1) << "	" << 0 << "	" << 0 << endl;
    }
    outfile << endl;
    outfile << "Types" << endl << endl;
    for(int n = 0; n < N; n++) outfile << n+1 << "	" << 1 << endl;
    outfile.close();

    sprintf(line,"molecule TemplateBacteria_%d template.bacteria_%d",k,k);
    lmp->input->one(line);
  }
}



void force_callback(void *ptr, bigint ntimestep, int nlocal, int *id, double **x, double **f) {
  double norm, norm2;
  int ns, ne, is, ie, state;
  int fx,fy;
  double dr[3];

  Info *info = (Info *) ptr;
  for(int q = 0; q < nlocal; q++) f[q][0] = f[q][1] = f[q][2] = 0.0;

//  bool check = false;

  char stationary[1000];
  char mobiles[1000];
  char buffer[1000];
  sprintf(stationary,"group new_stationary id ");
  sprintf(mobiles,"group new_mobiles id ");
  int n_stationary = 0;
  int n_mobiles = 0;
  float Inertia;

  for(int m = 0; m < Nmax_bacteria; m++){
    if(Bacterium[m].exist > 0){
      state = Bacterium[m].state;
      // get info about the versor of the molecule
      ns = Bacterium[m].nstart; // gloab id
      ne = Bacterium[m].nend;   // global id
      is = info->lmp->atom->map(ns); // local id
      ie = info->lmp->atom->map(ne);   // local id
      norm2 = 0;
      for(int k = 0; k < 3; k++){
        dr[k] = x[ie][k]-x[is][k]; if(dr[k] > 0.5*Lbox) dr[k] -= Lbox; if(dr[k] < -0.5*Lbox) dr[k] += Lbox;
        norm2 += dr[k]*dr[k];
      }
      norm = sqrt(norm2);
      // check mobile <-> mobile

      if(change_motility && (total_time > Bacterium[m].t_mobility)){ // need to change mobility
        if( (state == 0) || (state == 1) ){
          Bacterium[m].state = 2;
          sprintf(stationary,"%s %d:%d",stationary,ns,ne); n_stationary++;
        }else{
          if(Xrandom() < (t_run/(t_run+t_tumble))) Bacterium[m].state = 0; else Bacterium[m].state = 1;
          sprintf(mobiles,"%s %d:%d",mobiles,ns,ne); n_mobiles++;
        }
        fix_times(m);


        state = Bacterium[m].state;
      }

      // run
      if(state == 0){ // run
        f[ie][0] = Frun*dr[0]/norm;
        f[ie][1] = Frun*dr[1]/norm;
      }
      // tumble
      if(state == 1){ // tumble
        Inertia = powl(1+Bacterium[m].type*1.0/NumGrowMolecules,2); // we rescale the applied force by the momentum of inertia, so that rotations do not depend on the length (but for the effect of viscosity)
        fx = 0.5*Bacterium[m].Torque * dr[1]/norm2 * Inertia;
        fy = -0.5*Bacterium[m].Torque * dr[0]/norm2 * Inertia;
        f[is][0] = fx;
        f[is][1] = fy;
        f[ie][0] = -fx;
        f[ie][1] = -fy;
      }
      if( switch_run_tumble && ( (state == 0) || (state == 1) ) ){ // if mobile
        if(total_time > Bacterium[m].t_transition){ // need to change state from run to tumble
          if(Bacterium[m].state == 0){ // from running to tumbling
            Bacterium[m].state = 1;
            Bacterium[m].t_transition = total_time -1.0/t_tumble*log(Xrandom());
          }else if(Bacterium[m].state == 1){ // from tumbling to running
            Bacterium[m].state = 0;
            Bacterium[m].t_transition = total_time -1.0/t_run*log(Xrandom());
            if(Xrandom() < 0.5) Bacterium[m].Torque = Torque;
            else Bacterium[m].Torque = -Torque;
          }
        }
      };
    }
  }
  total_time = total_time+dt_integration;
  LAMMPS *lmp = info->lmp;
  if(n_stationary > 0){
    lmp->input->one(stationary);
    sprintf(buffer,"group stationary union stationary new_stationary"); lmp->input->one(buffer);
    sprintf(buffer,"group mobile subtract mobile new_stationary"); lmp->input->one(buffer);
    sprintf(buffer,"group new_stationary delete"); lmp->input->one(buffer);
  }
  if(n_mobiles > 0){
    lmp->input->one(mobiles);
    sprintf(buffer,"group mobile union mobile new_mobiles"); lmp->input->one(buffer);
    sprintf(buffer,"group stationary subtract stationary new_mobiles"); lmp->input->one(buffer);
    sprintf(buffer,"group new_mobiles delete"); lmp->input->one(buffer);
  }
}

void fix_times(int m){
  if(Bacterium[m].exist > 0){
    if(Bacterium[m].state == 0){ // the bacterium is running
       Bacterium[m].t_transition = total_time - 1.0/t_run*log(Xrandom());
       Bacterium[m].t_mobility = total_time - 1.0/t_mobile*log(Xrandom());
    };
    if(Bacterium[m].state == 1){
       Bacterium[m].t_transition = total_time - 1.0/t_tumble*log(Xrandom());
       Bacterium[m].t_mobility = total_time - 1.0/t_mobile*log(Xrandom());
       if(Xrandom() < 0.5) Bacterium[m].Torque = Torque;
       else Bacterium[m].Torque = -Torque;
    };
    if(Bacterium[m].state == 2){
       Bacterium[m].t_mobility = total_time - 1.0/t_immobile*log(Xrandom());
    }
  }
}


void run_and_tumble(){
  for(int m=0; m< Nmax_bacteria; m++) fix_times(m);
}

void compute_msd(LAMMPS* lmp, int lammps, int me){
  change_motility = false;
  reproduction = false;
  only_running = false;
  switch_run_tumble = true;
  int angle, type, state;
  create_template_molecule(lmp,me,NumGrowMolecules);
  int step_reproduction = int(t_reproduction/(dt_integration*NumGrowMolecules));
  for(int q = 0; q < number_initial; q++){
    angle = -10; // angle < 0; angle will be random;
    type = 0; // 0,1,2,.... ;
    state = 0; // used to fix run, sticke, tumble,...;
    if(only_running == true) state = 0;
    else{ if(Xrandom() < (t_run/(t_run+t_tumble))) state = 0; else state = 1;}
    create_bacterium(lmp, Xrandom()*Lbox, Xrandom()*Lbox, 0, angle, type, state,0);
    cout << q+1 << "	" << num_bacteria << endl;
  }
  run_and_tumble(); // fic the time of the transitions
  int index_center = int(NBeads_Molecule/2)+1;
  int natoms = NBeads_Molecule*number_initial;
  char buffer[100];
  sprintf(buffer,"pair_coeff      1 1 0 ${Diameter}"); lmp->input->one(buffer);
  sprintf(buffer,"group centers id %d:%d:%d",index_center,natoms,NBeads_Molecule); lmp->input->one(buffer);

  // thermalize
  for(int q = 0; q < 10; q++) update(lmp,me,lammps,step_reproduction);
  // msd
  sprintf(buffer,"compute MSD centers msd"); lmp->input->one(buffer);
  sprintf(buffer,"reset_timestep 0"); lmp->input->one(buffer);
  sprintf(buffer,"variable TIME_MSD equal step*dt"); lmp->input->one(buffer);
  sprintf(buffer,"variable cmsd equal c_MSD[1]+c_MSD[2]");lmp->input->one(buffer);
  sprintf(buffer,"thermo_style    custom step atoms temp pe ke press v_cmsd"); lmp->input->one(buffer);
  sprintf(buffer,"fix fixprint all print 100 '${TIME_MSD} ${cmsd}' file MSD.txt screen no"); lmp->input->one(buffer);

  for(int q = 0; q < 10000; q++){
    update(lmp,me,lammps,step_reproduction);
    if(reproduction) grow_molecules(lmp);
  }
}

void compute_velocity(LAMMPS* lmp){
  int count=0;
  double ** v = lmp->atom->v;
  double vsum=0.0;
  double vxsum=0.0;
  double vysum=0.0;
  double vzsum=0.0;
  int ne,ie;
  double Vx,Vy,Vz,v1,v_average,phi; //phi is a parameter

  for(int m = 0; m < Nmax_bacteria; m++){
    if(Bacterium[m].exist > 0){
      ne = Bacterium[m].nend;   // global id
      ie = lmp->atom->map(ne);   // local id
      Vx=v[ie][0];
      Vy=v[ie][1];
      Vz=v[ie][2];
      v1=sqrt(Vx*Vx+Vy*Vy);
      //v1=sqrt(v[ie][0]*v[ie][0]+v[ie][1]*v[ie][1]+v[ie][2]*v[ie][2]);
      vsum+=v1;
      vxsum+=Vx;
      vysum+=Vy;
      vzsum+=Vz;
      count+=1;
    }
  }
  v_average=vsum/count;
  phi=sqrt(abs(vxsum*vxsum+vysum*vysum+vzsum*vzsum))/(count*v_average);

  ofstream outfile("velocity.txt",ios::app);
  outfile<<v_average<<" "<<Vx<<" "<<Vy<<" "<<phi<<endl;
  outfile.close();
}

void visit_frequency(int Ncharge,int Nmesh){
  //no need for this for the moment, by use mesh dump output and read datas from python to realize this
  int Q[Ncharge];
  char buffer[1064];

  for(int k=0; k<Ncharge; k++){
    Q[k]=0;
  }

  if(psl){
    count2+=1;
  for(int i=0;i<Nmesh;i++){
    for(int j=0; j<Nmesh; j++){
      if(Trail[i][j].q!=0){
        for(int k=0; k<Ncharge; k++){
          if(Trail[i][j].q == k+1) Q[k]+=1;
        }
      }
    }
  }
}
  sprintf(buffer,"charge_%d.dat",count2);
  ofstream outfile(buffer);
  for(int k=0; k<Ncharge; k++){
    outfile<<k<<"	" <<Q[k]<<endl;
  }
  outfile.close();
}

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

void force_callback(void *, bigint, int, int *, double **, double **);

struct Info {
  int me;
  LAMMPS *lmp;
};

FixRigidNVESmall *fixrigid;
FixExternal *fixexternal;
//FixRigidNHSmall *fixrigid;

void simulate(LAMMPS* lmp);
double get_single_variable(LAMMPS* lmp, const char* pippo);
void setup(LAMMPS* lmp, int me);
void create_template_molecule(LAMMPS* lmp, int me, int ngrow);
void update(LAMMPS* lmp, int me, int lammps, int nsteps);
int const Nmax_bacteria = 10000;
int create_bacterium(LAMMPS* lmp, double x0, double y0, double z0, double angle, int type, int state, int tolerance, double vx = 0, double vy = 0, double vz = 0, bool check = true);
void grow_molecules(LAMMPS* lmp);
int touch(LAMMPS* lmp, int ns, int ne, double radius);
void changestate(LAMMPS* lmp);

int num_bacteria = 0;
double total_time = 0;
double t_run = 10;
double t_tumble = 2;
double t_reproduction = 30;
double t_immobile = 10;
//double p_move_again = 0.05;
double p_detach = 0.01;

// probabilites
double p_sibiling_bounded = 0.6;

double Lbox,Frun,Torque,Kn,DiameterPush,viscosity;
double Length_Molecule,Diameter,dt_integration;
int NBeads_Molecule, NumGrowMolecules;
void write_positions(LAMMPS* lmp, FILE* out, double tempo);
void open_dump(LAMMPS* lmp);
void close_dump(LAMMPS* lmp);
void delete_bacterium(LAMMPS* lmp, int m);
//void change_state(LAMMPS* lmp, int n);
void run_and_tumble(LAMMPS* lmp);
void immobile(LAMMPS* lmp);


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
  double Torque;
}Bacterium[Nmax_bacteria];

int main(int narg, char **arg){
  init_random(1,1);
  char run[1024];
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

  // open LAMMPS input script

  // run the input script thru LAMMPS one line at a time until end-of-file
  // driver proc 0 reads a line, Bcasts it to all procs
  // (could just send it to proc 0 of comm_lammps and let it Bcast)
  // all LAMMPS procs call input->one() on the line

  LAMMPS *lmp;
  if (lammps == 1) lmp = new LAMMPS(0,NULL,comm_lammps);

  setup(lmp, me);

  total_time = 0;

  int angle, type, state;
  create_template_molecule(lmp,me,NumGrowMolecules);
  for(int q = 0; q < 50; q++){
    angle = -10; // angle < 0; angle will be random;
    type = 0; // 0,1,2,.... ;
    state = 0; // used to fix run, sticke, tumble,...;
    if(Xrandom()<(t_immobile/(t_immobile+t_run+t_tumble))){
      state = 2;
    }else{
      if(Xrandom() < (t_run/(t_run+t_tumble))) state = 0; else state = 1;
    }
    create_bacterium(lmp, Xrandom()*Lbox, Xrandom()*Lbox, 0, angle, type, state,0);
    cout << q+1 << "	" << num_bacteria << endl;
  }
  run_and_tumble(lmp); // fic the time of the transitions


  FILE* out = fopen("start.txt","w");
  write_positions(lmp, out, 0);
  fclose(out);

  int ns = Bacterium[1].nstart;
  int ne = Bacterium[1].nend;
  sprintf(run,"group stationary id %d:%d", ns,ne);  lmp->input->one(run);
  sprintf(run,"group Delete id  %d:%d", ns,ne);  lmp->input->one(run);
  sprintf(run,"group stationary subtract stationary Delete"); lmp-> input->one(run);
  sprintf(run,"group Delete delete"); lmp->input->one(run);


  open_dump(lmp);

  // final time = 10000*0.002 = 20

  ofstream outfile("pippo.dat");

  // every t_reproduction a molecule split in two, if it can grow freely
  int step_reproduction = int(t_reproduction/(dt_integration*NumGrowMolecules));

//  for(int q = 0; q < 1000; q++){
//    update(lmp,me,lammps,step_reproduction);
//    grow_molecules(lmp);
//  }

  outfile.close();
  for(int q = 0; q < 100; q++){
    update(lmp,me,lammps,step_reproduction);
    grow_molecules(lmp);
//    create_bacterium(lmp, Xrandom()*Lbox, Xrandom()*Lbox, 0, angle, type, state,0);
  }


  close_dump(lmp);
  FILE* end = fopen("end.txt","w");
  write_positions(lmp, end, 0);
  fclose(end);

  MPI_Finalize();

  if (lammps == 1) delete lmp;
  return 0;
}

void open_dump(LAMMPS* lmp){
  char run[1024];
  sprintf(run,"group TOTAL union mobile stationary"); lmp-> input->one(run);

  sprintf(run,"reset_timestep 0");
  lmp->input->one(run);
//  sprintf(run,"dump WRITE mobile custom 100 dumps/dumpfile.*.txt x y mol type");
  sprintf(run,"dump WRITE TOTAL custom 200 dumps/dumpfile.*.txt x y q mol");
  lmp->input->one(run);
  sprintf(run,"dump_modify WRITE pad 10");
  lmp->input->one(run);
  sprintf(run,"dump_modify WRITE sort id");
  lmp->input->one(run);
}

void close_dump(LAMMPS* lmp){
  char run[1024];
  sprintf(run,"undump  WRITE");
  lmp->input->one(run);
}



void update(LAMMPS* lmp, int me, int lammps, int nsteps){

  char run[1024];
  int ifix;
//  int ns;
//  int ne;
//  sprintf(run,"echo none"); lmp->input->one(run);


  changestate(lmp);
/*
// the following code is meat for make some bacteria stop in the simulation
  for(int m = 0; m < Nmax_bacteria; m++)
  {   if(Bacterium[m].state == 0|| Bacterium[m].state == 1){
      if(Xrandom()< t_stationary/(t_stationary+t_tumble+t_run)) {
      Bacterium[m].state = 2;  // for the moment, I didn't set time for stationary period. maybe can write all in total
      //ns = Bacterium[m].nstart;
      //ne = Bacterium[m].nend;
      //sprintf(run,"group stationary id %d:%d", ns, ne); lmp->input->one(run);
      //sprintf(run,"group mobile subtract mobile stationary"); lmp-> input->one(run);
        }
    }
  }

  for(int m=0; m< Nmax_bacteria; m++){
    if(Bacterium[m].state == 2){
    if(Xrandom()< p_move_again) {
      if(Xrandom() < (t_run/(t_run+t_tumble))) Bacterium[m].state = 0; else Bacterium[m].state = 1;
      run_and_tumble(lmp, t_run, t_tumble);
      //ns = Bacterium[m].nstart;
      //ne = Bacterium[m].nend;
      //sprintf(run,"group moveagain id %d:%d", ns, ne); lmp->input->one(run);
      //sprintf(run,"group mobile union mobile moveagain"); lmp-> input->one(run);
    }
  }
}

*/

  //    delete some bacteria during the run, give a chance to detach from surface
  // not sure, maybe I should put this in update
  /*
  for(int m=0 ;m < Nmax_bacteria; m++){
      if(Xrandom()< p_detach){
      delete_bacterium(lmp,m);
     }
    }
// this part end here
*/
//  sprintf(run,"group TOTAL union mobile stationary"); lmp-> input->one(run);
  sprintf(run,"fix NVE mobile rigid/nve/small molecule"); lmp->input->one(run);

//  ifix = lmp->modify->find_fix("NVE");
//  fixrigid = (FixRigidNVESmall *) lmp->modify->fix[ifix];

  sprintf(run,"fix VIS     mobile     viscous %g", viscosity);
  lmp->input->one(run);
// the immobile bacteria has much higher viscosity
  double viscosity2= 50;
  sprintf(run,"fix VIS2     stationary     viscous %g", viscosity2);
  lmp->input->one(run);

//  sprintf(run,"neigh_modify exclude group all deleted"); lmp->input->one(run);
  sprintf(run,"neigh_modify exclude molecule all"); lmp->input->one(run);
  // make info avaiable to callback function
  // the purpose here is to addforce on previous forces per step
  Info info;
  info.me = me;
  info.lmp = lmp;
  // set callback to function inside fix external
  sprintf(run,"fix CALLBACK all external pf/callback 1 1"); // change all to mobile
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
  Bacterium[m].exist = -1;  // don't understand here, why k and k1
  int k1 = 0; for(int q = 0; q < Nmax_bacteria; q++) if(Bacterium[q].exist != -1) k1++;
  sprintf(line,"group deleted id %d:%d", ns, ne); lmp->input->one(line);
  sprintf(line,"neigh_modify exclude group all deleted"); lmp->input->one(line); // no need to compute forces between deleted molecules
  sprintf(line,"group molecule_del id %d:%d", ns, ne); lmp->input->one(line);
  if(Bacterium[m].state == 0||Bacterium[m].state == 1){
    sprintf(line,"group mobile subtract mobile molecule_del"); lmp->input->one(line);
  }else if(Bacterium[m].state == 2){
    sprintf(line,"group stationary subtract stationary molecule_del"); lmp->input->one(line);
  }
  sprintf(line,"group mobile subtract mobile molecule_del"); lmp->input->one(line);
//  sprintf(line,"delete_atoms group molecule_del mol yes compress no"); lmp->input->one(line);
//  we don't delete the molecule, but we shift it so that it does not interact with the others
  sprintf(line,"displace_atoms molecule_del move 0.0 0.0 10.0");lmp->input->one(line);
  sprintf(line,"group molecule_del delete"); lmp->input->one(line);
}

int create_bacterium(LAMMPS* lmp, double x0, double y0, double z0, double angle, int type, int state, int tolerance, double vx, double vy, double vz, bool check){
  int natoms = static_cast<int> (lmp->atom->natoms);
  int m = 0;
//  int n;
  int in_cylinder;
  char line[1024];
//  char line2[1024];
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
//    sprintf(line,"group molecule id %d:%d", ns, ne); lmp->input->one(line);
  }else{
    ns = Bacterium[m].nstart;
    ne = Bacterium[m].nend;
//    sprintf(line,"group molecule id %d:%d", ns, ne); lmp->input->one(line);
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
  // set the charge of the molecule; the charge is an index we use to understand that a macterium is the same when it grows
  sprintf(line,"set group molecule charge %d", num_bacteria+1); lmp->input->one(line);
  Bacterium[m].mol_id = num_bacteria+1;

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
  double* q = lmp->atom->q;
  imageint *image = lmp->atom->image;
  double unwrap_ns[3], unwrap_ne[3], norm, norm2, vel[3], dr[3];
  int new_ns, new_ns2, is, ie, ns, ne;
  double x1,y1,x2,y2, charge;
  bool stays;
  char line[1024];
  for(int m = 0; m < Nmax_bacteria; m++) grow[m] = Bacterium[m].exist;
  FILE* out;
  int p;
  bool fix_sibiling, ok_splitting;
  for(int m = 0; m < Nmax_bacteria; m++){
    if(grow[m] != -1){
      ns = Bacterium[m].nstart; // gloabal id
      ne = Bacterium[m].nend; // global id
      is = lmp->atom->map(ns); // local id
      ie = lmp->atom->map(ne);   // local id
      lmp->domain->unmap(x[is],image[is],unwrap_ns);
      lmp->domain->unmap(x[ie],image[ie],unwrap_ne);
      charge = q[ie]; // charge = index of the bacterium
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
           sprintf(line,"set group molecule_new charge %g", charge); lmp->input->one(line);
           sprintf(line,"group molecule_new delete"); lmp->input->one(line);
        }
      }else{ // split a molecule
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

/*
        if(Xrandom() < p_sibiling_bounded) stays = true; else stays = false; // the silbiling stays with some probability

        new_ns = new_ns2 = 0;
        new_ns = create_bacterium(lmp, x1, y1, 0, angle, 0, Bacterium[m].state, NBeads_Molecule, vel[0],vel[1],vel[2]);

        fix_sibiling = false;
        ok_splitting = false;
        if(new_ns > 0){
          ok_splitting = true;
          if(stays){
             if(Xrandom() < 0.5) angle += 180;
             new_ns2 = create_bacterium(lmp, x2, y2, 0, angle, 0, Bacterium[m].state, NBeads_Molecule, 0.0, 0.0,vel[2]);
             if(new_ns2 > 0) fix_sibiling = true; else ok_splitting = false;
          }
          if(ok_splitting){
            delete_bacterium(lmp,m);
            p = 0; while(Bacterium[p].nstart != new_ns) p++;
            Bacterium[p].t_transition = Bacterium[m].t_transition;
            Bacterium[p].Torque = Bacterium[m].Torque;
            sprintf(line,"group molecule_new id %d:%d", new_ns, new_ns+NBeads_Molecule-1); lmp->input->one(line);
            sprintf(line,"set group molecule_new charge %g", charge); lmp->input->one(line);
            sprintf(line,"group molecule_new delete"); lmp->input->one(line);
            if(fix_sibiling){
              p = 0; while(Bacterium[p].nstart != new_ns2) p++;
              Bacterium[p].Torque = Bacterium[m].Torque;
              if(Bacterium[p].state == 0) Bacterium[p].t_transition = total_time -1.0/t_run*log(Xrandom());
              else Bacterium[p].t_transition = total_time -1.0/t_tumble*log(Xrandom());
            }
          }else{
            p = 0; while(Bacterium[p].nstart != new_ns) p++;
            delete_bacterium(lmp,p);
          }
        };
*/

        // the splitting is surely successfull
//        out = fopen("prima.txt","w"); write_positions(lmp,out,0); fclose(out);


        delete_bacterium(lmp,m);

        if(Xrandom() < p_sibiling_bounded) stays = true; else stays = false; // the silbiling stays with some probability


        //father
        new_ns = create_bacterium(lmp, x1, y1, 0, angle, 0, Bacterium[m].state, 0, vel[0],vel[1],vel[2], false);
        p = 0; while(Bacterium[p].nstart != new_ns) p++;
        Bacterium[p].t_transition = Bacterium[m].t_transition;
        Bacterium[p].Torque = Bacterium[m].Torque;
        sprintf(line,"group molecule_new id %d:%d", new_ns, new_ns+NBeads_Molecule-1); lmp->input->one(line);
        sprintf(line,"set group molecule_new charge %g", charge); lmp->input->one(line);
        sprintf(line,"group molecule_new delete"); lmp->input->one(line);

        if(stays){
          if(Xrandom() < 0.5) angle += 180;
          new_ns2 = create_bacterium(lmp, x2, y2, 0, angle, 0, Bacterium[m].state, 0, 0.0, 0.0,vel[2], false);
          p = 0; while(Bacterium[p].nstart != new_ns2) p++;
          Bacterium[p].Torque = Bacterium[m].Torque;
          if(Bacterium[p].state == 0) Bacterium[p].t_transition = total_time -1.0/t_run*log(Xrandom());
          else if(Bacterium[p].state == 1) Bacterium[p].t_transition = total_time -1.0/t_tumble*log(Xrandom());
        };// might need to add bacterium[p].state == 2

//        out = fopen("dopo.txt","w"); write_positions(lmp,out,0); fclose(out);
//        cout << "CHECK " << endl; getchar();

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


double get_single_variable(LAMMPS* lmp, const char* pippo){
// char* l = pippo.c_str();
 char variable[100];
 strcpy(variable, pippo);
// strstr((char *)variable,pippo);
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
  int ns, ne, is, ie;
  int fx,fy;
  double dr[3];

  Info *info = (Info *) ptr;
  for(int q = 0; q < nlocal; q++) f[q][0] = f[q][1] = f[q][2] = 0.0;

  bool check = false;

//  cout << "--------------------------------------------------------------------- STATUS 88: " << ntimestep << " : "  << info->lmp->atom->map(88) << endl;

  for(int m = 0; m < Nmax_bacteria; m++){
    if(Bacterium[m].exist > 0){
//      if( (ntimestep >= 247460) && (m == 10)) check =true; else check = false;
      if(check) {cout << "START : " << ntimestep << "	" << m << endl;}
      // get info about the versor of the molecule
      ns = Bacterium[m].nstart; // gloab id
      ne = Bacterium[m].nend;   // global id
      is = info->lmp->atom->map(ns); // local id
      ie = info->lmp->atom->map(ne);   // local id
      norm2 = 0;
      if(check) {cout << "NS NE IS IE : " << ns << "	" << ne << "	" << is << "	" << ie << endl;}
      for(int k = 0; k < 3; k++){
        dr[k] = x[ie][k]-x[is][k]; if(dr[k] > 0.5*Lbox) dr[k] -= Lbox; if(dr[k] < -0.5*Lbox) dr[k] += Lbox;
        norm2 += dr[k]*dr[k];
      }
      norm = sqrt(norm2);
      if(check) {cout << "norm state : " << norm << "	" << Bacterium[m].state << endl;}
      if(Bacterium[m].state == 0){ // run
        f[ie][0] = Frun*dr[0]/norm;
        f[ie][1] = Frun*dr[1]/norm;
      }
      if(Bacterium[m].state == 1){ // tumble
        fx = 0.5*Bacterium[m].Torque * dr[1]/norm2;
        fy = -0.5*Bacterium[m].Torque * dr[0]/norm2;
        f[is][0] = fx;
        f[is][1] = fy;
        f[ie][0] = -fx;
        f[ie][1] = -fy;
      }
      if(Bacterium[m].state ==2){
        f[ie][0] = 0;
        f[ie][1] = 0;
      }

      if(check) {cout << "END : " << ntimestep << "	" << m << endl;}
    }
  }
  total_time = total_time+dt_integration;
  if(check) {cout << "END : " << ntimestep << endl;}
//  cout << "--------------------------------------------------------------------- FINAL 88: " << ntimestep << " : "  << info->lmp->atom->map(88) << endl;
}

void changestate(LAMMPS* lmp){
  int ns, ne;
  char run[1024];

  sprintf(run,"group mobile union mobile stationary"); lmp->input->one(run);
  for(int m = 0; m < Nmax_bacteria; m++){
  if(total_time > Bacterium[m].t_transition){ // need to change state
    if(Bacterium[m].state == 0){
      if(Xrandom()<(t_tumble/(t_tumble+t_immobile))){
        Bacterium[m].state = 1;
        Bacterium[m].t_transition = total_time -1.0/t_tumble*log(Xrandom());
      }else{
        Bacterium[m].state = 2;
        Bacterium[m].t_transition = total_time -1.0/t_immobile*log(Xrandom());
        ns = Bacterium[m].nstart;
        ne = Bacterium[m].nend;
        sprintf(run,"group stationary id %d:%d", ns, ne); lmp->input->one(run);
      }
    }else if(Bacterium[m].state == 1){
      if(Xrandom()<(t_run/(t_run+t_immobile))){
        Bacterium[m].state = 0;
        Bacterium[m].t_transition = total_time -1.0/t_run*log(Xrandom());
        if(Xrandom() < 0.5) Bacterium[m].Torque = Torque;
        else Bacterium[m].Torque = -Torque;
      }else{
        Bacterium[m].state = 2;
        Bacterium[m].t_transition = total_time -1.0/t_immobile*log(Xrandom());
        ns = Bacterium[m].nstart;
        ne = Bacterium[m].nend;
        sprintf(run,"group stationary id %d:%d", ns, ne); lmp->input->one(run);
      }
    }else if(Bacterium[m].state == 2){
      if(Xrandom()<(t_run/(t_run+t_tumble))){
        Bacterium[m].state = 0;
        Bacterium[m].t_transition = total_time -1.0/t_run*log(Xrandom());
      }else{
        Bacterium[m].state = 1;
        Bacterium[m].t_transition = total_time -1.0/t_tumble*log(Xrandom());
      }
      ns = Bacterium[m].nstart;
      ne = Bacterium[m].nend;
      sprintf(run,"group moveagain id %d:%d", ns, ne); lmp->input->one(run);
      sprintf(run,"group stationary subtract stationary moveagain"); lmp-> input->one(run);
      sprintf(run,"group moveagain delete"); lmp->input->one(run);
    }
  }
}
  sprintf(run,"group mobile subtract mobile stationary"); lmp-> input->one(run);
}

void run_and_tumble(LAMMPS* lmp){
  for(int m=0; m< Nmax_bacteria; m++){
    if(Bacterium[m].exist > 0){
      if(Bacterium[m].state == 0){ // the bacterium is running
         Bacterium[m].t_transition = total_time - 1.0/t_run*log(Xrandom());
      }else if(Bacterium[m].state == 1){
         Bacterium[m].t_transition = total_time - 1.0/t_tumble*log(Xrandom());
         if(Xrandom() < 0.5) Bacterium[m].Torque = Torque;
         else Bacterium[m].Torque = -Torque;
      }else if(Bacterium[m].state == 2){
        Bacterium[m].t_transition = total_time -1.0/t_immobile*log(Xrandom());
      }
    }
  }
}

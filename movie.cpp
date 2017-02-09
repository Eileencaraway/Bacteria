/***************************************************************************
                          main.cpp  -  description
                             -------------------
    begin                : Wed Nov 13 10:46:51 CST 2002
    copyright            : (C) 2002 by Massimo Pica Ciamarra
    email                : picaciam@na.infn.it
 ***************************************************************************/

/***************************************************************************
 *                                                                          *
 *   This program is free software; you can redistribute it and/or modify   *
 *   it under the terms of the GNU General Public License as published by   *
 *   the Free Software Foundation; either version 2 of the License, or      *
 *   (at your option) any later version.                                    *
 *                                                                          *
 ***************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include "xrand.h"

using namespace std;

int  images();

int every = 1;
int prepare_plot(float L, float D, int NumAtomMol, float total_time);
int extract_timestep(char buffer[100]);
int const nbin_visit = 1000;
float visit[nbin_visit][nbin_visit];
void print_visit(int timestep);
bool psl = false; 

int main(){
  char buffer[100];
  system("rm plot*png");
  images();
  system("./make_movie.x");
  sprintf(buffer,"mv output.avi movie.avi");
  system(buffer);
  return 0;
}


int images(){
  for(int i = 0; i < nbin_visit; i++){
    for(int j = 0; j < nbin_visit; j++) visit[i][j] = 0;
  }
  char buffer[100];
  char buffer2[100];  
  char buffer3[100];  
  sprintf(buffer,"ls ../dumps/dump* > name_plot.txt"); system(buffer); 
  sprintf(buffer,"ls ../dumps/s_dump* > name_plot2.txt"); system(buffer); 
  int n = 0;
  int NP, NumAtomMol, natoms, natoms2;
  float L,D, dt_integration, total_time;
  int q = 0;
  int num_read;
  system("cat ../in.setup | grep \"Lbox\" | grep \"variable\" | awk '{print $4}' > out");
  system("cat ../in.setup | grep -m 1 \"NumAtom\" | grep \"variable\" | awk '{print $4}' >> out");
  system("cat ../in.setup | grep -m 1 \"Diameter\" | grep \"variable\" | awk '{print $4}' >> out");
  system("cat ../in.setup | grep -m 1 \"dt_integration\" | grep \"variable\" | awk '{print $4}' >> out");
  FILE* in = fopen("out","r");
  fscanf(in,"%f",&L);
  fscanf(in,"%d",&NumAtomMol);
  fscanf(in,"%f",&D);
  fscanf(in,"%f",&dt_integration);
  fclose(in);
//  cout << L << "	" << NumAtomMol << "	" << D << "	" << dt_integration << endl; getchar();
  D /= L;
//  cout << "D " << D << endl; getchar();
  FILE* in_name = fopen("name_plot.txt","r");
  FILE* in_name2 = fopen("name_plot2.txt","r");
  ofstream outgrow("growth.dat");
  int timestep;
  while(fscanf(in_name,"%s",buffer) > 0){
    timestep = extract_timestep(buffer);
    total_time = dt_integration*timestep;
    fscanf(in_name2,"%s",buffer3);
    q++;
    if( (q == every) ){
      // read number of atoms
      sprintf(buffer2,"head -n 4 %s | tail -n 1 > out",buffer); system(buffer2);
      sprintf(buffer2,"head -n 4 %s | tail -n 1 >> out",buffer3); system(buffer2);
      in = fopen("out","r");
      fscanf(in,"%d",&natoms);
      fscanf(in,"%d",&natoms2);
      fclose(in);
      sprintf(buffer2,"tail -n %d %s > out",natoms,buffer); system(buffer2);
      sprintf(buffer2,"tail -n %d %s >> out",natoms2,buffer3); system(buffer2);
      num_read = prepare_plot(L,D,NumAtomMol,total_time);
      outgrow << total_time << "	" << num_read << "	" << natoms+natoms2 << endl;
  /*
      if(psl) system("gnuplot plot_plot_psl.plg");
      else system("gnuplot plot_plot.plg");
      system("cp plot.png final.png");
      if(n < 10000) sprintf(buffer,"mv plot.png plot_0%d.png",n);
      if(n < 1000) sprintf(buffer,"mv plot.png plot_00%d.png",n);
      if(n < 100) sprintf(buffer,"mv plot.png plot_000%d.png",n);
      if(n < 10) sprintf(buffer,"mv plot.png plot_0000%d.png",n);
      system(buffer);
*/
      n++;
      if(psl) print_visit(timestep);
      q = 0;
    }
  }
  outgrow.close();
  fclose(in_name);
  fclose(in_name2);
  return 0;
}

int prepare_plot(float L, float D, int NumAtomMol, float total_time){
  ofstream outfile("input_plot.txt");
  float x0,y0,n,x,y, type, x1, y1;
  float dx,dy, DX, DY;
  FILE* in; 

  int m = 1;
  char color[150][150];
  int ncolor = 0;
  in = fopen("colors.txt","r");
  while(fscanf(in,"%s",color[ncolor]) > 0) ncolor++;
  fclose(in);
  init_random(1,1);
//  sprintf(color[ncolor],"red"); ncolor++;
//  sprintf(color[ncolor],"blue"); ncolor++;
//  sprintf(color[ncolor],"cyan"); ncolor++;
//  sprintf(color[ncolor],"magenta"); ncolor++;
//  sprintf(color[ncolor],"yellow4"); ncolor++;

  int hours, minutes, seconds;
  seconds = total_time;
  hours = int(seconds/(60*60));
  seconds -= 60*60*hours;
  minutes = int(seconds/60);
  seconds -= 60*minutes;

  outfile << "set label 1 sprintf(\"%02d:%02d:%02d\"," << hours << "," << minutes << "," << seconds << ") at 1.02,0.9 font \"arial,14\"" << endl;
  outfile << "set style fill noborder" << endl;
  in = fopen("out","r");
  
  int nc = 0;
  int error = 0;
  int read = 0;
  int ibin, jbin;
  while(fscanf(in,"%f %f %f %f",&x, &y,&n, &type) > 0){
    x += 0.5*L; if(x > L) x -= L; y += 0.5*L; if(y > L) y -= L; //////////////// shift position
    nc = n; while(nc >= ncolor) nc -= ncolor;
    x0 = x; y0 = y;
    x /= L; y /= L; 
    outfile << "set object " << m << " circle at " << x << "," << y << " size " << D/2 << endl;
    outfile << "set object " << m << " fc rgb \"" << color[nc] << "\" fillstyle solid 1.0 noborder front" << endl;
    m++;
    for(int i = 0; i < NumAtomMol-1; i++){
      fscanf(in,"%f %f %f %f",&x, &y,&n, &type);
      x += 0.5*L; if(x > L) x -= L; y += 0.5*L; if(y > L) y -= L; ////////////////////////// shift position
      x1 = x, y1 = y;
      x /= L; y /= L;
      outfile << "set object " << m << " circle at " << x << "," << y << " size " << D/2 << endl;
      outfile << "set object " << m << " fc rgb \"" << color[nc] << "\" fillstyle solid 1.0 noborder front" << endl;
      m++;
      if(i == (1+int(NumAtomMol/2))){
        ibin = int(x*nbin_visit); if(ibin >= nbin_visit) ibin = nbin_visit-1; if(ibin < 0) ibin = 0;
        jbin = int(y*nbin_visit); if(jbin >= nbin_visit) jbin = nbin_visit-1; if(jbin < 0) jbin = 0;
        visit[ibin][jbin]++;
      } 
    }
    x = fabs(x1-x0); if(x > 0.5*L) x = x-L/2;
    y = fabs(y1-y0); if(y > 0.5*L) y = y-L/2;
    read++;
  }
  fclose(in);
  outfile.close();
  return read;
}

int extract_timestep(char buffer[100]){
  char buffer2[100];
  int k0 = 0;
  int n0 = 0;
  while(k0 != 3){ if(buffer[n0] == '.') k0++; n0++;}
  int n1 = n0+1;
  while(k0 != 4){ if(buffer[n1] == '.') k0++; n1++;}
  for(int s = n0; s < n1; s++){
    buffer2[s-n0] = buffer[s];
  }
  buffer2[n1-n0-1] = '\0';
  
//  cout << buffer << "	" << buffer2 << "	" << atoi(buffer2) << endl;
//  getchar();
  return atoi(buffer2);
}

void print_visit(int timestep){
  char buffer[100];
  sprintf(buffer,"visit/visit_%010d.dat", timestep);
  ofstream outfile(buffer);
  for(int i = 0; i < nbin_visit; i++){
    for(int j = 0; j < nbin_visit; j++){
      if(visit[i][j] > 0){
        outfile << (i*1.0)/nbin_visit << "	" << (j*1.0)/nbin_visit << endl;
        outfile << (1+i*1.0)/nbin_visit << "	" << (j*1.0)/nbin_visit << endl;
        outfile << (1+i*1.0)/nbin_visit << "	" << (1+j*1.0)/nbin_visit << endl;
        outfile << (i*1.0)/nbin_visit << "	" << (1+j*1.0)/nbin_visit << endl;
        outfile << (i*1.0)/nbin_visit << "	" << (j*1.0)/nbin_visit << endl;
        outfile << endl;
      }
    }
  }
  outfile.close();
  sprintf(buffer,"cp visit/visit_%010d.dat visit.dat",timestep);
  system(buffer);
}

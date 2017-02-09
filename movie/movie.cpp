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
int prepare_plot(float L, float D, int NumAtomMol);

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
  char buffer[100];
  char buffer2[100];
  char buffer3[100];
  sprintf(buffer,"ls ../dumps/dump* > name_plot.txt"); system(buffer);
  sprintf(buffer,"ls ../dumps/s_dump* > name_plot2.txt"); system(buffer);
  int n = 0;
  int NP, NumAtomMol, natoms, natoms2;
  float L,D;
  int q = 0;
  int num_read;
  system("cat ../in.setup | grep \"Lbox\" | grep \"variable\" | awk '{print $4}' > out");
  system("cat ../in.setup | grep \"NumAtom\" | grep \"variable\" | awk '{print $4}' >> out");

  getchar();
  FILE* in = fopen("out","r");
  fscanf(in,"%f",&L);
  fscanf(in,"%d",&NumAtomMol);
  fclose(in);
  D = 1;
  D /= L;
  FILE* in_name = fopen("name_plot.txt","r");
  FILE* in_name2 = fopen("name_plot2.txt","r");
  ofstream outgrow("growth.dat");
  while(fscanf(in_name,"%s",buffer) > 0){
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
      num_read = prepare_plot(L,D,NumAtomMol);
      outgrow << n << "	" << num_read << "	" << natoms+natoms2 << endl;
      system("gnuplot plot_plot.plg");
      system("cp plot.png final.png");
      if(n < 10000) sprintf(buffer,"mv plot.png plot_0%d.png",n);
      if(n < 1000) sprintf(buffer,"mv plot.png plot_00%d.png",n);
      if(n < 100) sprintf(buffer,"mv plot.png plot_000%d.png",n);
      if(n < 10) sprintf(buffer,"mv plot.png plot_0000%d.png",n);
      system(buffer);
      n++;
      q = 0;
    }
  }
  outgrow.close();
  fclose(in_name);
  fclose(in_name2);
  return 0;
}

int prepare_plot(float L, float D, int NumAtomMol){
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
  outfile << "set style fill noborder" << endl;
  in = fopen("out","r");

  int nc = 0;
  int error = 0;
  int read = 0;
  while(fscanf(in,"%f %f %f %f",&x, &y,&n, &type) > 0){
    x += 0.5*L; if(x > L) x -= L; y += 0.5*L; if(y > L) y -= L; //////////////// shift position
    nc = n; while(nc >= ncolor) nc -= ncolor;
    x0 = x; y0 = y;
    x /= L; y /= L;
    outfile << "set object " << m << " circle at " << x << "," << y << " size " << D/2 << endl;
    outfile << "set object " << m << " fc rgb \"" << color[nc] << "\" fillstyle solid 1.0 noborder" << endl;
    m++;
    for(int i = 0; i < NumAtomMol-1; i++){
      fscanf(in,"%f %f %f %f",&x, &y,&n, &type);
      x += 0.5*L; if(x > L) x -= L; y += 0.5*L; if(y > L) y -= L; ////////////////////////// shift position
      x1 = x, y1 = y;
      x /= L; y /= L;
      outfile << "set object " << m << " circle at " << x << "," << y << " size " << D/2 << endl;
      outfile << "set object " << m << " fc rgb \"" << color[nc] << "\" fillstyle solid 1.0 noborder" << endl;
      m++;
    }
    x = fabs(x1-x0); if(x > 0.5*L) x = x-L/2;
    y = fabs(y1-y0); if(y > 0.5*L) y = y-L/2;
    read++;
  }
  fclose(in);
  outfile.close();
  return read;
}

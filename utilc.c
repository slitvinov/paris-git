#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "inc.h"

extern void __module_solids_MOD_outfarray(int *carray);

/* check fortran indexes. 
   call comparef2c() in calling fortran90 program to check */

void comparef2c_()
{
  int carray[4];
  int nx=2,ny=2,nz=2;
#define ARRAY(i,j,k) (*(carray + (i) - 1 + ((j)-1) * nx + ((k)-1) * nx * ny))
#define CARRAY(i,j) (ARRAY(i,j,1))
  CARRAY(1,1) = 11;
  CARRAY(1,2) = 12;
  CARRAY(2,1) = 21;
  CARRAY(2,2) = 22;
  __module_solids_MOD_outfarray(carray);
}

/* sphere definition function */


double sphere_func(double x, double y, double z, double x0, double y0, double z0, double radius)
{
  return -((x - x0)*(x - x0) + (y - y0)*(y - y0) + (z - z0)*(z - z0)) + radius*radius;
}

/* example implicit solid definition function */

double solid_func_one_sphere(double x, double y, double z,double boxL)
{
  double x0=0.,y0=0.,z0=0.;
  double radius = 0.25*boxL;
  return sphere_func(x,y,z,x0,y0,z0,radius);
}

/* One basic CFC cell: one vertex + three faces */


double solid_func_CFC_scaled(double x, double y, double z)
{
  double radius = 0.25*sqrt(2.);
  /* 1 lower vertex and x=0 face */
  double a = MAX(sphere_func(x,y,z,0.,0.,0.,radius),sphere_func(x,y,z,0.,0.5,0.5,radius));
  double b = MAX(sphere_func(x,y,z,0.5,0.,0.5,radius),sphere_func(x,y,z,0.5,0.5,0.,radius));
   return MAX(a,b);
}

/* One basic cube with CFC array and all vertices and faces */

double solid_func_CFC_(double x, double y, double z, double boxL)
{
  /* rescale by boxL: no need */
  x /= boxL;   y /= boxL;   z /= boxL;
  //  printf("HAHAHAHA %g %g %g  %g \n",x,y,z,boxL);
  // exit(1);
  /* shift to lower left corner: no need ? (add shift variable later ? ) 
  x += 0.5; y += 0.5; z += 0.5; */
  /* cells at root vertex and three first neighbors */
  double a = MAX(solid_func_CFC_scaled(x,y,z),solid_func_CFC_scaled(x-1.,y,z));
  double b = MAX(solid_func_CFC_scaled(x,y-1.,z),solid_func_CFC_scaled(x,y,z-1.));
  a = MAX(a,b);
  /* three next lower vertices */
  double c = MAX(solid_func_CFC_scaled(x-1.,y-1.,z-1.),solid_func_CFC_scaled(x-1.,y-1.,z));
  b = MAX(solid_func_CFC_scaled(x,y-1.,z-1.),solid_func_CFC_scaled(x-1.,y,z-1.));
  a = MAX(c,a);
  return MAX(b,a);
}

FILE * fd = NULL;
int np;

void make_visit_file_()
{
  int err;
  if((err =  MPI_Comm_size(MPI_COMM_WORLD,&np) ) != MPI_SUCCESS ) 
    {
      fprintf(stderr,"MPI error %d, aborting\n",err);
      exit(1);
    }
  fd = fopen("parallel.visit","w");
  fprintf(fd,"!NBLOCKS %d\n",np);
}

void append_visit_file_(char * rootname, int * padding)
{
#define STOPCHAR '-'

  int pprank, np, err;
  
  if((err =  MPI_Comm_rank(MPI_COMM_WORLD,&pprank) ) != MPI_SUCCESS ) 
    {
      fprintf(stderr,"MPI error %d, aborting\n",err);
      exit(1);
    }
  if(pprank != 0) return;
  if((err =  MPI_Comm_size(MPI_COMM_WORLD,&np) ) != MPI_SUCCESS ) 
    {
      fprintf(stderr,"MPI error %d, aborting\n",err);
      exit(1);
    }
  if( fd == NULL)  make_visit_file_();

  /* Dirty trick to terminate fortran-generated string */
  char * prootname = rootname;
  while( *(prootname++) != STOPCHAR );
  *prootname = '\0';

  for(pprank=0;pprank<np;pprank++)
    fprintf(fd,"%s%0*d.vtk\n",rootname,*padding,pprank);  // 3 is the padding in ftc
  fflush(fd);
}
 
void close_visit_file_()
{
  fclose(fd);
}



	   

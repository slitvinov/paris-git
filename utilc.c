#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ftcinc.h"

extern void outfarray_(int *carray);
extern void __module_solids_MOD_outfarray(int *carray);
/* check fortran indexes. 
   call comparef2c() in calling fortran program to check */

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

double solid_func_CFC(double x, double y, double z, double boxL)
{
  /* rescale by boxL */
  x /= boxL;   y /= boxL;   z /= boxL;
  /* shift to lower left corner */
  x += 0.5; y += 0.5; z += 0.5;
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

/* initialize one solid */

void inisolids_(double *solids,int *nx, int *ny, int *nz,double boxL)
{
#define SOLIDS(i,j,k) (*(solids + (i) - 1 + ((j)-1) * *nx + ((k)-1) * *nx * *ny))

  /* dx uses boxL but h elsewhere in the code is 1/(nx-2) */

  double dx = boxL/((double) *nx-2.0);  

  int i,j,k;
  for(i=1;i<=*nx;i++) 
    for(j=1;j<=*ny;j++) 
      for(k=1;k<=*nz;k++) 
	if(solid_func_CFC (COORD(i,*nx,dx),COORD(j,*ny,dx),COORD(k,*nz,dx),boxL)>0.)
	    SOLIDS(i,j,k)=1.;
}

  

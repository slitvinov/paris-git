#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define ABS(a)   ((a) >  0  ? (a) : -(a))

#define MAXLINES 100000

int main (int argc, char * argv[])
{
  if(argc < 3) 
    {
      printf("%s : error: not enough command line arguments.\n\n"
             "Usage: %s n, phi \n"
             "n is the number of spheres in a volume of unspecified size.\n"
	     "(the dimensionless result should be independent of volume size.)\n"
	     "phi is the porosity\n"
	     "Example %s 1600 0.19\n"
	     "\n",argv[0],argv[0],argv[0]);
      exit(1);
    }
  double size=64e0;  
  double K, s;
  float n,phi;
  double pi=M_PI;
  double R=0.0625*size;  /* gives R=4 as in Cancelliere et al. */
  sscanf(argv[1],"%g",&n); 
  sscanf(argv[2],"%g",&phi); 
  n = n / (size*size*size);
  s = 4*pi*R*R*n*phi;
  K = phi*phi*phi/(6.0*s*s);
  printf("Dimensionless Kozeny-Carman estimate K/R^2 = %g\n",K/(R*R));
  //  printf("%g %g %g  %g\n",pi,n,s,phi);

  double phipred;
  phipred = exp ( - 4*pi*R*R*R*n/3. );
  double npred = - log(phi)*(3./(4*pi*R*R*R));
    printf("predicted phi = %g predicted n*volume = %g \n",phipred,npred*size*size*size);
}



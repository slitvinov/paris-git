#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define ABS(a)   ((a) >  0  ? (a) : -(a))

#define MAXLINES 100000

int main (int argc, char * argv[])
{
  if(argc < 1) 
    {
      printf("%s : error: not enough command line arguments.\n\n"
             "Usage: %s FILE \n"
             "FILE should contain two columns.\n"
             "the log of the time derivative is commputed.\n"
	     "\n",argv[0],argv[0]);
      exit(1);
    }
  FILE * fd1;
  float * t = ( float *) malloc(MAXLINES*sizeof(float));
  float * u = ( float *) malloc(MAXLINES*sizeof(float));
  int nlines=0;
  if((fd1 = fopen(argv[1],"r")) == NULL) 
    {
      fprintf(stderr,"argv[0]: error: could not open file %s\n",argv[1]);
      exit(3);
    }
  int returnscan1=0;
  double du=0.;
  double small=1e-7;
  int delay=10;
  while((returnscan1 = fscanf(fd1,"%g %g",t,u)) != EOF)
    { 
      nlines++;
      if(returnscan1 < 2) 
	{ 
	  fprintf(stderr,"%s: error on line %d counted %d items\n",argv[0],nlines,returnscan1);
	  exit(2); 
	}
      //      printf("%g %g %g\n",*t,dum,*u);
      if(nlines>delay)
	{
	  du = ABS( (*u - *(u-delay))/(*t - *(t-delay)));
	  du = MAX(du,small);
	  du = log(du);
	  printf("%g %g\n",*t,du);
	}
      t++;u++;
    }
      return 0;
}



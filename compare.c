#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAXLINES 10000

int main (int argc, char * argv[])
{
  if(argc < 4) 
    {
      printf("compare: error: not enough command line arguments.\n\n"
             "Usage: compare FILE1 FILE2 TOLERANCE [RELATIVE=0] [OUTPUT=0]\n"
             "FILE1 and FILE2 should contain two columns each.\n"
             "The norm of the difference of the second columns is computed.\n"
	     "If RELATIVE=1 the relative error is computed.\n"
	     "If OUTPUT=1 the relative error is printed even if the test passes.\n"
	     "\n");
      exit(1);
    }
  float tolerance;
  sscanf(argv[3],"%g",&tolerance);
  int relative=0;
  if(argc >= 5)  sscanf(argv[4],"%d",&relative);
  int output=0;
  if(argc >= 6)  sscanf(argv[5],"%d",&output);

  FILE * fd1;FILE * fd2;
  float * x1 = ( float *) malloc(MAXLINES*sizeof(float));
  float * x2 = ( float *) malloc(MAXLINES*sizeof(float));
  float * y1 = ( float *) malloc(MAXLINES*sizeof(float));
  float * y2 = ( float *) malloc(MAXLINES*sizeof(float));
  int nlines=0;
  if((fd1 = fopen(argv[1],"r")) == NULL) 
    {
      fprintf(stderr,"argv[0]: error: could not open file %s\n",argv[1]);
      exit(3);
    }
  if((fd2 = fopen(argv[2],"r")) == NULL) 
    {
      fprintf(stderr,"%s: error: could not open file %s\n",argv[0],argv[2]);
      exit(3);
    }
  int returnscan1=0;
  int returnscan2=0;
  double diff=0.;
  double scaling;
  while((returnscan1 = fscanf(fd1,"%g %g",x1,y1)) != EOF && nlines < MAXLINES && (returnscan2 = fscanf(fd2,"%g %g",x2,y2)) )
    { 
      nlines++;
      if(returnscan1 < 2 ) 
	{ 
	  fprintf(stderr,"In directory "); system("pwd");
	  fprintf(stderr,"%s: error on line %d \n"
		  "counted only %d item instead of at least 2 in file %s\n", 
		  argv[0],nlines,returnscan1,argv[1]);
	  exit(2); 
	}
      if(returnscan2 < 2) 
	{ 
	  fprintf(stderr,"In directory "); system("pwd");
	  fprintf(stderr,"%s: error on line %d \n"
		  "counted only %d item instead of at least 2 in file %s\n", 
		  argv[0],nlines,returnscan2,argv[2]);
	  exit(2); 
	}
      //      printf("%g %g %g %g\n",*x1,*y1,*x2,*y2);
      scaling = (1. - relative) + relative * *y2 * *y2;
      diff += (*y1-*y2)*(*y1-*y2)/scaling;
      x1++;
      y1++;
      x2++;
      y2++;
    }
  diff = sqrt(diff/nlines);
  if(diff < (double) tolerance) 
    {
      if(output < 2) 
	{
	  printf("\033[32;1m PASS\033[0m");
	  if(output) printf(" L2 relative error norm = %g",diff);
	  printf("\n");
	}
      else
	printf("second: %g\n",diff);
    }
  else
    {
      printf("\033[31;1m FAIL\033[0m error=%g \n",diff);
    }
  return 0;
}



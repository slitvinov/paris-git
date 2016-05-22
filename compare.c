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
             "The L2 AND L_\\infty norms of the difference of the second columns are computed.\n"
	     "If RELATIVE=1 the relative error is computed.\n"
	     "If OUTPUT=1 the relative error is printed even if the test passes.\n"
	     "The relative error is defined as || y_1 - y_2 ||/|| y_2 ||\n"
	     "\n");
      exit(1);
    }
  double tolerance;
  sscanf(argv[3],"%lf",&tolerance);
  int relative=0;
  if(argc >= 5)  sscanf(argv[4],"%d",&relative);
  int output=0;
  if(argc >= 6)  sscanf(argv[5],"%d",&output);

  double errmax = 0.;
  double  error, error2, y2norm2=0.;
  
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
  int ierr;
  while((returnscan1 = fscanf(fd1,"%g %g",x1,y1)) != EOF && nlines < MAXLINES && (returnscan2 = fscanf(fd2,"%g %g",x2,y2)) )
    { 
      nlines++;
      if(returnscan1 < 2 ) 
	{ 
	  fprintf(stderr,"In directory "); ierr=system("pwd");
	  if(ierr != 0) exit(3);
	  fprintf(stderr,"%s: error on line %d \n"
		  "counted only %d item instead of at least 2 in file %s\n", 
		  argv[0],nlines,returnscan1,argv[1]);
	  exit(2); 
	}
      if(returnscan2 < 2) 
	{ 
	  fprintf(stderr,"In directory "); ierr=system("pwd");
	  if(ierr != 0) exit(3);
	  fprintf(stderr,"%s: error on line %d \n"
		  "counted only %d item instead of at least 2 in file %s\n", 
		  argv[0],nlines,returnscan2,argv[2]);
	  exit(2); 
	}
      //      printf("%g %g %g %g\n",*x1,*y1,*x2,*y2);
      y2norm2 +=  *y2 * *y2;
      error2 = (*y1-*y2)*(*y1-*y2);
      error = sqrt(error2);
      if(error > errmax) errmax = error;
      diff += error2;
      x1++;
      y1++;
      x2++;
      y2++;
    }
  // compute ||y2||_2  (L2 norm of second column) 
  y2norm2 /= nlines;
  scaling = (1. - relative) + relative * y2norm2;
  // rescale for relative errors
  diff /= scaling;
  errmax /= sqrt(scaling);
  // get L2 of error
  diff = sqrt(diff/nlines);
  if(diff < (double) tolerance) 
    {
      if(output < 2) 
	{
	  printf("\033[32;1m PASS\033[0m");
	  if(output)
	    {
	      if(relative==1) 
		printf(" L2 relative error norm = %g  L_\\infty relative error norm = %g",diff,errmax);
	      else
		printf(" L2 error norm = %g  L_\\infty error norm = %g",diff,errmax);
	    }
	  printf("\n");
	}
      else
	printf(" %g %g\n",diff,errmax);
    }
  else
    {
      if(output < 2) 
	{
	  if(relative==1)
	    printf("\033[31;1m FAIL\033[0m L2 relative error norm = %g L_\\infty relative error norm = %g \n",diff,errmax);
	  else
	    printf("\033[31;1m FAIL\033[0m L2 error norm = %g L_\\infty error norm = %g \n",diff,errmax);
	}
      else
	  printf("%g %g \n",diff,errmax);
    }
  return 0;
}



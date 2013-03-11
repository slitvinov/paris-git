#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAXLINES 100
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define ABS(a) ((a) > 0. ? (a) : (-(a)))

float minflt_(void);

int main (int argc, char * argv[])
{
  if(argc < 4) 
    {
      printf("%s: error: not enough command line argument\n",argv[0]);
      exit(1);
    }
    //    strcpy(filename[1],argv[1]);
  //  double tolerance=sscanf(argv[3],"%d",&tolerance);
  float tolerance;
  sscanf(argv[3],"%g",&tolerance);
  //  printf(" %g ",tolerance);
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
      fprintf(stderr,"argv[0]: error: could not open file %s\n",argv[2]);
      exit(3);
    }
  int returnscan1=0;
  int returnscan2=0;
  double diff=0.;
  while((returnscan1 = fscanf(fd1,"%g %g",x1,y1)) != EOF && nlines < MAXLINES && (returnscan2 = fscanf(fd2,"%g %g",x2,y2)) )
    { 
      nlines++;
      if(returnscan1 < 2 || returnscan2 < 2) 
	{ 
	  fprintf(stderr,"argv[0]: error on line %d counted %d %d items\n",nlines,returnscan1,returnscan2);
	  exit(2); 
	}
      //      printf("%g %g %g %g\n",*x1,*y1,*x2,*y2);
      diff += (*y1-*y2)*(*y1-*y2);
      x1++;
      y1++;
      x2++;
      y2++;
    }
  diff = sqrt(diff/nlines);
  if(diff < (double) tolerance) 
    {
      printf("PASS\n");
    }
  else
    {
      printf("FAIL error=%g \n",diff);
    }
  exit(0);
}



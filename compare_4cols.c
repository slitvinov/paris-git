/* Program to compare two files of discrete data in 3D. 
   Files has to be in format: x y z var. */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAXLINES 98777

void write_diffs_file(double x, double y, double z, double out, int read_flag, char param[20]);

int main(int argc, char * argv[]) // argc: argument count, # of arguments read from command line. argv: argument vector
{
  if(argc < 4) 
    {
      printf("compare: error: not enough command line arguments.\n\n"
             "Usage: pariscompare4 PAR_HOLDING_FILE SERIAL_FILE TOLERANCE VAR_DISPLAY\n"
             "PAR_HOLDING_FILE contains a list of the parallel files to be compared. The first line should be NPROCS N, N an int.\n"
	     "SERIAL_FILE and those referenced in PAR_HOLDING_FILE should contain four columns each, x, y, z, var.\n"
	     "TOLERANCE specifies the maximum difference between individual coordinate reference points when comparing two data files.\n"
	     "VAR_DISPLAY gives the name, max 10 characters, that will be displayed for the variable in ViSit as well as the 3D file name.\n"
	     "\n");
      exit(1);
    }
  FILE *par_holder; FILE *mono;
  int read_holder_head=0;
  int n_pfiles=0;
  double cell_tolerance;
  char var[20];
  double error2=0., errmax=0.;
  int total_found=0;
  par_holder = fopen(argv[1],"r");
  if(par_holder == NULL)
    {
      fprintf(stderr,"Reading error: could not open file %s\n",argv[1]);
      exit(1);
    }
  else
    { 
      printf("Reading parallel holding file: %s \n",argv[1]);
    }
  mono = fopen(argv[2],"r");
  if(mono == NULL)
    {
      fprintf(stderr,"Reading error: could not open serial file %s\n",argv[2]);
      exit(1);
    }
  else
    { 
      printf("Reading serial file: %s \n",argv[2]);
    }
  read_holder_head = fscanf(par_holder,"NPROCS %d",&n_pfiles);
  if (read_holder_head == EOF) {
    printf("Error in reading holder header");
    exit(1);
  }
  else {
    printf("Number of blocks written at top: %d \n",n_pfiles);
  }
  sscanf(argv[3],"%lg",&cell_tolerance);
  printf("Tolerance in cell locations: %g\n",cell_tolerance);
  sscanf(argv[4],"%s",var);
  printf("Variable display name: %s\n",var);
  int i;
  int head = 0;
  for (i=1; i<=n_pfiles; i++) {
    FILE *par_file;
    int read_par_file, read_par_line, n_lines;
    int serial_line;
    char filename[50];
    double ref[4];
    read_par_file = fscanf(par_holder,"%s",filename);
    if (read_par_file != EOF) {
      printf("File name read: %s \n",filename);
    }
    else {
      printf("Error in reading file name in holder file");
      exit(1);
    }
    n_lines = 0; 
    par_file = fopen(filename,"r");
    read_par_line = fscanf(par_file,"%lg %lg %lg %lg",&ref[0],&ref[1],&ref[2],&ref[3]);
    if (read_par_line == EOF) {
      printf("Error in reading line from par file: %s\n",filename);
      exit(1);
    }
    rewind(mono); // Starts at first line of single file for every new file  

    while ((read_par_line != EOF) && (n_lines < MAXLINES)) { 
      double comp[4];
      short int found;
  
      n_lines++;
      // printf("%lg %lg %lg %lg \n",ref[0],ref[1],ref[2],ref[3]);
      // Now match line in serial file and calculate difference
      found = 0;
      serial_line = fscanf(mono,"%lg %lg %lg %lg",&comp[0],&comp[1],&comp[2],&comp[3]);
      double ser_start[3] = {comp[0],comp[1],comp[2]};
      
      while ((serial_line != EOF) && (found == 0)) {
      	double diff;
      	//test if found
      	if ((fabs(ref[0]-comp[0])<cell_tolerance) && (fabs(ref[1]-comp[1])<cell_tolerance) && (fabs(ref[2]-comp[2])<cell_tolerance)) {
      	  //calc diff
      	  diff = ref[3]-comp[3];
	  error2+= diff*diff;
	  if (fabs(diff) > errmax) { 
	    errmax = fabs(diff); 
	  }
      	  //write diff to file
      	  found = 1;
	  total_found++;
      	  //call writing file
      	  write_diffs_file(comp[0],comp[1],comp[2],diff,head,var);
      	  head = 1;
      	}
      	else {
      	  serial_line = fscanf(mono,"%lg %lg %lg %lg",&comp[0],&comp[1],&comp[2],&comp[3]);
      	  if ((fabs(ser_start[0]-comp[0])<cell_tolerance) && (fabs(ser_start[1]-comp[1])<cell_tolerance) && (fabs(ser_start[2]-comp[2])<cell_tolerance)) {
      	    printf("Error! Discrete point unmatched\n");
      	    exit(1);
      	  }
      	  if ((serial_line == EOF) && (found == 0)) {
      	    rewind(mono);
      	    serial_line = fscanf(mono,"%lg %lg %lg %lg",&comp[0],&comp[1],&comp[2],&comp[3]);
      	  }
      	}
      } //serial line loop
      
      read_par_line = fscanf(par_file,"%lg %lg %lg %lg",&ref[0],&ref[1],&ref[2],&ref[3]);
    } //parallel line loop
    printf("Number of lines read in file %s: %d \n",filename,n_lines);
    fclose(par_file);
  }
  if (total_found >= 1) {error2 = sqrt(error2/total_found);}
  printf("Total points compared: %10d \n",total_found);
  printf("L2 norm: %14.9f L_inf norm: %14.9f \n",error2,errmax);
  fclose(mono);
  fclose(par_holder);
  return 0;
}

void write_diffs_file(double x, double y, double z, double out, int read_flag, char param[20]) {
  FILE *diff_file;
  if (read_flag==0) {
    diff_file = fopen("diff.3D","w");
    fprintf(diff_file,"X Y Z %s\n",param);
    fprintf(diff_file,"%14.9f %14.9f %14.9f %14.9f\n",x,y,z,out);
  }
  else {
    diff_file = fopen("diff.3D","a");
    fprintf(diff_file,"%14.9f %14.9f %14.9f %14.9f\n",x,y,z,out);
  }
  fclose(diff_file);
}

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define COLORED 2
#define IS_COLOREDP(x) ( *(x) == COLORED ) 
#define max_size_for_print 6

typedef enum
{
  LEFT = 0,
  RIGHT,
  BOTTOM,
  TOP,
  BACK,
  FRONT
} Directions;

int * bigrock;

int cubesize;
int size;
int max;

struct point {
  int xcoord[3];
  int * paddress;
  };


typedef struct point point_t;

FILE * fd;

int  * makept                  (point_t * currentp,int x,int y,int z);
int    pt_is_on_any_face       (point_t * currentp);
int    any_neighbor_is_colored (point_t * currentp);
int    pt_is_in_face           (point_t * currentp,Directions face);
int  * get_neighbor_address    (point_t * currentp,Directions dir);
void   printrock               (int n);
void   error_exit              (char * message);
double read_data               ();
int    max_X                   (double * fraction);
void   ghost_face              (Directions face);
void   split_print             (int nx, int npx);

int main (int argc, char * argv[])
{
  //  if(1) printf("1\n"); if(0) printf("0\n"); exit(0);
	  
  if(argc < 2) 
    {
      printf("%s: error: not enough command line arguments.\n\n"
	     "Usage: %s size [nPx] \n"
             "size is the integer length of a cube edge\n"
	     "nPx is the number of processes per direction (optional)\n"
             "various rock statistics are computed\n"
	     "\n",argv[0],argv[0]);
      exit(1);
    }
  sscanf(argv[1],"%d",&cubesize);
  //  printf(" argc = %d, argv[1] = %s\n",argc,argv[1]);
  // exit(1);
  point_t currentp;
  int x,y,z;

  size=cubesize+2;
  max=cubesize+1;

  bigrock = (int *) malloc( size*size*size*sizeof(int));

  //if((fd = fopen("roche.txt", "r")) == NULL)
  //if((fd = fopen("test3cube.txt", "r")) == NULL)
  //    error_exit("Can't open file.\n");

  fd=stdin;
  double totalporosity=read_data();

  if(argc==3)
    {
      int nx=cubesize;
      int npx;
      sscanf(argv[2],"%d",&npx);
      split_print(nx,npx);
    }

  printf("\n\nFile read and written.\nPorosity=%f %% \n\n",totalporosity*100.);
  printrock(0);
  printf("\n------\n");
  x=max-1;
    printf("\n");
  for(y=1;y<max;y++)
    {
    for(z=1;z<max;z++)
      {
	printf("%d",*makept(&currentp,x,y,z));
      }
    printf("\n");
    }
    printf("\n");


  // fill the ghost layers

  Directions face;
  for(face=0;face<6;face++)
    ghost_face(face);

  printf("Periodized:\n");
  printrock(1);
  printf("\n------\n");

  // initialize face by coloring the fluid voxels
  // and print

  Directions flow_dir=LEFT;
  int comp_flow=flow_dir/2;

  for(x=1;x<max;x++)
    for(y=1;y<max;y++)
      for(z=1;z<max;z++)
	{
	  int * pcolor = makept(&currentp,x,y,z);
	  if( pt_is_in_face(&currentp,flow_dir)) 
	      if ( *pcolor ) 
		  *pcolor = COLORED;
	}

  printf("Face initialized:\n");
  printrock(0);
  printf("\n------\n");

  int iter;
  int itermax=size;

  int colcount=0;
  for(iter=0;iter<itermax;iter++)
    {
      colcount=0;
       for(x=1;x<max;x++)
	for(y=1;y<max;y++)
	  for(z=1;z<max;z++)
	    {
	      int * pcolor = makept(&currentp,x,y,z);
	      if(*pcolor != 0) 
		if(any_neighbor_is_colored(&currentp)) 
		  {
		    *pcolor = COLORED;
		    if(pt_is_in_face(&currentp,RIGHT)) 
		      colcount++;
		  }
	    } 
      double fraction;
      int maxi=max_X(&fraction);
      if(size <= max_size_for_print) printf("********\n");
      printf("iter=%d: max_x=%d colored fraction %f%% colored cells on end face %d\n",iter+1,maxi,fraction*100,colcount);
      printrock(0);
       for(face=0;face<6;face++)
	 {
	   int component=face/2;
	   if(component != comp_flow) ghost_face(face);
	 }
     }

  colcount=0;
  face=RIGHT;
  for(x=1;x<max;x++)
    for(y=1;y<max;y++)
      for(z=1;z<max;z++)
 	{
	  int * pcolor = makept(&currentp,x,y,z);
	  if( pt_is_in_face(&currentp,face)) 
	    if ( * pcolor == COLORED) colcount++;
	}

  printf("\n\nNumber of colored points on RIGHT face is %d.\n",colcount);
}

int * makept(point_t * currentp,int x,int y,int z)
{
  int i;
  currentp->xcoord[0] = x;
  currentp->xcoord[1] = y;
  currentp->xcoord[2] = z;
  for(i=0;i<3;i++)
    if(currentp->xcoord[i] < 0 || currentp->xcoord[i] > max) error_exit("out of bounds");
  currentp->paddress = bigrock + x*size*size + y*size + z; 
  return currentp->paddress; 
}

int pt_is_on_any_face(point_t * currentp)
 {
   int is_in=0;
   Directions face;
   for(face=0;face<6;face++)
     {
       is_in += pt_is_in_face(currentp,face);
     }
   return is_in;
 }

int any_neighbor_is_colored(point_t * currentp)
{
  Directions dir;
  int sum_neighbors=0;
  int neigh_color;
  for(dir=0;dir<6;dir++)
    {
      int * p = get_neighbor_address(currentp,dir); 
      if(p == NULL) 
	neigh_color=0;
      else
	neigh_color = IS_COLOREDP(p);
      sum_neighbors += neigh_color;
      //      printf("dir %d neigh %d\n",dir,neigh_color);
    }
  return sum_neighbors;
}

int pt_is_in_face(point_t * currentp,Directions dir)
 {
   int component = dir/2;
   int orientation = dir - component*2;  
   if(!orientation)
     {
       if(currentp->xcoord[component] == 1) return 1;
     }
   else
     {
       if(currentp->xcoord[component] == cubesize) return 1;
     }
   return 0;
 }

int * get_neighbor_address(point_t * currentp,Directions dir)
 {
   int component = dir/2;
   int orientation = dir - component*2;
   int xn[3];
   point_t neighbor_pt;
   int * preturn;

   xn[0] = currentp->xcoord[0];
   xn[1] = currentp->xcoord[1];
   xn[2] = currentp->xcoord[2];
   
   if(!orientation)
     {
       xn[component]--;
     }
   else
     {
       xn[component]++;
     }
     return makept(&neighbor_pt,xn[0],xn[1],xn[2]);
 }
 
void printrock(int m)
{
  int x, y, z;
  point_t currentp;
  int n = 1 - m;
  if(size > max_size_for_print) return;
  printf("------\n");
  for(x=n;x<size-n;x++)
    {
    for(y=n;y<size-n;y++)
      {
	for(z=n;z<size-n;z++)
	  {
	      printf("%d ",*makept(&currentp,x,y,z));
	  }
	printf("\n");
      }
    if(x<size-n-1) printf("\n");
    else printf("------\n");
    }
}


void   error_exit(char * message)
{
  printf("error exit:");
  printf("%s\n",message);
  printf("\n");
  exit(1);
}

double read_data()
{
  double totalporosity=0.;
  int x,y,z,readitems;
  int *p;
  point_t currentp;
  long int data_size=0;
  long int n_voids=0;
  for(x=1;x<max;x++)
    {
      int porosity =0;
      for(y=1;y<max;y++)
	for(z=1;z<max;z++)
	  {
	    int * p = makept(&currentp,x,y,z);
	    readitems = fscanf(fd, "%d", p );
	    data_size += readitems;
	    if(readitems != 1)
	      {
		printf("FAILED TO READ AT x y z %d %d %d\n",x,y,z);
		exit(1);
	      }
	    // if(*p == 2) *p=1;// asume label 2 is a rock also (perhaps it is water on grains ?)
	    if(*p == 2) 
	      {
		fprintf(stderr,"\n wrong value: *p == 2\n");
		exit(0);
	      }
	    if(*p>1 || *p<0) 
	      {
		printf("wrong value x y z %d %d %d *p=%d",x,y,z,*p);
		exit(1);
	      }
	    // comment out the following because of reversal made by Bertrand
	    //	    *p = 1 - *p;  
	    porosity += *p;
	    n_voids += *p;
	  }
	    double layerporosity = porosity/((double) cubesize*cubesize);
	    if(layerporosity > 1e-4) printf("."); else printf("@");
	    totalporosity += layerporosity;
	    //	    printf("\n number of void pixels :%ld \n",n_voids);
    }
  if(data_size != cubesize*cubesize*cubesize) 
    {
      printf("\nread items = %ld.\n",data_size);
      error_exit("read items differs from given size value.");
    }
  if(stdin != fd) fclose(fd);
  printf("\n number of void pixels :%ld \n",n_voids);
  return totalporosity /= cubesize;
}

int max_X(double *fraction)
{
  double totalporosity=0.;
  int x,y,z;
  int *p;
  point_t currentp;
  int x_max=0;
  for(x=1;x<max;x++)
    {
      int porosity =0;
      for(y=1;y<max;y++)
	for(z=1;z<max;z++)
	  {
	    int * p = makept(&currentp,x,y,z);
	    porosity += (*p > 1 ? 1 : 0);
	    if( *p == COLORED) x_max = x;
	  }
	    double layerporosity = porosity/((double) cubesize*cubesize);
	    totalporosity += layerporosity;
    }
  *fraction = totalporosity / cubesize;
  return x_max;
}



void ghost_face(Directions dir)
{
       int component = dir/2;
       int orientation = dir - component*2;
       int c,i,x[3],xb[3],xg[3];
       point_t bulkpt,ghostpt;
       
       if(!orientation)
	 {
	   xg[component]=0;
	   xb[component]=cubesize;
	 }
       else
	 {
	   xg[component]=max;
	   xb[component]=1;
	 }

       for(x[1]=0;x[1]<size;x[1]++)
	 for(x[2]=0;x[2]<size;x[2]++)
	   {
	     c=0;
	     i=1;
	     while(c<3)
	       {
		 if(c!=component)
		   {
		     xg[c]=xb[c]=x[i];
		     i++;
		   }
		 c++;
	       }
	     *makept(&ghostpt,xg[0],xg[1],xg[2]) = *makept(&bulkpt,xb[0],xb[1],xb[2]);
	   }
}       
void split_print(int nx, int npx)
{
  int i,j,k,rank;
  int is,js,ks,ie,je,ke;
  int coords1,coords2,coords3;
  int ny,nz;
  int npy,npz;
  int mx,my,mz;
  int Ng=2;
  ny=nz=nx;
  npy=npz=npx;
  rank=0;
  char filename[17];
  // xyz style of processor ordering
  mx=nx/npx;my=ny/npy;mz=nz/npz;
  for(coords1=0;coords1<npx;coords1++)
    for(coords2=0;coords2<npy;coords2++)
      for(coords3=0;coords3<npz;coords3++)
	{
	  is=coords1*mx+1; ie=coords1*mx+mx;
	  js=coords2*my+1; je=coords2*my+my;
	  ks=coords3*mz+1; ke=coords3*mz+mz;
	  //if(npx*npy*npz>1) 
	  snprintf(filename,17,"bitmap-%05d.txt",rank);
	    //else
	    //strcpy(filename,"bitmap.txt");
	  //printf("Writing to %s\n",filename);
	  FILE * bitmapfd = fopen(filename,"w");
	  int value;
	  for(i=is;i<=ie;i++)
	    for(j=js;j<=je;j++)
	      for(k=ks;k<=ke;k++)
		{
		  value = *(bigrock + i*size*size + j*size + k); 
		  fprintf(bitmapfd,"%d\n",value);
		}
	  rank++;
	  fclose(bitmapfd);
	}
}
	

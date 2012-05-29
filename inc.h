/*
#ifdef HAVE_MPI
int PPRANK;
#endif
*/

#define COORD(i,n,dx)  (((double) ( (i) - (n)/2) - 0.5)*(dx))
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))


// 1-D heat transfer
// Each site take the value of the mean of left and right neighbors.
// First iteration, mpi send and receive for two processors
/******************************************************************************/
// Usage: convection [options]
//
// Options:
//    -o <file>         Sends output to a file             (default = out.txt)
//    -n <grid_pts>     Sets the number of gris_pts to use (default = 10)
//    -t <time_step>    Sets number of time steps to use   (default = 10)
/******************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
char* program_name;

void arrayPrint(double ** array_ptr, int len, int height, FILE* f){
      //double* next = (U_t+(t_index)*N);
      // Each worker will create output for their own array
   int row;
   int col;
   double * array = *array_ptr;
   double value;
   for(row=0; row < height; row++){
     for(col = 0; col < len; col ++){
       value = *(array+(row*len+col));
       fprintf(f,"%4f ", value);
     }
   fprintf(f,"\n");
   }
}

void timeStep(double** previous_ptr, double** next_ptr, int len, int world_rank,int world_size){
  /*
  1. all workers calculate on inner elements (ignore edge effects)
  2. Calculate the worldLeft
    - (left most on worker 0) // periodic boundary conditions
  3. rank 1 calculates its localLeft (edge of local array)
  4. rank 0 calculates its localRight (edge of local array)
  5 calculate the worldRight
    - (right most on last worker)  periodic boundary conditions
  */
  // Step 1
  // Calculate next step on inner elements (ignore edge effects)
  //printf("Calculating inner elements\n");

  double* previous = *previous_ptr;
  double* next = *next_ptr;

  int i;
  for(i=1; i < len-1; i++){
     // Parallelizable portion
     // but step depend on previous
     next[i] = (previous[i-1] + previous[i+1])/2;
     //printf("Rank %d: next[%d] = %f\n", world_rank, i,next[i]);
     //printf("Rank %d, next[%d] = %f\n",world_rank,i,next[i]);
  }
  //double worldLeft,    worldRight;
  double leftNeighbor, rightNeighbor;
  int leftRank,     rightRank;

  double localLeft = previous[0];
  double localRight = previous[len-1];
  //printf("Setting left and right ranks\n");
// Set left Rank
  if (world_rank == 0){
    // Set periodic boundary conditions
    leftRank = world_size-1;
  }else{
		leftRank = world_rank-1;
	}

// Set right rank
  if (world_rank == world_size-1){
    // Set periodic boundary conditions
    rightRank = 0;
  }else{
		rightRank = world_rank+1;
	}

  //printf("Rank %d: Left rank is %d Right rank is %d\n",world_rank,leftRank,rightRank);
  //printf("Rank %d computing local left\n",world_size);
  // Step 2 Compute local left
// The edge ranks
  if( world_rank == 0 || world_rank == (world_size-1) ){
    // Send localRight to rightRank if you're at the edge of world
		MPI_Send(&localRight,1,MPI_DOUBLE,rightRank,0,MPI_COMM_WORLD);
  }
  // Everyone should receive leftNeighbor from leftRank
	MPI_Recv(&leftNeighbor,1, MPI_DOUBLE,leftRank,0,
						 MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  // The middle ranks
	if ( world_rank != 0 && world_rank!=(world_size-1) ){
		MPI_Send(&localRight,1,MPI_DOUBLE,rightRank,0,MPI_COMM_WORLD);
	}

// Everyone do this calculation
	next[0] = (leftNeighbor + previous[1])/2;
  //printf("Rank %d: next[0] = %f\n", world_rank,next[0]);

  // Step 3 Calculate local right
  //printf("Rank %d computing local righi\n",world_size);
  if( world_rank == 0 || world_rank == (world_size-1) ){
  // Send localLeft to leftRank if you're at the edge of the world
    MPI_Send(&localLeft,1,MPI_DOUBLE,leftRank,0,MPI_COMM_WORLD);
  }
  // receive rightNeighbor from rightRank
  MPI_Recv(&rightNeighbor,1, MPI_DOUBLE,rightRank,0,
             MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  // The middle ranks
  if ( world_rank != 0 && world_rank!=(world_size-1) ){
    MPI_Send(&localLeft,1,MPI_DOUBLE,leftRank,0,MPI_COMM_WORLD);
  }
  next[len-1] = (rightNeighbor + previous[len-1-1])/2;
  //printf("Rank %d: next[%d] = %f\n\n", world_rank,len-1-1,next[len-1-1]);

}
void usage(){
   printf("Usage is %s [options]\n",program_name);
   printf("Options:\n");
   printf("-o <file>       Sends output to a file            (default = out.txt)\n");
   printf("-n <grid_pts>   Sets number of grid points to use (default=10)\n");
   printf("-t <time_steps> Sets number of grid points to use (default=10)\n");
   exit(8);
}

void set_args(int argc, char* argv[],int* N, int* T, char* fileName){
  // Parse input arguments and flags
   int  tIsSet=0;
   int  nIsSet=0;
   int fnIsSet=0;

   while ((argc >1) && (argv[1][0] == '-')){
   // argv[1][1] is the option character
      switch(argv[1][1]){
         case 'o':
            //set output file
            strcpy(fileName,argv[2]);
            break;
         case 'n':
            //set grid points
            *N=atoi(argv[2]);
				nIsSet=1;
            break;
         case 't':
            //set time steps
            *T=atoi(argv[2]);
				tIsSet=1;
            break;
         default:
            printf("Bad option %s\n",argv[1]);
            usage();
      }
      argv+=2;
      argc-=2;
   }
	if(!tIsSet){
		*T=10;
	}
	if(!nIsSet){
		*N=10;
	}
  if (!fnIsSet){
    // The 'x' is for debugging purposes for now
    // each process will create it's out outfile
    char out_file[] = "outfile.txt_x";
    //fileName=&out_file;
    strcpy(fileName,out_file);
  }
}

void initialize(int argc, char* argv[], int* N, int* T,int* n_part,
                int world_rank, int world_size,
                FILE** f,char* fileName, double** U_t_Ptr) {
                  //TODO update for >2 workers
// Initialize constants and starting heat data

  //
  set_args(argc, argv, N, T, fileName);
  // Replace the _x in the outfile name with the world_rank
  // so each process creates a unique output file
  int max_n = *N;
  int max_t = *T;
  fileName[12]='0'+world_rank;

  double dx = 1/(double)max_n;

// Set output file
  *f = fopen(fileName,"w");
  if (f == NULL){
     printf("Error opening file!\n");
     exit(1);
  }

  // Num of grid points per worker

  *n_part=max_n/world_size;
  // if there are an odd number of points, & we are the the last portion,
  // we'll take the extra grid point
  // This only works for two workers!
  if (max_n%world_size!=0 & world_rank==(world_size-1) ){
    *n_part+=(max_n%world_size);
  }

  // Set initial conditions
  // Pre allocate memory for U_t[t,x] data for heat as function of
  // position and time.

  double* U_t = calloc(max_n*max_t,sizeof(double));
  *U_t_Ptr = U_t;
  //printf("Setting initial conditions, U_t size = %dX%d = %lu\n", max_n,max_t,sizeof(*U_t));
  int i;
  int buffer=(max_n/world_size)*world_rank;
  for(i=0; i<*n_part;i++){
    // Initial heat is some crazy function, don't worry too much about it

     *(U_t+i) = 10*exp( - (
               pow((double)(i+buffer)*(dx)-0.5,2)
               /0.05 ) );


     //*(U_t+i) = (double)(i+buffer);
     //printf("In rank %d, U_t_Ptr[%d] = %f\n",world_rank, i,U_t[i]);
  }


}

int main(int argc, char * argv[]) {
   program_name = argv[0];

  //Initialize MPI environment
   MPI_Init(&argc,&argv);

   // Find out rank, size
   int world_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);

   printf("Starting rank %d/%d\n", world_rank,world_size);

   // We are assuming at least 2 processes for this task
   if (world_size < 2) {
     fprintf(stderr, "World size must at least 2 for %s\n", argv[0]);
     MPI_Abort(MPI_COMM_WORLD, 1);
   }

   // Set constants
   int N; // Number of total grid points
   int T; // Number of time steps
   int n_part; // The number of grid points for this rank
   int partner_rank;
   FILE* f;
   double* U_t;
   char out_file[] = "outfile.txt_x";

   initialize(argc, argv, &N, &T,&n_part,world_rank,world_size, &f,out_file, &U_t);

     int t_index;
// For each time step
   for(t_index = 1; t_index < T; t_index++){
     //  Do a forward step
      // Pointer magic, *U_t+X means
      // go to where U_t begins and advance +X elements.
      // A 2D array of dimension nxm starts a new row after m
      // elements.

      double* previous = ( U_t+( (t_index-1)*n_part) );
      double* next = ( U_t+(t_index*n_part) );
    //  timeStep(double* previous, double* next, int len, int world_rank,int partner_rank,int world_size)
      timeStep(&previous, &next,n_part,world_rank,world_size);
      //arrayPrint(next,n_part,f);
   }
  arrayPrint(&U_t,n_part,T,f);
  printf("Ending rank %d/%d\n", world_rank,world_size);
  fclose(f);
  MPI_Finalize();
  return 0;
}

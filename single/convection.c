// 1-D heat transfer
// Each site take the value of the mean of left and right neighbors.
/******************************************************************************/
// Usage: convection [options]
//
// Options:
//    -o <file>         Sends output to a file             (default = out.txt)
//    -n <grid_pts>     Sets the number of gris_pts to use (default = 5)
//    -t <time_step>    Sets number of time steps to use   (default = 10)
/******************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
char* program_name;

void arrayPrint(double * array, int len, FILE* f){
   int index = 0;
   for(index = 0; index < len; index ++){
     double value = array[index];
     fprintf(f,"%4f ", array[index]);
   }
   fprintf(f,"\n");
 }

void timeStep(double* previous, double* next, int len){
   
   double leftEdge  = previous[0];
   double rightEdge = previous[len-1];

   int i;
   for(i=1; i < len-1; i++){
      // Parallelizable portion
      // Next step depend on previous
      next[i] = (previous[i-1] + previous[i+1])/2;
   }
   // Enforce periodic boundary conditions

   next[0] = (previous[1]+rightEdge)/2;
   // len-1 is end of array, len-1-1 is left neighbor of end
   next[len-1] = (leftEdge+previous[len-1-1])/2;
}
void usage(){
   printf("Usage is %s [options]\n",program_name);
   printf("Options:\n");
   printf("-o <file>       Sends output to a file            (default = out.txt)\n");
   printf("-n <grid_pts>   Sets number of grid points to use (default=5)\n");
   printf("-t <time_steps> Sets number of grid points to use (default=10)\n");
   exit(8);
}

void set_args(int argc, char* argv[],int* N, int* T, char* fileName){
   int newT;
   int newN;
   if (argc ==1){
   // There are no arguments
   newT = 10;
   *T = newT;
   newN = 5;
   *N = newN;
   

   } 
   while ((argc >1) && (argv[1][0] == '-')){
   // argv[1][1] is the option character
      switch(argv[1][1]){
         case 'o':
            //set output file
            fileName=&argv[2];
            break;
         case 'n':
            //set grid points
            newN=atoi(argv[2]);
            *N=newN;
            break;
         case 't':
            //set time steps
            newT=atoi(argv[2]);
            *T=newT;
            break;
         default:
            printf("Bad option %s\n",argv[1]);
            usage();
      }
      argv+=2;
      argc-=2;
   }
}

int main(int argc, char *argv[]) {
   program_name=argv[0];
   // Set constants
   int N; // Number of grid points
   int T; // Number of time steps
   char out_file[] = "outfile.txt";
   // Override defaults for options
   set_args(argc, argv, &N, &T, out_file);

   double dx = 1/(double)N;
 
 
   // Pre allocate memory for U_t[t,x]
   double  U_t[T][N];
   //memset(U_t,0,sizeof(double)*N*T);

   // Set output file  
   FILE* f = fopen(out_file,"w");
   if (f == NULL){
      printf("Error opening file!\n");
      exit(1);
   }
   
   // Set initial conditions
   int i;
   for(i=0; i<N;i++){
      U_t[0][i] = 10*exp( - ( pow((double)i*dx-0.5,2)/0.05 ) );
   }
   arrayPrint(*U_t,N,f);

   int t_index;
// For each time step
   for(t_index = 1; t_index < T; t_index++){
//  Do a forward step 
      
      // Pointer magic, *U_t+X means
      // go to where U_t begins and advance +X elements.
      // A 2D array of dimension nxm starts a new row after m
      // elements.
      
      double* previous = (*U_t+(t_index-1)*N); 
      double* next = (*U_t+(t_index)*N);
      timeStep(previous, next,N);
      
      arrayPrint(next,N,f);
   }

  return 0;
}

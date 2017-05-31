// 1-D heat transfer
// Each site take the value of the mean of left and right neighbors.
//
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

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

int main(int argc, char const *argv[]) {
   // Set constants
   int N = 5e0; // Number of grid points
   int T = 1e1; // Number of time steps
 
   double dx = 1/(double)N;
 
 
   // Pre allocate memory for U_t[t,x]
   double  U_t[T][N];
   //memset(U_t,0,sizeof(double)*N*T);

   // Set output file  
   FILE* f = fopen("out.txt","w");
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

// 1-D convection equation
// du/dt + v du/dx = 0
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

int main(int argc, char const *argv[]) {
   // Set constants
   int N = 5; // Number of grid points
   int T = 10;
 
   double dx = 1/(double)N;
   double dt = 1/(double)T;
   double tfinal = 1.0;
   double t = 0;
 
   // heat transfer
   double v = 1.0;
 
   double b = dt*dx/v;
   double a = (1+b);
   
   // Pre allocate memory for U_t[t,x]
   double  U_t[T][N];
   memset(U_t,0,sizeof(double)*N*T);
   // Set output file  
   FILE* f = fopen("out.txt","w");
   if (f == NULL){
      printf("Error opening file!\n");
      exit(1);
   }
   
   // Set initial conditions
   int i;
   for(i=0; i<N;i++){
      U_t[0][i] = 10*exp( - ( pow(i-(double)N/3,2) ) );
   }
   arrayPrint(*U_t,N,f);

   int t_index;
// For each time step
   for(t_index = 1; t_index < T; t_index++){
//  Do a forward step (This part parallizable)
      for(i=1; i < N-1; i++){
         U_t[t_index][i] = (U_t[t_index-1][i-1]+U_t[t_index-1][i+1])/2;
      }
    //Enforce Periodic boundary conditions
   U_t[t_index][N-1] = (U_t[t_index-1][N-2]+U_t[t_index-1][0])/2;
   U_t[t_index][0] = (U_t[t_index-1][1]+U_t[t_index-1][N-1])/2;

   double * array_to_print = (*U_t+t_index*N);
   arrayPrint(array_to_print,N,f);
   }

  return 0;
}

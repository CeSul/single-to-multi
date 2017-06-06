#include <stdio.h>
#include <stdlib.h>

void initialize(int** A_Ptr, int* N){
  N[0]=10;
  int len=N[0];
  int* A = calloc(len,sizeof(int));
  *A_Ptr = A;

  printf("Size of A is %lu\n", sizeof(A));
  int i;
  for(i=0; i<(len); i++){
    A[i]=2*i+1;
  }
}


int main(int argc, char const *argv[]) {
  int* A;
  int N;
  initialize(&A,&N);
  int i;
  for(i=0; i<N; i++){
    printf("A[%d] = %d\n",i,A[i] );
  }
  return 0;
}

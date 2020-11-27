/*
* file:   sequential.c
* Iplemantation of Reverse Cuthill Mckee algorithm (sequential version)
*
* author: Charalampos Papadakis (9128)
* email: papadakic@ece.auth.gr
*/


#include <stdio.h>
#include <stdlib.h>
#include "queue.h"
#include <time.h>
#include <sys/time.h>


#define N 35000 // adjacency matrix's size NxN
#define SPARSITY 0.7 //the percentage of zeros
#define QUEUE_SIZE 10000

//Use it only if I want to see the final reordered table after finding R
int* reorder(int* matrix, int* R, int n) {
  int* newMatrix =(int*) calloc(n*n, sizeof(int));
  if(newMatrix == NULL){
    printf("Error at memory Initialization \n");
    exit(0);
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if(*(matrix+n*R[i]+j) == 1) {
        for (int k = 0; k < n; k++) {
          if(R[k] == j) {
            *(newMatrix+n*i+k) = 1;
          }
        }
      }
    }
  }

  return newMatrix;
}

// //pass the results in a text to visualize them
// void results(int* matrix, int rows, int col, char* file_path)
// {
//   FILE* file = fopen(file_path, "w");
//   if(file == NULL)
//     exit(0);
//
//   for (size_t i = 0; i < rows; i++) {
//     for (size_t j = 0; j < col; j++) {
//       fprintf(file, "%d, ", *(matrix+i*col+j));
//     }
//     fprintf(file, "\n");
//   }
//   fclose(file);
// }


void swapElement(int* a, int* b) {
       int temp = *a;
       *a = *b;
       *b = temp;
     }


//Matrix Initialization function
void sparseSymmetric (int * matrix , int n , double sparsity){
//sparsityLimit the percentage of the zeros
  int notZero = (n*n) - (n*n*sparsity);
  int  x=0, y=0;
  srand(time(NULL));
  for (int i = 0; i < notZero; i+=2) {
    do {
      x = rand() % n;
      y = rand() % n;
    } while(x == y);

    *(matrix+x*n+y) = 1;
    *(matrix+y*n+x) = 1;
  }

  // //passing results to a text file to use it for visualisation of the problem
  // FILE* file = fopen("generatedMatrix.txt", "w");
  //     if(file == NULL)
  //       exit(0);
  //
  //     for (int i = 0; i < n; i++) {
  //       for (int j = 0; j < n; j++) {
  //         fprintf(file, "%d, ", *(matrix+i*n+j));
  //       }
  //       fprintf(file, "\n");
  //     }
  //     fclose(file);
  }

  //sorts the connected nodes with the one we are adding to R, in asceding order depending on their degrees
  void sortByDegree(int *connections ,int  *degrees , int connectedNo){

  int *tempArray = (int *)malloc(connectedNo*sizeof(int));
  if(tempArray == NULL){
    printf("Error at memory Initialization \n");
    exit(0);
  }

  for(int i = 0; i < connectedNo ; i++){
    tempArray[i] = degrees[connections[i]];
  }
  //Sort the array of the connections in increasing order by their degree
  for (int i = 0; i < connectedNo-1; i++){
      for (int j = 0; j < connectedNo-i-1; j++){
        if (tempArray[j] > tempArray[j+1]){
            swapElement(&tempArray[j], &tempArray[j+1]);
            swapElement(&connections[j], &connections[j+1]);
        }
      }
    }
  }


void ReverseCuthilMckee(int *matrix , int *degrees ,int *  R){

Queue *Q  = createQueue(QUEUE_SIZE);
int Rcounter = 0; //counts how many nodes are inside R , we need to reach N nodes
int *notAdded = (int*)malloc(N*sizeof(int));

int currentIndex ;

  for (int  i = 0; i < N; i++) {
    notAdded[i] = 1;
  }

//Iterate for all the rows
while(Rcounter != N){

  int minDegIndex = 0 ;
  for (int i = 0; i< N ; i++){
    if(degrees[i] < degrees[minDegIndex] && notAdded[i]==1){ //if the node is not already visited and has smaller degree we make it minimum
      minDegIndex = i;
    }
  }

  enqueue( Q, minDegIndex); //we add to the queue the index in which is the minimum degree
  notAdded[minDegIndex] = 0 ; //we define that it is now visited

// until the queue empties we dequeue nodes and do the operation with the neighbors again
  while(!(isEmpty(Q))){

   currentIndex = dequeue(Q);
   //Find the connections of the current Index
   int *connections = (int *)malloc(degrees[currentIndex]*sizeof(int));
   if(connections == NULL){
     printf("Error at memory Initialization \n");
     exit(0);
   }
   int connectedNo = 0;


   for (int i = 0; i < N ; i++){
     if( i!= currentIndex && *(matrix+(currentIndex)*N + i) == 1 && notAdded[i] == 1 ){
       connections[connectedNo++] = i;
       notAdded[i] = 0;
     }
   }

   //sort the neighbors found by their degree value in order to get the increasing order
   sortByDegree(connections , degrees , connectedNo);

//we add to the queue the connections in ascending order
   for(int i = 0; i < connectedNo ; i++){
       enqueue( Q, connections[i]);
   }
   R[Rcounter++] = currentIndex;

   free(connections);
   }
  }
   free(Q);

   //Reversing R Matrix
   int nSize  = N;

   if (nSize % 2 == 0)
       nSize -= 1;

      nSize= nSize / 2;

       for (int i = 0; i <= nSize; i++) {
           int j = R[N - 1 - i];
           swapElement(&R[N - 1 - i] , &R[i]);
           R[i] = j;
         }

   }


//we find the degree of all the nodes (how many connections each one has)
   void degFinder(int *matrix,int *  degrees){
     int sum = 0 ;
     for(int i = 0; i< N; i++){
       for(int j = 0; j < N ; j++){
         sum +=*(matrix + i*N + j);
       }
       degrees[i]=sum;
       sum = 0;
     }
   }

int main(int  argc, const char* argv[]){

int *  matrix = (int * )calloc(N*N,sizeof(int)); //the sparse symmetric matrix
int * degrees = (int *)calloc(N,sizeof(int)); //the matrix with all the degrees
int  * R = (int *)calloc(N,sizeof(int)); //the matrix that will show how we should reorder

if(matrix == NULL || degrees == NULL || R == NULL){
  printf("Error at memory Initialization \n");
  exit(0);
}

struct timeval start, end;

sparseSymmetric(matrix, N , SPARSITY);


gettimeofday(&start, NULL);
degFinder(matrix , degrees) ;
ReverseCuthilMckee(matrix , degrees , R);
gettimeofday(&end, NULL);

double time = ((double)((end.tv_sec*1e6 + end.tv_usec) - (start.tv_sec*1e6 + start.tv_usec)))*1e-6;
printf(" >>> Time it took to execute the Reverse Cuthill Mckee algorithm to the given matrix was: %lf sec\n", time);

// for(int i=0; i<N; i++){
//   printf("%d\n",R[i]);
// }
// int* reorderedMatrix = reorder( matrix, R, N) ;
// results(reorderedMatrix, N, N, "output_matrix.txt");
}

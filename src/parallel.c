/*
* file:   parallel.c
* Iplemantation of Reverse Cuthill Mckee algorithm (parallel-openMP version)
*
* author: Charalampos Papadakis (9128)
* email: papadakic@ece.auth.gr
*/


#include <stdio.h>
#include <stdlib.h>
#include "queue.h"
#include <time.h>
#include <sys/time.h>
#include <omp.h>



#define N 35000 // adjacency matrix's size NxN
#define SPARSITY 0.7 //the percentage of zeros
#define QUEUE_SIZE 10000
#define MAX_THREADS 2 //the threads we will use for parallel tasks

int CHUNKSIZE = N / MAX_THREADS;

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


void swapElement(int* a, int* b) {
       int temp = *a;
       *a = *b;
       *b = temp;
     }


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

//Matrix Initialization function
void sparseSymmetric (int * matrix , int n , double sparsity){
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

  int minIndex(int* degrees, int* notVisited, int n, int t) {
    int deg[t];
    int idx[t];
    int tid;

    // find the min degree
    int i;
    int CHUNKSIZE = n / t;
    #pragma omp parallel num_threads(MAX_THREADS) private(i, tid)
    {
      tid = omp_get_thread_num();
      int minIndex = -1;  // The pos of min degree node in matrix
      int min = n+5; // A node can not have degree > SIZE
      #pragma omp for schedule(static, CHUNKSIZE)
        for (i = 0; i < n; i++) {
          if(degrees[i] < min && notVisited[i] == 1) {
            minIndex = i;
            min = degrees[i];
          }
        }
        deg[tid] = min;
        idx[tid] = minIndex;
    }

    int totalIdx = idx[0];
    int totalMin = deg[0];

    for (size_t i = 1; i < t; i++) {
      if(deg[i] < totalMin) {
        totalMin = deg[i];
        totalIdx = idx[i];
      }
    }

    return totalIdx;
  }

void ReverseCuthilMckee(int *matrix , int *degrees ,int *  R){

Queue *Q  = createQueue(QUEUE_SIZE);
int Rcounter = 0;
int *notAdded = (int*)malloc(N*sizeof(int));
if( notAdded == NULL){
  printf("Error at memory Initialization \n");
  exit(0);
}
int currentIndex ;
int i;

omp_lock_t writelock;
omp_init_lock(&writelock);

#pragma omp parallel for  schedule(static,CHUNKSIZE) num_threads (MAX_THREADS)  private(i) shared(notAdded)
  for (i = 0; i < N; i++) {
    notAdded[i] = 1;
  }

//Iterate for all the rows
while(Rcounter != N){


  int minDegreeIndex  = minIndex(degrees, notAdded , N, MAX_THREADS);

  enqueue( Q, minDegreeIndex);
  notAdded[minDegreeIndex] = 0 ;

  while(!(isEmpty(Q))){

   currentIndex = dequeue(Q);
   //Find the connections of the current Index
   int *connections = (int *)malloc(degrees[currentIndex]*sizeof(int));
   if(connections == NULL){
     printf("Error at memory Initialization \n");
     exit(0);
   }
   int connectedNo = 0;

   #pragma omp parallel for  schedule(static,CHUNKSIZE) num_threads (MAX_THREADS)  private(i) shared(notAdded , matrix)
   for (int i = 0; i < N ; i++){
     if( i!= currentIndex && *(matrix+(currentIndex)*N + i) == 1 && notAdded[i] == 1 ){
       omp_set_lock(&writelock);
       connections[connectedNo++] = i;
       notAdded[i] = 0;
       omp_unset_lock(&writelock);
     }
   }

   //sort the neighbors found by their degree value in order to get the increasing order
   sortByDegree(connections , degrees , connectedNo);

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

int *  matrix = (int * )calloc(N*N,sizeof(int));
int * degrees = (int *)calloc(N,sizeof(int));
int  * R = (int *)calloc(N,sizeof(int));

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
//

// int* reorderedMatrix = reorder( matrix, R, N) ;
// results(reorderedMatrix, N, N, "output_matrix.txt");

}

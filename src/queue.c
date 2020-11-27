#include <stdio.h>
#include <stdlib.h>
#include "queue.h"
// function to create a queue
// of given capacity.
// It initializes size of queue as 0
 Queue * createQueue(unsigned capacity)
{
   Queue * queue = ( Queue*)malloc(
        sizeof(Queue));
    queue->capacity = capacity;
    queue->first = queue->size = 0;

    // This is important, see the enqueue
    queue->last = capacity - 1;
    queue->array = (int*)malloc(
        queue->capacity * sizeof(int));
    return queue;
}

// Queue is full when size becomes
// equal to the capacity
int isFull( Queue * queue)
{
    return (queue->size == queue->capacity);
}

// Queue is empty when size is 0
int isEmpty( Queue * queue)
{
    return (queue->size == 0);
}

// Function to add an item to the queue.
// It changes last and size
void enqueue( Queue * queue, int item)
{
    if (isFull(queue))
        return;
    queue->last = (queue->last + 1)
                  % queue->capacity;
    queue->array[queue->last] = item;
    queue->size = queue->size + 1;
  //  printf("%d enqueued to queue\n", item);
}

// Function to remove an item from queue.
// It changes first and size
int dequeue(Queue * queue)
{
    if (isEmpty(queue))
        return 0;
    int item = queue->array[queue->first];
    queue->first = (queue->first + 1)
                   % queue->capacity;
    queue->size = queue->size - 1;
    return item;
}

// Function to get first element of queue
int first(Queue * queue)
{
    if (isEmpty(queue))
        return 0;
    return queue->array[queue->first];
}

// Function to get last element of queue
int last( Queue * queue)
{
    if (isEmpty(queue))
        return 0;
    return queue->array[queue->last];
}

#ifndef PRIORITYQ_H
#define PRIORITYQ_H
// this is a specialized heap based priority Q designed for use in dijkstras
// you must know in advance how big the capacity will need to be
// moroever, the 'key's (vertices) must go from 0 -> capicity-1
// (perhaps 'key' & value are swapped

// data structures
typedef struct
{
    double value;
    int key;      // not really necessary since can figure it out from position in array
    int heapIdx;
} entry;


typedef struct
{
    int capacity; // how many entries there (not mutable!)
    int nextOpen; // the lowest idx free slot
    entry** heap; // actual 'q'
    entry* values;// raw data that gets orginized in the heap
} PriorityQ;

// interface
PriorityQ* newPriorityQ(int cap);
void freePriorityQ(PriorityQ* pq);
int pqIsEmpty(PriorityQ* pq);
void pqInsert(PriorityQ* pq, int key, double value);
int pqTop(PriorityQ* pq);
void pqPop(PriorityQ* pq);
void pqDecreaseKey(PriorityQ* pq, int key, double newValue);

#endif // PRIORITYQ_H

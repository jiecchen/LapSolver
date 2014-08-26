#pragma once
#include <utility>
#include <algorithm>
#include "Heap.h"

template <typename Key, typename Value>
class BinaryHeap : public Heap<Key, Value>
{
public:
    BinaryHeap(int capacity)
        : capacity(capacity),
          size(0),
          heap(new Value[capacity+1]),
          heap_index(new int[capacity+1]),
          priority(new Key[capacity+1])
    {
        heap_index[0:capacity] = -1;
        priority[0:capacity] = std::numeric_limits<Key>::max();
    }

    virtual ~BinaryHeap()
    {
        delete[] heap;
        delete[] heap_index;
        delete[] priority;
    }

    inline void pop(Key *oKey, Value *oValue);
    inline void push(Key key, Value pri);
    inline void decreaseKey(Key key, Value pri);
    inline bool isEmpty()
    {
        return size == 0;
    }

    void print()
    {
        for (unsigned i = 0; i < size; ++i)
        {
            printf("%d:%lf ", heap[i], priority[heap[i]]);
        }
        printf("\n");
    }

private:
    const int capacity;
    int size;

    Value *heap;
    Key *priority;
    int *heap_index;

    inline void upHeap(Value value);
    inline void downHeap(Value value);
    inline void swapKeys(Value value1, Value key2);
};


template <typename Key, typename Value>
void BinaryHeap<Key, Value>::pop(Key *oKey, Value *oValue)
{
    *oKey = priority[heap[0]];
    *oValue = heap[0];

    size--;
    swapKeys(heap[0], heap[size]);
    heap_index[heap[size]] = -1;

    downHeap(heap[0]);
}


template <typename Key, typename Value>
void BinaryHeap<Key, Value>::push(Key pri, Value key)
{
    heap[size] = key;
    heap_index[key] = size++;
    priority[key] = pri;

    upHeap(key);
}

template <typename Key, typename Value>
void BinaryHeap<Key, Value>::decreaseKey(Key pri, Value key)
{
    priority[key] = pri;
    upHeap(key);
    downHeap(key);
}

template <typename Key, typename Value>
void BinaryHeap<Key, Value>::upHeap(Value value)
{
    if (heap_index[value] == 0) return;
    int parentIdx = (heap_index[value] - 1) / 2;
    Value parent = heap[parentIdx];

    if (priority[value] < priority[parent])
    {
        swapKeys(value, parent);
        upHeap(value);
    }
}

template <typename Key, typename Value>
void BinaryHeap<Key, Value>::downHeap(Value value)
{
    if(value == -1) return;
    int leftIdx = (heap_index[value] * 2) + 1;
    int rightIdx = (heap_index[value] * 2) + 2;
    int to_swap = value;

    if (leftIdx < size && priority[heap[leftIdx]] < priority[to_swap])
        to_swap = heap[leftIdx];
    else if (rightIdx < size && priority[heap[rightIdx]] < priority[to_swap])
        to_swap = heap[rightIdx];

    if(to_swap != value) {
        swapKeys(value, to_swap);
        downHeap(value);
    }
}

template <typename Key, typename Value>
void BinaryHeap<Key, Value>::swapKeys(Value value1, Value value2)
{
    std::swap(heap[heap_index[value1]], heap[heap_index[value2]]);
    std::swap(heap_index[value1], heap_index[value2]);
}

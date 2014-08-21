#pragma once

#include <utility>

class PriorityQueue {
  public:
    PriorityQueue (int capacity);
    ~PriorityQueue ();

    void Push (int key, double pri);
    void DecreaseKey (int key, double pri);
    std::pair<int, double> Pop ();

  private:
    void UpHeap (int key);
    void DownHeap (int key);
    void SwapKeys (int key1, int key2);

    int* heap;
    int* heap_index;
    double* priority;
    int size;
};

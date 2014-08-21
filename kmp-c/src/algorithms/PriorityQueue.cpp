#include <utility>
#include <algorithm>
#include <set>

#include "PriorityQueue.h"

using namespace std;

PriorityQueue::PriorityQueue (int capacity) {
    heap = new int[capacity];
    heap_index = new int[capacity];
    priority = new double[capacity];

    heap_index[0:capacity] = -1;

    size = 0;
}

PriorityQueue::~PriorityQueue () {
    delete[] heap;
    delete[] heap_index;
    delete[] priority;
}

void PriorityQueue::Push (int key, double pri) {
    heap[size] = key;
    heap_index[key] = size;
    priority[key] = pri;
    size++;

    UpHeap(key);
}

void PriorityQueue::DecreaseKey (int key, double pri) {
    priority[key] = pri;
    UpHeap(key);
    DownHeap(key);
}

pair<int, double> PriorityQueue::Pop () {
    pair<int, double> ret(heap[0], priority[heap[0]]);

    SwapKeys(0, size-1);
    heap_index[heap[size-1]] = -1;
    size--;

    DownHeap(heap[0]);

    return ret;
}

void PriorityQueue::UpHeap (int key) {
    if (heap_index[key] == 0) return;

    int parent = (heap_index[key]-1) / 2;
    if (priority[heap[parent]] > priority[key]) {
        SwapKeys (key, heap[parent]);
        UpHeap (key);
    }
}

void PriorityQueue::DownHeap (int key) {
    int left_child = 2*(heap_index[key]+1)-1, right_child = 2*(heap_index[key]+1);
    int to_swap = key;

    if (left_child < size && priority[heap[left_child]] < priority[heap[to_swap]]) {
        to_swap = left_child;
    }

    if (right_child < size && priority[heap[right_child]] < priority[heap[to_swap]]) {
        to_swap = right_child;
    }

    if (to_swap != key) {
        SwapKeys (key, heap[to_swap]);
        DownHeap (key);
    }
}

void PriorityQueue::SwapKeys (int key1, int key2) {
    swap(heap[key1], heap[key2]);
    swap(priority[heap[key1]], priority[heap[key2]]);
    swap(heap_index[heap[key1]], heap_index[heap[key2]]);
}

#pragma once
#include <utility>
#include <algorithm>
#include <vector>
#include "Heap.h"

template <typename Key, typename Value>
class BinaryHeap : public Heap<Key, Value>
{
    typedef HeapElement<Key, Value> Element;
    static const Key maxValue = std::numeric_limits<Key>::max();
    static const Key minValue = std::numeric_limits<Key>::min();

public:
    BinaryHeap(int capacity)
        : capacity(capacity),
          size(0),
          data_v(std::vector<Element>(capacity + 2)),
          data(data_v.data())
    {
        Value defVal = Value();
        data[0].key = minValue;
        data[1:capacity + 1] = (Element) { maxValue, defVal };
    }

    virtual ~BinaryHeap() {}

    inline void pop(Key *oKey, Value *oValue);
    inline void push(Key key, Value pri);
    inline void decreaseKey(Key key, Value pri);
    inline bool isEmpty()
    {
        return size == 0;
    }

private:
    const int capacity;
    int size;

    std::vector<Element> data_v;
    Element *data;
};

template <typename Key, typename Value>
void BinaryHeap<Key, Value>::pop(Key *oKey, Value *oValue)
{
    *oKey = data[1].key;
    *oValue = data[1].value;

    data[1] = data[size];
    data[size] = (Element) { maxValue, Value() };
    size--;

    int cur = 1;
    int succ = 2;
    while (succ <= size)
    {
        int left = cur << 1;
        int right = left + 1;
        succ = (data[left].key < data[right].key) ? left : right;

        if (data[cur].key > data[succ].key)
            std::swap(data[cur], data[succ]);
        else break;
        cur = succ;
        succ = cur << 1;
    }
}


template <typename Key, typename Value>
void BinaryHeap<Key, Value>::push(Key pri, Value val)
{
    int cur = ++size;
    int pred = size >> 1;
    data[size] = (Element) { pri, val };

    Key predKey = data[pred].key;

    while (predKey > pri)
    {
        std::swap(data[cur], data[pred]);
        cur = pred;
        predKey = data[pred >>= 1].key;
    }
}

template <typename Key, typename Value>
void BinaryHeap<Key, Value>::decreaseKey(Key pri, Value key)
{
    fprintf(stderr, "BinaryHeap::decreaseKey not implemented!\n");
    exit(0);
}

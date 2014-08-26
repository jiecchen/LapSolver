#include <limits>
#include "Heap.h"
#include "util/aligned.h"

template <typename Key, typename Value>
class Aligned4aryHeap : public Heap<Key, Value>
{
public:
    Aligned4aryHeap(int capacity)
        : capacity(capacity),
          nElements(0),
          key_v(aligned_vector<Key, AlignCache>(capacity + 4)),
          value_v(aligned_vector<Value, AlignCache>(capacity + 4)),
          keys(key_v.data()),
          values(value_v.data())
    {
        keys[0] = minValue;
        keys[1:capacity + 3] = maxValue;
    }

    virtual ~Aligned4aryHeap() {}

    inline void pop(Key *oKey, Value *oValue);
    inline void push(Key key, Value pri);
    inline void decreaseKey(Key key, Value pri);
    inline bool isEmpty()
    {
        return nElements == 0;
    }

private:
    static const Key maxValue = std::numeric_limits<Key>::max();
    static const Key minValue = std::numeric_limits<Key>::min();

    int capacity;
    int nElements;

    aligned_vector<Key, AlignCache> key_v;
    aligned_vector<Value, AlignCache> value_v;
    Key *keys;
    Value *values;
};

template <typename Key, typename Value>
void Aligned4aryHeap<Key, Value>::pop(Key *oKey, Value *oValue)
{
    *oKey = keys[1];
    *oValue = values[1];

    Key minKey, childKey;
    int child;

    int hole = 1;
    int succ = 2;
    int sz = nElements--;

    // Starting at the root
    while (succ < sz)
    {
        // Find the smallest child
        minKey = keys[succ]; child = 0;

        childKey = keys[succ + 1];
        if (childKey < minKey)
        {
            minKey = childKey;
            child = 1;
        }

        childKey = keys[succ + 2];
        if (childKey < minKey)
        {
            minKey = childKey;
            child = 2;
        }

        childKey = keys[succ + 3];
        if (childKey < minKey)
        {
            minKey = childKey;
            child = 3;
        }

        // Move it up
        succ += child;
        keys[hole] = minKey;
        values[hole] = values[succ];

        hole = succ;
        succ = (hole << 2) - 2;
    }

    // Now there's a hole where the end of the min-chain was
    Key bubble = keys[sz];
    int pred = (hole + 2) >> 2;

    // Find the new position for the rightmost element
    while (keys[pred] > bubble)
    {
        keys[hole] = keys[pred];
        values[hole] = values[pred];
        hole = pred;
        pred = (hole + 2) >> 2;
    }

    // Move in the data
    keys[hole] = bubble;
    values[hole] = values[sz];

    // Clear it out
    keys[sz] = maxValue;
}

template <typename Key, typename Value>
void Aligned4aryHeap<Key, Value>::push(Key key, Value pri)
{
    int hole = ++nElements;
    int pred = (hole + 2) >> 2;
    Key predKey = keys[pred];

    while (predKey > key)
    {
        keys[hole] = predKey;
        values[hole] = values[pred];

        hole = pred;
        pred = (hole + 2) >> 2;
        predKey = keys[pred];
    }

    keys[hole] = key;
    values[hole] = pri;
}

template <typename Key, typename Value>
void Aligned4aryHeap<Key, Value>::decreaseKey(Key key, Value pri)
{
    fprintf(stderr, "Aligned4aryHeap::decreaseKey not implemented!\n");
    exit(0);
}

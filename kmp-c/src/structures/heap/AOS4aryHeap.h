#include <limits>
#include "Heap.h"
#include "util/aligned.h"

template <typename Key, typename Value>
class AOS4aryHeap : public Heap<Key, Value>
{
    typedef HeapElement<Key, Value> Element;

public:
    AOS4aryHeap(int capacity)
        : capacity(capacity),
          nElements(0),
          data_v(aligned_vector<Element, AlignCache>(capacity + 4)),
          data(data_v.data())
    {
        Value defVal = Value();
        data[0].key = minValue;
        data[1:capacity+3] = (Element) { maxValue, defVal };
    }

    virtual ~AOS4aryHeap() {}

    inline void pop(Key *oKey, Value *oValue);
    inline void push(Key key, Value pri);
    inline void decreaseKey(Key key, Value pri);
    inline bool isEmpty() { return nElements == 0; }

private:
    static const Key maxValue = std::numeric_limits<Key>::max();
    static const Key minValue = std::numeric_limits<Key>::min();

    int capacity;
    int nElements;

    aligned_vector<Element, AlignCache> data_v;
    Element* data;
};

template <typename Key, typename Value>
void AOS4aryHeap<Key, Value>::pop(Key *oKey, Value *oValue)
{
    *oKey = data[1].key;
    *oValue = data[1].value;

    Key minKey, childKey;
    int child;

    int hole = 1;
    int succ = 2;
    int sz = nElements--;

    // Starting at the root
    while (succ < sz) {
        // Find the smallest child
        minKey = data[succ].key; child = 0;

        childKey = data[succ + 1].key;
        if(childKey < minKey) { minKey = childKey; child = 1; }

        childKey = data[succ + 2].key;
        if(childKey < minKey) { minKey = childKey; child = 2; }

        childKey = data[succ + 3].key;
        if(childKey < minKey) { minKey = childKey; child = 3; }

        // Move it up
        succ += child;
        data[hole].key = minKey;
        data[hole].value = data[succ].value;

        hole = succ;
        succ = (hole << 2) - 2;
    }

    // Now there's a hole where the end of the min-chain was
    Key bubble = data[sz].key;
    int pred = (hole + 2) >> 2;

    // Find the new position for the rightmost element
    while (data[pred].key > bubble) {
        data[hole] = data[pred];
        hole = pred;
        pred = (hole + 2) >> 2;
    }

    // Move in the data
    data[hole].key = bubble;
    data[hole].value = data[sz].value;

    // Clear it out
    data[sz].key = maxValue;
}

template <typename Key, typename Value>
void AOS4aryHeap<Key, Value>::push(Key key, Value pri)
{
    int hole = ++nElements;
    int pred = (hole + 2) >> 2;
    Key predKey = data[pred].key;

    while (predKey > key) {
        data[hole].key = predKey;
        data[hole].value = data[pred].value;

        hole = pred;
        pred = (hole + 2) >> 2;
        predKey = data[pred].key;
    }

    data[hole].key = key;
    data[hole].value = pri;
}

template <typename Key, typename Value>
void AOS4aryHeap<Key, Value>::decreaseKey(Key key, Value pri)
{
    fprintf(stderr, "AOS4aryHeap::decreaseKey not implemented!\n");
    exit(0);
}

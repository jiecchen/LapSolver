#include <limits>
#include "Heap.h"
#include "util/aligned.h"

template <typename Key, typename Value>
class AOS4aryHeap : public Heap<Key, Value>
{
    typedef HeapElement<Key, Value> Element;
    static const int offset = (AlignCache - 1) / sizeof(Element);

public:
    AOS4aryHeap(int capacity)
        : capacity(capacity),
          size(0),
          data_v(aligned_vector<Element, AlignCache>(capacity + 4 + offset)),
          data(data_v.data() + offset - 1)
    {
        Value defVal = Value();
        data[0].key = minValue;
        data[1:capacity + 3] = (Element) { maxValue, defVal };
    }

    virtual ~AOS4aryHeap() {}

    inline void pop(Key *oKey, Value *oValue);
    inline void push(Key key, Value pri);
    inline void decreaseKey(Key key, Value pri);
    inline bool isEmpty()
    {
        return size == 0;
    }

private:
    static const Key maxValue = std::numeric_limits<Key>::max();
    static const Key minValue = std::numeric_limits<Key>::min();

    int capacity;
    int size;

    aligned_vector<Element, AlignCache> data_v;
    Element *data;
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
    int sz = size--;

    Element *_data = data;

    // Starting at the root
    while (succ < sz)
    {
        // Find the smallest child
        minKey = _data[succ].key; child = 0; // all the read misses happen here

        childKey = _data[succ + 1].key;
        if (childKey < minKey)
        {
            minKey = childKey;
            child = 1;
        }

        childKey = _data[succ + 2].key;
        if (childKey < minKey)
        {
            minKey = childKey;
            child = 2;
        }

        childKey = _data[succ + 3].key;
        if (childKey < minKey)
        {
            minKey = childKey;
            child = 3;
        }

        // Move it up
        succ += child;
        _data[hole].key = minKey;
        _data[hole].value = _data[succ].value;

        hole = succ;
        succ = (hole << 2) - 2;
    }

    // Now there's a hole where the end of the min-chain was
    Key bubble = _data[sz].key;
    int pred = (hole + 2) >> 2;

    // Find the new position for the rightmost element
    while (_data[pred].key > bubble)
    {
        _data[hole] = _data[pred];
        hole = pred;
        pred = (hole + 2) >> 2;
    }

    // Move in the data
    _data[hole].key = bubble;
    _data[hole].value = _data[sz].value;

    // Clear it out
    _data[sz].key = maxValue;
}

template <typename Key, typename Value>
void AOS4aryHeap<Key, Value>::push(Key key, Value pri)
{
    int hole = ++size;
    int pred = (hole + 2) >> 2;
    Key predKey = data[pred].key;

    while (predKey > key)
    {
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

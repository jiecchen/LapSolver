#include <limits>
#include <vector>
#include <set>
#include <cstdio>
#include "Heap.h"
using namespace std;

class ArrayCompare
{
public:
    ArrayCompare(double *p): p(p) { };

    bool operator() (const int &x, const int &y) const
    {
        if (p[x] == p[y])
            return x < y;
        return p[x] < p[y];
    }

private:
    double *p;
};

template <typename Key, typename Value>
class TreeHeap : public Heap<Key, Value>
{
    typedef HeapElement<Key, Value> Element;

public:
    TreeHeap(int capacity) :
        keys(std::vector<Key>(capacity)),
        tree(std::set<Value, ArrayCompare>(ArrayCompare(keys.data())))
    {
    }

    virtual ~TreeHeap() {}

    inline void pop(Key *oPri, Value *oVal)
    {
        *oVal = *tree.begin();
        *oPri = keys[*oVal];
        tree.erase(*oVal);
    }

    inline void push(Key pri, Value val)
    {
        keys[val] = pri;
        tree.insert(val);
    }

    inline void decreaseKey(Key pri, Value val)
    {
        tree.erase(val);
        push(pri, val);
    }

    inline bool isEmpty()
    {
        return tree.empty();
    }

private:
    std::vector<Key> keys;
    std::set<Value, ArrayCompare> tree;
};

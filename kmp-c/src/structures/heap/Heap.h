#pragma once

template<typename Key, typename Value>
class Heap {
public:
    virtual ~Heap() {}

    virtual void pop(Key *oKey, Value *oValue) = 0;
    virtual void push(Key key, Value pri) = 0;
    virtual void decreaseKey(Key key, Value pri) = 0;
    virtual bool isEmpty() = 0;
};

template <typename Key, typename Value>
struct HeapElement { Key key; Value value; };

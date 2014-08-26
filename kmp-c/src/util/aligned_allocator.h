#pragma once
#include <memory>
#include "mkl.h"

enum Alignment : size_t
{
    AlignSIMD = 16,
    AlignAVX = 32,
    AlignCache = 64
};

template <typename T, size_t TALIGN>
class aligned_allocator : public std::allocator<T>
{
    typedef typename std::allocator<T>::pointer pointer;
    typedef typename std::allocator<T>::size_type size_type;
public:
    pointer allocate(size_type n, const void *hint = nullptr)
    {
        size_t count = sizeof(T) * n;
        return reinterpret_cast<pointer>
               ((!hint)
                ? mkl_calloc(sizeof(T), n, TALIGN)
                : mkl_realloc((void *) hint, count));
    }

    void deallocate(pointer p, size_type n)
    {
        mkl_free(p);
    }

    template <typename U> struct rebind
    {
        typedef aligned_allocator<U, TALIGN> other;
    };
};

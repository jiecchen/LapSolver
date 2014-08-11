#pragma once
#include <memory>
#include "mkl.h"

enum Alignment : size_t
{
    AlignSIMD = 16,
    AlignAVX = 32
};

template <typename T, size_t TALIGN = AlignSIMD, size_t TBLOCK = 4>
class aligned_allocator : public std::allocator<T>
{
    typedef typename std::allocator<T>::pointer pointer;
    typedef typename std::allocator<T>::size_type size_type;
public:
    pointer allocate(size_type n, const void *hint = nullptr)
    {
        size_t count = sizeof(T) * n;
        size_t count_left = count % TBLOCK;
        if (count_left != 0)
            count += TBLOCK - count_left;
        return reinterpret_cast<pointer>
               ((!hint)
                ? mkl_malloc(count, TALIGN)
                : mkl_realloc((void *) hint, count));
    }

    void deallocate(pointer p, size_type n)
    {
        mkl_free(p);
    }
};

#include <algorithm>
#include <type_traits>
#include <cilk/cilk.h>
#include "edgelist.h"

template <typename T>
static T arrayMax(T *arr, unsigned int n)
{
    return __sec_reduce_max(arr[0:n]);
}

EdgeList::EdgeList(int ne)
    : ne(ne),
      nv(0)
{
    u = aligned_vector<int>(ne, 0);
    v = aligned_vector<int>(ne, 0);
    w = aligned_vector<double>(ne, 0.0);
}

EdgeList::EdgeList(int ne, int *u, int *v)
    : ne(ne),
      nv(1 + std::max(arrayMax(u, ne),
                      arrayMax(v, ne)))
{
    this->u = aligned_vector<int>(u, u + ne);
    this->v = aligned_vector<int>(v, v + ne);
    this->w = aligned_vector<double>(ne, 0.0);
}

EdgeList::EdgeList(int ne, int *u, int *v, double *w)
    : ne(ne),
      nv(1 + std::max(arrayMax(u, ne),
                      arrayMax(v, ne)))
{
    this->u = aligned_vector<int>(u, u + ne);
    this->v = aligned_vector<int>(v, v + ne);
    this->w = aligned_vector<double>(w, w + ne);
}

EdgeList::EdgeList(const aligned_vector<int> &u, const aligned_vector<int> &v)
    : ne(u.size()),
      nv(1 + std::max(arrayMax(u.data(), ne),
                      arrayMax(v.data(), ne)))
{
    this->u = aligned_vector<int>(u);
    this->v = aligned_vector<int>(v);
    this->w = aligned_vector<double>(ne, 0.0);
}

EdgeList::EdgeList(const aligned_vector<int> &u,
                   const aligned_vector<int> &v,
                   const aligned_vector<double> &w)
    : ne(u.size()),
      nv(1 + std::max(arrayMax(u.data(), ne),
                      arrayMax(v.data(), ne)))
{
    this->u = aligned_vector<int>(u);
    this->v = aligned_vector<int>(v);
    this->w = aligned_vector<double>(w);
}

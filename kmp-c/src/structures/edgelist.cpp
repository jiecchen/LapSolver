#include <algorithm>
#include <cilk/cilk.h>
#include <cilk/reducer_opand.h>
#include "edgelist.h"

EdgeList::EdgeList(int ne): ne(ne)
{
    u = aligned_vector<int>(ne, 0);
    v = aligned_vector<int>(ne, 0);
    w = aligned_vector<double>(ne, 0.0);
}

EdgeList::EdgeList(int ne, int *u, int *v): ne(ne)
{
    this->u = aligned_vector<int>(u, u + ne);
    this->v = aligned_vector<int>(v, v + ne);
    this->w = aligned_vector<double>(ne, 0.0);
}

EdgeList::EdgeList(int ne, int *u, int *v, double *w): ne(ne)
{
    this->u = aligned_vector<int>(u, u + ne);
    this->v = aligned_vector<int>(v, v + ne);
    this->w = aligned_vector<double>(w, w + ne);
}

EdgeList::EdgeList(const aligned_vector<int> &u, const aligned_vector<int> &v) : ne(u.size())
{
    this->u = aligned_vector<int>(u);
    this->v = aligned_vector<int>(v);
    this->w = aligned_vector<double>(ne, 0.0);
}

EdgeList::EdgeList(const aligned_vector<int> &u,
                   const aligned_vector<int> &v,
                   const aligned_vector<double> &w)
    : ne(u.size())
{
    this->u = aligned_vector<int>(u);
    this->v = aligned_vector<int>(v);
    this->w = aligned_vector<double>(w);
}

bool EdgeList::isSymmetric()
{
    for(int i = 0; i < ne; ++i)
        if(u[i] > v[i])
        	return false;
    return true;
}

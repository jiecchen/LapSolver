#pragma once

class UnionFind {
  public:
    UnionFind(int n);
    ~UnionFind();
    void link(int u, int v);
    int find_set(int u);

  private:
    int* parent;
    int* to_compress;
};

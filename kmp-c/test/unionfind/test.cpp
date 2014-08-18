// Test for UnionFind

#include <cstdlib>
#include <iostream>
#include <string>
#include <algorithms/UnionFind.h>

using namespace std;

int main(int argc, char *argv[])
{
    // read n and number of links
    int n, q;
    cin >> n >> q;

    UnionFind uf(n);

    // perform links
    for (int i = 0; i <= q; i++) {
        int u, v;
        cin >> u >> v;
        uf.link(u, v);
        
        // report state after each link
        for (int j = 0; j < n; j++) {
            cout << uf.find_set(j) << ' ';
        } cout << '\n';
    }

    return 0;
}

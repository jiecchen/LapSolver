// Test for GraphLoader

#include <cstdlib>
#include <string>
#include <structures/Graph.h>
#include <structures/GraphLoader.h>

int main(int argc, char const *argv[])
{
    atexit(mkl_free_buffers);
    Graph g = GraphLoader::fromStdin();

    printf("%d vertices\n", g.nv);
    for (int i = 0; i < g.nv; i++) {
        printf("deg(%d) = %d\n", i, g.getDegree(i));
    }

    return 0;
}

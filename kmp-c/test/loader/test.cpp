// Test for GraphLoader

#include <cstdlib>
#include <string>
#include <structures/Graph.h>
#include <structures/GraphLoader.h>

int main(int argc, char const *argv[])
{
    Graph g;
    std::string inFileName;
    try
    {
        if (argc < 2)
            inFileName = "inG.ijv";
        else
            inFileName = argv[1];
        g = GraphLoader::fromFile(inFileName);
    }
    catch (int e)
    {
        if (e > 0)
            fprintf(stderr, "Parse error in '%s' on line %d\n", inFileName.c_str(), e + 1);
        else
            fprintf(stderr, "Failed to open file '%s'\n", inFileName.c_str());
        return 0;
    }

    printf("%d vertices\n", g.nv);
    for (int i = 0; i < g.nv; i++) {
        printf("deg(%d) = %d\n", i, g.getDegree(i));
    }

    return 0;
}

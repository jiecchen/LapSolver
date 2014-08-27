// Test for Priority Queue

#include <cstdlib>
#include <string>
#include <memory>
#include <structures/GraphLoader.h>
#include <algorithms/PriorityQueue.h>
#include <util/Benchmark.h>

using namespace std;

int main(int argc, char *argv[])
{
    atexit(mkl_free_buffers);

    int N, M;
    cin >> N >> M;

    std::shared_ptr<PriorityQueue> pq;

    auto bench = make_benchmark(argc, argv, [&] () {
        pq.reset(new PriorityQueue(N));

        for (int i = 0; i < M; i++) {
            string OP;
            cin >> OP;

            if (OP == "PUSH") {
                int x;
                double v;
                cin >> x >> v;

                pq->Push(x,v);
            }

            if (OP == "POP") {
                pair <int, double> p;
                p = pq->Pop();

                cout << p.first << " " << p.second << endl;
            }

            if (OP == "DEC") {
                int x;
                double v;
                cin >> x >> v;

                pq->DecreaseKey(x,v);    
            }
        }
    });

    

    return 0;
}

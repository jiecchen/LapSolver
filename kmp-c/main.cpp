#include <iostream>
#include <vector>
#include <cstdlib>
#include <sys/time.h>
#include <ctime>
#include <numeric>
#include "mkl.h"
#include "util/aligned.h"
using namespace std;

#define RANDOM (((double) rand()) / RAND_MAX)

long getTimeMS()
{
    struct timeval tv;

    gettimeofday(&tv, NULL);

    auto ret = tv.tv_usec;
    ret /= 1000;
    ret += (tv.tv_sec * 1000);

    return ret;
}

int main(int argc, char const *argv[])
{
    vmlSetMode(VML_EP);
    srand(time(NULL));

    const int SIZE = 1000000;
    aligned_vector<double> x(SIZE);
    aligned_vector<double> y(SIZE);
    for (int i = 0; i < SIZE; ++i) {
        x[i] = RANDOM;
        y[i] = RANDOM;
    }

    auto start = getTimeMS();
    auto z = x + y;
    auto stop = getTimeMS();
    cout <<  stop-start << "ms" << endl;

    return 0;
}
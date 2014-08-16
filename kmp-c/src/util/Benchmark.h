#pragma once
#include <string>
#include <chrono>
#include <stdexcept>
#include <cstdio>

template<typename Func>
class Benchmark
{
public:
    Benchmark(int argc, char **argv, Func f)
        : name(argv[0])
    {
        for (int i = 1; i < argc; ++i)
        {
            std::string flag(argv[i]);
            if (flag == "-v")
                echo = true;
            else if (flag == "-n" && i + 1 < argc)
            {
                try
                {
                    nTrials = std::stoi(argv[i + 1]);
                }
                catch (...)
                {
                    fprintf(stderr, "bad value for flag -n '%s'\n", argv[i + 1]);
                    break;
                }
            }
        }
        runBenchmark(f);
    }

    Benchmark(std::string name, Func f, bool echo = false, int nTrials = 1)
        : name(name), echo(echo), nTrials(nTrials)
    {
        runBenchmark(f);
    }

    ~Benchmark()
    {
        if (echo)
            fprintf(stderr, "%s completed in %dms (on average over %d trials)\n", name.c_str(), avgTime.count(), nTrials);
    }

    void runBenchmark(Func f)
    {
        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < nTrials; ++i)
            f();
        auto stop = std::chrono::high_resolution_clock::now();
        avgTime = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start) / nTrials;
    }

private:
    std::string name;
    bool echo = false;
    int nTrials = 1;
    std::chrono::milliseconds avgTime;
};

template<typename Func>
Benchmark<Func> make_benchmark(int argc, char **argv, Func f)
{
    return Benchmark<Func>(argc, argv, f);
}


#include <benchmark/benchmark.h>
#include <vector>
#include <list>
#include <random>
#include <numeric>

const long fromRange = 8;

const long toRange = 1 << 20;

static void BM_AccumulateList(benchmark::State &state)
{
    const size_t size = state.range(0);
    std::list<float> v1;
    std::default_random_engine generator;
    std::uniform_real_distribution<float> distribution(-100.0f, 100.0f);
    for (size_t i = 0; i < size; i++)
    {
        v1.push_back(distribution(generator));

    }
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(std::accumulate(v1.cbegin(), v1.cend(), 0.0f));
    }
}

BENCHMARK(BM_AccumulateList)->Range(fromRange, toRange);

static void BM_AccumulateVector(benchmark::State &state)
{
    const size_t size = state.range(0);
    std::vector<float> v1(size);
    std::default_random_engine generator;
    std::uniform_real_distribution<float> distribution(-100.0f, 100.0f);
    std::generate(v1.begin(), v1.end(), [&distribution, &generator]()
    { return distribution(generator); });
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(std::accumulate(v1.cbegin(), v1.cend(), 0.0f));
    }
}

BENCHMARK(BM_AccumulateVector)->Range(fromRange, toRange);

BENCHMARK_MAIN ();
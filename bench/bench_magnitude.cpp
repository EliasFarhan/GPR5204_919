//
// Created by efarhan on 19/02/2021.
//


#include <benchmark/benchmark.h>
#include <vector>
#include <random>
#include "vector3.h"
#include "bench_utils.h"

const long fromRange = 8;

const long toRange = 1 << 20;

static void BM_MagnitudeNaive(benchmark::State &state)
{
    std::vector<maths::Vec3f> v1(state.range(0));
    std::for_each(v1.begin(), v1.end(),[](maths::Vec3f& v){FillRandom(v);});
    for(auto _ : state)
    {
        for(auto& v : v1)
        {
            benchmark::DoNotOptimize(v.Magnitude());
        }
    }
}

BENCHMARK(BM_MagnitudeNaive)->Range(fromRange, toRange);

static void BM_MagnitudeAoSoA(benchmark::State &state)
{
    std::vector<maths::FourVec3f> v1(state.range(0));
    std::for_each(v1.begin(), v1.end(),[](maths::FourVec3f& v){FillRandom(v);});
    for(auto _ : state)
    {
        for(std::size_t i = 0; i < state.range(0)/4; i++)
        {
            benchmark::DoNotOptimize(v1[i].Magnitude());
        }
    }
}
BENCHMARK(BM_MagnitudeAoSoA)->Range(fromRange, toRange);

BENCHMARK_MAIN ();
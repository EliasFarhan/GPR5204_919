//
// Created by efarhan on 19/02/2021.
//

#include "bench_utils.h"

#include <random>

void FillRandom(maths::Vec3f& v)
{
    static std::random_device rd;  //Will be used to obtain a seed for the random number engine
    static std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    static std::uniform_real_distribution<float> dis(-100.0f, 100.0f);
    v[0] = dis(gen);
    v[1] = dis(gen);
    v[2] = dis(gen);
}

void FillRandom(maths::FourVec3f& v)
{
    static std::random_device rd;  //Will be used to obtain a seed for the random number engine
    static std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    static std::uniform_real_distribution<float> dis(-100.0f, 100.0f);
    for(auto& xs : v.Xs())
    {
        xs = dis(gen);
    }
    for(auto& ys : v.Ys())
    {
        ys = dis(gen);
    }
    for(auto& zs : v.Zs())
    {
        zs = dis(gen);
    }
}

void FillRandom(maths::EightVec3f& v)
{
    static std::random_device rd;  //Will be used to obtain a seed for the random number engine
    static std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    static std::uniform_real_distribution<float> dis(-100.0f, 100.0f);
    for(auto& xs : v.Xs())
    {
        xs = dis(gen);
    }
    for(auto& ys : v.Ys())
    {
        ys = dis(gen);
    }
    for(auto& zs : v.Zs())
    {
        zs = dis(gen);
    }
}
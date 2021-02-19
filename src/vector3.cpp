#include <vector3.h>
#include <instrinsics.h>

#include <cmath>

namespace maths
{

Vec3f::Vec3f(const float *ptr)
{
    values[0] = ptr[0];
    values[1] = ptr[1];
    values[2] = ptr[2];
}
float Vec3f::Magnitude() const
{
    return std::sqrt(x*x+y*y+z*z);
}

float Dot(const Vec3f &v1, const Vec3f &v2)
{
    return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
}

FourVec3f::FourVec3f(const Vec3f *ptr)
{
    for(int i = 0; i < 4; i++)
    {
        xs[i] = ptr[i][0];
        ys[i] = ptr[i][1];
        zs[i] = ptr[i][2];
    }
}


#if defined(__SSE__)
std::array<float, 4> FourVec3f::Dot(const FourVec3f &v1, const FourVec3f &v2)
{
    auto x1 = _mm_loadu_ps(v1.xs.data());
    auto y1 = _mm_loadu_ps(v1.ys.data());
    auto z1 = _mm_loadu_ps(v1.zs.data());

    auto x2 = _mm_loadu_ps(v2.xs.data());
    auto y2 = _mm_loadu_ps(v2.ys.data());
    auto z2 = _mm_loadu_ps(v2.zs.data());

    x1 = _mm_mul_ps(x1, x2);
    y1 = _mm_mul_ps(y1, y2);
    z1 = _mm_mul_ps(z1, z2);

    x1 = _mm_add_ps(x1, y1);
    x1 = _mm_add_ps(x1, z1);

    alignas(4 * sizeof(float))
    std::array<float, 4> result;
    _mm_store_ps(result.data(), x1);
    return result;
}
std::array<float, 4> FourVec3f::Magnitude() const
{
    auto x1 = _mm_loadu_ps(xs.data());
    auto y1 = _mm_loadu_ps(ys.data());
    auto z1 = _mm_loadu_ps(zs.data());

    x1 = _mm_mul_ps(x1, x1);
    y1 = _mm_mul_ps(y1, y1);
    z1 = _mm_mul_ps(z1, z1);

    x1 = _mm_add_ps(x1, y1);
    x1 = _mm_add_ps(x1, z1);
    x1 = _mm_sqrt_ps(x1);

    alignas(4 * sizeof(float))
    std::array<float, 4> result;
    _mm_store_ps(result.data(), x1);
    return result;
}
#else
std::array<float, 4> FourVec3f::Dot(const FourVec3f &v1, const FourVec3f &v2)
{
    alignas(4 * sizeof(float))
    std::array<float, 4> result;

    v4sf x1 = {v1.xs[0], v1.xs[1], v1.xs[2], v1.xs[3]};
    v4sf y1 = {v1.ys[0], v1.ys[1], v1.ys[2], v1.ys[3]};
    v4sf z1 = {v1.zs[0], v1.zs[1], v1.zs[2], v1.zs[3]};

    v4sf x2 = {v2.xs[0], v2.xs[1], v2.xs[2], v2.xs[3]};
    v4sf y2 = {v2.ys[0], v2.ys[1], v2.ys[2], v2.ys[3]};
    v4sf z2 = {v2.zs[0], v2.zs[1], v2.zs[2], v2.zs[3]};

    x1 = x1 * x2;
    y1 = y1 * y2;
    z1 = z1 * z2;

    x1 = x1 + y1;
    x1 = x1 + z1;

    result[0] = x1[0];
    result[1] = x1[1];
    result[2] = x1[2];
    result[3] = x1[3];

    return result;
}
#endif


EightVec3f::EightVec3f(const Vec3f *ptr)
{
    for(int i = 0; i < 8; i++)
    {
        xs[i] = ptr[i][0];
        ys[i] = ptr[i][1];
        zs[i] = ptr[i][2];
    }
}


#if defined(__AVX2__)
std::array<float, 8> EightVec3f::Dot(const EightVec3f &v1, const EightVec3f &v2)
{
    auto x1 = _mm256_loadu_ps(v1.xs.data());
    auto y1 = _mm256_loadu_ps(v1.ys.data());
    auto z1 = _mm256_loadu_ps(v1.zs.data());

    auto x2 = _mm256_loadu_ps(v2.xs.data());
    auto y2 = _mm256_loadu_ps(v2.ys.data());
    auto z2 = _mm256_loadu_ps(v2.zs.data());

    x1 = _mm256_mul_ps(x1, x2);
    y1 = _mm256_mul_ps(y1, y2);
    z1 = _mm256_mul_ps(z1, z2);

    x1 = _mm256_add_ps(x1, y1);
    x1 = _mm256_add_ps(x1, z1);

    alignas(8 * sizeof(float))
    std::array<float, 8> result;
    _mm256_store_ps(result.data(), x1);
    return result;
}
std::array<float, 8> EightVec3f::Magnitude() const
{
    auto x1 = _mm256_loadu_ps(xs.data());
    auto y1 = _mm256_loadu_ps(ys.data());
    auto z1 = _mm256_loadu_ps(zs.data());

    x1 = _mm256_mul_ps(x1, x1);
    y1 = _mm256_mul_ps(y1, y1);
    z1 = _mm256_mul_ps(z1, z1);

    x1 = _mm256_add_ps(x1, y1);
    x1 = _mm256_add_ps(x1, z1);
    x1 = _mm256_sqrt_ps(x1);

    alignas(8 * sizeof(float))
    std::array<float, 8> result;
    _mm256_store_ps(result.data(), x1);
    return result;
}
#endif
}
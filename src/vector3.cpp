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

Vec3f Vec3f::operator+(const Vec3f &v) const
{
    Vec3f result(x+v.x, y+v.y, z+v.z);
    return result;
}

Vec3f Vec3f::operator-(const Vec3f &v) const
{
    Vec3f result(x-v.x, y-v.y, z-v.z);
    return result;
}

Vec3f Vec3f::operator*(float f) const
{
    Vec3f result(x*f, y*f, z*f);
    return result;
}

Vec3f Vec3f::operator/(float f) const
{
    Vec3f result(x/f, y/f, z/f);
    return result;
}
Vec3f Vec3f::Normalized() const
{
    return *this/Magnitude();
}
float Vec3f::SqrMagnitude() const
{
    return Dot(*this, *this);
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


std::array<Vec3f, 4> FourVec3f::vectors() const
{
    std::array<Vec3f, 4> v;
    for(std::size_t i = 0; i < 4; i++)
    {
        v[i][0] = xs[i];
        v[i][1] = ys[i];
        v[i][2] = zs[i];
    }
    return v;
}
FourVec3f::FourVec3f(Vec3f v)
{
    for(std::size_t i = 0; i < 4; i++)
    {
        xs[i] = v.x;
        ys[i] = v.y;
        zs[i] = v.z;
    }
}

FourVec3f FourVec3f::operator+(const FourVec3f &v) const
{
    FourVec3f fv3f;
    for(std::size_t i = 0; i < 4; i++)
    {
        fv3f.xs[i] = xs[i] + v.xs[i];
        fv3f.ys[i] = ys[i] + v.ys[i];
        fv3f.zs[i] = zs[i] + v.zs[i];
    }
    return fv3f;
}

FourVec3f FourVec3f::operator-(const FourVec3f &v) const
{
    FourVec3f fv3f;
    for(std::size_t i = 0; i < 4; i++)
    {
        fv3f.xs[i] = xs[i] - v.xs[i];
        fv3f.ys[i] = ys[i] - v.ys[i];
        fv3f.zs[i] = zs[i] - v.zs[i];
    }
    return fv3f;
}

std::array<float, 4> FourVec3f::SqrMagnitude() const
{
    return Dot(*this, *this);
}


FourVec3f::FourVec3f(
    const std::array<float, 4>& xs,
    const std::array<float, 4>& ys,
    const std::array<float, 4>& zs) :
    xs(xs), ys(ys), zs(zs)
{

}

FourVec3f FourVec3f::Normalized() const
{
    return *this/Magnitude();
}


FourVec3f FourVec3f::operator*(const std::array<float, 4>& values) const
{
    FourVec3f result;
    for(std::size_t i = 0; i < 4; i++)
    {
        result.xs[i] = xs[i]*values[i];
        result.ys[i] = ys[i]*values[i];
        result.zs[i] = zs[i]*values[i];
    }
    return result;

}
FourVec3f FourVec3f::operator/(const std::array<float, 4>& values) const
{
    FourVec3f result;
    for(std::size_t i = 0; i < 4; i++)
    {
        result.xs[i] = xs[i]/values[i];
        result.ys[i] = ys[i]/values[i];
        result.zs[i] = zs[i]/values[i];
    }
    return result;
}
FourVec3f FourVec3f::operator*(float value) const
{
    FourVec3f result;
    for(std::size_t i = 0; i < 4; i++)
    {
        result.xs[i] = xs[i]*value;
        result.ys[i] = ys[i]*value;
        result.zs[i] = zs[i]*value;
    }
    return result;
}

#if defined(__SSE__)
std::array<float, 4> Sqrt(const std::array<float, 4> &v)
{
    auto vs = _mm_loadu_ps(v.data());
    vs = _mm_sqrt_ps(vs);
    alignas(4 * sizeof(float))
    std::array<float, 4> result;
    _mm_store_ps(result.data(), vs);
    return result;
}
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

EightVec3f::EightVec3f(Vec3f v)
{
    for(int i = 0; i < 8; i++)
    {
        xs[i] = v.x;
        ys[i] = v.y;
        zs[i] = v.z;
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
std::array<Vec3f, 8> EightVec3f::vectors() const
{
    std::array<Vec3f, 8> v;
    for(std::size_t i = 0; i < 8; i++)
    {
        v[i][0] = xs[i];
        v[i][1] = ys[i];
        v[i][2] = zs[i];
    }
    return v;
}
#endif
}
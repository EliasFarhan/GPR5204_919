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
    const Vec3f result(x+v.x, y+v.y, z+v.z);
    return result;
}

Vec3f Vec3f::operator-(const Vec3f &v) const
{
    const Vec3f result(x-v.x, y-v.y, z-v.z);
    return result;
}

Vec3f Vec3f::operator*(float f) const
{
    const Vec3f result(x*f, y*f, z*f);
    return result;
}

Vec3f Vec3f::operator/(float f) const
{
    const Vec3f result(x/f, y/f, z/f);
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




#if defined(__SSE__)


FourVec3f FourVec3f::operator+(const FourVec3f& v) const
{
    FourVec3f fv3f;
    auto x1 = _mm_loadu_ps(xs.data());
    auto y1 = _mm_loadu_ps(ys.data());
    auto z1 = _mm_loadu_ps(zs.data());

    const auto x2 = _mm_loadu_ps(v.xs.data());
    const auto y2 = _mm_loadu_ps(v.ys.data());
    const auto z2 = _mm_loadu_ps(v.zs.data());

    x1 = _mm_add_ps(x1, x2);
    y1 = _mm_add_ps(y1, y2);
    z1 = _mm_add_ps(z1, z2);

    _mm_storeu_ps(fv3f.xs.data(), x1);
    _mm_storeu_ps(fv3f.ys.data(), y1);
    _mm_storeu_ps(fv3f.zs.data(), z1);
    return fv3f;
}

FourVec3f FourVec3f::operator-(const FourVec3f& v) const
{
    FourVec3f fv3f;
    auto x1 = _mm_loadu_ps(xs.data());
    auto y1 = _mm_loadu_ps(ys.data());
    auto z1 = _mm_loadu_ps(zs.data());

    const auto x2 = _mm_loadu_ps(v.xs.data());
    const auto y2 = _mm_loadu_ps(v.ys.data());
    const auto z2 = _mm_loadu_ps(v.zs.data());

    x1 = _mm_sub_ps(x1, x2);
    y1 = _mm_sub_ps(y1, y2);
    z1 = _mm_sub_ps(z1, z2);

    _mm_storeu_ps(fv3f.xs.data(), x1);
    _mm_storeu_ps(fv3f.ys.data(), y1);
    _mm_storeu_ps(fv3f.zs.data(), z1);
    return fv3f;
}

FourVec3f FourVec3f::operator*(const std::array<float, 4>& values) const
{
    FourVec3f result;
    auto x = _mm_loadu_ps(xs.data());
    auto y = _mm_loadu_ps(ys.data());
    auto z = _mm_loadu_ps(zs.data());
    const auto v = _mm_loadu_ps(values.data());

    x = _mm_mul_ps(x, v);
    y = _mm_mul_ps(y, v);
    z = _mm_mul_ps(z, v);

    _mm_storeu_ps(result.Xs().data(), x);
    _mm_storeu_ps(result.Ys().data(), y);
    _mm_storeu_ps(result.Zs().data(), z);
    return result;

}
FourVec3f FourVec3f::operator/(const std::array<float, 4>& values) const
{
    FourVec3f result;
    auto x = _mm_loadu_ps(xs.data());
    auto y = _mm_loadu_ps(ys.data());
    auto z = _mm_loadu_ps(zs.data());
    const auto v = _mm_loadu_ps(values.data());

    x = _mm_div_ps(x, v);
    y = _mm_div_ps(y, v);
    z = _mm_div_ps(z, v);

    _mm_storeu_ps(result.Xs().data(), x);
    _mm_storeu_ps(result.Ys().data(), y);
    _mm_storeu_ps(result.Zs().data(), z);
    return result;
}
FourVec3f FourVec3f::operator*(float value) const
{
    FourVec3f result;
    auto x = _mm_loadu_ps(xs.data());
    auto y = _mm_loadu_ps(ys.data());
    auto z = _mm_loadu_ps(zs.data());
    const auto v = _mm_load1_ps(&value);

    x = _mm_mul_ps(x, v);
    y = _mm_mul_ps(y, v);
    z = _mm_mul_ps(z, v);

    _mm_storeu_ps(result.Xs().data(), x);
    _mm_storeu_ps(result.Ys().data(), y);
    _mm_storeu_ps(result.Zs().data(), z);
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
    
    std::array<float, 4> result;
    _mm_storeu_ps(result.data(), x1);
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
    
    std::array<float, 4> result;
    _mm_storeu_ps(result.data(), x1);
    return result;
}
template<>
std::array<float, 4> Multiply(const std::array<float, 4> &v1, const std::array<float, 4> &v2)
{
    auto v1s = _mm_loadu_ps(v1.data());
    auto v2s = _mm_loadu_ps(v2.data());
    v1s = _mm_mul_ps(v1s, v2s);

    std::array<float, 4> result;
    _mm_storeu_ps(result.data(), v1s);
    return result;
}
template<>
std::array<float, 4> Multiply(const std::array<float, 4> &v1, float value)
{
    auto v1s = _mm_loadu_ps(v1.data());
    auto v2 = _mm_load1_ps(&value);
    v1s = _mm_mul_ps(v1s, v2);

    std::array<float, 4> result;
    _mm_storeu_ps(result.data(), v1s);
    return result;
}
template<>
std::array<float, 4> Sqrt(const std::array<float, 4> &v)
{
    auto vs = _mm_loadu_ps(v.data());
    vs = _mm_sqrt_ps(vs);

    std::array<float, 4> result;
    _mm_storeu_ps(result.data(), vs);
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

EightVec3f::EightVec3f(const std::array<float, 8> &xs, const std::array<float, 8> &ys, const std::array<float, 8> &zs) :
    xs(xs), ys(ys), zs(zs)
{

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


EightVec3f EightVec3f::Normalized() const
{
    return *this/Magnitude();
}

#if defined(__AVX2__)


EightVec3f EightVec3f::operator+(const EightVec3f& v) const
{
    EightVec3f result;

    auto x1 = _mm256_loadu_ps(xs.data());
    auto y1 = _mm256_loadu_ps(xs.data());
    auto z1 = _mm256_loadu_ps(xs.data());

    const auto x2 = _mm256_loadu_ps(v.xs.data());
    const auto y2 = _mm256_loadu_ps(v.xs.data());
    const auto z2 = _mm256_loadu_ps(v.xs.data());

    x1 = _mm256_add_ps(x1, x2);
    y1 = _mm256_add_ps(y1, y2);
    z1 = _mm256_add_ps(z1, z2);

    _mm256_storeu_ps(result.xs.data(), x1);
    _mm256_storeu_ps(result.ys.data(), y1);
    _mm256_storeu_ps(result.zs.data(), z1);

    return result;
}

EightVec3f EightVec3f::operator-(const EightVec3f& v) const
{
    EightVec3f result;

    auto x1 = _mm256_loadu_ps(xs.data());
    auto y1 = _mm256_loadu_ps(xs.data());
    auto z1 = _mm256_loadu_ps(xs.data());

    const auto x2 = _mm256_loadu_ps(v.xs.data());
    const auto y2 = _mm256_loadu_ps(v.xs.data());
    const auto z2 = _mm256_loadu_ps(v.xs.data());

    x1 = _mm256_sub_ps(x1, x2);
    y1 = _mm256_sub_ps(y1, y2);
    z1 = _mm256_sub_ps(z1, z2);

    _mm256_storeu_ps(result.xs.data(), x1);
    _mm256_storeu_ps(result.ys.data(), y1);
    _mm256_storeu_ps(result.zs.data(), z1);

    return result;
}

EightVec3f EightVec3f::operator*(const std::array<float, 8>& values) const
{
    EightVec3f result;
    auto x = _mm256_loadu_ps(xs.data());
    auto y = _mm256_loadu_ps(ys.data());
    auto z = _mm256_loadu_ps(zs.data());
    const auto v = _mm256_loadu_ps(values.data());

    x = _mm256_mul_ps(x, v);
    y = _mm256_mul_ps(y, v);
    z = _mm256_mul_ps(z, v);

    _mm256_storeu_ps(result.Xs().data(), x);
    _mm256_storeu_ps(result.Ys().data(), y);
    _mm256_storeu_ps(result.Zs().data(), z);
    return result;
}

EightVec3f EightVec3f::operator*(float value) const
{
    EightVec3f result;
    auto x = _mm256_loadu_ps(xs.data());
    auto y = _mm256_loadu_ps(ys.data());
    auto z = _mm256_loadu_ps(zs.data());
    const auto v = _mm256_broadcast_ss(&value);

    x = _mm256_mul_ps(x, v);
    y = _mm256_mul_ps(y, v);
    z = _mm256_mul_ps(z, v);

    _mm256_storeu_ps(result.Xs().data(), x);
    _mm256_storeu_ps(result.Ys().data(), y);
    _mm256_storeu_ps(result.Zs().data(), z);
    return result;
}

EightVec3f EightVec3f::operator/(const std::array<float, 8>& values) const
{
    EightVec3f result;
    auto x = _mm256_loadu_ps(xs.data());
    auto y = _mm256_loadu_ps(ys.data());
    auto z = _mm256_loadu_ps(zs.data());
    const auto v = _mm256_loadu_ps(values.data());

    x = _mm256_div_ps(x, v);
    y = _mm256_div_ps(y, v);
    z = _mm256_div_ps(z, v);

    _mm256_storeu_ps(result.Xs().data(), x);
    _mm256_storeu_ps(result.Ys().data(), y);
    _mm256_storeu_ps(result.Zs().data(), z);
    return result;
}

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
    
    std::array<float, 8> result;
    _mm256_storeu_ps(result.data(), x1);
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
    
    std::array<float, 8> result;
    _mm256_storeu_ps(result.data(), x1);
    return result;
}

template<>
std::array<float, 8> maths::Multiply(const std::array<float, 8> &v1, const std::array<float, 8> &v2)
{
    auto v1s = _mm256_loadu_ps(v1.data());
    auto v2s = _mm256_loadu_ps(v2.data());
    v1s = _mm256_mul_ps(v1s, v2s);

    std::array<float, 8> result;
    _mm256_storeu_ps(result.data(), v1s);
    return result;
}

template<>
std::array<float, 8> Multiply(const std::array<float, 8>& v1, float value)
{
    auto v1s = _mm256_loadu_ps(v1.data());
    auto v2 = _mm256_broadcast_ss(&value);
    v1s = _mm256_mul_ps(v1s, v2);

    std::array<float, 8> result;
    _mm256_storeu_ps(result.data(), v1s);
    return result;
}

template<>
std::array<float, 8> Sqrt(const std::array<float, 8> &v)
{
    auto vs = _mm256_loadu_ps(v.data());
    vs = _mm256_sqrt_ps(vs);

    std::array<float, 8> result;
    _mm256_storeu_ps(result.data(), vs);
    return result;
}

#endif
}
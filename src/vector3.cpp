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
#endif

#if defined(__aarch64__)
FourVec3f FourVec3f::operator+(const FourVec3f& v) const
{
    FourVec3f fv3f;
    auto x1 = vld1q_f32(xs.data());
    auto y1 = vld1q_f32(ys.data());
    auto z1 = vld1q_f32(zs.data());

    const auto x2 = vld1q_f32(v.xs.data());
    const auto y2 = vld1q_f32(v.ys.data());
    const auto z2 = vld1q_f32(v.zs.data());

    x1 = vaddq_f32(x1, x2);
    y1 = vaddq_f32(y1, y2);
    z1 = vaddq_f32(z1, z2);

    vst1q_f32(fv3f.xs.data(), x1);
    vst1q_f32(fv3f.ys.data(), y1);
    vst1q_f32(fv3f.zs.data(), z1);
    return fv3f;
}

FourVec3f FourVec3f::operator-(const FourVec3f& v) const
{
    FourVec3f fv3f;
    auto x1 = vld1q_f32(xs.data());
    auto y1 = vld1q_f32(ys.data());
    auto z1 = vld1q_f32(zs.data());

    const auto x2 = vld1q_f32(v.xs.data());
    const auto y2 = vld1q_f32(v.ys.data());
    const auto z2 = vld1q_f32(v.zs.data());

    x1 = vsubq_f32(x1, x2);
    y1 = vsubq_f32(y1, y2);
    z1 = vsubq_f32(z1, z2);

    vst1q_f32(fv3f.xs.data(), x1);
    vst1q_f32(fv3f.ys.data(), y1);
    vst1q_f32(fv3f.zs.data(), z1);
    return fv3f;
}

FourVec3f FourVec3f::operator*(const std::array<float, 4>& values) const
{
    FourVec3f fv3f;
    auto x1 = vld1q_f32(xs.data());
    auto y1 = vld1q_f32(ys.data());
    auto z1 = vld1q_f32(zs.data());

    const auto x2 = vld1q_f32(values.data());
    const auto y2 = vld1q_f32(values.data());
    const auto z2 = vld1q_f32(values.data());

    x1 = vmulq_f32(x1, x2);
    y1 = vmulq_f32(y1, y2);
    z1 = vmulq_f32(z1, z2);

    vst1q_f32(fv3f.xs.data(), x1);
    vst1q_f32(fv3f.ys.data(), y1);
    vst1q_f32(fv3f.zs.data(), z1);
    return fv3f;

}
FourVec3f FourVec3f::operator/(const std::array<float, 4>& values) const
{
    FourVec3f fv3f;
    auto x1 = vld1q_f32(xs.data());
    auto y1 = vld1q_f32(ys.data());
    auto z1 = vld1q_f32(zs.data());

    const auto x2 = vld1q_f32(values.data());
    const auto y2 = vld1q_f32(values.data());
    const auto z2 = vld1q_f32(values.data());

    x1 = vdivq_f32(x1, x2);
    y1 = vdivq_f32(y1, y2);
    z1 = vdivq_f32(z1, z2);

    vst1q_f32(fv3f.xs.data(), x1);
    vst1q_f32(fv3f.ys.data(), y1);
    vst1q_f32(fv3f.zs.data(), z1);
    return fv3f;
}
FourVec3f FourVec3f::operator*(float value) const
{
	FourVec3f fv3f;
    auto x1 = vld1q_f32(xs.data());
    auto y1 = vld1q_f32(ys.data());
    auto z1 = vld1q_f32(zs.data());

    const auto v = vld1q_dup_f32(&value);

    x1 = vmulq_f32(x1, v);
    y1 = vmulq_f32(y1, v);
    z1 = vmulq_f32(z1, v);

    vst1q_f32(fv3f.xs.data(), x1);
    vst1q_f32(fv3f.ys.data(), y1);
    vst1q_f32(fv3f.zs.data(), z1);
    return fv3f;
}

std::array<float, 4> FourVec3f::Dot(const FourVec3f &v1, const FourVec3f &v2)
{
	FourVec3f fv3f;
    auto x1 = vld1q_f32(v1.xs.data());
    auto y1 = vld1q_f32(v1.ys.data());
    auto z1 = vld1q_f32(v1.zs.data());

    const auto x2 = vld1q_f32(v2.xs.data());
    const auto y2 = vld1q_f32(v2.ys.data());
    const auto z2 = vld1q_f32(v2.zs.data());
    

    x1 = vmulq_f32(x1, x2);
    y1 = vmulq_f32(y1, y2);
    z1 = vmulq_f32(z1, z2);

    x1 = vaddq_f32(x1, y1);
    x1 = vaddq_f32(x1, z1);
    
    std::array<float, 4> result;
    vst1q_f32(result.data(), x1);
    return result;
}
std::array<float, 4> FourVec3f::Magnitude() const
{
    auto x1 = vld1q_f32(xs.data());
    auto y1 = vld1q_f32(ys.data());
    auto z1 = vld1q_f32(zs.data());

    x1 = vmulq_f32(x1, x1);
    y1 = vmulq_f32(y1, y1);
    z1 = vmulq_f32(z1, z1);

    x1 = vaddq_f32(x1, y1);
    x1 = vaddq_f32(x1, z1);
    x1 = vsqrtq_f32(x1);
    
    std::array<float, 4> result;
    vst1q_f32(result.data(), x1);
    return result;
}
template<>
std::array<float, 4> Multiply(const std::array<float, 4> &v1, const std::array<float, 4> &v2)
{
    auto v1s = vld1q_f32(v1.data());
    auto v2s = vld1q_f32(v2.data());
    v1s = vmulq_f32(v1s, v2s);

    std::array<float, 4> result;
    vst1q_f32(result.data(), v1s);
    return result;
}
template<>
std::array<float, 4> Multiply(const std::array<float, 4> &v1, float value)
{
    auto v1s = vld1q_f32(v1.data());
    auto v2 = vld1q_dup_f32(&value);
    v1s = vmulq_f32(v1s, v2);

    std::array<float, 4> result;
    vst1q_f32(result.data(), v1s);
    return result;
}
template<>
std::array<float, 4> Sqrt(const std::array<float, 4> &v)
{
    auto vs = vld1q_f32(v.data());
    vs = vsqrtq_f32(vs);

    std::array<float, 4> result;
    vst1q_f32(result.data(), vs);
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

#if defined(__aarch64__)

EightVec3f EightVec3f::operator+(const EightVec3f& v) const
{
    EightVec3f result;
	for(int i = 0; i < 8; i++)
	{
		result.xs[i] = xs[i]+v.xs[i];
		result.ys[i] = ys[i]+v.ys[i];
		result.zs[i] = zs[i]+v.zs[i];
	}
    return result;
}

EightVec3f EightVec3f::operator-(const EightVec3f& v) const
{
    EightVec3f result;
	for(int i = 0; i < 8; i++)
	{
		result.xs[i] = xs[i]-v.xs[i];
		result.ys[i] = ys[i]-v.ys[i];
		result.zs[i] = zs[i]-v.zs[i];
	}
    return result;
}

EightVec3f EightVec3f::operator*(const std::array<float, 8>& values) const
{
    EightVec3f result;
	for(int i = 0; i < 8; i++)
	{
		result.xs[i] = xs[i]*values[i];
		result.ys[i] = ys[i]*values[i];
		result.zs[i] = zs[i]*values[i];
	}
    return result;
}

EightVec3f EightVec3f::operator*(float value) const
{
    EightVec3f result;
	for(int i = 0; i < 8; i++)
	{
		result.xs[i] = xs[i]*value;
		result.ys[i] = ys[i]*value;
		result.zs[i] = zs[i]*value;
	}
    return result;
}

EightVec3f EightVec3f::operator/(const std::array<float, 8>& values) const
{
    EightVec3f result;
	for(int i = 0; i < 8; i++)
	{
		result.xs[i] = xs[i]/values[i];
		result.ys[i] = ys[i]/values[i];
		result.zs[i] = zs[i]/values[i];
	}
    return result;
}

std::array<float, 8> EightVec3f::Dot(const EightVec3f &v1, const EightVec3f &v2)
{
	std::array<float, 8> result;
	for(int i = 0; i < 8; i++)
	{
		result[i] = v1.xs[i]*v2.xs[i]+v1.ys[i]*v2.ys[i]+v1.zs[i]*v2.zs[i];
	}
    return result;
}
std::array<float, 8> EightVec3f::Magnitude() const
{
    std::array<float, 8> result;
	for(int i = 0; i < 8; i++)
	{
		result[i] = std::sqrt(
		xs[i]*xs[i]+ys[i]*ys[i]+zs[i]*zs[i]);
	}
    return result;
}

template<>
std::array<float, 8> maths::Multiply(const std::array<float, 8> &v1, const std::array<float, 8> &v2)
{
	std::array<float, 8> result;
	for(int i = 0; i < 8; i++)
	{
		result[i] = v1[i]*v2[i];
	}
    return result;
}

template<>
std::array<float, 8> Multiply(const std::array<float, 8>& v1, float value)
{
    std::array<float, 8> result;
	for(int i = 0; i < 8; i++)
	{
		result[i] = v1[i]*value;
	}
    return result;
}

template<>
std::array<float, 8> Sqrt(const std::array<float, 8> &v)
{
    std::array<float, 8> result;
	for(int i = 0; i < 8; i++)
	{
		result[i] = std::sqrt(v[i]);
	}
    return result;
}

#endif
}

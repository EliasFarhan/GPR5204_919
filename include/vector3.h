#pragma once
#include <array>

namespace maths
{
union Vec3f
{
    struct
    {
        float x,y,z;
    };
    std::array<float, 3> values{};

    constexpr Vec3f() : x(0), y(0), z(0){}
    Vec3f(float x, float y, float z) : x(x), y(y), z(z){}
    Vec3f(const float* ptr);

    Vec3f operator+(const Vec3f& v) const;
    Vec3f operator-(const Vec3f& v) const;
    Vec3f operator*(float f) const;
    Vec3f operator/(float f) const;

    float& operator[](std::size_t idx) { return values[idx]; }
    float operator[](std::size_t idx) const { return values[idx]; }

    [[nodiscard]] float Magnitude() const;
    [[nodiscard]] float SqrMagnitude() const;
    Vec3f Normalized() const;
    static const maths::Vec3f zero;
};

inline Vec3f const Vec3f::zero = Vec3f();

float Dot(const Vec3f& v1, const Vec3f& v2);

class FourVec3f
{
public:
    FourVec3f() = default;
    FourVec3f(const Vec3f *ptr);
    FourVec3f(Vec3f v);
    FourVec3f(const std::array<float, 4>& xs, const std::array<float, 4>& ys, const std::array<float, 4>& zs);

    [[nodiscard]] std::array<float, 4> Magnitude() const;
    [[nodiscard]] std::array<float, 4> SqrMagnitude() const;
    static std::array<float, 4> Dot(const FourVec3f &v1, const FourVec3f &v2);
    FourVec3f Normalized() const;

    [[nodiscard]] const std::array<float, 4> & Xs() const {return xs;}
    std::array<float, 4> & Xs() {return xs;}
    std::array<float, 4>& Ys(){return ys;}
    [[nodiscard]] const std::array<float, 4> & Zs() const {return zs;}
    std::array<float, 4> & Zs() {return zs;}

    FourVec3f operator+(const FourVec3f& v) const;
    FourVec3f operator-(const FourVec3f& v) const;
    FourVec3f operator*(const std::array<float, 4>& values) const;
    FourVec3f operator*(float value) const;
    FourVec3f operator/(const std::array<float, 4>& values) const;

    std::array<Vec3f, 4> vectors() const;
private:
    std::array<float, 4> xs{};
    std::array<float, 4> ys{};
    std::array<float, 4> zs{};
};

class EightVec3f
{
public:
    EightVec3f() = default;
    EightVec3f(const Vec3f *ptr);
    EightVec3f(Vec3f v);

    [[nodiscard]] std::array<float, 8> Magnitude() const;
    static std::array<float, 8> Dot(const EightVec3f &v1, const EightVec3f &v2);
    std::array<float, 8>& Xs(){return xs;}
    std::array<float, 8>& Ys(){return ys;}
    std::array<float, 8>& Zs(){return zs;}

    std::array<Vec3f, 8> vectors() const;
private:
    std::array<float, 8> xs{};
    std::array<float, 8> ys{};
    std::array<float, 8> zs{};
};

template<std::size_t N>
std::array<float, N> Multiply(const std::array<float, N>& v1, const std::array<float, N>& v2)
{
    std::array<float, N> result;
    for(std::size_t i = 0; i < N; i++)
    {
        result[i] = v1[i] * v2[i];
    }
    return result;
}

template<std::size_t N>
std::array<float, N> Multiply(const std::array<float, N>& v1, float value)
{
    std::array<float, N> result;
    for(std::size_t i = 0; i < N; i++)
    {
        result[i] = v1[i] * value;
    }
    return result;
}

template<std::size_t N>
std::array<float, N> Inverse(const std::array<float,N>& v)
{
    std::array<float, N> result;
    for(std::size_t i = 0; i < N; i++)
    {
        result[i] = 1.0f / v[i];
    }
    return result;
}

template<std::size_t N>
std::array<float, N> Negative(const std::array<float,N>& v)
{
    std::array<float, N> result;
    for(std::size_t i = 0; i < N; i++)
    {
        result[i] = -v[i];
    }
    return result;
}

std::array<float, 4> Sqrt(const std::array<float,4>& v);
}
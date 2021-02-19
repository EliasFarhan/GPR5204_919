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

    Vec3f() : x(0), y(0), z(0){}
    Vec3f(float x, float y, float z) : x(x), y(y), z(z){}
    Vec3f(const float* ptr);

    Vec3f operator+(const Vec3f& v);
    Vec3f operator-(const Vec3f& v);

    float& operator[](std::size_t idx) { return values[idx]; }
    float operator[](std::size_t idx) const { return values[idx]; }

    float Magnitude() const;

};

float Dot(const Vec3f& v1, const Vec3f& v2);

class FourVec3f
{
public:
    FourVec3f() = default;
    FourVec3f(const Vec3f *ptr);

    std::array<float, 4> Magnitude() const;
    static std::array<float, 4> Dot(const FourVec3f &v1, const FourVec3f &v2);
    std::array<float, 4>& Xs(){return xs;}
    std::array<float, 4>& Ys(){return ys;}
    std::array<float, 4>& Zs(){return zs;}
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

    std::array<float, 8> Magnitude() const;
    static std::array<float, 8> Dot(const EightVec3f &v1, const EightVec3f &v2);
    std::array<float, 8>& Xs(){return xs;}
    std::array<float, 8>& Ys(){return ys;}
    std::array<float, 8>& Zs(){return zs;}
private:
    std::array<float, 8> xs{};
    std::array<float, 8> ys{};
    std::array<float, 8> zs{};
};
}
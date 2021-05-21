//
// Created by efarhan on 23/02/2021.
//

#pragma once

#include <array>
#include <cstdlib>
#include <cmath>

#include "instrinsics.h"

namespace maths
{
template<std::size_t N>
class FloatArray
{
public:
    FloatArray() noexcept = default;
    FloatArray(const std::array<float, N>& arg)  noexcept : values_(arg){}

    float& operator[](std::size_t idx) noexcept { return values_[idx]; }
    const float& operator[](std::size_t idx) const noexcept { return values_[idx]; }

    FloatArray operator+(const FloatArray& rhs) const noexcept;
    FloatArray operator-(const FloatArray& rhs) const noexcept;
    FloatArray operator-() const noexcept;
    FloatArray operator*(const FloatArray& rhs) const noexcept;
    FloatArray operator*(float rhs) const noexcept;
    static FloatArray Sqrt(const FloatArray& rhs) noexcept;
    static FloatArray Reciprocal(const FloatArray& rhs) noexcept;

    [[nodiscard]] const std::array<float, N>& array() const noexcept {return values_;}
private:
    std::array<float, N> values_;
};

template<std::size_t N>
FloatArray<N> FloatArray<N>::operator+(const FloatArray &rhs) const noexcept
{
    FloatArray<N> result;
    for(std::size_t i = 0; i < N; i++)
    {
        result[i] = values_[i]+rhs[i];
    }
    return result;
}

template<std::size_t N>
FloatArray<N> FloatArray<N>::operator-(const FloatArray &rhs) const noexcept
{
    FloatArray<N> result;
    for(std::size_t i = 0; i < N; i++)
    {
        result[i] = values_[i]-rhs[i];
    }
    return result;
}

template<std::size_t N>
FloatArray<N> FloatArray<N>::operator*(const FloatArray &rhs) const noexcept
{
    FloatArray<N> result;
    for(std::size_t i = 0; i < N; i++)
    {
        result[i] = values_[i]*rhs[i];
    }
    return result;
}

template<std::size_t N>
FloatArray<N> FloatArray<N>::operator*(float rhs) const noexcept
{
    FloatArray<N> result;
    for(std::size_t i = 0; i < N; i++)
    {
        result[i] = values_[i]*rhs;
    }
    return result;
}

template<std::size_t N>
FloatArray<N> FloatArray<N>::Sqrt(const FloatArray &rhs) noexcept
{
    FloatArray<N> result;
    for(std::size_t i = 0; i < N; i++)
    {
        result[i] = std::sqrt(rhs[i]);
    }
    return result;
}
template<std::size_t N>
FloatArray<N> FloatArray<N>::Reciprocal(const FloatArray &rhs) noexcept
{
    FloatArray<N> result;
    for(std::size_t i = 0; i < N; i++)
    {
        result[i] = 1.0f / rhs[i];
    }
    return result;
}

template<std::size_t N>
FloatArray<N> FloatArray<N>::operator-() const noexcept
{
    FloatArray<N> result;
    for(std::size_t i = 0; i < N; i++)
    {
        result[i] = -values_[i];
    }
    return result;
}

using FourFloat = FloatArray<4>;
using EightFloat = FloatArray<8>;

#if defined(__SSE__)
template<>
inline FourFloat FourFloat::Sqrt(const FloatArray<4>& rhs) noexcept
{
    auto vs = _mm_loadu_ps(&rhs[0]);
    vs = _mm_sqrt_ps(vs);

    FourFloat result;
    _mm_storeu_ps(&result[0], vs);
    return result;
}

template<>
inline FourFloat FourFloat::operator*(const FloatArray<4>& rhs) const noexcept
{
    auto v1s = _mm_loadu_ps(values_.data());
    auto v2s = _mm_loadu_ps(rhs.values_.data());
    v1s = _mm_mul_ps(v1s, v2s);

    FourFloat result;
    _mm_storeu_ps(result.values_.data(), v1s);
    return result;
}

template<>
inline FourFloat FourFloat::operator*(float rhs) const noexcept
{
    auto v1s = _mm_loadu_ps(values_.data());
    auto v2 = _mm_load1_ps(&rhs);
    v1s = _mm_mul_ps(v1s, v2);

    FourFloat result;
    _mm_storeu_ps(result.values_.data(), v1s);
    return result;
}

template<>
inline FourFloat FourFloat::Reciprocal(const FloatArray<4>& rhs) noexcept
{
    auto vs = _mm_loadu_ps(&rhs[0]);
    vs = _mm_rcp_ps(vs);

    FourFloat result;
    _mm_storeu_ps(&result[0], vs);
    return result;
}
#endif

#if defined(__AVX2__)
template<>
inline EightFloat EightFloat::Sqrt(const EightFloat& rhs) noexcept
{
    auto vs = _mm256_loadu_ps(rhs.values_.data());
    vs = _mm256_sqrt_ps(vs);

    EightFloat result;
    _mm256_storeu_ps(result.values_.data(), vs);
    return result;
}

template<>
inline EightFloat EightFloat::operator*(const EightFloat& rhs) const noexcept
{
    auto v1s = _mm256_loadu_ps(values_.data());
    auto v2s = _mm256_loadu_ps(rhs.values_.data());
    v1s = _mm256_mul_ps(v1s, v2s);

    EightFloat result;
    _mm256_storeu_ps(result.values_.data(), v1s);
    return result;
}

template<>
inline EightFloat EightFloat::operator*(float rhs) const noexcept
{
    auto v1s = _mm256_loadu_ps(values_.data());
    auto v2 = _mm256_broadcast_ss(&rhs);
    v1s = _mm256_mul_ps(v1s, v2);

    EightFloat result;
    _mm256_storeu_ps(result.values_.data(), v1s);
    return result;
}

template<>
inline EightFloat EightFloat::Reciprocal(const EightFloat& rhs) noexcept
{
    auto vs = _mm256_loadu_ps(rhs.values_.data());
    vs = _mm256_rcp_ps(vs);

    EightFloat result;
    _mm256_storeu_ps(result.values_.data(), vs);
    return result;
}

#endif

#if defined(__aarch64__)
template<>
inline FourFloat FourFloat::Sqrt(const FloatArray<4>& rhs) noexcept
{
    auto vs = vld1q_f32(rhs.values_.data());
    vs = vsqrtq_f32(vs);

    FourFloat result;
    vst1q_f32(&result[0], vs);
    return result;
}

template<>
inline FourFloat FourFloat::operator*(const FloatArray<4>& rhs) const noexcept
{
    auto v1s = vld1q_f32(values_.data());
    auto v2s = vld1q_f32(rhs.values_.data());
    v1s = vmulq_f32(v1s, v2s);

    FourFloat result;
    vst1q_f32(&result[0], v1s);
    return result;
}

template<>
inline FourFloat FourFloat::operator*(float rhs) const noexcept
{
    auto v1s = vld1q_f32(values_.data());
    auto v2 = vld1q_dup_f32(&rhs);
    v1s = vmulq_f32(v1s, v2);

    FourFloat result;
    vst1q_f32(&result[0], v1s);
    return result;
}
#endif
}
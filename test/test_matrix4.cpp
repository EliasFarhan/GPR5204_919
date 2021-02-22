//
// Created by efarhan on 2/22/21.
//

#include <gtest/gtest.h>
#include <math_utils.h>

#include "matrix4.h"

TEST(Matrix, Inverse)
{
    maths::Mat4f identity = maths::Mat4f::Identity;
    const auto result = identity.Inverse();
    const auto resultOpti = identity.InverseOpti();
    EXPECT_FLOAT_EQ(0.0f, maths::Mat4f::MatrixDiff(identity, result));
    EXPECT_FLOAT_EQ(0.0f, maths::Mat4f::MatrixDiff(identity, resultOpti));

    const maths::Mat4f m = maths::Mat4f(std::array<maths::Vec4f, 4>
                                              {
                                                      maths::Vec4f(1.0f,5.0f,9.0f,13.0f),
                                                      maths::Vec4f(2.0f,6.0f,-10.0f,14.0f),
                                                      maths::Vec4f(3.0f,-7.0f,11.0f,15.0f),
                                                      maths::Vec4f(4.0f,8.0f,12.0f,16.0f)
                                              });

    const auto mInvCalculus = m.Inverse();
    const auto mInvCalculusOpti = m.InverseOpti();

    const maths::Mat4f mInv = maths::Mat4f(std::array<maths::Vec4f, 4>
                                                 {
                                                         maths::Vec4f(-198.0f,7.0f,20.0f,136.0f),
                                                         maths::Vec4f(10.0f,0.0f,-30.0f,20.0f),
                                                         maths::Vec4f(14.0f,-21.0f,0.0f,7.0f),
                                                         maths::Vec4f(34.0f,14.0f,10.0f,-23.0f)
                                                 })*(1.0f/420.0f);


    EXPECT_TRUE(maths::Equal(maths::Mat4f::MatrixDiff(mInvCalculus, mInv), 0.0f));
    EXPECT_TRUE(maths::Equal(maths::Mat4f::MatrixDiff(mInvCalculusOpti, mInv), 0.0f));
}
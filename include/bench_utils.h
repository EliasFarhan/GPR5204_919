#pragma once

#include <vector>

#include "vector3.h"
#include "matrix4.h"


void FillRandom(maths::Vec3f& v);
void FillRandom(maths::FourVec3f& v);
void FillRandom(maths::EightVec3f& v);
void FillRandom(maths::Mat4f& m);
void FillRandom(std::vector<int>& v, int low, int high);
#pragma once

#include <vector>

#include "vector3.h"


void FillRandom(maths::Vec3f& v);
void FillRandom(maths::FourVec3f& v);
void FillRandom(maths::EightVec3f& v);
void FillRandom(std::vector<int>& v, int low, int high);
#include "RGBImageStudent.hpp"
#include "IntensityImageStudent.hpp"

IntensityImageStudent RGBToIntensity(RGBImageStudent& before) {
	auto after = IntensityImageStudent(before.getWidth(), before.getHeight());
	for (int i = 0; i < (before.getWidth() * before.getHeight()); i++) {
		after.imageVector[i] = (before.imageVector[i].r + before.imageVector[i].g + before.imageVector[i].b)/3;
	}
}

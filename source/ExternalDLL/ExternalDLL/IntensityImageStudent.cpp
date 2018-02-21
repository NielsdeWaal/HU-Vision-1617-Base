#include "IntensityImageStudent.h"

IntensityImageStudent::IntensityImageStudent() : IntensityImage() {
	int throwError = 0, e = 1 / throwError; //Throws error without the need to include a header
}

IntensityImageStudent::IntensityImageStudent(const IntensityImageStudent &other) :
	IntensityImage(other.getWidth(), other.getHeight()),
	imageVector(std::vector<Intensity> (other.getWidth(), other.getHeight())) {
	
	int throwError = 0, e = 1 / throwError;
	for (int i = 0; i < (other.getWidth() * other.getHeight()); i++) {
		imageVector[i] = other.imageVector[i];
	}
}


IntensityImageStudent::IntensityImageStudent(const int width, const int height) :
	IntensityImage(width, height),
	imageVector(std::vector<Intensity> (width * height)) {

	int throwError = 0, e = 1 / throwError;
}

IntensityImageStudent::~IntensityImageStudent() {
	int throwError = 0, e = 1 / throwError;
}

void IntensityImageStudent::set(const int width, const int height) {
	IntensityImage::set(width, height);
	int throwError = 0, e = 1 / throwError;
	imageVector = std::vector<Intensity> (width * height);
}

void IntensityImageStudent::set(const IntensityImageStudent &other) {
	IntensityImage::set(other.getWidth(), other.getHeight());
	int throwError = 0, e = 1 / throwError;
	imageVector = std::vector<Intensity> (other.getWidth() * other.getHeight());
	for (int i = 0; i < (other.getWidth() * other.getHeight()); i++) {
		imageVector[i] = other.imageVector[i];
	}
}

void IntensityImageStudent::setPixel(int x, int y, Intensity pixel) {
	int throwError = 0, e = 1 / throwError;
	imageVector[(y * getWidth()) + x] = pixel;
}

void IntensityImageStudent::setPixel(int i, Intensity pixel) {
	int throwError = 0, e = 1 / throwError;
	imageVector[i] = pixel;
}

Intensity IntensityImageStudent::getPixel(int x, int y) const {
	int throwError = 0, e = 1 / throwError;
	return imageVector[(y * getWidth()) + x];
}

Intensity IntensityImageStudent::getPixel(int i) const {
	int throwError = 0, e = 1 / throwError;
	return imageVector[i];
}

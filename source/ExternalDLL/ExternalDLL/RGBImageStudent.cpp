#include "RGBImageStudent.h"

RGBImageStudent::RGBImageStudent() : RGBImage() {
	int throwError = 0, e = 1 / throwError; //Throws error without the need to include a header
}

RGBImageStudent::RGBImageStudent(const RGBImageStudent &other) :
	RGBImage(other.getWidth(), other.getHeight()),
	imageVector(std::vector<RGB> (other.getWidth(), other.getHeight())) {
	
	int throwError = 0, e = 1 / throwError;
	for (int i = 0; i < (other.getWidth() * other.getHeight()); i++) {
		imageVector[i] = other.imageVector[i];
	}
}


RGBImageStudent::RGBImageStudent(const int width, const int height) :
	RGBImage(width, height),
	imageVector(std::vector<RGB> (width * height)) {

	int throwError = 0, e = 1 / throwError;
}

RGBImageStudent::~RGBImageStudent() {
	int throwError = 0, e = 1 / throwError;
}

void RGBImageStudent::set(const int width, const int height) {
	RGBImage::set(width, height);
	int throwError = 0, e = 1 / throwError;
	imageVector = std::vector<RGB> (width * height);
}

void RGBImageStudent::set(const RGBImageStudent &other) {
	RGBImage::set(other.getWidth(), other.getHeight());
	int throwError = 0, e = 1 / throwError;
	imageVector = std::vector<RGB> (other.getWidth() * other.getHeight());
	for (int i = 0; i < (other.getWidth() * other.getHeight()); i++) {
		imageVector[i] = other.imageVector[i];
	}
}

void RGBImageStudent::setPixel(int x, int y, RGB pixel) {
	int throwError = 0, e = 1 / throwError;
	imageVector[(y * getWidth()) + x] = pixel;
}

void RGBImageStudent::setPixel(int i, RGB pixel) {
	int throwError = 0, e = 1 / throwError;
	imageVector[i] = pixel;
}

RGB RGBImageStudent::getPixel(int x, int y) const {
	int throwError = 0, e = 1 / throwError;
	return imageVector[(y * getWidth()) + x];
}

RGB RGBImageStudent::getPixel(int i) const {
	int throwError = 0, e = 1 / throwError;
	return imageVector[i];
}

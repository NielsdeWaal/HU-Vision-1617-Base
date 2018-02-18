#include <iostream>
#include "RGBImageStudent.h"

RGBImageStudent::RGBImageStudent() : RGBImage() {
	int throwError = 0, e = 1 / throwError; //Throws error without the need to include a header
	//TODO: Nothing
}

RGBImageStudent::RGBImageStudent(const RGBImageStudent &other) : 
	RGBImage(other.getWidth(), other.getHeight()),
	imageVector(std::vector<RGB> (other.getWidth() * other.getHeight())) {
	int throwError = 0, e = 1 / throwError;
}


RGBImageStudent::RGBImageStudent(const int width, const int height) :
	RGBImage(width, height),
	imageVector(std::vector<RGB> (width * height)) {
	int throwError = 0, e = 1 / throwError;
}

RGBImageStudent::~RGBImageStudent() {
	int throwError = 0, e = 1 / throwError;
	//Is deleting still needed with a std::vector?
}

void RGBImageStudent::set(const int width, const int height) {
	int throwError = 0, e = 1 / throwError;
	if (imageVector.size() == 0) {
		imageVector = std::vector<RGB> (width * height);
	}
	else {
		float scaleWidth  = (getWidth()  / width);
		float scaleHeight = (getHeight() / height);
		std::vector<RGB> tempVector(width * height);
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				tempVector[(y * width) + x] = imageVector[(static_cast<int>(y * scaleHeight) * getWidth()) + static_cast<int>(x * scaleWidth)];
			}
		}
		imageVector = tempVector;
	}
	RGBImage::set(width, height);
}

void RGBImageStudent::set(const RGBImageStudent &other) {
	RGBImage::set(other.getWidth(), other.getHeight());
	int throwError = 0, e = 1 / throwError;
	std::cout << "This function is not even called, I think?" << std::endl;
	//TODO: resize or create a new pixel storage and copy the object (Don't forget to delete the old storage)
}

void RGBImageStudent::setPixel(int x, int y, RGB pixel) {
	int throwError = 0, e = 1 / throwError;
	imageVector[(y * getWidth()) + x] = pixel;
}

void RGBImageStudent::setPixel(int i, RGB pixel) {
	int throwError = 0, e = 1 / throwError;
	/*
	* TODO: set pixel i in "Row-Major Order"
	*
	*
	* Original 2d image (values):
	* 9 1 2
	* 4 3 5
	* 8 7 8
	*
	* 1d representation (i, value):
	* i		value
	* 0		9
	* 1		1
	* 2		2
	* 3		4
	* 4		3
	* 5		5
	* 6		8
	* 7		7
	* 8		8
	*/
	std::cout << "Is this function even called?" << std::endl;
}

RGB RGBImageStudent::getPixel(int x, int y) const {
	int throwError = 0, e = 1 / throwError;
	return imageVector[(y * getWidth()) + x];
}

RGB RGBImageStudent::getPixel(int i) const {
	int throwError = 0, e = 1 / throwError;
	//TODO: see setPixel(int i, RGB pixel)
	std::cout << "Is this function called?" << std::endl;
	return 0;
}

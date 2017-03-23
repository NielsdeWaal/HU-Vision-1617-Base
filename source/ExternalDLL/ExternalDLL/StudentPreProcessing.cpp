#include "StudentPreProcessing.h"
#include "DefaultPreProcessing.h"

#include <cstdint>
#include <cmath>
#include "ImageFactory.h"

#include <iostream>

IntensityImage * StudentPreProcessing::stepToIntensityImage(const RGBImage &image) const {
    return nullptr;
}

std::tuple<uint,uint,double> getNewDimensions(const IntensityImage &image, size_t desiredPixelCount = 40000) {
    
    double scale = sqrt(40000 / ((double)image.getWidth() * (double)image.getHeight()));

    return std::make_tuple(scale*image.getWidth(), scale*image.getHeight(), scale);
}

IntensityImage *scaleNearestNeighbor(const IntensityImage &image) {

    auto dim = getNewDimensions(image);
	IntensityImage *out = ImageFactory::newIntensityImage(std::get<0>(dim), std::get<1>(dim));

    // Size of other image relative to ours.
    double otherScale = 1 / std::get<2>(dim);
    
    for (uint y = 0; y < std::get<1>(dim); ++y) {
        uint otherY = round(otherScale * y);
        for (uint x = 0; x < std::get<0>(dim); ++x) {
            uint otherX = round(otherScale * x);
            out->setPixel(x, y, image.getPixel(otherX, otherY));
            std::cout << "other: (" << otherX << "," << otherY << ") << " << otherScale << "\n";
        }
    }
    
    return out;
}

IntensityImage * StudentPreProcessing::stepScaleImage(const IntensityImage &image) const {

    if (image.getWidth() * image.getHeight() > 40000)
        return scaleNearestNeighbor(image);
    else
        return ImageFactory::newIntensityImage(image);

    // DefaultPreProcessing bla;
    // return bla.stepScaleImage(image);
}

IntensityImage * StudentPreProcessing::stepEdgeDetection(const IntensityImage &image) const {
    return nullptr;
}

IntensityImage * StudentPreProcessing::stepThresholding(const IntensityImage &image) const {
    return nullptr;
}

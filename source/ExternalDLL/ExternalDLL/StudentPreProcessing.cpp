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
    double origScale = 1 / std::get<2>(dim);
    
    for (uint y = 0; y < std::get<1>(dim); ++y) {

        uint origY = round(origScale * y);

        for (uint x = 0; x < std::get<0>(dim); ++x) {
            uint origX = round(origScale * x);
            out->setPixel(x, y, image.getPixel(origX, origY));
            //std::cout << "orig: (" << origX << "," << origY << ") << " << origScale << "\n";
        }
    }
    
    return out;
}

#include <chrono>

IntensityImage *StudentPreProcessing::stepScaleImage(const IntensityImage &image) const {

    // std::this_thread::sleep_for(std::chrono::milliseconds(4200));

    IntensityImage *ret = nullptr;

    using Clock = std::chrono::steady_clock;
    auto start  = Clock::now();

    // Switching here so we can measure the performance of the default
    // implementation as well.
    constexpr bool USE_STUDENT_SCALING = true;

    if (USE_STUDENT_SCALING) {
        if (image.getWidth() * image.getHeight() > 40000)
            ret = scaleNearestNeighbor(image);
        else
            ret = ImageFactory::newIntensityImage(image);
    } else {
        DefaultPreProcessing bla;
        ret = bla.stepScaleImage(image);
    }

    // DefaultPreProcessing bla;
    // return bla.stepScaleImage(image);

    auto end = Clock::now();

    auto duration = end - start;
    double scale = (double)Clock::period::num / Clock::period::den;
    std::cout << "Duration: " << ((double)duration.count() * (scale * 1000)) << "ms\n";

    return ret;

    // if (image.getWidth() * image.getHeight() > 40000)
    //     return scaleNearestNeighbor(image);
    // else
    //     return ImageFactory::newIntensityImage(image);

    // DefaultPreProcessing bla;
    // return bla.stepScaleImage(image);
}

IntensityImage * StudentPreProcessing::stepEdgeDetection(const IntensityImage &image) const {
    return nullptr;
}

IntensityImage * StudentPreProcessing::stepThresholding(const IntensityImage &image) const {
    return nullptr;
}

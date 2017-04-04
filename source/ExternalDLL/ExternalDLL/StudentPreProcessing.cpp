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

static IntensityImage *scaleNearestNeighbor(const IntensityImage &image) {

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

#include <thread>
#include <vector>

template<class F, class ...Args>
static void draadificeer(unsigned jobs, F fun, Args ...args) {
    auto cores = std::thread::hardware_concurrency();
    if (!cores)
        cores = 2;
    std::cerr << "got " << cores << " threads\n";

    size_t jobsPerThread = jobs / cores;

    std::vector<std::thread> pool;

    for (size_t i = 0; i < cores; i++) {
        // i * jobsPerThread,
        //     jobsPerThread + (i == cores-1 ? std::get<1>(dim) % jobsPerThread : 0))
        pool.emplace_back(std::bind(fun,
                                    i * jobsPerThread,
                                    jobsPerThread + (i == cores-1 ? jobs % jobsPerThread : 0),
                                    args...));
    }

    for (auto &t : pool)
        t.join();
}

static IntensityImage *scaleNearestNeighborMt(const IntensityImage &image,
                                              std::tuple<unsigned,unsigned,double> dim,
                                              IntensityImage *out) {

    const double origScale = 1 / std::get<2>(dim);


    draadificeer(std::get<1>(dim), [&image, &dim, origScale, out](size_t rowStart, size_t rowCount) {
            for (uint y = rowStart; y < rowStart + rowCount; ++y) {

                uint origY = round(origScale * y);

                for (uint x = 0; x < std::get<0>(dim); ++x) {
                    uint origX = round(origScale * x);
                    out->setPixel(x, y, image.getPixel(origX, origY));
                }
            }
        });

    return out;
}

#include <chrono>

IntensityImage *StudentPreProcessing::stepScaleImage(const IntensityImage &image) const {

    // std::this_thread::sleep_for(std::chrono::milliseconds(4200));

    IntensityImage *ret = nullptr;

    std::unique_ptr<DefaultPreProcessing> bla;

    // Switching here so we can measure the performance of the default
    // implementation as well.
    constexpr bool USE_STUDENT_SCALING = true;

    if (!USE_STUDENT_SCALING)
        bla = std::make_unique<DefaultPreProcessing>();

    using Clock = std::chrono::steady_clock;
    auto start  = Clock::now();

    if (USE_STUDENT_SCALING) {
        auto dim = getNewDimensions(image);
        IntensityImage *out = ImageFactory::newIntensityImage(std::get<0>(dim), std::get<1>(dim));

        if (image.getWidth() * image.getHeight() > 40000)
            ret = scaleNearestNeighborMt(image, dim, out);
        else
            ret = ImageFactory::newIntensityImage(image);
    } else {
        ret = bla->stepScaleImage(image);
    }

    auto end = Clock::now();

    auto duration = end - start;
    double scale = (double)Clock::period::num / Clock::period::den;
    std::cout << "Duration: " << ((double)duration.count() * (scale * 1000)) << "ms\n";

    return ret;
}

IntensityImage * StudentPreProcessing::stepEdgeDetection(const IntensityImage &image) const {
    return nullptr;
}

IntensityImage * StudentPreProcessing::stepThresholding(const IntensityImage &image) const {
    return nullptr;
}

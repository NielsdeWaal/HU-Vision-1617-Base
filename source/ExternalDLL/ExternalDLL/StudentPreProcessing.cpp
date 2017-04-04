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

constexpr double cubicInterpolate(const double p[4], const double x){
    return p[1] + 0.5 * x * (p[2] - p[0] + x * (2.0 * p[0] - 5.0 * p[1] + 4.0 * p[2] - p[3] + x * (3.0 * (p[1] - p[2]) + p[3] - p[0])));
}

double bicubicInterpolate(const double p[4][4], const double x, const double y){
    double h[4];
    for(uint8_t i = 0; i < 4; ++i){
        h[i] = cubicInterpolate(p[i], y);
    }
    return cubicInterpolate(h, x);
}

static IntensityImage *scaleBicubic(const IntensityImage &image,
                                              std::tuple<unsigned,unsigned,double> dim,
                                              IntensityImage *out) {

    const double origScale = 1 / std::get<2>(dim);


    draadificeer(std::get<1>(dim), [&image, &dim, origScale, out](size_t rowStart, size_t rowCount) {
            for (uint y = rowStart; y < rowStart + rowCount; ++y) {

                double origY = origScale * y;

                for (uint x = 0; x < std::get<0>(dim); ++x) {
                    double origX = origScale * x;

                    double p[4][4];

                    for(int row = 0; row < 4; ++row){
                        for(int col = 0; col < 4; ++col){
                            p[col][row] = image.getPixel((unsigned)origX - 1 + col, (unsigned)origY - 1 + row);
                        }
                    }

                    out->setPixel(x, y, bicubicInterpolate(p, origX-floor(origX), origY-floor(origY)));
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

        // Intensity testpuntjes[4][4] = {{0, 30, 50, 20}, {10, 40, 70, 20}, {10, 2, 50, 30}, {0, 20, 80, 30}};
        // std::cout << "Testbicubic: " << bicubicInterpolate(testpuntjes, 0.1, 0.1) << std::endl;
        // std::cout << "Testcubic: " << cubicInterpolate(testpuntjes[0], 0.1) << std::endl;


        if (image.getWidth() * image.getHeight() > 40000)
            ret = scaleBicubic(image, dim, out);
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

#include "StudentPreProcessing.h"
#include "DefaultPreProcessing.h"

#include <cstdint>
#include <cmath>
#include "ImageFactory.h"
#include "PixelType.h"
#include "RGBImage.h"

#include <iostream>
#include <functional>
#include <thread>
#include <vector>
#include <algorithm>

//#include <boost/compute.hpp>

//#define BOOST_COMPUTE_DEBUG_KERNEL_COMPILATION = 1

//namespace compute = boost::compute;

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


IntensityImage * StudentPreProcessing::stepToIntensityImage(const RGBImage &image) const {
	auto* intensityImage = new IntensityImageStudent(image.getWidth(), image.getHeight());

	for (int i = 0; i < (image.getWidth() * image.getHeight()); i++) {

		//Average
		//intensityImage->setPixel(i, ((image.getPixel(i).r + image.getPixel(i).g + image.getPixel(i).b) / 3));

		//Luma
		//intensityImage->setPixel(i, (image.getPixel(i).r * 0.2126 + image.getPixel(i).g * 0.7152 + image.getPixel(i).b * 0.0722));

		//Decomposition; max
		//intensityImage->setPixel(i, std::max({image.getPixel(i).r, image.getPixel(i).g, image.getPixel(i).b}));

		//Decomposition; min
		//intensityImage->setPixel(i, std::min({image.getPixel(i).r, image.getPixel(i).g, image.getPixel(i).b}));

		//Only red channel
		//intensityImage->setPixel(i, image.getPixel(i).r);

		//Only green channel
		//intensityImage->setPixel(i, image.getPixel(i).g);

		//Only blue channel
		//intensityImage->setPixel(i, image.getPixel(i).b);
	}

        /*
        // get the default device
        auto device = compute::system::default_device();
        auto context = compute::context(device);
        auto queue = compute::command_queue(context, device);

        std::vector<RGB> data;
        for (int i = 0; i < (image.getWidth() * image.getHeight()); i++) {
            data.push_back(image.getPixel(i));
        }

        compute::vector<RGB> pixels(data.size(), context);
        compute::copy(data.begin(), data.end(), pixels.begin(), queue);

        std::cout << "test" << std::endl;

        std::string source = BOOST_COMPUTE_STRINGIZE_SOURCE(
            __kernel void custom_kernel(__global const RGB *pixels,
                                        __global unsigned char *greyValues)
            {
                const uint i = get_global_id(0);
                const __global RGB *pixel = &pixels[i];

                //const float4 center = { 0, 0, 0, 0 };
                //const float4 position = { atom->x, atom->y, atom->z, 0 };
                const char4 resultVal = {(pixel->r + pixel->g + pixel->b) / 3};

                //distances[i] = distance(position, center);
                greyValues[i] = resultVal;
            }
        );

        std::cout << "test" << std::endl;

        // add type definition for Atom to the start of the program source
        source = compute::type_definition<RGB>() + "\n" + source;

        std::cout << "test" << std::endl;

        compute::program program = compute::program::build_with_source(source, context);

        std::cout << "test" << std::endl;

        compute::vector<char> greyResultValues(pixels.size(), context);

        std::cout << "test" << std::endl;

        compute::kernel custom_kernel = program.create_kernel("custom_kernel");
        custom_kernel.set_arg(0, pixels);
        custom_kernel.set_arg(1, greyResultValues);

        queue.enqueue_1d_range_kernel(custom_kernel, 0, pixels.size(), 1);

        for(int i = 0; i < greyResultValues.size() ; i++) {
            intensityImage->setPixel(i, greyResultValues[i]);
        }

	return intensityImage;

        //const int imageSize = image.getWidth() * image.getHeight();

        //std::vector<unsigned char> host_r;
        //std::vector<unsigned char> host_g;
        //std::vector<unsigned char> host_b;

        //std::vector<unsigned char> host_dest;

        //std::cout << "test" << std::endl;

        //for (int i = 0; i < (image.getWidth() * image.getHeight()); i++) {
        //    host_r.push_back(image.getPixel(i).r);
        //    host_g.push_back(image.getPixel(i).g);
        //    host_b.push_back(image.getPixel(i).b);
        //}

        //compute::vector<unsigned char> vector_r(host_r.size(), context);
        //compute::vector<unsigned char> vector_g(host_g.size(), context);
        //compute::vector<unsigned char> vector_b(host_b.size(), context);

        //size_t vecSize = vector_r.size();

        //compute::vector<unsigned char> vector_dest(vecSize);

        //assert(vector_r.size() == vector_dest.size());

        //auto R = compute::lambda::get<0>(compute::_1);
        //auto G = compute::lambda::get<1>(compute::_1);
        //auto B = compute::lambda::get<2>(compute::_1);

        //auto I = compute::lambda::get<3>(compute::_1);

        //std::cout << "test" << std::endl;

        ////std::cout << vector_r[0] << "\n";
        ////std::cout << vector_r[1] << "\n";
        ////std::cout << vector_r[2] << "\n";
        ////std::cout << vector_r[3] << std::endl;;

        //compute::for_each(
        //    compute::make_zip_iterator(
        //        boost::make_tuple(
        //            vector_r.begin(), vector_g.begin(), vector_b.begin(), vector_dest.begin()
        //        )
        //    ),
        //    compute::make_zip_iterator(
        //        boost::make_tuple(
        //            vector_r.end(), vector_g.end(), vector_b.end(), vector_dest.end()
        //        )
        //    ),
        //    compute::lambda::make_tuple(
        //        I = ((R + G + B)/ 3)
        //    ),
        //    queue
        //);

        //std::cout << "test" << std::endl;

        //compute::copy(vector_dest.begin(), vector_dest.end(), host_dest.begin(), queue);  
	*/
}

std::tuple<uint,uint,double> getNewDimensions(const IntensityImage &image, size_t desiredPixelCount = 40000) {
    
    double scale = sqrt(40000 / ((double)image.getWidth() * (double)image.getHeight()));

    return std::make_tuple(scale*image.getWidth(), scale*image.getHeight(), scale);
}

static IntensityImage *scaleNearestNeighbor(const IntensityImage &image,
                                            std::tuple<unsigned,unsigned,double> dim,
                                            IntensityImage *out) {

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

static IntensityImage *scaleBilinearMt(const IntensityImage &image,
                                       std::tuple<unsigned,unsigned,double> dim,
                                       IntensityImage *out) {

    const double origScale = 1 / std::get<2>(dim);

    draadificeer(std::get<1>(dim), [&image, &dim, origScale, out](size_t rowStart, size_t rowCount) {

            for (uint y = rowStart; y < rowStart + rowCount; ++y) {

                double origY = origScale * y;
                unsigned y1  = floor(origY);
                unsigned y2  =  ceil(origY);

                double yFrac = origY - (long)origY;

                for (uint x = 0; x < std::get<0>(dim); ++x) {
                    double origX = origScale * x;
                    unsigned x1  = floor(origX);
                    unsigned x2  =  ceil(origX);

                    double xFrac = origX - (long)origX;

                    auto px1y1 = image.getPixel(x1, y1);
                    auto px1y2 = image.getPixel(x1, y2);
                    auto px2y1 = image.getPixel(x2, y1);
                    auto px2y2 = image.getPixel(x2, y2);

                    auto h1 = px1y1 + xFrac * (px2y1 - px1y1);
                    auto h2 = px1y2 + xFrac * (px2y2 - px1y2);

                    auto v = h1 + yFrac * (h2 - h1);

                    out->setPixel(x, y, v);
                }
            }
        });

    return out;
}

double cubicInterpolate(const double p[4], double x){
    return p[1] + 0.5 * x * (p[2] - p[0] + x * (2.0 * p[0] - 5.0 * p[1] + 4.0 * p[2] - p[3] + x * (3.0 * (p[1] - p[2]) + p[3] - p[0])));
}

double bicubicInterpolate(const double p[4][4], double x, double y){
    double h[4];
    for(int i = 0; i < 4; ++i){
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

                    // NB: Here we assume that the type of Intensity is an unsigned char.
                    out->setPixel(x, y, std::max(0, std::min(255, (int)bicubicInterpolate(p, origX-(long)origX, origY-(long)origY))));
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

        if (image.getWidth() * image.getHeight() > 40000) {
            ret = scaleBicubic(image, dim, out);
            // ret = scaleNearestNeighborMt(image, dim, out);
            // ret = scaleBilinearMt(image, dim, out);
        } else {
            ret = ImageFactory::newIntensityImage(image);
        }
    } else {
        ret = bla->stepScaleImage(image);
    }

    auto end = Clock::now();

    auto duration = end - start;
    double scale = (double)Clock::period::num / Clock::period::den;

    // Verify that the clock period is short enough. We require 1µs precision.
    static_assert(Clock::period::num == 1 && Clock::period::den >= std::micro::den,
                  "Clock resolution too low for accurate performance testing");

    std::cout << "Duration: " << ((double)duration.count() * (scale * 1000)) << "ms\n";

    return ret;
}

IntensityImage * StudentPreProcessing::stepEdgeDetection(const IntensityImage &image) const {
    return nullptr;
}

IntensityImage * StudentPreProcessing::stepThresholding(const IntensityImage &image) const {
    return nullptr;
}

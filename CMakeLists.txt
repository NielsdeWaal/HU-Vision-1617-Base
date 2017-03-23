cmake_minimum_required(VERSION 2.6)

project(vissen)

enable_testing()

set(warnings
    "-Wall"
    "-Wextra"
    "-pedantic-errors"
    "-Wfatal-errors"
    "-Wcast-align"
    "-Wmissing-declarations"
    "-Wredundant-decls"
    "-Wuninitialized"
    "-Wno-unused-parameter"
    "-Wno-missing-field-initializers")

set(cxxflags
    "-std=c++14"
    "-g3"
    "-ggdb"
    "-O3"
    #"-flto"
    ${warnings})

string(REGEX REPLACE ";" " " cxxflags "${cxxflags}")

#set(ldflags
#    "-flto")

string(REGEX REPLACE ";" " " ldflags "${ldflags}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${cxxflags}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ldflags}")

set (sources
    source/ExternalDLL/ExternalDLL/export.cpp
    source/ExternalDLL/ExternalDLL/PostProcessing.cpp
    source/ExternalDLL/ExternalDLL/DefaultLocalization.cpp
    source/ExternalDLL/ExternalDLL/DefaultExtraction.cpp
    source/ExternalDLL/ExternalDLL/Feature.cpp
    source/ExternalDLL/ExternalDLL/GrayscaleAlgorithm.cpp
    source/ExternalDLL/ExternalDLL/SonnetXVIII.cpp
    source/ExternalDLL/ExternalDLL/IntensityImage.cpp
    source/ExternalDLL/ExternalDLL/FeatureMap.cpp
    source/ExternalDLL/ExternalDLL/ImageIO.cpp
    source/ExternalDLL/ExternalDLL/StudentExtraction.cpp
    source/ExternalDLL/ExternalDLL/RGBImageStudent.cpp
    source/ExternalDLL/ExternalDLL/DLLExecution.cpp
    source/ExternalDLL/ExternalDLL/DefaultPreProcessing.cpp
    source/ExternalDLL/ExternalDLL/HereBeDragons.cpp
    source/ExternalDLL/ExternalDLL/IntensityImageStudent.cpp
    source/ExternalDLL/ExternalDLL/StudentPreProcessing.cpp
    source/ExternalDLL/ExternalDLL/RGBImage.cpp
    source/ExternalDLL/ExternalDLL/Main.cpp
    source/ExternalDLL/ExternalDLL/StudentLocalization.cpp
    source/ExternalDLL/ExternalDLL/IntensityImagePrivate.cpp
    source/ExternalDLL/ExternalDLL/Image.cpp
    source/ExternalDLL/ExternalDLL/ImageFactory.cpp
    source/ExternalDLL/ExternalDLL/Representation.cpp
    source/ExternalDLL/ExternalDLL/RGBImagePrivate.cpp)

find_package(OpenCV REQUIRED)
add_executable(vissen ${sources})
target_link_libraries(vissen ${OpenCV_LIBS})

// Copyright (c) 2016 Giorgio Marcias
//
// This file is part of distance_transform, a C++11 implementation of the
// algorithm in "Distance Transforms of Sampled Functions"
// Pedro F. Felzenszwalb, Daniel P. Huttenlocher
// Theory of Computing, Vol. 8, No. 19, September 2012
//
// This source code is subject to Apache 2.0 License.
//
// Author: Giorgio Marcias
// email: marcias.giorgio@gmail.com

#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>

#include <distance_transform.hpp>
using namespace dt;

int main(int argc, char *argv[])
{
    std::vector<std::size_t> size({10, 10});
    MMArray<float, 2> f(size.data());
    MMArray<std::size_t, 2> indices(size.data());
    for (std::size_t i = 0; i < size[0]; ++i)
        for (std::size_t j = 0; j < size[1]; ++j) {
            if (i == j)
                f[i][j] = 0.0f;
            else
                f[i][j] = std::numeric_limits<float>::max();
        }

    DistanceTransform::initializeIndices(indices);  // this is not necessary, since distanceTransformL2() already does it

    std::cout << "indices:" << std::endl;
    for (std::size_t i = 0; i < size[0]; ++i) {
        for (std::size_t j = 0; j < size[1]; ++j)
            std::cout << std::setw(7) << indices[i][j] << ' ';
        std::cout << std::endl;
    }

    std::cout << std::endl << "Window [2:4;3:6]:" << std::endl;
    std::vector<std::size_t> winStart({2, 3});
    std::vector<std::size_t> winSize({3, 4});
    MArray<std::size_t, 2> win = indices.window(winStart.data(), winSize.data());
    for (std::size_t i = 0; i < win.size(); ++i) {
        for (std::size_t j = 0; j < win[i].size(); ++j)
            std::cout << std::setw(7) << win[i][j] << ' ';
        std::cout << std::endl;
    }

    std::cout << std::endl << "Slice 2 at dimension 0:" << std::endl;
    MArray<std::size_t, 1> sl = indices.slice(0, 2);
    for (std::size_t j = 0; j < sl.size(); ++j)
        std::cout << std::setw(7) << sl[j] << ' ';
    std::cout << std::endl << "Slice 2 at dimension 1:" << std::endl;
    indices.slice(1, 2, sl);
    for (std::size_t j = 0; j < sl.size(); ++j)
        std::cout << std::setw(7) << sl[j] << ' ';
    std::cout << std::endl;

    std::cout << std::endl << "Window [4:9] of slide 2 at dimension 1:" << std::endl;
    sl = sl.window(4, 6);
    for (std::size_t j = 0; j < sl.size(); ++j)
        std::cout << std::setw(7) << sl[j] << ' ';
    std::cout << std::endl << std::endl;

    std::cout << "f:" << std::endl;
    for (std::size_t i = 0; i < size[0]; ++i) {
        for (std::size_t j = 0; j < size[1]; ++j)
            std::cout << std::setw(4) << std::setprecision(1) << std::scientific << f[i][j] << ' ';
        std::cout << std::endl;
    }


    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    DistanceTransform::distanceTransformL2(f, f, indices, true);    // true for keeping squared distances, false for square roots

    std::cout << std::endl << "2D distance function computed in: " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count() << " ns." << std::endl;


    std::cout << std::endl << "D (squared):" << std::endl;
    for (std::size_t i = 0; i < size[0]; ++i) {
        for (std::size_t j = 0; j < size[1]; ++j)
            std::cout << std::setw(4) << std::setprecision(1) << std::fixed << f[i][j] << ' ';
        std::cout << std::endl;
    }

    std::cout << "indices:" << std::endl;
    for (std::size_t i = 0; i < size[0]; ++i) {
        for (std::size_t j = 0; j < size[1]; ++j)
            std::cout << std::setw(7) << indices[i][j] << ' ';
        std::cout << std::endl;
    }

    // 2D
    size = {320, 240};
    MMArray<float, 2> f2D(size.data());
    for (std::size_t i = 0; i < size[0]; ++i)
        for (std::size_t j = 0; j < size[1]; ++j)
            f2D[i][j] = std::numeric_limits<float>::max();
    f2D[0][0] = 0.0f;
    start = std::chrono::steady_clock::now();
    DistanceTransform::distanceTransformL2(f2D, f2D);
    std::cout << std::endl << size[0] << 'x' << size[1] << " distance function computed in: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count() << " ms." << std::endl;

    // 3D
    size = {400, 320, 240};
    MMArray<float, 3> f3D(size.data());
    for (std::size_t i = 0; i < size[0]; ++i)
        for (std::size_t j = 0; j < size[1]; ++j)
            for (std::size_t k = 0; k < size[2]; ++k)
                f3D[i][j][k] = std::numeric_limits<float>::max();
    f3D[0][0][0] = 0.0f;
    start = std::chrono::steady_clock::now();
    DistanceTransform::distanceTransformL2(f3D, f3D);
    std::cout << std::endl << size[0] << 'x' << size[1] << 'x' << size[2] << " distance function computed in: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count() << " ms." << std::endl;
    std::cout << "f3D[0][0][0] == 0.0 : " << std::boolalpha << (f3D[0][0][0] == 0.0f) << " -  f3D[399][319][239] == "  << std::sqrt(399.0f*399.0f + 319.0f*319.0f + 239.0f*239.0f) << " : " << std::boolalpha << (f3D[399][319][239] == std::sqrt(399.0f*399.0f + 319.0f*319.0f + 239.0f*239.0f)) << std::endl << std::endl;

    // 6D
    size = {5, 5, 5, 5, 5, 5};
    MMArray<float, 6> f6D(size.data());
    for (std::size_t i = 0; i < size[0]; ++i)
        for (std::size_t j = 0; j < size[1]; ++j)
            for (std::size_t k = 0; k < size[2]; ++k)
                for (std::size_t l = 0; l < size[3]; ++l)
                    for (std::size_t m = 0; m < size[4]; ++m)
                        for (std::size_t n = 0; n < size[5]; ++n)
                            f6D[i][j][k][l][m][n] = std::numeric_limits<float>::max();
    f6D[0][0][0][0][0][0] = 0.0f;
    start = std::chrono::steady_clock::now();
    DistanceTransform::distanceTransformL2(f6D, f6D);
    std::cout << std::endl << size[0] << 'x' << size[1] << 'x' << size[2] << 'x' << size[3] << 'x' << size[4] << 'x' << size[5] << " distance function computed in: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count() << " ms." << std::endl;
    std::cout << "f6D[0][0][0][0][0][0] == 0.0 : " << std::boolalpha << (f6D[0][0][0][0][0][0] == 0.0f) << " - f6D[4][4][4][4][4] == "  << std::sqrt(4.0f*4.0f * 6) << " : " << std::boolalpha << (f6D[4][4][4][4][4][4] == std::sqrt(4.0f*4.0f * 6)) << " - f6D[0][0][0][0][4] == "  << std::sqrt(4.0f*4.0f) << " : " << std::boolalpha << (f6D[0][0][0][0][0][4] == std::sqrt(4.0f*4.0f)) << std::endl << std::endl;

    return 0;
}

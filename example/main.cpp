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
    std::cout << "Slice 2:" << std::endl;
    MArray<std::size_t, 1> sl = indices.slice(1, 2);
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

    return 0;
}

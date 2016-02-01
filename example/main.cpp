// Copyright (c) 2015 Giorgio Marcias
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
#include <distance_transform.hpp>


int main(int argc, char *argv[])
{
    std::size_t size[2] = {5, 5};
    MMArray<float, 2> f(size);
    MMArray<std::size_t, 2> indices(size);
    for (std::size_t i = 0; i < size[0]; ++i)
        for (std::size_t j = 0; j < size[1]; ++j) {
            if (i == 1 || j == 3)
                f[i][j] = 0.0f;
            else
                f[i][j] = std::numeric_limits<float>::max();
            indices[i][j] = i * size[1] + j;
        }
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

    DistanceTransform::distanceTransformL2(f, f, indices, true);

    std::cout << std::endl << "D:" << std::endl;
    for (std::size_t i = 0; i < size[0]; ++i) {
        for (std::size_t j = 0; j < size[1]; ++j)
            std::cout << std::setw(4) << std::setprecision(1) << std::scientific << f[i][j] << ' ';
        std::cout << std::endl;
    }
    std::cout << "indices:" << std::endl;
    for (std::size_t i = 0; i < size[0]; ++i) {
        for (std::size_t j = 0; j < size[1]; ++j)
            std::cout << std::setw(7) << indices[i][j] << ' ';
        std::cout << std::endl;
    }

    return 0;
}

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

#ifndef distance_transform_h
#define distance_transform_h

#include <cmath>
#include <thread>
#include <vector>
#include "multiple_array.hpp"

namespace dt {

/// The DistanceTransform class provides static method to compute a distance field over any multi-dimensional regularly sampled function.
/// The dimension is fixed at compile time.
/// It is also possible to compute the index of the nearest minimum of each sample.
class DistanceTransform {
public:
    /**
     *    @brief Compute the L2-norm distance field D of a DIM-dimensional sampled function f. D gets the distance from the local minima of f.
     *    @param f              A DIM-dimensional, regularly sampled function.
     *    @param D              The resulting distance field of f.
     *    @param squared        Compute squared distances (L2)^2 - avoiding to compute square roots - (true) or keep them normal (false - default).
     *    @param nThreads       The number of threads for parallel computation. If <= 1, the computation will be sequential.
     *    @note Arrays f and D can also be the same.
     */
    template < typename Scalar, std::size_t DIM >
    inline static void distanceTransformL2(const MArray<Scalar, DIM> &f, MArray<Scalar, DIM> &D, const bool squared = false, const std::size_t nThreads = std::thread::hardware_concurrency())
    {
        std::size_t fSize[DIM], DSize[DIM];
        f.size(fSize);
        D.size(DSize);
        for (std::size_t d = 0; d < DIM; ++d)
            if (DSize[d] != fSize[d])
                throw std::out_of_range("Matrixes do not have same size.");

        MMArray<Scalar, DIM> fCopy(fSize);
        fCopy.import(f);
        MMArray<Scalar, DIM> DCopy(DSize);

        MArray<Scalar, DIM> tmpF(fCopy);
        MArray<Scalar, DIM> tmpD(DCopy);

        std::size_t order[DIM];

        for (std::size_t d = 0; d < DIM; ++d) {
            // permute rotate
            for (std::size_t o = 0; o < DIM; ++o)
                order[o] = (d + o) % DIM;
            MArray<Scalar, DIM> tmpF_rotated = tmpF.permute(order);
            MArray<Scalar, DIM> tmpD_rotated = tmpD.permute(order);


            std::size_t winStart[DIM] = {0}, winSize[DIM];
            tmpF_rotated.size(winSize);

            std::size_t range = winSize[0];
            if (nThreads < range) {
                range += range % nThreads;
                range /= nThreads;
            }
            std::size_t nWindows = winSize[0] / range + (winSize[0] % range != 0 ? 1 : 0);

            if (nWindows > 1) {
                std::vector<MArray<Scalar, DIM>> tmpWindowsF(nWindows);
                std::vector<MArray<Scalar, DIM>> tmpWindowsD(nWindows);
                std::vector<std::thread> threads(nWindows);

                for (std::size_t i = 0; i < nWindows; ++i) {
                    winStart[0] = i * range;
                    winSize[0] = std::min(range, tmpF_rotated.sizeAt(0) - winStart[0]);
                    tmpWindowsF.at(i) = tmpF_rotated.window(winStart, winSize);
                    tmpWindowsD.at(i) = tmpD_rotated.window(winStart, winSize);
                    winStart[0] = 0;
                    winSize[0] = tmpF_rotated.sizeAt(0);
                    threads.at(i) = std::thread(static_cast<void(*)(const MArray<Scalar, DIM> &, MArray<Scalar, DIM> &)>(&distanceL2Helper), std::cref(tmpWindowsF.at(i)), std::ref(tmpWindowsD.at(i)));
                }
                for (std::size_t i = 0; i < nWindows; ++i)
                    threads.at(i).join();
            } else {
                distanceL2Helper(tmpF_rotated, tmpD_rotated);
            }

            std::swap(tmpD, tmpF);
        }

        if (DIM % 2 == 0)
            DCopy = std::move(fCopy);

        D.import(DCopy);

        if (!squared)
            element_wiseSquareRoot(D);
    }

    /**
     *    @brief Compute the L2-norm distance field D of a 1-dimensional sampled function f. D gets the distance from the local minima of f.
     *    @param f              A 1-dimensional, regularly sampled function.
     *    @param D              The resulting distance field of f.
     *    @param squared        Compute squared distances (L2)^2 - avoiding to compute square roots - (true) or keep them normal (false - default).
     *    @param nThreads       The number of threads for parallel computation. Actually NOT used, since it's not easy to run a single row computation in parallel.
     *    @note Arrays f and D can also be the same.
     */
    template < typename Scalar >
    inline static void distanceTransformL2(const MArray<Scalar, 1> &f, MArray<Scalar, 1> &D, const bool squared = false, const std::size_t nThreads = std::thread::hardware_concurrency())
    {
        (void)nThreads; // unused
        std::size_t fSize[1], DSize[1];
        f.size(fSize);
        D.size(DSize);
        if (DSize[0] != fSize[0])
            throw std::out_of_range("Matrixes do not have same size.");

        distanceL2(f, D);

        if (!squared)
            element_wiseSquareRoot(D);
    }

    /**
     *    @brief Compute the L2-norm distance field D of a DIM-dimensional sampled function f. D gets the distance from the local minima of f.
     *           Compute also the (index of the) nearest local minimum of each sample.
     *    @param f              A DIM-dimensional, regularly sampled function.
     *    @param D              The resulting distance field of f.
     *    @param I              Resulting array containing the (index of the) local minimum for each sample.
     *    @param squared        Compute squared distances (L2)^2 - avoiding to compute square roots - (true) or keep them normal (false - default).
     *    @param nThreads       The number of threads for parallel computation. If <= 1, the computation will be sequential.
     *    @note Arrays f and D can also be the same.
     */
    template < typename Scalar, std::size_t DIM >
    inline static void distanceTransformL2(const MArray<Scalar, DIM> &f, MArray<Scalar, DIM> &D, MArray<std::size_t, DIM> &I, const bool squared = false, const std::size_t nThreads = std::thread::hardware_concurrency())
    {
        std::size_t fSize[DIM], DSize[DIM], ISize[DIM];
        f.size(fSize);
        D.size(DSize);
        I.size(ISize);
        for (std::size_t d = 0; d < DIM; ++d)
            if (DSize[d] != fSize[d] || ISize[d] != fSize[d])
                throw std::out_of_range("Matrixes do not have same size.");

        // initialize I
        initializeIndices(I);

        MMArray<Scalar, DIM> fCopy(fSize);
        fCopy.import(f);
        MMArray<Scalar, DIM> DCopy(DSize);
        MMArray<std::size_t, DIM> ICopyPre(ISize), ICopyPost(ISize);
        ICopyPre.import(I);

        MArray<Scalar, DIM> tmpF(fCopy);
        MArray<Scalar, DIM> tmpD(D);
        MArray<std::size_t, DIM> Ipre(ICopyPre);
        MArray<std::size_t, DIM> Ipost(ICopyPost);

        std::size_t order[DIM];

        for (std::size_t d = 0; d < DIM; ++d) {
            // permute rotate
            for (std::size_t o = 0; o < DIM; ++o)
                order[o] = (d + o) % DIM;
            MArray<Scalar, DIM> tmpF_rotated = tmpF.permute(order);
            MArray<Scalar, DIM> tmpD_rotated = tmpD.permute(order);
            MArray<std::size_t, DIM> Ipre_rotated = Ipre.permute(order);
            MArray<std::size_t, DIM> Ipost_rotated = Ipost.permute(order);

            std::size_t winStart[DIM] = {0}, winSize[DIM];
            tmpF_rotated.size(winSize);

            std::size_t range = winSize[0];
            if (nThreads < range) {
                range += range % nThreads;
                range /= nThreads;
            }
            std::size_t nWindows = winSize[0] / range + (winSize[0] % range != 0 ? 1 : 0);

            if (nWindows > 1) {
                std::vector<MArray<Scalar, DIM>> tmpWindowsF(nWindows);
                std::vector<MArray<Scalar, DIM>> tmpWindowsD(nWindows);
                std::vector<MArray<std::size_t, DIM>> tmpWindowsIPre(nWindows);
                std::vector<MArray<std::size_t, DIM>> tmpWindowsIPost(nWindows);
                std::vector<std::thread> threads(nWindows);

                for (std::size_t i = 0; i < nWindows; ++i) {
                    winStart[0] = i * range;
                    winSize[0] = std::min(range, tmpF_rotated.sizeAt(0) - winStart[0]);
                    tmpWindowsF.at(i) = tmpF_rotated.window(winStart, winSize);
                    tmpWindowsD.at(i) = tmpD_rotated.window(winStart, winSize);
                    tmpWindowsIPre.at(i) = Ipre_rotated.window(winStart, winSize);
                    tmpWindowsIPost.at(i) = Ipost_rotated.window(winStart, winSize);
                    winStart[0] = 0;
                    winSize[0] = tmpF_rotated.sizeAt(0);
                    threads.at(i) = std::thread(static_cast<void(*)(const MArray<Scalar, DIM> &, MArray<Scalar, DIM> &, const MArray<std::size_t, DIM> &, MArray<std::size_t, DIM> &)>(&distanceL2Helper), std::cref(tmpWindowsF.at(i)), std::ref(tmpWindowsD.at(i)), std::cref(tmpWindowsIPre.at(i)), std::ref(tmpWindowsIPost.at(i)));
                }
                for (std::size_t i = 0; i < nWindows; ++i)
                    threads.at(i).join();
            } else {
                distanceL2Helper(tmpF_rotated, tmpD_rotated, Ipre_rotated, Ipost_rotated);
            }

            std::swap(tmpD, tmpF);
            std::swap(Ipost, Ipre);
        }

        if (DIM % 2 == 0) {
            DCopy = std::move(fCopy);
            ICopyPost = std::move(ICopyPre);
        }

        D.import(DCopy);
        I.import(ICopyPost);

        if (!squared)
            element_wiseSquareRoot(D);
    }

    /**
     *    @brief Compute the L2-norm distance field D of a 1-dimensional sampled function f. D gets the distance from the local minima of f.
     *           Compute also the (index of the) nearest local minimum of each sample.
     *    @param f              A 1-dimensional, regularly sampled function.
     *    @param D              The resulting distance field of f.
     *    @param I              Resulting array containing the (index of the) local minimum for each sample.
     *    @param squared        Compute squared distances (L2)^2 - avoiding to compute square roots - (true) or keep them normal (false - default).
     *    @param nThreads       The number of threads for parallel computation. Actually NOT used, since it's not easy to run a single row computation in parallel.
     *    @note Arrays f and D can also be the same.
     */
    template < typename Scalar >
    inline static void distanceTransformL2(const MArray<Scalar, 1> &f, MArray<Scalar, 1> &D, MArray<std::size_t, 1> &I, const bool squared = false, const std::size_t nThreads = std::thread::hardware_concurrency())
    {
        (void)nThreads; // unused
        std::size_t fSize[1], DSize[1], ISize[1];
        f.size(fSize);
        D.size(DSize);
        I.size(ISize);
        if (DSize[0] != fSize[0] || ISize[0] != fSize[0])
            throw std::out_of_range("Matrixes do not have same size.");

        // initialize I
        initializeIndices(I);

        distanceL2(f, D, I);

        if (!squared)
            element_wiseSquareRoot(D);
    }

    /**
     *    @brief Set up the initial indices of a DIM-dimensional sampled function.
     *    @param I              Resulting array containing the initial index for each sample.
     */
    template < std::size_t DIM >
    inline static void initializeIndices(MArray<std::size_t, DIM> &I)
    {
        MArray<std::size_t, DIM-1> I_q;
        for (std::size_t q = 0; q < I.sizeAt(0); ++q) {
            I.at(q, I_q);
            initializeIndices(I_q);
        }
    }

    /**
     *    @brief Set up the initial indices of a 1-dimensional sampled function.
     *    @param I              Resulting array containing the initial index for each sample.
     */
    inline static void initializeIndices(MArray<std::size_t, 1> &I)
    {
        for (std::size_t q = 0; q < I.sizeAt(0); ++q)
            I[q] = I.accumulatedOffset(q);
    }



private:
    /**
     *    @brief The loop iteration process that can be executed sequentially and concurrently as well.
     *    @param f             A DIM-dimensional, regularly sampled function (a window, in multi-threading).
     *    @param D             The resulting distance field of f (a window, in multi-threading).
     *    @param d             The dimension where to slice.
     *    @param order         The order in which to permute the slices.
     */
    template < typename Scalar, std::size_t DIM >
    inline static void distanceL2Helper(const MArray<Scalar, DIM> &f, MArray<Scalar, DIM> &D)
    {
        MArray<Scalar, DIM-1> f_dq;
        MArray<Scalar, DIM-1> D_dq;

        for (std::size_t q = 0; q < f.sizeAt(0); ++q) {
            f.slice(0, q, f_dq);
            D.slice(0, q, D_dq);
            distanceL2(f_dq, D_dq);
        }
    }

    /**
     *    @brief The actual distance field computation is done by recursive calls of this method, in lower dimenional sub-functions.
     *    @param f              A DIM-dimensional, regularly sampled function.
     *    @param D              The resulting distance field of f.
     */
    template < typename Scalar, std::size_t DIM >
    inline static void distanceL2(const MArray<Scalar, DIM> &f, MArray<Scalar, DIM> &D)
    {
        MArray<Scalar, DIM-1> f_q, D_q;
        // compute distance at lower dimensions for each hyperplane
        for (std::size_t q = 0; q < f.sizeAt(0); ++q) {
            f.at(q, f_q);
            D.at(q, D_q);
            distanceL2(f_q, D_q);
        }
    }

    /**
     *    @brief The actual distance field computation as in the "Distance Transforms of Sampled Functions" paper, performed in a mono-dimensional function.
     *    @param f              A 1-dimensional, regularly sampled function.
     *    @param D              The resulting distance field of f.
     */
    template < typename Scalar >
    inline static void distanceL2(const MArray<Scalar, 1> &f, MArray<Scalar, 1> &D)
    {
        if (f.sizeAt(0) == 0 || f.sizeAt(0) > D.sizeAt(0))
            return;
        if (f.sizeAt(0) == 1) {
            D[0] = f[0];
            return;
        }
        std::size_t k = 0;                          // index of rightmost parabola in lower envelope
        std::size_t *v = new std::size_t[f.sizeAt(0)]; // locations of parabolas in lower envelope
        double *z = new double[f.sizeAt(0) + 1];       // locations of boundaries between parabolas
        double s = double(0);
        // initialization
        v[0] = 0;
        z[0] = -std::numeric_limits<double>::max();
        z[1] = std::numeric_limits<double>::max();
        // compute lowest envelope:
        for (std::size_t q = 1; q < f.sizeAt(0); ++q) {
            ++k;    // this compensates for first line of next do-while block
            do {
                --k;
                // compute horizontal position of intersection between the parabola from q and the current lowest parabola
                s = ((f[q] + q*q) - static_cast<double>(f[v[k]] + v[k]*v[k])) / (2*q - static_cast<double>(2*v[k]));
            } while (s <= z[k]);
            ++k;
            v[k] = q;
            z[k] = s;
            z[k+1] = std::numeric_limits<double>::max();
        }
        // fill in values of distance transform
        k = 0;
        for (std::size_t q = 0; q < f.sizeAt(0); ++q) {
            while(z[k+1] < static_cast<double>(q))
                ++k;
            D[q] = f[v[k]] + (q - static_cast<Scalar>(v[k]))*(q - static_cast<Scalar>(v[k]));
        }
        // delete allocated memory
        delete[] z;
        delete[] v;
    }

    /**
     *    @brief The loop iteration process that can be executed sequentially and concurrently as well.
     *    @param f             A DIM-dimensional, regularly sampled function (a window, in multi-threading).
     *    @param D             The resulting distance field of f (a window, in multi-threading).
     *    @param Ipre          Array containing the initial inedx of local minimum for each sample.
     *    @param Ipost         Array containing the resulting index of local minimum for each sample.
     */
    template < typename Scalar, std::size_t DIM >
    inline static void distanceL2Helper(const MArray<Scalar, DIM> &f, MArray<Scalar, DIM> &D, const MArray<std::size_t, DIM> &Ipre, MArray<std::size_t, DIM> &Ipost)
    {
        MArray<Scalar, DIM-1> f_dq;
        MArray<Scalar, DIM-1> D_dq;
        MArray<std::size_t, DIM-1> Ipre_dq;
        MArray<std::size_t, DIM-1> Ipost_dq;

        for (std::size_t q = 0; q < f.sizeAt(0); ++q) {
            f.slice(0, q, f_dq);
            D.slice(0, q, D_dq);
            Ipre.slice(0, q, Ipre_dq);
            Ipost.slice(0, q, Ipost_dq);
            distanceL2(f_dq, D_dq, Ipre_dq, Ipost_dq);
        }
    }

    /**
     *    @brief The actual distance field computation is done by recursive calls of this method, in lower dimenional sub-functions.
     *    @param f              A DIM-dimensional, regularly sampled function.
     *    @param D              The resulting distance field of f.
     *    @param I              Resulting array containing the index of the local minimum for each sample.
     */
    template < typename Scalar, std::size_t DIM >
    inline static void distanceL2(const MArray<Scalar, DIM> &f, MArray<Scalar, DIM> &D, const MArray<std::size_t, DIM> &Ipre, MArray<std::size_t, DIM> &Ipost)
    {
        MArray<Scalar, DIM-1> f_q, D_q;
        MArray<std::size_t, DIM-1> Ipre_q, Ipost_q;
        // compute distance at lower dimensions for each hyperplane
        for (std::size_t q = 0; q < f.sizeAt(0); ++q) {
            f.at(q, f_q);
            D.at(q, D_q);
            Ipre.at(q, Ipre_q);
            Ipost.at(q, Ipost_q);
            distanceL2(f_q, D_q, Ipre_q, Ipost_q);
        }
    }

    /**
     *    @brief The actual distance field computation as in the "Distance Transforms of Sampled Functions" paper, performed in a mono-dimensional function.
     *    @param f              A 1-dimensional, regularly sampled function.
     *    @param D              The resulting distance field of f.
     *    @param I              Resulting array containing the (index of the) local minimum for each sample.
     */
    template < typename Scalar >
    inline static void distanceL2(const MArray<Scalar, 1> &f, MArray<Scalar, 1> &D, const MArray<std::size_t, 1> &Ipre, MArray<std::size_t, 1> &Ipost)
    {
        if (f.sizeAt(0) == 0 || f.sizeAt(0) > D.sizeAt(0))
            return;
        if (f.sizeAt(0) == 1) {
            D[0] = f[0];
            Ipost[0] = Ipre[0];
            return;
        }
        std::size_t k = 0;                          // index of rightmost parabola in lower envelope
        std::size_t *v = new std::size_t[f.sizeAt(0)]; // locations of parabolas in lower envelope
        double *z = new double[f.sizeAt(0) + 1];       // locations of boundaries between parabolas
        double s = double(0);
        // initialization
        v[0] = 0;
        z[0] = -std::numeric_limits<double>::max();
        z[1] = std::numeric_limits<double>::max();
        // compute lowest envelope:
        for (std::size_t q = 1; q < f.sizeAt(0); ++q) {
            ++k;    // this compensates for first line of next do-while block
            do {
                --k;
                // compute horizontal position of intersection between the parabola from q and the current lowest parabola
                s = ((f[q] + q*q) - static_cast<double>(f[v[k]] + v[k]*v[k])) / (2*q - static_cast<double>(2*v[k]));
            } while (s <= z[k]);
            ++k;
            v[k] = q;
            z[k] = s;
            z[k+1] = std::numeric_limits<double>::max();
        }
        // fill in values of distance transform
        k = 0;
        for (std::size_t q = 0; q < f.sizeAt(0); ++q) {
            while(z[k+1] < static_cast<double>(q))
                ++k;
            D[q] = f[v[k]] + (q - static_cast<Scalar>(v[k]))*(q - static_cast<Scalar>(v[k]));
            Ipost[q] = Ipre[v[k]];
        }
        // delete allocated memory
        delete[] z;
        delete[] v;
    }

public:
    /**
     *    @brief Compute the square root of each individual element of a DIM-dimensional array.
     *    @param m              A DIM-dimensioanl array whose element have to be square rooted.
     */
    template < typename Scalar, std::size_t DIM >
    inline static void element_wiseSquareRoot(MArray<Scalar, DIM> &m)
    {
        MArray<Scalar, DIM-1> mm;
        for (std::size_t q = 0; q < m.sizeAt(0); ++q) {
            m.at(q, mm);
            element_wiseSquareRoot(mm);
        }
    }

    /**
     *    @brief Compute the square root of each individual element of a 1-dimensional array.
     *    @param m              A 1-dimensioanl array whose element have to be square rooted.
     */
    template < typename Scalar >
    inline static void element_wiseSquareRoot(MArray<Scalar, 1> &m)
    {
        for (std::size_t q = 0; q < m.sizeAt(0); ++q)
            m[q] = static_cast<Scalar>(std::sqrt(m[q]));
    }
};

}

#endif /* distance_transform_h */

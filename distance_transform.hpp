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

#ifndef distance_transform_h
#define distance_transform_h

#include <cmath>
#include "multiple_array.hpp"

class DistanceTransform{
public:
    template < typename Scalar = float, std::size_t DIM = 2 >
    inline static void distanceTransformL2(const MMArray<Scalar, DIM> &f, MMArray<Scalar, DIM> &D, const bool squared = false)
    {
        MMArray<Scalar, DIM> fCopy(f);
        // compute for each slice
        for (std::size_t d = 0; d < DIM; ++d)
            for (std::size_t q = 0; q < fCopy.size(d); ++q)
                distanceL2(fCopy.slice(d, q), fCopy.slice(d, q));
        D = std::move(fCopy);
        if (!squared)
            element_wiseSquareRoot(D);
    }
    
    template < typename Scalar = float >
    inline static void distanceTransformL2(const MArray<Scalar, 1> &f, MArray<Scalar, 1> &D, const bool squared = false)
    {
        distanceL2(f, D);
        if (!squared) {
            element_wiseSquareRoot(D);
        }
    }

    template < typename Scalar = float, std::size_t DIM = 2 >
    inline static void distanceTransformL2(const MMArray<Scalar, DIM> &f, MMArray<Scalar, DIM> &D, MArray<std::size_t, DIM> &I, const bool squared = false)
    {
        MMArray<Scalar, DIM> fCopy(f);
        MArray<Scalar, DIM-1> fCopy_dq;
        MArray<std::size_t, DIM-1> I_dq;
        // compute for each slice
        for (std::size_t d = 0; d < DIM; ++d)
            for (std::size_t q = 0; q < fCopy.size(d); ++q) {
                fCopy_dq = fCopy.slice(d, q);
                I_dq = I.slice(d, q);
                distanceL2(fCopy_dq, fCopy_dq, I_dq);
            }
        D = std::move(fCopy);
        if (!squared)
            element_wiseSquareRoot(D);
    }
    
    template < typename Scalar = float >
    inline static void distanceTransformL2(const MArray<Scalar, 1> &f, MArray<Scalar, 1> &D, MArray<std::size_t, 1> &I, const bool squared = false)
    {
        distanceL2(f, D, I);
        if (!squared) {
            element_wiseSquareRoot(D);
        }
    }

private:
    template < typename Scalar = float, std::size_t DIM >
    inline static void distanceL2(const MArray<Scalar, DIM> &f, MArray<Scalar, DIM> &D)
    {
        // compute distance at lower dimensions for each hyperplane
        for (std::size_t q = 0; q < f.size(); ++q)
            distanceL2(f[q], D[q]);
    }

    template < typename Scalar = float >
    inline static void distanceL2(const MArray<Scalar, 1> &f, MArray<Scalar, 1> &D)
    {
        if (f.size() == 0 || f.size() > D.size())
            return;
        if (f.size() == 1) {
            D[0] = f[0];
            return;
        }
        std::size_t k = 0;                          // index of rightmost parabola in lower envelope
        std::size_t *v = new std::size_t[f.size()]; // locations of parabolas in lower envelope
        double *z = new double[f.size() + 1];       // locations of boundaries between parabolas
        double s = double(0);
        // initialization
        v[0] = 0;
        z[0] = std::numeric_limits<double>::lowest();
        z[1] = std::numeric_limits<double>::max();
        // compute lowest envelope:
        for (std::size_t q = 1; q < f.size(); ++q) {
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
        for (std::size_t q = 0; q < f.size(); ++q) {
            while(z[k+1] < static_cast<double>(q))
                ++k;
            D[q] = f[v[k]] + (q - static_cast<Scalar>(v[k]))*(q - static_cast<Scalar>(v[k]));
        }
        // delete allocated memory
        delete[] z;
        delete[] v;
    }

    template < typename Scalar = float, std::size_t DIM >
    inline static void distanceL2(const MArray<Scalar, DIM> &f, MArray<Scalar, DIM> &D, MArray<std::size_t, DIM> &I)
    {
        // compute distance at lower dimensions for each hyperplane
        for (std::size_t q = 0; q < f.size(); ++q)
            distanceL2(f[q], D[q], I[q]);
    }

    template < typename Scalar = float >
    inline static void distanceL2(const MArray<Scalar, 1> &f, MArray<Scalar, 1> &D, MArray<std::size_t, 1> &I)
    {
        if (f.size() == 0 || f.size() > D.size())
            return;
        if (f.size() == 1) {
            D[0] = f[0];
            return;
        }
        std::size_t k = 0;                          // index of rightmost parabola in lower envelope
        std::size_t *v = new std::size_t[f.size()]; // locations of parabolas in lower envelope
        double *z = new double[f.size() + 1];       // locations of boundaries between parabolas
        double s = double(0);
        // initialization
        v[0] = 0;
        z[0] = std::numeric_limits<double>::lowest();
        z[1] = std::numeric_limits<double>::max();
        // compute lowest envelope:
        for (std::size_t q = 1; q < f.size(); ++q) {
            ++k;    // this compensates for first line of next do-while block
            do {
                --k;
                // compute horizontal position of intersection between the parabola from q and the current lowest parabola
                double a = f[q] + q*q;
                double b = f[v[k]] + v[k]*v[k];
                double den = 2 * (q - static_cast<double>(v[k]));
                s = (a - b) / den;
                s = ((f[q] + q*q) - static_cast<double>(f[v[k]] + v[k]*v[k])) / (2*q - static_cast<double>(2*v[k]));
            } while (s <= z[k]);
            ++k;
            v[k] = q;
            z[k] = s;
            z[k+1] = std::numeric_limits<double>::max();
        }
        // fill in values of distance transform
        k = 0;
        for (std::size_t q = 0; q < f.size(); ++q) {
            while(z[k+1] < static_cast<double>(q))
                ++k;
            D[q] = f[v[k]] + (q - static_cast<Scalar>(v[k]))*(q - static_cast<Scalar>(v[k]));
            I[q] = f.accumulatedOffset(v[k]);
        }
        // delete allocated memory
        delete[] z;
        delete[] v;
    }

public:
    template < typename Scalar = float, std::size_t DIM >
    inline static void element_wiseSquareRoot(MArray<Scalar, DIM> &m)
    {
        for (std::size_t q = 0; q < m.size(); ++q) {
            MArray<Scalar, DIM-1> mm = m[q];
            element_wiseSquareRoot(mm);
        }
    }

    template < typename Scalar = float >
    inline static void element_wiseSquareRoot(MArray<Scalar, 1> &m)
    {
        for (std::size_t q = 0; q < m.size(); ++q)
            m[q] = static_cast<Scalar>(std::sqrt(m[q]));
    }
};

#endif /* distance_transform_h */

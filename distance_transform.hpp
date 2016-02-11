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
#include "multiple_array.hpp"

class DistanceTransform{
public:
    template < typename Scalar = float, std::size_t DIM = 2 >
    inline static void distanceTransformL2(const MMArray<Scalar, DIM> &f, MMArray<Scalar, DIM> &D, const bool squared = false)
    {
        if (&D != &f) {
            std::size_t size[DIM];
            f.size(size);
            D.resize(size);
        }
        MMArray<Scalar, DIM> fCopy(f);
        MArray<Scalar, DIM> tmpF(fCopy);
        MArray<Scalar, DIM> tmpD(D);
        MArray<Scalar, DIM-1> f_dq;
        MArray<Scalar, DIM-1> D_dq;
        // compute for each slice
        for (std::size_t d = 0; d < DIM; ++d) {
            for (std::size_t q = 0; q < tmpF.size(d); ++q) {
                tmpF.slice(d, q, f_dq);
                tmpD.slice(d, q, D_dq);
                distanceL2(f_dq, D_dq);
            }
            std::swap(tmpD, tmpF);
        }
        if (DIM % 2 == 0)
            D = std::move(fCopy);
        if (!squared)
            element_wiseSquareRoot(D);
    }
    
    template < typename Scalar = float >
    inline static void distanceTransformL2(const MMArray<Scalar, 1> &f, MMArray<Scalar, 1> &D, const bool squared = false)
    {
        std::size_t size[1];
        f.size(size);
        if (&D != &f)
            D.resize(size);
        distanceL2(f, D);
        if (!squared) {
            element_wiseSquareRoot(D);
        }
    }

    template < typename Scalar = float, std::size_t DIM = 2 >
    inline static void distanceTransformL2(const MMArray<Scalar, DIM> &f, MMArray<Scalar, DIM> &D, MMArray<std::size_t, DIM> &I, const bool squared = false)
    {
        std::size_t size[DIM];
        f.size(size);
        if (&D != &f)
            D.resize(size);
        I.resize(size);
        MMArray<Scalar, DIM> fCopy(f);          // make a safe copy of f
        MArray<Scalar, DIM> tmpF(fCopy);
        MArray<Scalar, DIM> tmpD(D);
        MMArray<std::size_t, DIM> ICopy(I);     // make a safe copy of I
        MArray<std::size_t, DIM> Ipre(ICopy);
        MArray<std::size_t, DIM> Ipost(I);
        MArray<Scalar, DIM-1> f_dq;
        MArray<Scalar, DIM-1> D_dq;
        MArray<std::size_t, DIM-1> Ipre_dq;
        MArray<std::size_t, DIM-1> Ipost_dq;
        // initialize I
        initializeIndices(Ipre);
        // compute for each slice
        for (std::size_t d = 0; d < DIM; ++d) {
            for (std::size_t q = 0; q < tmpF.size(d); ++q) {
                tmpF.slice(d, q, f_dq);
                tmpD.slice(d, q, D_dq);
                Ipre.slice(d, q, Ipre_dq);
                Ipost.slice(d, q, Ipost_dq);
                distanceL2(f_dq, D_dq, Ipre_dq, Ipost_dq);
            }
            std::swap(tmpD, tmpF);
            std::swap(Ipost, Ipre);
        }
        if (DIM % 2 == 0) {
            D = std::move(fCopy);
            I = std::move(ICopy);
        }
        if (!squared)
            element_wiseSquareRoot(D);
    }
    
    template < typename Scalar = float >
    inline static void distanceTransformL2(const MMArray<Scalar, 1> &f, MMArray<Scalar, 1> &D, MMArray<std::size_t, 1> &I, const bool squared = false)
    {
        std::size_t size[1];
        f.size(size);
        if (&D != &f)
            D.resize(size);
        I.resize(size);
        distanceL2(f, D, I);
        if (!squared) {
            element_wiseSquareRoot(D);
        }
    }

    template < std::size_t DIM >
    inline static void initializeIndices(MArray<std::size_t, DIM> &I)
    {
        MArray<std::size_t, DIM-1> I_q;
        for (std::size_t q = 0; q < I.size(); ++q) {
            I.at(q, I_q);
            initializeIndices(I_q);
        }
    }
    
    inline static void initializeIndices(MArray<std::size_t, 1> &I)
    {
        for (std::size_t q = 0; q < I.size(); ++q)
            I[q] = I.accumulatedOffset(q);
    }
    
private:
    template < typename Scalar = float, std::size_t DIM >
    inline static void distanceL2(const MArray<Scalar, DIM> &f, MArray<Scalar, DIM> &D)
    {
        MArray<Scalar, DIM-1> f_q, D_q;
        // compute distance at lower dimensions for each hyperplane
        for (std::size_t q = 0; q < f.size(); ++q) {
            f.at(q, f_q);
            D.at(q, D_q);
            distanceL2(f_q, D_q);
        }
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
        z[0] = -std::numeric_limits<double>::max();
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
    inline static void distanceL2(const MArray<Scalar, DIM> &f, MArray<Scalar, DIM> &D, const MArray<std::size_t, DIM> &Ipre, MArray<std::size_t, DIM> &Ipost)
    {
        MArray<Scalar, DIM-1> f_q, D_q;
        MArray<std::size_t, DIM-1> Ipre_q, Ipost_q;
        // compute distance at lower dimensions for each hyperplane
        for (std::size_t q = 0; q < f.size(); ++q) {
            f.at(q, f_q);
            D.at(q, D_q);
            Ipre.at(q, Ipre_q);
            Ipost.at(q, Ipost_q);
            distanceL2(f_q, D_q, Ipre_q, Ipost_q);
        }
    }

    template < typename Scalar = float >
    inline static void distanceL2(const MArray<Scalar, 1> &f, MArray<Scalar, 1> &D, const MArray<std::size_t, 1> &Ipre, MArray<std::size_t, 1> &Ipost)
    {
        if (f.size() == 0 || f.size() > D.size())
            return;
        if (f.size() == 1) {
            D[0] = f[0];
            Ipost[0] = Ipre[0];
            return;
        }
        std::size_t k = 0;                          // index of rightmost parabola in lower envelope
        std::size_t *v = new std::size_t[f.size()]; // locations of parabolas in lower envelope
        double *z = new double[f.size() + 1];       // locations of boundaries between parabolas
        double s = double(0);
        // initialization
        v[0] = 0;
        z[0] = -std::numeric_limits<double>::max();
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
            Ipost[q] = Ipre[v[k]];
        }
        // delete allocated memory
        delete[] z;
        delete[] v;
    }

public:
    template < typename Scalar = float, std::size_t DIM >
    inline static void element_wiseSquareRoot(MArray<Scalar, DIM> &m)
    {
        MArray<Scalar, DIM-1> mm;
        for (std::size_t q = 0; q < m.size(); ++q) {
            m.at(q, mm);
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

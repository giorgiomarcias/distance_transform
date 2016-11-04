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

#include <cmath>
#include <distance_transform/distance_transform.hpp>

namespace dt {

	using namespace dope;

	template < typename Scalar, SizeType DIM >
	inline void DistanceTransform::distanceTransformL2(const DopeVector<Scalar, DIM> &f, DopeVector<Scalar, DIM> &D, const bool squared, const std::size_t nThreads)
	{
		Index<DIM> fSize, DSize;
		fSize = f.allSizes();
		DSize = D.allSizes();
		if (DSize != fSize)
			throw std::out_of_range("Matrixes do not have same size.");

		Grid<Scalar, DIM> fCopy(fSize);
		fCopy.import(f);
		Grid<Scalar, DIM> DCopy(DSize);

		DopeVector<Scalar, DIM> tmpF(fCopy);
		DopeVector<Scalar, DIM> tmpD(DCopy);

		Index<DIM> order;

		for (SizeType d = static_cast<SizeType>(0); d < DIM; ++d) {
			// permute rotate
			for (SizeType o = static_cast<SizeType>(0); o < DIM; ++o)
				order[o] = (d + o) % DIM;
			DopeVector<Scalar, DIM> tmpF_rotated = tmpF.permute(order);
			DopeVector<Scalar, DIM> tmpD_rotated = tmpD.permute(order);


			Index<DIM> winStart = Index<DIM>::Zero(), winSize;
			tmpF_rotated.allSizes(winSize);

			SizeType range = winSize[0];
			if (nThreads < range) {
				range += range % nThreads;
				range /= nThreads;
			}
			std::size_t nWindows = winSize[0] / range + (winSize[0] % range != 0 ? 1 : 0);

			if (nWindows > 1) {
				std::vector<DopeVector<Scalar, DIM>> tmpWindowsF(nWindows);
				std::vector<DopeVector<Scalar, DIM>> tmpWindowsD(nWindows);
				std::vector<std::thread> threads(nWindows);

				for (std::size_t i = 0; i < nWindows; ++i) {
					winStart[0] = i * range;
					winSize[0] = std::min(range, tmpF_rotated.sizeAt(0) - winStart[0]);
					tmpWindowsF.at(i) = tmpF_rotated.window(winStart, winSize);
					tmpWindowsD.at(i) = tmpD_rotated.window(winStart, winSize);
					winStart[0] = 0;
					winSize[0] = tmpF_rotated.sizeAt(0);
					threads.at(i) = std::thread(static_cast<void(*)(const DopeVector<Scalar, DIM> &, DopeVector<Scalar, DIM> &)>(&distanceL2Helper), std::cref(tmpWindowsF.at(i)), std::ref(tmpWindowsD.at(i)));
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



	template < typename Scalar >
	inline void DistanceTransform::distanceTransformL2(const DopeVector<Scalar, 1> &f, DopeVector<Scalar, 1> &D, const bool squared, const std::size_t)
	{
		Index1 fSize, DSize;
		fSize = f.allSizes();
		DSize = D.allSizes();
		if (DSize != fSize)
			throw std::out_of_range("Matrixes do not have same size.");

		distanceL2(f, D);

		if (!squared)
			element_wiseSquareRoot(D);
	}



	template < typename Scalar, SizeType DIM >
	inline void DistanceTransform::distanceTransformL2(const DopeVector<Scalar, DIM> &f, DopeVector<Scalar, DIM> &D, DopeVector<SizeType, DIM> &I, const bool squared, const std::size_t nThreads)
	{
		Index<DIM> fSize, DSize, ISize;
		fSize = f.allSizes();
		DSize = D.allSizes();
		ISize = I.allSizes();
		if (DSize != fSize || ISize != fSize)
			throw std::out_of_range("Matrixes do not have same size.");

		// initialize I
		initializeIndices(I);

		Grid<Scalar, DIM> fCopy(fSize);
		fCopy.import(f);
		Grid<Scalar, DIM> DCopy(DSize);
		Grid<std::size_t, DIM> ICopyPre(ISize), ICopyPost(ISize);
		ICopyPre.import(I);

		DopeVector<Scalar, DIM> tmpF(fCopy);
		DopeVector<Scalar, DIM> tmpD(D);
		DopeVector<std::size_t, DIM> Ipre(ICopyPre);
		DopeVector<std::size_t, DIM> Ipost(ICopyPost);

		Index<DIM> order;

		for (SizeType d = static_cast<SizeType>(0); d < DIM; ++d) {
			// permute rotate
			for (SizeType o = static_cast<SizeType>(0); o < DIM; ++o)
				order[o] = (d + o) % DIM;
			DopeVector<Scalar, DIM> tmpF_rotated = tmpF.permute(order);
			DopeVector<Scalar, DIM> tmpD_rotated = tmpD.permute(order);
			DopeVector<std::size_t, DIM> Ipre_rotated = Ipre.permute(order);
			DopeVector<std::size_t, DIM> Ipost_rotated = Ipost.permute(order);

			Index<DIM> winStart = Index<DIM>::Zero(), winSize;
			tmpF_rotated.allSizes(winSize);

			std::size_t range = winSize[0];
			if (nThreads < range) {
				range += range % nThreads;
				range /= nThreads;
			}
			std::size_t nWindows = winSize[0] / range + (winSize[0] % range != 0 ? 1 : 0);

			if (nWindows > 1) {
				std::vector<DopeVector<Scalar, DIM>> tmpWindowsF(nWindows);
				std::vector<DopeVector<Scalar, DIM>> tmpWindowsD(nWindows);
				std::vector<DopeVector<SizeType, DIM>> tmpWindowsIPre(nWindows);
				std::vector<DopeVector<SizeType, DIM>> tmpWindowsIPost(nWindows);
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
					threads.at(i) = std::thread(static_cast<void(*)(const DopeVector<Scalar, DIM> &, DopeVector<Scalar, DIM> &, const DopeVector<SizeType, DIM> &, DopeVector<SizeType, DIM> &)>(&distanceL2Helper), std::cref(tmpWindowsF.at(i)), std::ref(tmpWindowsD.at(i)), std::cref(tmpWindowsIPre.at(i)), std::ref(tmpWindowsIPost.at(i)));
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




	template < typename Scalar >
	inline void DistanceTransform::distanceTransformL2(const DopeVector<Scalar, 1> &f, DopeVector<Scalar, 1> &D, DopeVector<SizeType, 1> &I, const bool squared, const std::size_t)
	{
		Index1 fSize, DSize, ISize;
		fSize = f.allSizes();
		DSize = D.allSizes();
		ISize = I.allSizes();
		if (DSize != fSize || ISize != fSize)
			throw std::out_of_range("Matrixes do not have same size.");

		// initialize I
		initializeIndices(I);

		distanceL2(f, D, I);

		if (!squared)
			element_wiseSquareRoot(D);
	}



	template < SizeType DIM >
	inline void DistanceTransform::initializeIndices(DopeVector<SizeType, DIM> &I)
	{
		DopeVector<SizeType, DIM-1> I_q;
		for (SizeType q = static_cast<SizeType>(0); q < I.sizeAt(0); ++q) {
			I.at(q, I_q);
			initializeIndices(I_q);
		}
	}



	inline void DistanceTransform::initializeIndices(DopeVector<SizeType, 1> &I)
	{
		for (SizeType q = static_cast<SizeType>(0); q < I.sizeAt(0); ++q)
			I[q] = I.accumulatedOffset(q);
	}



	template < typename Scalar, SizeType DIM >
	inline void DistanceTransform::distanceL2Helper(const DopeVector<Scalar, DIM> &f, DopeVector<Scalar, DIM> &D)
	{
		DopeVector<Scalar, DIM-1> f_dq;
		DopeVector<Scalar, DIM-1> D_dq;

		for (SizeType q = static_cast<SizeType>(0); q < f.sizeAt(0); ++q) {
			f.slice(0, q, f_dq);
			D.slice(0, q, D_dq);
			distanceL2(f_dq, D_dq);
		}
	}



	template < typename Scalar, SizeType DIM >
	inline void DistanceTransform::distanceL2(const DopeVector<Scalar, DIM> &f, DopeVector<Scalar, DIM> &D)
	{
		DopeVector<Scalar, DIM-1> f_q, D_q;
		// compute distance at lower dimensions for each hyperplane
		for (SizeType q = static_cast<SizeType>(0); q < f.sizeAt(0); ++q) {
			f.at(q, f_q);
			D.at(q, D_q);
			distanceL2(f_q, D_q);
		}
	}



	template < typename Scalar >
	inline void DistanceTransform::distanceL2(const DopeVector<Scalar, 1> &f, DopeVector<Scalar, 1> &D)
	{
		if (f.sizeAt(0) == static_cast<SizeType>(0) || f.sizeAt(0) > D.sizeAt(0))
			return;
		if (f.sizeAt(0) == static_cast<SizeType>(1)) {
			D[0] = f[0];
			return;
		}
		SizeType k = static_cast<SizeType>(0);          // index of rightmost parabola in lower envelope
		SizeType *v = new SizeType[f.sizeAt(0)];	    // locations of parabolas in lower envelope
		double *z = new double[f.sizeAt(0) + 1];        // locations of boundaries between parabolas
		double s = double(0);
		// initialization
		v[0] = static_cast<SizeType>(0);
		z[0] = -std::numeric_limits<double>::max();
		z[1] = std::numeric_limits<double>::max();
		// compute lowest envelope:
		for (SizeType q = static_cast<SizeType>(1); q < f.sizeAt(0); ++q) {
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
		k = static_cast<SizeType>(0);
		for (SizeType q = static_cast<SizeType>(0); q < f.sizeAt(0); ++q) {
			while(z[k+1] < static_cast<double>(q))
				++k;
			D[q] = f[v[k]] + (q - static_cast<Scalar>(v[k]))*(q - static_cast<Scalar>(v[k]));
		}
		// delete allocated memory
		delete[] z;
		delete[] v;
	}



	template < typename Scalar, SizeType DIM >
	inline void DistanceTransform::distanceL2Helper(const DopeVector<Scalar, DIM> &f, DopeVector<Scalar, DIM> &D, const DopeVector<SizeType, DIM> &Ipre, DopeVector<SizeType, DIM> &Ipost)
	{
		DopeVector<Scalar, DIM-1> f_dq;
		DopeVector<Scalar, DIM-1> D_dq;
		DopeVector<SizeType, DIM-1> Ipre_dq;
		DopeVector<SizeType, DIM-1> Ipost_dq;

		for (SizeType q = static_cast<SizeType>(0); q < f.sizeAt(0); ++q) {
			f.slice(0, q, f_dq);
			D.slice(0, q, D_dq);
			Ipre.slice(0, q, Ipre_dq);
			Ipost.slice(0, q, Ipost_dq);
			distanceL2(f_dq, D_dq, Ipre_dq, Ipost_dq);
		}
	}



	template < typename Scalar, SizeType DIM >
	inline void DistanceTransform::distanceL2(const DopeVector<Scalar, DIM> &f, DopeVector<Scalar, DIM> &D, const DopeVector<SizeType, DIM> &Ipre, DopeVector<SizeType, DIM> &Ipost)
	{
		DopeVector<Scalar, DIM-1> f_q, D_q;
		DopeVector<SizeType, DIM-1> Ipre_q, Ipost_q;
		// compute distance at lower dimensions for each hyperplane
		for (SizeType q = static_cast<SizeType>(0); q < f.sizeAt(0); ++q) {
			f.at(q, f_q);
			D.at(q, D_q);
			Ipre.at(q, Ipre_q);
			Ipost.at(q, Ipost_q);
			distanceL2(f_q, D_q, Ipre_q, Ipost_q);
		}
	}



	template < typename Scalar >
	inline void DistanceTransform::distanceL2(const DopeVector<Scalar, 1> &f, DopeVector<Scalar, 1> &D, const DopeVector<SizeType, 1> &Ipre, DopeVector<SizeType, 1> &Ipost)
	{
		if (f.sizeAt(0) == static_cast<SizeType>(0) || f.sizeAt(0) > D.sizeAt(0))
			return;
		if (f.sizeAt(0) == static_cast<SizeType>(1)) {
			D[0] = f[0];
			Ipost[0] = Ipre[0];
			return;
		}
		SizeType k = static_cast<SizeType>(0);	        // index of rightmost parabola in lower envelope
		SizeType *v = new SizeType[f.sizeAt(0)];        // locations of parabolas in lower envelope
		double *z = new double[f.sizeAt(0) + 1];        // locations of boundaries between parabolas
		double s = double(0);
		// initialization
		v[0] = static_cast<SizeType>(0);
		z[0] = -std::numeric_limits<double>::max();
		z[1] = std::numeric_limits<double>::max();
		// compute lowest envelope:
		for (SizeType q = static_cast<SizeType>(1); q < f.sizeAt(0); ++q) {
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
		k = static_cast<SizeType>(0);
		for (SizeType q = static_cast<SizeType>(0); q < f.sizeAt(0); ++q) {
			while(z[k+1] < static_cast<double>(q))
				++k;
			D[q] = f[v[k]] + (q - static_cast<Scalar>(v[k]))*(q - static_cast<Scalar>(v[k]));
			Ipost[q] = Ipre[v[k]];
		}
		// delete allocated memory
		delete[] z;
		delete[] v;
	}



	template < typename Scalar, SizeType DIM >
	inline void DistanceTransform::element_wiseSquareRoot(DopeVector<Scalar, DIM> &m)
	{
		DopeVector<Scalar, DIM-1> mm;
		for (SizeType q = static_cast<SizeType>(0); q < m.sizeAt(0); ++q) {
			m.at(q, mm);
			element_wiseSquareRoot(mm);
		}
	}



	template < typename Scalar >
	inline void DistanceTransform::element_wiseSquareRoot(DopeVector<Scalar, 1> &m)
	{
		for (SizeType q = static_cast<SizeType>(0); q < m.sizeAt(0); ++q)
			m[q] = static_cast<Scalar>(std::sqrt(m[q]));
	}

}

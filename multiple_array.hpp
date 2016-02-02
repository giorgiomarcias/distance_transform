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

#ifndef multiple_array_hpp
#define multiple_array_hpp

#include <cstddef>
#include <cstring>
#include <memory>
#include <type_traits>
#include <exception>
#include <sstream>
#include <utility>

template < typename T, std::size_t D >
class MArray {
public:
    MArray() : _array(nullptr), _accumulatedOffset(0)
    {
        for (std::size_t i = 0; i < D; ++i) {
            _size[i] = 0;
            _offset[i] = 0;
        }
    }

    MArray(const T *array, const std::size_t accumulatedOffset, const std::size_t size[D])
        : _array(const_cast<T*>(array))
        , _accumulatedOffset(_array ? accumulatedOffset : 0)
    {
        for (std::size_t i = 0; i < D; ++i)
            _size[i] = _array ? size[i] : 0;
        _offset[D-1] = 1;
        for (std::size_t j = D-1; j > 0; --j)
            _offset[j-1] = _size[j-1] * _offset[j];
    }

    MArray(const T *array, const std::size_t accumulatedOffset, const std::size_t size[D], const std::size_t offset[D])
        : _array(const_cast<T*>(array))
        , _accumulatedOffset(_array ? accumulatedOffset : 0)
    {
        for (std::size_t i = 0; i < D; ++i) {
            _size[i] = _array ? size[i] : 0;
            _offset[i] = _array ? offset[i] : 0;
        }
    }

    MArray(const MArray &ma)
        : _array(ma._array)
        , _accumulatedOffset(ma._accumulatedOffset)
    {
        for (std::size_t i = 0; i < D; ++i)
            _size[i] = ma._size[i];
        _offset[D-1] = 1;
        for (std::size_t j = D-1; j > 0; --j)
            _offset[j-1] = _size[j-1] * _offset[j];
    }

    MArray(MArray &&ma)
        : _array(ma._array)
        , _accumulatedOffset(ma._accumulatedOffset)
    {
        ma._array = nullptr;
        for (std::size_t i = 0; i < D; ++i) {
            _size[i] = ma._size[i];
            ma._size[i] = 0;
            _offset[i] = ma._offset[i];
            ma._offset[i] = 0;
        }
    }

    inline MArray & operator=(const MArray &ma)
    {
        if (&ma != this) {
            _array = ma._array;
            _accumulatedOffset = ma._accumulatedOffset;
            for (std::size_t i = 0; i < D; ++i) {
                _size[i] = ma._size[i];
                _offset[i] = ma._offset[i];
            }
        }
        return *this;
    }

    inline MArray & operator=(MArray &&ma)
    {
        if (&ma != this) {
            _array = ma._array;
            ma._array = nullptr;
            _accumulatedOffset = ma._accumulatedOffset;
            ma._accumulatedOffset = 0;
            for (std::size_t i = 0; i < D; ++i) {
                _size[i] = ma._size[i];
                ma._size[i] = 0;
                _offset[i] = ma._offset[i];
                ma._offset[i] = 0;
            }
        }
        return *this;
    }

    inline MArray<T, D-1> operator[](const std::size_t i) const
    {
        if (i >= _size[0]) {
            std::stringstream stream;
            stream << "Index " << i << " is out of range [0, " << _size[0]-1 << ']';
            throw std::out_of_range(stream.str());
        }
        return MArray<T, D-1>(_array + _offset[0] * i, accumulatedOffset(i), _size + 1, _offset + 1);
    }

    inline MArray<T, D-1> slice(const std::size_t d, const std::size_t i) const
    {
        if (d >= D) {
            std::stringstream stream;
            stream << "Index " << d << " is out of range [0, " << D-1 << ']';
            throw std::out_of_range(stream.str());
        }
        if (i >= _size[d]) {
            std::stringstream stream;
            stream << "Index " << i << " is out of range [0, " << _size[d]-1 << ']';
            throw std::out_of_range(stream.str());
        }
        std::size_t size[D-1];
        std::size_t offset[D-1];
        std::size_t k = 0;
        for (std::size_t j = 0; j < d; ++j, ++k) {
            size[k] = _size[j];
            offset[k] = _offset[j];
        }
        for (std::size_t j = d+1; j < D; ++j, ++k) {
            size[k] = _size[j];
            offset[k] = _offset[j];
        }
        return MArray<T, D-1>(_array + _offset[d] * i, accumulatedOffset(i, d), size, offset);
    }

    inline std::size_t size(const std::size_t d = 0) const
    {
        return _size[d];
    }

    inline void size(std::size_t s[D]) const
    {
        for (std::size_t i = 0; i < D; ++i)
            s[i] = _size[i];
    }

    inline std::size_t totalSize() const
    {
        std::size_t total = _size[0];
        for (std::size_t i = 1; i < D; ++i)
            total *= _size[i];
        return total;
    }

    inline std::size_t accumulatedOffset(const std::size_t i, const std::size_t d = 0) const
    {
        if (d >= D) {
            std::stringstream stream;
            stream << "Index " << d << " is out of range [0, " << D-1 << ']';
            throw std::out_of_range(stream.str());
        }
        if (i >= _size[d]) {
            std::stringstream stream;
            stream << "Index " << i << " is out of range [0, " << _size[d]-1 << ']';
            throw std::out_of_range(stream.str());
        }
        return  _accumulatedOffset + _offset[d] * i;
    }

private:
    T                  *_array;
    std::size_t         _accumulatedOffset;

protected:
    std::size_t         _size[D];

private:
    std::size_t         _offset[D];
};

template < typename T >
class MArray<T, 1> {
public:
    MArray()
        : _array(nullptr)
        , _accumulatedOffset(0)
        , _size(0)
        , _offset(0)
    { }

    MArray(const T *array, const std::size_t accumulatedOffset, const std::size_t size)
        : _array(array)
        , _accumulatedOffset(_array ? accumulatedOffset : 0)
        , _size(_array ? size : 0)
        , _offset(0)
    { }

    MArray(const T *array, const std::size_t accumulatedOffset, const std::size_t size, const std::size_t offset)
        : _array(array)
        , _accumulatedOffset(_array ? accumulatedOffset : 0)
        , _size(_array ? size : 0)
        , _offset(_array ? offset : 0)
    { }

    MArray(const T *array, const std::size_t accumulatedOffset, const std::size_t size[1])
        : _array(const_cast<T*>(array))
        , _accumulatedOffset(_array ? accumulatedOffset : 0)
        , _size(_array ? size[0] : 0)
        , _offset(0)
    { }

    MArray(const T *array, const std::size_t accumulatedOffset, const std::size_t size[1], const std::size_t offset[1])
        : _array(const_cast<T*>(array))
        , _accumulatedOffset(_array ? accumulatedOffset : 0)
        , _size(_array ? size[0] : 0)
        , _offset(_array ? offset[0] : 0)
    { }

    MArray(const MArray &ma)
        : _array(ma._array)
        , _accumulatedOffset(ma._accumulatedOffset)
        , _size(ma._size)
        , _offset(ma._offset)
    { }

    MArray(MArray &&ma)
        : _array(ma._array)
        , _accumulatedOffset(ma._accumulatedOffset)
        , _size(ma._size)
        , _offset(ma._offset)
    {
        ma._array = nullptr;
        ma._accumulatedOffset = 0;
        ma._size = 0;
        ma._offset = 0;
    }

    MArray & operator=(const MArray &ma)
    {
        if (&ma != this) {
            _array = ma._array;
            _size = ma._size;
            _offset = ma._offset;
        }
        return *this;
    }

    MArray & operator=(MArray &&ma)
    {
        if (&ma != this) {
            _array = ma._array;
            _size = ma._size;
            _offset = ma._offset;
            ma._array = nullptr;
            ma._size = 0;
            ma._offset = 0;
        }
        return *this;
    }

    inline const T & operator[](const std::size_t i) const
    {
        if (i >= _size) {
            std::stringstream stream;
            stream << "Index " << i << " is out of range [0, " << _size-1 << ']';
            throw std::out_of_range(stream.str());
        }
        return *(_array + i * _offset);
    }

    inline T & operator[](const std::size_t i)
    {
        if (i >= _size) {
            std::stringstream stream;
            stream << "Index " << i << " is out of range [0, " << _size-1 << ']';
            throw std::out_of_range(stream.str());
        }
        return *(_array + i * _offset);
    }

    inline const T & slice(const std::size_t i) const
    {
        return *this[i];
    }

    inline T & slice(const std::size_t i)
    {
        return *this[i];
    }

    inline std::size_t size() const
    {
        return _size;
    }

    inline std::size_t accumulatedOffset(const std::size_t i = 0) const
    {
        if (i >= _size) {
            std::stringstream stream;
            stream << "Index " << i << " is out of range [0, " << _size-1 << ']';
            throw std::out_of_range(stream.str());
        }
        return _accumulatedOffset + _offset * i;
    }

private:
    T                  *_array;
    std::size_t         _accumulatedOffset;

protected:
    std::size_t         _size;

private:
    std::size_t         _offset;
};

template < typename T, std::size_t D >
class MMArray : public MArray<T, D> {
public:
    MMArray()
        : MArray<T, D>()
        , _array(nullptr)
    { }

    MMArray(const std::size_t size[D])
    {
        resize(size);
    }

    MMArray(const MMArray<T, D> &mma)
    {
        *this = mma;
    }

    MMArray(MMArray &&mma)
        : MArray<T, D>(std::forward<MMArray<T, D>>(mma))
        , _array(std::move(mma._array))
    { }

    MMArray & operator=(const MMArray &mma)
    {
        if (&mma != this) {
            _array.reset(new T[mma.totalSize()]);
            std::memcpy(_array.get(), mma._array.get(), mma.totalSize() * sizeof(T));
            std::size_t size[D];
            mma.size(size);
            MArray<T, D>::operator=(MArray<T, D>(_array.get(), 0, size));
        }
        return *this;
    }

    MMArray & operator=(MMArray &&mma)
    {
        if (&mma != this) {
            MArray<T, D>::operator=(std::forward<MMArray<T, D>>(mma));
            _array = std::move(mma._array);
        }
        return *this;
    }
    
    inline void resize(const std::size_t size[D])
    {
        std::size_t total = size[0];
        for (std::size_t i = 1; i < D; ++i)
            total *= size[i];
        _array.reset(new T[total]);                                     // Be aware: data is LOST!
        MArray<T, D>::operator=(MArray<T, D>(_array.get(), 0, size));
    }

private:
    std::unique_ptr<T[]>    _array;
};

#endif /* multiple_array_hpp */

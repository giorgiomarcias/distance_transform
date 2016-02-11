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
            _offset[j-1] = _size[j] * _offset[j];
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
        for (std::size_t i = 0; i < D; ++i) {
            _size[i] = ma._size[i];
            _offset[i] = ma._offset[i];
        }
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
    
    inline void at(const std::size_t i, MArray<T, D-1> &s) const
    {
        if (i >= _size[0]) {
            std::stringstream stream;
            stream << "Index " << i << " is out of range [0, " << _size[0]-1 << ']';
            throw std::out_of_range(stream.str());
        }
        s._array = _array + _offset[0] * i;
        s._accumulatedOffset = accumulatedOffset(i);
        for (std::size_t j = 1; j < D; ++j) {
            s._size[j-1] = _size[j];
            s._offset[j-1] = _offset[j];
        }
    }

    inline MArray<T, D-1> operator[](const std::size_t i) const
    {
        MArray<T, D-1> s;
        at(i, s);
        return s;
    }
    
    inline MArray<T, D-1> slice(const std::size_t d, const std::size_t i) const
    {
        MArray<T, D-1> s;
        slice(d, i, s);
        return s;
    }

    inline void slice(const std::size_t d, const std::size_t i, MArray<T, D-1> &s) const
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
        s._array = _array + _offset[d] * i;
        s._accumulatedOffset = accumulatedOffset(i, d);
        std::size_t k = 0;
        for (std::size_t j = 0; j < d; ++j, ++k) {
            s._size[k] = _size[j];
            s._offset[k] = _offset[j];
        }
        for (std::size_t j = d+1; j < D; ++j, ++k) {
            s._size[k] = _size[j];
            s._offset[k] = _offset[j];
        }
    }
    
    inline MArray<T, D> permute(const std::size_t order[D])
    {
        bool included[D];
        for (std::size_t d = 0; d < D; ++d)
            included[d] = false;
        for (std::size_t d = 0; d < D; ++d) {
            if (order[d] >= D) {
                std::stringstream stream;
                stream << "Index " << order[d] << " is out of range [0, " << D-1 << ']';
                throw std::out_of_range(stream.str());
            }
            if (included[order[d]]) {
                std::stringstream stream;
                stream << "Dimension " << order[d] << " duplicated.";
                throw std::out_of_range(stream.str());
            }
            included[order[d]] = true;
        }
        std::size_t size[D];
        std::size_t offset[D];
        for (std::size_t d = 0; d < D; ++d) {
            size[d] = _size[order[d]];
            offset[d] = _offset[order[d]];
        }
        return MArray<T, D>(_array, _accumulatedOffset, size, offset);
    }

    inline std::size_t size(const std::size_t d = 0) const
    {
        if (d >= D) {
            std::stringstream stream;
            stream << "Index " << d << " is out of range [0, " << D-1 << ']';
            throw std::out_of_range(stream.str());
        }
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
    friend class MArray<T, D+1>;
    T                  *_array;
    std::size_t         _accumulatedOffset;
    std::size_t         _size[D];
    std::size_t         _offset[D];
};

template < typename T >
class MArray<T, 1> {
public:
    MArray()
        : _array(nullptr)
        , _accumulatedOffset(0)
    {
        _size[0] = 0;
        _offset[0] = 0;
    }

    MArray(const T *array, const std::size_t accumulatedOffset, const std::size_t size)
        : _array(array)
        , _accumulatedOffset(_array ? accumulatedOffset : 0)
    {
        _size[0] = _array ? size : 0;
        _offset[0] = _array ? 1 : 0;
    }

    MArray(const T *array, const std::size_t accumulatedOffset, const std::size_t size, const std::size_t offset)
        : _array(array)
        , _accumulatedOffset(_array ? accumulatedOffset : 0)
    {
        _size[0] = _array ? size : 0;
        _offset[0] = _array ? offset : 0;
    }

    MArray(const T *array, const std::size_t accumulatedOffset, const std::size_t size[1])
        : _array(const_cast<T*>(array))
        , _accumulatedOffset(_array ? accumulatedOffset : 0)
    {
        _size[0] = _array ? size[0] : 0;
        _offset[0] = _array ? 1 : 0;
    }

    MArray(const T *array, const std::size_t accumulatedOffset, const std::size_t size[1], const std::size_t offset[1])
        : _array(const_cast<T*>(array))
        , _accumulatedOffset(_array ? accumulatedOffset : 0)
    {
        _size[0] = _array ? size[0] : 0;
        _offset[0] = _array ? offset[0] : 0;
    }

    MArray(const MArray &ma)
        : _array(ma._array)
        , _accumulatedOffset(ma._accumulatedOffset)
    {
        _size[0] = ma._size[0];
        _offset[0] = ma._offset[0];
    }

    MArray(MArray &&ma)
        : _array(ma._array)
        , _accumulatedOffset(ma._accumulatedOffset)
    {
        _size[0] = ma._size[0];
        _offset[0] = ma._offset[0];
        ma._array = nullptr;
        ma._accumulatedOffset = 0;
        ma._size[0] = 0;
        ma._offset[0] = 0;
    }

    MArray & operator=(const MArray &ma)
    {
        if (&ma != this) {
            _array = ma._array;
            _accumulatedOffset = ma._accumulatedOffset;
            _size[0] = ma._size[0];
            _offset[0] = ma._offset[0];
        }
        return *this;
    }

    MArray & operator=(MArray &&ma)
    {
        if (&ma != this) {
            _array = ma._array;
            _accumulatedOffset = ma._accumulatedOffset;
            _size[0] = ma._size[0];
            _offset[0] = ma._offset[0];
            ma._array = nullptr;
            ma._accumulatedOffset = 0;
            ma._size[0] = 0;
            ma._offset[0] = 0;
        }
        return *this;
    }

    inline const T & operator[](const std::size_t i) const
    {
        if (i >= _size[0]) {
            std::stringstream stream;
            stream << "Index " << i << " is out of range [0, " << _size[0]-1 << ']';
            throw std::out_of_range(stream.str());
        }
        return *(_array + i * _offset[0]);
    }

    inline T & operator[](const std::size_t i)
    {
        if (i >= _size[0]) {
            std::stringstream stream;
            stream << "Index " << i << " is out of range [0, " << _size[0]-1 << ']';
            throw std::out_of_range(stream.str());
        }
        return *(_array + i * _offset[0]);
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
        return _size[0];
    }

    inline std::size_t accumulatedOffset(const std::size_t i = 0) const
    {
        if (i >= _size[0]) {
            std::stringstream stream;
            stream << "Index " << i << " is out of range [0, " << _size[0]-1 << ']';
            throw std::out_of_range(stream.str());
        }
        return _accumulatedOffset + _offset[0] * i;
    }

private:
    friend class MArray<T, 2>;
    T                  *_array;
    std::size_t         _accumulatedOffset;
    std::size_t         _size[1];
    std::size_t         _offset[1];
};

template < typename T, std::size_t D >
class MMArray : public MArray<T, D> {
public:
    MMArray()
        : MArray<T, D>()
        , _arrayPtr(nullptr)
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
        , _arrayPtr(std::move(mma._arrayPtr))
    { }

    MMArray & operator=(const MMArray &mma)
    {
        if (&mma != this) {
            std::size_t size[D];
            mma.size(size);
            resize(size);
            std::memcpy(_arrayPtr.get(), mma._arrayPtr.get(), MArray<T, D>::totalSize() * sizeof(T));
        }
        return *this;
    }

    MMArray & operator=(MMArray &&mma)
    {
        if (&mma != this) {
            MArray<T, D>::operator=(std::forward<MMArray<T, D>>(mma));
            _arrayPtr = std::move(mma._arrayPtr);
        }
        return *this;
    }
    
    inline void resize(const std::size_t size[D])
    {
        std::size_t total = size[0];
        for (std::size_t i = 1; i < D; ++i)
            total *= size[i];
        if (total > 0) {
            _arrayPtr.reset(new T[total]);                                     // Be aware: data is LOST!
            MArray<T, D>::operator=(MArray<T, D>(_arrayPtr.get(), 0, size));
        } else {
            clear();
        }
    }
    
    inline void clear()
    {
        _arrayPtr.reset(nullptr);
        MArray<T, D>::operator=(MArray<T, D>());
    }

private:
    std::unique_ptr<T[]>    _arrayPtr;
};

#endif /* multiple_array_hpp */

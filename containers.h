#if !defined(H_CONTAINERS)
#define H_CONTAINERS

#include "nocopy_integer.h"
#include <array>
#include <cstdint>
#include <functional>
#include <linbox/randiter/random-prime.h>
#include <ostream>

using namespace std;

template <class T>
class PtrVector : public vector<T*> {
public:
    PtrVector(size_t len)
        : vector<T*>(len)
    {
    }
    T& val(size_t i) const { return *(this->operator[](i)); }
    T*& ptr(size_t i) { return this->operator[](i); }
    ~PtrVector()
    {
        this->erase(this->begin());
    }
};

template <class T>
class NumPtrVector : public PtrVector<T> {
public:
    NumPtrVector(size_t len)
        : PtrVector<T>(len)
    {
    }

    T& operator[](const size_t i) const
    {
        return this->val(i);
    }

    friend ostream& operator<<(ostream& out, const NumPtrVector<T>& arr)
    {
        size_t N = arr.size();
        out << "[";
        for (size_t i = 0; i < N; i++) {
            out << arr[i];
            if (i != N - 1) {
                out << " , ";
            }
        }
        out << "]";
        return out;
    }
};

template <class T>
using PtrAllocator = std::function<T*(size_t index)>;

// used to store a list of large primes as pointers
template <class T, size_t N>
class PtrArray : public array<T*, N> {
public:
    PtrArray() = default;

    PtrArray(const PtrAllocator<T>& allocator)
    {
        for (size_t i = 0; i < N; i++) {
            this->operator[](i) = allocator(i);
        }
    }

    PtrVector<T>* EXPENSIVE_NEW_PTR_VECTOR() const
    {
        auto vec = new PtrVector<T>(N);
        for (size_t i = 0; i < N; i++) {
            vec->ptr(i) = this->ptr(i);
        }
        return vec;
    }

    // default C++ array<T,N> destructor frees the memory
    ~PtrArray() = default;

    PtrArray(const PtrArray&) = delete;
    const PtrArray& operator=(const PtrArray&) = delete;

    T& val(size_t i) const { return *(this->operator[](i)); }
    T* ptr(size_t i) const { return this->operator[](i); }

    size_t count() const { return N; }

    friend ostream& operator<<(ostream& out, const PtrArray<T, N>& arr)
    {
        out << "[";
        for (size_t i = 0; i < N; i++) {
            out << arr[i];
            if (i != N - 1) {
                out << " , ";
            }
        }
        out << "]";
        return out;
    };
};

template <class T, size_t N>
class NumPtrArray : public PtrArray<T, N> {
public:
    NumPtrArray() = default;

    NumPtrArray(const PtrAllocator<T>& allocator)
        : PtrArray<T, N>(allocator)
    {
    }

    T& operator[](const size_t i) const
    {
        return this->val(i);
    }

    NumPtrVector<T>* EXPENSIVE_NEW_NUM_PTR_VECTOR() const
    {
        auto vec = new NumPtrVector<T>(N);
        for (size_t i = 0; i < N; i++) {
            vec->ptr(i) = this->ptr(i);
        }
        return vec;
    }

    vector<Givaro::Integer>* EXPENSIVE_NEW_INTEGER_VECTOR() const
    {
        auto vec = new vector<Givaro::Integer>(N);
        for (size_t i = 0; i < N; i++) {
            vec->operator[](i) = this->val(i);
        }
        return vec;
    }

    friend ostream& operator<<(ostream& out, const NumPtrArray<T, N>& arr)
    {
        out << "[";
        for (size_t i = 0; i < N; i++) {
            out << arr[i];
            if (i != N - 1) {
                out << " , ";
            }
        }
        out << "]";
        return out;
    }
    ~NumPtrArray() = default;

    NumPtrArray(const NumPtrArray&) = delete;
    const NumPtrArray& operator=(const NumPtrArray&) = delete;
};

template <class T, size_t N>
class ConstNumPtrArray : public NumPtrArray<const T, N> {
private:
    Givaro::Integer _product;
    Givaro::Integer _sum;

public:
    const Givaro::Integer& product;
    const Givaro::Integer& sum;

    ConstNumPtrArray(const PtrAllocator<const T>& allocator)
        : NumPtrArray<const T, N>(allocator)
        , _product(N > 0 ? 1 : 0)
        , _sum(0)
        , product(_product)
        , sum(_sum)
    {
        for (size_t i = 0; i < N; i++) {
            _product *= this->operator[](i);
            _sum += this->operator[](i);
        }
    }

    const T* ptr(const size_t& i) const
    {
        return PtrArray<const T, N>::operator[](i);
    }

    const T& operator[](const size_t& i) const
    {
        return *(PtrArray<const T, N>::operator[](i));
    }

    NumPtrVector<T>* EXPENSIVE_NEW_NUM_PTR_VECTOR() const
    {
        auto vec = new NumPtrVector<T>(N);
        for (size_t i = 0; i < N; i++) {
            vec->ptr(i) = new T{ this->val(i) };
        }
        return vec;
    }

    friend ostream& operator<<(ostream& out,
        const ConstNumPtrArray<T, N>& arr)
    {
        out << "[";
        for (size_t i = 0; i < N; i++) {
            out << arr[i];
            if (i != N - 1) {
                out << " , ";
            }
        }
        out << "]";
        return out;
    }

    ~ConstNumPtrArray() = default;

    ConstNumPtrArray(const ConstNumPtrArray&) = delete;
    const ConstNumPtrArray& operator=(const ConstNumPtrArray&) = delete;
};

#endif // H_CONTAINERS

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
    // default C++ array<T,N> destructor frees the memory
    ~PtrArray() = default;

    PtrArray(const PtrArray&) = delete;
    const PtrArray& operator=(const PtrArray&) = delete;

    const T& ref(size_t i) const { return *(this->operator[](i)); }

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

    T* ptr(const size_t& i) const
    {
        return PtrArray<T, N>::operator[](i);
    }

    T& operator[](const size_t& i) const
    {
        return *(PtrArray<T, N>::operator[](i));
    }

    friend ostream& operator<<(ostream& out,
        const NumPtrArray<T, N>& arr)
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

    const vector<Integer> EXPENSIVE_TO_VECTOR() const
    {
        vector<Integer> vec(N);
        for (size_t i = 0; i < N; i++) {
            vec[i] = this->operator[](i);
        }
        return vec;
    }

    const T* ptr(const size_t& i) const
    {
        return PtrArray<const T, N>::operator[](i);
    }

    const T& operator[](const size_t& i) const
    {
        return *(PtrArray<const T, N>::operator[](i));
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

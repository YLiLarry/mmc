#if !defined(H_PRIME_GEN)
#define H_PRIME_GEN

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

    friend ostream& operator<<(ostream& out, const PtrArray<T, N>& arr)
    {
        out << "[";
        for (size_t i = 0; i < N; i++) {
            out << arr[i];
            if (i != N - 1) {
                out << ",";
            }
        }
        out << "]";
        return out;
    };
};

template <class T, size_t N>
class ConstNumPtrArray : public PtrArray<const T, N> {
private:
    T _product;
    T _sum;

public:
    const T& product;
    const T& sum;

    ConstNumPtrArray(const PtrAllocator<const T>& allocator)
        : PtrArray<const T, N>(allocator)
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

    friend ostream& operator<<(ostream& out,
        const ConstNumPtrArray<T, N>& arr)
    {
        out << "[";
        for (size_t i = 0; i < N; i++) {
            out << arr[i];
            if (i != N - 1) {
                out << ",";
            }
        }
        out << "]";
        return out;
    }

    ~ConstNumPtrArray() = default;

    ConstNumPtrArray(const ConstNumPtrArray&) = delete;
    const ConstNumPtrArray& operator=(const ConstNumPtrArray&) = delete;
};

// class Primeallocator {
//     T _allocator;

// public:
//     Primeallocator()
//         : _allocator(LinBox::RandomPrimeIter(B))
//     {
//     }
//     T* next()
//     {
//         return _allocator.random
//     }
// }

// generate N random primes each of B-bits length
template <typename T, size_t N, uint_fast64_t B>
class PrimeGen : public ConstNumPtrArray<T, N> {
protected:
    static LinBox::RandomPrimeIter _primeItr;

public:
    PrimeGen(const PtrAllocator<T>& allocator)
        : ConstNumPtrArray<T, N>(allocator)
    {
    }
    ~PrimeGen() = default;
    PrimeGen(const PrimeGen&) = delete;
    PrimeGen& operator&=(const PrimeGen&) = delete;
};

template <typename T, size_t N, uint_fast64_t B>
LinBox::RandomPrimeIter PrimeGen<T, N, B>::_primeItr{ B };
// template <class T, uint_fast8_t B>
// class Primeallocator {
//     T _allocator;

// public:
//     Primeallocator()
//         : _allocator(LinBox::RandomPrimeIter(B))
//     {
//     }
//     T* next()
//     {
//         return _allocator.random
//     }
// }

// return an array containing primes in type T, eg. T = int64_t
// all primes are exactly B-b   its long, ie, between [2^(B-1), 2^B-1]
// eg. PrimeGen<uint64_t, 10, 64> p{};
//     generates 10 primes each exactly 64 bits
// access primes using p[index]
template <typename T, size_t N, uint_fast64_t B>
class PrimeGenExact : public PrimeGen<T, N, B> {
public:
    PrimeGenExact()
        : PrimeGen<T, N, B>([&](size_t index) {
            T* t = new T;
            PrimeGen<T, N, B>::_primeItr.random_exact(*t);
            return t;
        })
    {
    }

    ~PrimeGenExact() = default;

    PrimeGenExact(const PrimeGenExact&) = delete;
    PrimeGenExact& operator&=(const PrimeGenExact&) = delete;
};

// return an array containing primes in type T, eg. T = int64_t
// all primes are at most B-bits long, ie, between [2, 2^B-1]
// eg. PrimeGen<uint64_t, 10, 64> p{};
//     generates 10 primes each at most 64 bits
// access primes using p[index]
template <typename T, size_t N, uint_fast64_t B>
class PrimeGenMost : public PrimeGen<T, N, B> {
public:
    PrimeGenMost()
        : PrimeGen<T, N, B>([this](size_t index) {
            T* t = new T;
            PrimeGen<T, N, B>::_primeItr.random(*t);
            return t;
        })
    {
    }
    ~PrimeGenMost() = default;
    PrimeGenMost(const PrimeGenMost&) = delete;
    PrimeGenMost& operator&=(const PrimeGenMost&) = delete;
};

#endif // H_PRIME_GEN

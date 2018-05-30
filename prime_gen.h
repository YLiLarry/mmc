#if !defined(H_PRIME_GEN)
#define H_PRIME_GEN

#include "nocopy_integer.h"
#include <array>
#include <cstdint>
#include <linbox/randiter/random-prime.h>
#include <ostream>

using namespace std;

// template <typename T, size_t N>
// ostream&
// operator<<(ostream& out, array<T, N>& arr)
// {
//     out << "[";
//     for (size_t i = 0; i < N; i++) {
//         out << arr[i];
//         if (i != N - 1) {
//             out << ",";
//         }
//     }
//     out << "]";
//     return out;
// };

// used to store a list of large primes as pointers
template <class T, size_t N>
class ConstPtrArray : public array<T const*, N> {
public:
    ConstPtrArray(std::function<const T*(size_t index)> generator)
    {
        for (size_t i = 0; i < N; i++) {
            (*this)[i] = generator(i);
        }
    }
    // default C++ array<T,N> destructor frees the memory
    ~ConstPtrArray() = default;

    ConstPtrArray(const ConstPtrArray&) = delete;
    const ConstPtrArray& operator=(const ConstPtrArray&) = delete;

    const T& ref(size_t i) const { return *((*this)[i]); }

    friend ostream& operator<<(ostream& out, const ConstPtrArray<T, N>& arr)
    {
        out << "[";
        for (size_t i = 0; i < N; i++) {
            out << arr.ref(i);
            if (i != N - 1) {
                out << ",";
            }
        }
        out << "]";
        return out;
    };
};

// helper class to compute the product in initialization list
template <class T, size_t N>
class ConstNumPtrArray : public ConstPtrArray<T, N> {
public:
    T product;
    T sum;
    ConstNumPtrArray(std::function<const T*(size_t index)> generator)
        : ConstPtrArray<T, N>(generator)
        , product(N > 0 ? 1 : 0)
        , sum(0)
    {
        for (int i = 0; i < N; i++) {
            product *= this->ref(i);
            sum += this->ref(i);
        }
    }
};

template <class T, size_t N>
class ConstNumArray {
protected:
    const ConstNumPtrArray<T, N> _ptrarr;

public:
    const T& product;
    const T& sum;
    ConstNumArray(std::function<const T*(size_t index)> generator)
        : _ptrarr(generator)
        , product(_ptrarr.product)
        , sum(_ptrarr.sum)
    {
    }
    const T& operator[](const size_t& i) const { return _ptrarr.ref(i); };
    friend ostream& operator<<(ostream& out, ConstNumArray<T, N>& _ptrarr)
    {
        out << _ptrarr._ptrarr;
        return out;
    }
};

// class PrimeGenerator {
//     T _generator;

// public:
//     PrimeGenerator()
//         : _generator(LinBox::RandomPrimeIter(B))
//     {
//     }
//     T* next()
//     {
//         return _generator.random
//     }
// }

// generate N random primes each of B-bits length
template <typename T, size_t N, uint_fast64_t B>
class PrimeGen : public ConstNumArray<T, N> {
protected:
    static LinBox::RandomPrimeIter _primeItr;

public:
    PrimeGen(std::function<const T*(size_t index)> generator)
        : ConstNumArray<T, N>(generator)
    {
    }
    ~PrimeGen() = default;
    PrimeGen(const PrimeGen&) = delete;
    PrimeGen& operator&=(const PrimeGen&) = delete;
};

template <typename T, size_t N, uint_fast64_t B>
LinBox::RandomPrimeIter PrimeGen<T, N, B>::_primeItr{ B };
// template <class T, uint_fast8_t B>
// class PrimeGenerator {
//     T _generator;

// public:
//     PrimeGenerator()
//         : _generator(LinBox::RandomPrimeIter(B))
//     {
//     }
//     T* next()
//     {
//         return _generator.random
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
        : PrimeGen<T, N, B>{ [&](size_t index) {
            T* t = new T;
            PrimeGen<T, N, B>::_primeItr.random_exact(*t);
            return t;
        } }
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

#if !defined(H_PRIME_GEN)
#define H_PRIME_GEN

#include "containers.h"
#include "nocopy_integer.h"
#include <array>
#include <cstdint>
#include <functional>
#include <givaro/givintprime.h>
#include <linbox/randiter/random-prime.h>
#include <ostream>
using namespace std;

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
LinBox::RandomPrimeIter PrimeGen<T, N, B>::_primeItr{ B, (uint64_t)BaseTimer::seed() };

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
            Integer p;
            PrimeGen<T, N, B>::_primeItr.random_exact(p);
            T* t = new T(p);
            return t;
        })
    {
#ifdef ASSERT_MMC
        Givaro::IntPrimeDom primeDom;
        for (size_t i = 0; i < N; i++) {
            assert(primeDom.isprime(this->val(i)));
        }
#endif
    }

    friend ostream& operator<<(ostream& out, const PrimeGenExact<T, N, B>& arr)
    {
        out << "PrimeGenExact [";
        for (size_t i = 0; i < N; i++) {
            out << arr[i];
            if (i != N - 1) {
                out << ",";
            }
        }
        out << "]" << endl
            << " - count: " << arr.count() << endl
            << " - max bit length " << B << endl
            << " - product bit length: " << arr.product.bitsize() << endl;
        return out;
    }

    ~PrimeGenExact() = default;

    PrimeGenExact(const PrimeGenExact&) = delete;
    PrimeGenExact& operator&=(const PrimeGenExact&) = delete;
};

// template <size_t N, uint_fast64_t B>
// class PrimeGenExact<LInteger, N, B> : public PrimeGen<LInteger, N, B> {
// public:
//     PrimeGenExact()
//         : PrimeGen<LInteger, N, B>([&](size_t index) {
//             LInteger* p = new LInteger();
//             PrimeGen<LInteger, N, B>::_primeItr.random_exact(*p);
//             return p;
//         })
//     {
//     }

//     ~PrimeGenExact() = default;

//     PrimeGenExact(const PrimeGenExact&) = delete;
//     PrimeGenExact& operator&=(const PrimeGenExact&) = delete;
// };

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

#if !defined(H_PRIME_GEN)
#define H_PRIME_GEN

#include "nocopy_integer.h"
#include <array>
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
protected:
    ConstPtrArray* _self = this;

public:
    ConstPtrArray(std::function<const T*()> generator)
    {
        for (int i = 0; i < N; i++) {
            (*this)[i] = generator();
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

// helper class to compute the product in intialization list
template <class T, size_t N>
class PtrProduct {
public:
    T value;
    PtrProduct(const ConstPtrArray<T, N>& arr)
        : value(N > 0 ? 1 : 0)
    {
        for (int i = 0; i < N; i++) {
            value *= arr.ref(i);
        }
    }
};

template <class T, size_t N>
class ConstNumArray : public ConstPtrArray<T, N> {
private:
    const PtrProduct<T, N> _product;

public:
    const T& product;
    ConstNumArray(std::function<const T*()> generator)
        : ConstPtrArray<T, N>(generator)
        , _product(*this)
        , product(_product.value)
    {
    }
};

// generate N random primes each of B-bits length
template <typename T, size_t N, uint_fast8_t B>
class PrimeGen {
protected:
    LinBox::RandomPrimeIter _primeItr{ B };
    const ConstNumArray<T, N> _primes;

public:
    const T& product;
    PrimeGen(std::function<const T*()> generator)
        : _primes(generator)
        , product(_primes.product)
    {
    }
    ~PrimeGen() = default;
    PrimeGen(const PrimeGen&) = delete;
    PrimeGen& operator&=(const PrimeGen&) = delete;

    friend ostream& operator<<(ostream& out, PrimeGen<T, N, B>& pg)
    {
        out << pg._primes;
        return out;
    }

    const T& operator[](const size_t& i) const { return _primes[i]; };
};

// return an array containing primes in type T, eg. T = int64_t
// all primes are exactly B-b   its long, ie, between [2^(B-1), 2^B-1]
// eg. PrimeGen<uint64_t, 10, 64> p{};
//     generates 10 primes each exactly 64 bits
// access primes using p[index]
template <typename T, size_t N, uint_fast8_t B>
class PrimeGenExact : public PrimeGen<T, N, B> {
public:
    PrimeGenExact()
        : PrimeGen<T, N, B>([=]() {
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
template <typename T, size_t N, uint_fast8_t B>
class PrimeGenMost : public PrimeGen<T, N, B> {
public:
    PrimeGenMost()
        : PrimeGen<T, N, B>([=]() {
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

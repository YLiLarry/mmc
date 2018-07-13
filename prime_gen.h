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
#include <iterator>
#include "coprime_gen_abstract.h"

// return an array containing unique primes in type T, eg. T = int64_t
// all primes are exactly B-b   its long, ie, between [2^(B-1), 2^B-1]
// eg. PrimeGenMost<uint64_t, 10, 64> p{};
//     generates 10 primes each exactly 64 bits
// access primes using p[index]
template <typename T>
class PrimeGenMost : public CoprimeGenAbstract<T>
{
  protected:
    Givaro::IntPrimeDom _int_prime_domain;
    Givaro::Integer _product = 1;
    uint_fast64_t _max_bitsize;
    T _max;

  public:
    inline virtual Givaro::Integer product() const override { return _product; }
    inline virtual uint_fast64_t product_bitsize() const override { return _product.bitsize(); }
    inline virtual uint_fast64_t max_bitsize() const override { return _max_bitsize; }
    inline virtual const T &max() const override { return _max; }

  public:
    PrimeGenMost(uint_fast64_t product_bound, uint_fast64_t max_bound)
    {
        Givaro::Integer prime_bound = 1;
        // linbox bug: uint64 could be undefined
        prime_bound <<= max_bound;
        _int_prime_domain.prevprimein(prime_bound);
        _max = prime_bound;
        _max_bitsize = prime_bound.bitsize();
        while (product_bitsize() <= product_bound)
        {
            if (prime_bound < 2)
            {
                cerr << "we ran out of primes, consider increasing max_bound." << endl;
                abort();
            }
            vector<T>::push_back(prime_bound);
            _product *= prime_bound;
            _int_prime_domain.prevprimein(prime_bound);
        }

#ifdef TEST_MMC
        size_t N = this->count();
        Givaro::IntPrimeDom primeDom;
        for (size_t i = 0; i < N; i++)
        {
            if (!primeDom.isprime(this->val(i)))
            {
                cerr << "assertion failed " << this->val(i) << " is not a prime" << endl;
            }
            assert(primeDom.isprime(this->val(i)));
        }
#endif
        assert(_max >= 2);
        assert(_product >= 2);
    }

    ~PrimeGenMost() = default;

    PrimeGenMost(const PrimeGenMost &) = delete;
    PrimeGenMost &operator&=(const PrimeGenMost &) = delete;
};

#endif // H_PRIME_GEN

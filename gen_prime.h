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
#include "gen_coprime_abstract.h"

// return an array containing unique primes in type T, eg. T = int64_t
// all primes are exactly B-b   its long, ie, between [2^(B-1), 2^B-1]
// eg. GenPrimeMost<uint64_t, 10, 64> p{};
//     generates 10 primes each exactly 64 bits
// access primes using p[index]
template <typename T>
class GenPrimeMost : public GenCoprimeAbstract<T>
{
  protected:
    Givaro::IntPrimeDom _int_prime_domain;
    Givaro::Integer _product = 1;

  public:
    inline virtual Givaro::Integer product() const override { return _product; }
    inline virtual uint_fast64_t product_bitsize() const override { return product().bitsize(); }
    inline virtual uint_fast64_t max_bitsize() const override { return Givaro::Integer(max()).bitsize(); }
    inline virtual const T &max() const override { return this->operator[](0); }

  public:
    GenPrimeMost(uint_fast64_t product_bound, uint_fast64_t max_bound)
    {
        assert(product_bound > max_bound);
        assert(max_bound > 1);
        Givaro::Integer prime_bound = 1;
        // linbox bug: uint64 could be undefined
        prime_bound <<= max_bound;
        _int_prime_domain.prevprimein(prime_bound);
        while (product_bitsize() <= product_bound)
        {
            if (prime_bound < 2)
            {
                cerr << "We ran out of primes, consider increasing max_bound." << endl;
                abort();
            }
            vector<T>::push_back(prime_bound);
            _product *= prime_bound;
            _int_prime_domain.prevprimein(prime_bound);
        }

#if CHECK_MMC
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
        assert(max() >= 2);
        assert(product() >= 2);
    }

    ~GenPrimeMost() = default;

    GenPrimeMost(const GenPrimeMost &) = delete;
    GenPrimeMost &operator&=(const GenPrimeMost &) = delete;
};

#endif // H_PRIME_GEN

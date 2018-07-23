#if !defined(H_GEN_MARGE_LEAST)
#define H_GEN_MARGE_LEAST

#include <gmp++/gmp++.h>
#include "gen_prime.h"
#include <cmath>
#include <iterator>

// generates N co-primes of form 2^n - 1,
// starting from the largets one 2^max_bound, and count down
// until the product is larger than 2^product_bound
class GenMargeLeast : public GenCoprimeAbstract<Givaro::Integer>
{
  protected:
    Givaro::IntPrimeDom _int_prime_domain;
    Givaro::Integer _product = 1;

  public:
    inline virtual Givaro::Integer product() const override { return _product; }
    inline virtual uint_fast64_t product_bitsize() const override { return product().bitsize(); }
    inline virtual uint_fast64_t max_bitsize() const override { return max().bitsize(); }
    inline virtual const Givaro::Integer &max() const override { return this->operator[](this->count() - 1); }

  public:
    GenMargeLeast(uint_fast64_t product_bound, uint_fast64_t min_bound)
    {
        assert(product_bound > min_bound && "The product of moduli must be greater than any moduli.");
        assert(min_bound >= 0);
        Givaro::Integer prime_expo = min_bound;
        _int_prime_domain.nextprimein(prime_expo);
        while (product_bitsize() < product_bound)
        {
            Givaro::Integer marge(1);
            marge <<= (prime_expo % product_bound);
            marge--;
            // marge = 2^prime_expo - 1
            this->push_back(marge);
            _product *= marge;
            _int_prime_domain.nextprimein(prime_expo);
        }
#if CHECK_MMC
        size_t N = this->count();
        for (size_t i = 0; i < N; i++)
        {
            for (size_t j = 0; j < N; j++)
            {
                if (i != j && gcd(this->val(i), this->val(j)) != 1)
                {
                    cerr << "GenMargeLeast is not generating coprimes." << endl
                         << " - got: " << this->val(i) << " and " << this->val(j) << endl;
                    abort();
                }
            }
        }
#endif
        assert(max() >= 3);
        assert(product() >= 3);
    }

    ~GenMargeLeast() = default;
    GenMargeLeast(GenMargeLeast &) = delete;
};

#endif // H_GEN_MARGE_LEAST

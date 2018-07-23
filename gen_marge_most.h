#if !defined(H_GEN_MARGE_MOST)
#define H_GEN_MARGE_MOST

#include <gmp++/gmp++.h>
#include "gen_prime.h"
#include <cmath>
#include <iterator>

// generates N co-primes of form 2^n - 1,
// starting from the largets one 2^max_bound, and count down
// until the product is larger than 2^product_bound
class GenMargeMost : public GenCoprimeAbstract<Givaro::Integer>
{
  protected:
    Givaro::IntPrimeDom _int_prime_domain;
    Givaro::Integer _product = 1;

  public:
    inline virtual Givaro::Integer product() const override { return _product; }
    inline virtual uint_fast64_t product_bitsize() const override { return product().bitsize(); }
    inline virtual uint_fast64_t max_bitsize() const override { return max().bitsize(); }
    inline virtual const Givaro::Integer &max() const override { return this->operator[](0); }

  public:
    GenMargeMost(uint_fast64_t product_bound, uint_fast64_t max_bound)
    {
        assert(product_bound > max_bound && "The product of moduli must be greater than any moduli.");
        assert(max_bound > 1 && "Any moduli must be at least 2 bits");
        Givaro::Integer prime_expo = max_bound;
        _int_prime_domain.prevprimein(prime_expo);
        while (product_bitsize() < product_bound)
        {
            if (prime_expo < 2)
            {
                cerr << "Failure to generate coprimes: We ran out of marge coprimes. Consider increasing max bitsize bound." << endl;
                abort();
            }
            Givaro::Integer marge(1);
            marge <<= (prime_expo % max_bound);
            marge--;
            // marge = 2^prime_expo - 1
            this->push_back(marge);
            _product *= marge;
            _int_prime_domain.prevprimein(prime_expo);
        }
#if CHECK_MMC
        size_t N = this->count();
        for (size_t i = 0; i < N; i++)
        {
            for (size_t j = 0; j < N; j++)
            {
                if (i != j && gcd(this->val(i), this->val(j)) != 1)
                {
                    cerr << "GenMargeMost is not generating coprimes." << endl
                         << " - got: " << this->val(i) << " and " << this->val(j) << endl;
                    abort();
                }
            }
        }
#endif
        assert(max() >= 3);
        assert(product() >= 3);
    }

    ~GenMargeMost() = default;
    GenMargeMost(GenMargeMost &) = delete;
};

#endif // H_GEN_MARGE_MOST

#if !defined(H_FERMAT_GEN)
#define H_FERMAT_GEN

#include <gmp++/gmp++.h>
#include "gen_prime.h"
#include <cmath>
#include <vector>

// generates N co-primes of form 2^(2^n) + 1 where
// the largest one is at most most_bits, and
// the product of them are at most 2^product_bits bits long,
// since there are not many fermat numbers we don't generate randomly
// we just start from the largest and count down until the product is large enough
class GenFermatMost : public GenCoprimeAbstract<Givaro::Integer>
{
  protected:
    Givaro::IntPrimeDom _int_prime_domain;
    Givaro::Integer _product = 1;

  public:
    inline virtual Givaro::Integer product() const override { return _product; }
    inline virtual uint_fast64_t product_bitsize() const override { return product().bitsize(); }
    inline virtual uint_fast64_t max_bitsize() const override { return max().bitsize(); }
    inline virtual const Givaro::Integer &max() const override { return this->operator[](0); }

    // because it never makes sense to use a single Fermat number as the level 1 moduli,
    // and Fermat numbers have the property, F{n} = ! F{n-1}, we favor using (n-1) Fermat 
    // numbers under the max_bound rather than F{n}.
  public:
    GenFermatMost(uint_fast64_t product_bound, uint_fast64_t max_bound)
    {
        assert(product_bound > max_bound);
        assert(max_bound > 1);
        assert(std::floor(Givaro::logtwo(max_bound)) == std::ceil(Givaro::logtwo(max_bound)) && 
                "For Fermat numbers the max_bound (max exponent) must be a power of two");

        max_bound <<= 1; // find the previous power of 2
        while (product_bitsize() <= product_bound)
        {
            if (max_bound < 1)
            {
                std::cerr << "We ran out of Fermat numbers, consider increasing max_bound" << std::endl;
                std::cerr << " - currently generated: " << std::endl
                          << *this << std::endl;
            }
            assert(max_bound >= 1);
            Givaro::Integer t = 1;
            t <<= max_bound;
            t++;
            this->push_back(t); // t = 2^(2^max_bound) + 1
            _product *= t;
            max_bound <<= 1; // find the previous power of 2
        }
#if CHECK_MMC
        size_t N = this->size();
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
    }

    ~GenFermatMost() = default;
    GenFermatMost(GenFermatMost &) = delete;
};

#endif // H_FERMAT_GEN

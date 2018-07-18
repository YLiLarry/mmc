#if !defined(H_FERMAT_GEN)
#define H_FERMAT_GEN

#include <gmp++/gmp++.h>
#include "prime_gen.h"
#include <cmath>

// generates N co-primes of form 2^(2^n) + 1 where
// the largest one is at most most_bits, and
// the product of them are at most 2^product_bits bits long,
// since there are not many fermat numbers we don't generate randomly
// we just start from the largest and count down until the product is large enough
class FermatGenMost : public CoprimeGenAbstract<Givaro::Integer>
{
  protected:
    Givaro::IntPrimeDom _int_prime_domain;
    Givaro::Integer _product = 1;
    Givaro::Integer _max = 1;

  public:
    inline virtual Givaro::Integer product() const override { return _product; }
    inline virtual uint_fast64_t product_bitsize() const override { return _product.bitsize(); }
    inline virtual uint_fast64_t max_bitsize() const override { return _max.bitsize(); }
    inline virtual const Givaro::Integer &max() const override { return this->operator[](0); }

  public:
    FermatGenMost(uint_fast64_t product_bound, uint_fast64_t max_bound)
    {
        assert(product_bound > 1);
        assert(max_bound > 1);
        uint_fast64_t n = 1;
        max_bound >>= 1; // find the cloest power of 2
        while (n < max_bound)
        {
            n <<= 1;
        }
        while (product() <= product_bound)
        {
            Givaro::Integer t = 1;
            t <<= n;
            t++;
            this->push_back(t); // t = 2^(2^n) + 1
            _product *= t;
            n >>= 1;
        }
    }

    ~FermatGenMost() = default;
    FermatGenMost(FermatGenMost &) = delete;
};

#endif // H_FERMAT_GEN

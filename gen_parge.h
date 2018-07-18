#if !defined(H_GEN_PARGE)
#define H_GEN_PARGE

#include <gmp++/gmp++.h>
#include "gen_prime.h"
#include <cmath>
#include <vector>

// generates N co-primes of form 2^(2^n) + 1 where
// the largest one is at most most_bits, and
// the product of them are at most 2^product_bits bits long,
// since there are not many Parge numbers we don't generate randomly
// we just start from the largest and count down until the product is large enough
class GenPargeMost : public GenCoprimeAbstract<Givaro::Integer>
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
    GenPargeMost(uint_fast64_t product_bound, uint_fast64_t max_bound)
    {
        assert(product_bound > 1);
        assert(max_bound > 1);
        uint_fast64_t n = 1;
        max_bound >>= 1; // find the cloest power of 2
        while (n < max_bound)
        {
            n <<= 1;
        }
        while (product_bitsize() <= product_bound)
        {
            if (n < 1)
            {
                std::cerr << "We ran out of Parge numbers, consider increasing max_bound" << std::endl;
                std::cerr << " - currently generated: " << std::endl
                          << *this << std::endl;
            }
            assert(n >= 1);
            Givaro::Integer t = 1;
            t <<= n;
            t++;
            this->push_back(t); // t = 2^(2^n) + 1
            _product *= t;
            n >>= 1;
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

    ~GenPargeMost() = default;
    GenPargeMost(GenPargeMost &) = delete;
};

#endif // H_GEN_PARGE

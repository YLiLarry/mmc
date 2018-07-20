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
class GenPargeShift : public GenCoprimeAbstract<Givaro::Integer>
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
    GenPargeShift(uint_fast64_t product_bound, uint_fast64_t max_bound, uint_fast64_t coefficient)
    {
        assert(product_bound > max_bound && "The product of moduli must be greater than any moduli.");
        assert(max_bound > 1 && "Any moduli must be at least 2 bits");

        assert(max_bound % coefficient == 0 && "For Parge 2^(n*c) + 1 numbers the max bitsize (n*c) must be divisiable by the coefficient c.");
        assert(std::floor(Givaro::logtwo(max_bound)) == std::ceil(Givaro::logtwo(max_bound)) &&
               "For Parge numbers the max bitsize must be a power of two");

        uint_fast64_t n = max_bound / coefficient;
        while (product_bitsize() < product_bound)
        {
            if (n < 1)
            {
                std::cerr << "Failure to generate coprimes: We ran out of Parge numbers. Consider increasing max bitsize bound." << std::endl;
                std::cerr << " - currently generated: " << std::endl
                          << *this << std::endl;
                abort();
            }
            assert(n >= 1);
            Givaro::Integer t = 1;
            t <<= (coefficient * n);
            t++;
            this->push_back(t); // t = 2^(2^n * c) + 1
            _product *= t;
            n >>= 1; // max_bound = previous 2^n
        }
#if CHECK_MMC
        size_t N = this->size();
        for (size_t i = 0; i < N; i++)
        {
            for (size_t j = 0; j < N; j++)
            {
                if (i != j && gcd(this->val(i), this->val(j)) != 1)
                {
                    cerr << "GenPargeShift is not generating coprimes." << endl
                         << " - got: " << this->val(i) << " and " << this->val(j) << endl;
                    abort();
                }
            }
        }
#endif
    }

    ~GenPargeShift() = default;
    GenPargeShift(GenPargeShift &) = delete;
};

#endif // H_FERMAT_GEN

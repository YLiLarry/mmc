#if !defined(H_GEN_PARGE)
#define H_GEN_PARGE

#include <gmp++/gmp++.h>
#include "gen_prime.h"
#include <cmath>
#include <vector>
#include <bitset>

// generate parge numbers using the block generation scheme
class GenPargeBlock : public GenCoprimeAbstract<Givaro::Integer>
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
    GenPargeBlock(uint_fast64_t product_bound, uint_fast64_t max_bound)
    {
        assert(product_bound > max_bound && "The product of moduli must be greater than any moduli.");
        assert(max_bound > 1 && "Any moduli must be at least 2 bits");
        uint_fast64_t prev_pow_2 = 1;
        while (prev_pow_2 <= max_bound)
        {
            prev_pow_2 <<= 1;
        }
        prev_pow_2 >>= 1;
        uint_fast64_t mask = prev_pow_2 - 1; // mask = 0xffff... (contains max_bound 1s)
        int_fast64_t n = -1;                 // n = 0xffff... (contains 64 bits of 1s)
        while (product_bitsize() < product_bound)
        {
            uint_fast64_t expo = n & mask;
            std::cerr << std::bitset<64>(expo) << endl;
            if (expo <= 0)
            {
                std::cerr << "Failure to generate coprimes: We ran out of parge numbers. Consider increasing max bitsize bound." << std::endl;
                std::cerr << " - currently generated: " << std::endl
                          << *this << std::endl;
            }
            assert(expo > 0);
            Givaro::Integer t = 1;
            t <<= expo;
            t++;
            this->push_back(t); // t = 2^(2^n) + 1
            _product *= t;
            n <<= 1; // n = 0xfffffffe
        }
#if CHECK_MMC
        size_t N = this->size();
        for (size_t i = 0; i < N; i++)
        {
            for (size_t j = 0; j < N; j++)
            {
                if (i != j && gcd(this->val(i), this->val(j)) != 1)
                {
                    cerr << "GenPargeBlock is not generating coprimes." << endl
                         << " - got: " << this->val(i) << " and " << this->val(j) << endl;
                    abort();
                }
            }
        }
#endif
    }

    ~GenPargeBlock() = default;
    GenPargeBlock(GenPargeBlock &) = delete;
};

#endif // H_GEN_PARGE

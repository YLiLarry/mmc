#if !defined(H_GEN_PARGE_BLOCK)
#define H_GEN_PARGE_BLOCK

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
    Givaro::Integer _max = -1;

  public:
    inline virtual Givaro::Integer product() const override { return _product; }
    inline virtual uint_fast64_t product_bitsize() const override { return product().bitsize(); }
    inline virtual uint_fast64_t max_bitsize() const override { return max().bitsize(); }
    inline virtual const Givaro::Integer &max() const override { return _max; }

  public:
    GenPargeBlock(uint_fast64_t product_bound, uint_fast64_t max_bound, uint_fast64_t block_size)
    {
        assert(block_size > 0);
        assert(block_size < 64);
        uint_fast64_t mask = (1 << block_size) - 1; // mask = 0xffff... (contains block_size 1s)
        int_fast64_t n = -1;                        // n = 0xffff... (contains 64 bits of 1s)
        size_t i_block = 1;
        bool last_block = false;
        while (product_bitsize() < product_bound)
        {
            uint_fast64_t expo = n & mask;
            // std::cerr << "mask: " << std::bitset<64>(mask) << std::endl;
            // std::cerr << "n: " << std::bitset<64>(n) << std::endl;
            if (expo <= 0)
            {
                if (last_block)
                {
                    std::cerr << "We ran out of parge numbers, consider increasing block size or max bound." << std::endl;
                    std::cerr << "Currently generated: " << std::endl
                              << " - max_bound: " << max_bound << std::endl
                              << *this << std::endl;
                    //abort();
                    // reset 
                    block_size++;
                    mask = (1 << block_size) - 1;
                    n = -1;
                    i_block = 1;
                    last_block = false;
                    _product = 1;
                    _max = -1;
                    this->clear();
                    std::cerr << "Trying block size " << block_size << std::endl;
                    continue;
                }
                i_block++;
                mask = (1 << i_block * block_size) - 1;
                if (mask >= max_bound)
                {
                    last_block = true;
                    mask = (1 << max_bound) - 1;
                }
                continue;
            }
            std::cerr << "expo: " << std::bitset<64>(expo) << std::endl;
            assert(expo > 0);
            Givaro::Integer t = 1;
            t <<= expo;
            t++;
            this->push_back(t); // t = 2^(2^n) + 1
            _product *= t;
            if (t > _max)
            {
                _max = t;
            }
            n <<= 1; // n = 0xfffffffe
        }
#if CHECK_MMC
        std::cerr << "CHECK_MMC: checking correctness..." << std::endl;
        size_t N = this->size();
        for (size_t i = 0; i < N; i++)
        {
            for (size_t j = 0; j < N; j++)
            {
                if (i != j && gcd(this->val(i), this->val(j)) != 1)
                {
                    std::cerr << "GenPargeBlock is not generating coprimes." << std::endl
                              << " - got: " << this->val(i) << " and " << this->val(j) << std::endl;
                    abort();
                }
            }
        }
        std::cerr << "CHECK_MMC: passed." << std::endl;
#endif
    }

    ~GenPargeBlock() = default;
    GenPargeBlock(GenPargeBlock &) = delete;
};

#endif // H_GEN_PARGE_BLOCK

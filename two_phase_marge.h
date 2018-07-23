#if !defined(H_TWO_STEPS_ALGO)
#define H_TWO_STEPS_ALGO

#include "containers.h"
#include "gen_marge.h"
#include "gen_prime.h"
#include "sim_rns.h"
#include <iostream>
#include <gmp++/gmp++.h>
#include <fflas-ffpack/field/rns-double.h>
#include "two_phase_abstract.h"

namespace CNMA
{
extern "C"
{
#include "cnma/marge_num.h"
#include "cnma/reconstruct_marge.h"
}
} // namespace CNMA
// Phase 1:
// level_1_moduli_count is the number of co-primes moduli, each of bit length 2^B_F.
// co-primes are stored as T_F type in memory.
// Phase 2:
// N_M is the number of prime moduli, each of bit length B_M.
// primes are stored as T_M type in memory.
// The following relation must be true: N_M * B_M > 2^B_F
class TwoPhaseMarge : public TwoPhaseAbstract
{
  public:
    TwoPhaseMarge(uint_fast64_t level_1_product_bitsize,
                  uint_fast64_t level_1_moduli_bitsize,
                  uint_fast64_t level_2_product_bitsize,
                  uint_fast64_t level_2_moduli_bitsize)
        : TwoPhaseAbstract(new GenMargeMost(level_1_product_bitsize, level_1_moduli_bitsize),
                           new GenPrimeMost<double>(level_2_product_bitsize, level_2_moduli_bitsize))
    {
    }

    ~TwoPhaseMarge()
    {
        delete this->level_1_moduli;
        delete this->level_2_moduli;
    }

    TwoPhaseMarge(const TwoPhaseMarge &) = delete;
    TwoPhaseMarge &operator=(const TwoPhaseMarge &) = delete;

    ///////////////////////////////////////////////////////////////////////////////////////////

    /*
        this helper method is used by matrix_product(...)
    */
  protected:
    virtual const vector<Phase1_Int> matrix_reduce_phase_1(const vector<Givaro::Integer> &inputs) const override
    {
        size_t len_inputs = inputs.size();
        // phase 1 begins
        // p1_reduced stores multi-moduli representation of each input
        vector<Phase1_Int> p1_reduced(len_inputs * level_1_moduli_count);
        for (size_t i = 0; i < len_inputs; i++)
        {
            for (size_t f = 0; f < level_1_moduli_count; f++)
            {
                Phase1_Int &t = p1_reduced[i * level_1_moduli_count + f];
                t = inputs[i];
                CNMA::dc_reduce_minus(t.get_mpz(), (level_1_moduli->val(f) + 1).bitsize() - 1);
            }
#if TIME_MMC
            // print a dot for every 100 entries
            if (i % 100 == 0)
            {
                cerr << ".";
            }
#endif
        }
#if TIME_MMC
        cerr << endl;
#endif
        return p1_reduced;
    }

  protected:
    /* 
        use this method to recover from a single reduced matrix to phase 1 representations
    */
    virtual const vector<Givaro::Integer> matrix_recover_phase_1(const vector<Phase1_Int> &phase2_recovered) const override
    {
        // phase 1 recovery begins
        size_t out_len = phase2_recovered.size() / level_1_moduli_count;
        vector<Givaro::Integer> phase1_recovered(out_len);

        // initialization
        mpz_t _Mi[level_1_moduli_count]; // tmp
        mpz_t _f[level_1_moduli_count];
        mpz_t _r[level_1_moduli_count];
        for (size_t i = 0; i < level_1_moduli_count; i++)
        {
            mpz_init(_Mi[i]);
            mpz_init(_r[i]);
            mpz_init_set(_f[i], level_1_moduli->val(i).get_mpz());
        }
        CNMA::precompute_Mi(_Mi, _f, level_1_moduli_count);
        // recover
        for (size_t i = 0; i < out_len; i++)
        {
            Givaro::Integer &t = phase1_recovered[i];
            for (size_t f = 0; f < level_1_moduli_count; f++)
            {
                const Phase1_Int &in = phase2_recovered[i * level_1_moduli_count + f];
                if (in >= level_1_moduli->product())
                {
                    cerr << "Computation overflows. Recovered an integer that is greater than the product of level 1 moduli." << endl;
                }
                assert(in < level_1_moduli->product());
                mpz_mod(_r[f], in.get_mpz(), level_1_moduli->val(f).get_mpz());
            }
            CNMA::garner(t.get_mpz(), level_1_moduli_count, _r, _f, _Mi);
            mpz_mod(t.get_mpz(), t.get_mpz(), level_1_moduli->product().get_mpz());
#if TIME_MMC
            if (i % 100 == 0)
            {
                cerr << ".";
            }
#endif
        }
        for (size_t i = 0; i < level_1_moduli_count; i++)
        {
            mpz_clear(_Mi[i]);
            mpz_clear(_r[i]);
            mpz_clear(_f[i]);
        }
#if TIME_MMC
        cerr << endl;
#endif
        return phase1_recovered;
    }
};

#endif // H_TWO_STEPS_ALGO

#ifndef H_TWO_PHASE_PARGE_SHIFT
#define H_TWO_PHASE_PARGE_SHIFT

#include "two_phase_parge_abstract.h"
#include "gen_parge_shift.h"

class TwoPhasePargeShift : public TwoPhasePargeAbstract
{
    uint_fast64_t _level_1_moduli_bitsize_coefficient;

  public:
    TwoPhasePargeShift(uint_fast64_t level_1_product_bitsize,
                       uint_fast64_t level_1_moduli_bitsize,
                       uint_fast64_t level_1_moduli_bitsize_coefficient,
                       uint_fast64_t level_2_moduli_bitsize)
        : TwoPhasePargeAbstract(new GenPargeShift(level_1_product_bitsize, level_1_moduli_bitsize, level_1_moduli_bitsize_coefficient),
                           new GenPrimeMost<double>(level_1_moduli_bitsize, level_2_moduli_bitsize))
                           , _level_1_moduli_bitsize_coefficient(level_1_moduli_bitsize_coefficient)
    {
    }

    ~TwoPhasePargeShift()
    {
        delete this->level_1_moduli;
        delete this->level_2_moduli;
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
        mpz_t _f[level_1_moduli_count];
        mpz_t _r[level_1_moduli_count];
        for (size_t i = 0; i < level_1_moduli_count; i++)
        {
            mpz_init(_r[i]);
            mpz_init_set(_f[i], level_1_moduli->val(i).get_mpz());
        }
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
            CNMA::garner_parge(t.get_mpz(), level_1_moduli_count, _r, _f, _level_1_moduli_bitsize_coefficient);
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
            mpz_clear(_r[i]);
            mpz_clear(_f[i]);
        }
#if TIME_MMC
        cerr << endl;
#endif
        return phase1_recovered;
    }
};

#endif

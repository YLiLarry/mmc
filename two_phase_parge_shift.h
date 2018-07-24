#ifndef H_TWO_PHASE_PARGE_SHIFT
#define H_TWO_PHASE_PARGE_SHIFT

#include "two_phase_parge_abstract.h"
#include "gen_parge_shift.h"

namespace CNMA
{
extern "C"
{
#include "cnma/parge_num.h"
#include "cnma/reconstruct_parge_shift.h"
}
} // namespace CNMA

class TwoPhasePargeShift : public TwoPhasePargeAbstract
{
    uint_fast64_t m_level_1_moduli_bitsize_coefficient;

  public:
    TwoPhasePargeShift(uint_fast64_t level_1_product_bitsize,
                       uint_fast64_t level_1_moduli_bitsize,
                       uint_fast64_t level_1_moduli_bitsize_coefficient,
                       uint_fast64_t level_2_moduli_bitsize)
        : TwoPhasePargeAbstract(new GenPargeShift(level_1_product_bitsize, level_1_moduli_bitsize, level_1_moduli_bitsize_coefficient),
                                new GenPrimeMost<double>(level_1_moduli_bitsize, level_2_moduli_bitsize)),
          m_level_1_moduli_bitsize_coefficient(level_1_moduli_bitsize_coefficient)
    {
    }

  protected:
    /* 
        use this method to recover from a single reduced matrix to phase 1 representations
    */
    virtual const vector<Givaro::Integer> matrix_recover_phase_1(const vector<Phase1_Int> &phase2_recovered) const override
    {
        // phase 1 recovery begins
        size_t out_len = phase2_recovered.size() / m_level_1_moduli_count;
        vector<Givaro::Integer> phase1_recovered(out_len);

        // initialization
        mpz_t _f[m_level_1_moduli_count];
        mpz_t _r[m_level_1_moduli_count];
        for (size_t i = 0; i < m_level_1_moduli_count; i++)
        {
            mpz_init(_r[i]);
            mpz_init_set(_f[i], m_level_1_moduli->val(i).get_mpz());
        }
        // recover
        for (size_t i = 0; i < out_len; i++)
        {
            Givaro::Integer &t = phase1_recovered[i];
            for (size_t f = 0; f < m_level_1_moduli_count; f++)
            {
                const Phase1_Int &in = phase2_recovered[i * m_level_1_moduli_count + f];
                mpz_mod(_r[f], in.get_mpz(), m_level_1_moduli->val(f).get_mpz());
            }
            CNMA::garner_parge_shift(t.get_mpz(), m_level_1_moduli_count, _r, _f, m_level_1_moduli_bitsize_coefficient);
            mpz_mod(t.get_mpz(), t.get_mpz(), m_level_1_moduli->product().get_mpz());
#if TIME_MMC
            if (i % 100 == 0)
            {
                cerr << ".";
            }
#endif
        }
        for (size_t i = 0; i < m_level_1_moduli_count; i++)
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

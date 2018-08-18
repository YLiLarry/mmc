#ifndef H_TWO_PHASE_PARGE_BLOCK
#define H_TWO_PHASE_PARGE_BLOCK

#include "two_phase_parge_abstract.h"
#include "gen_parge_block.h"

#include "cnma/parge_num.h"
#include "cnma/reconstruct_parge_block.h"


class TwoPhasePargeBlock : public TwoPhasePargeAbstract
{
  public:
    TwoPhasePargeBlock(uint_fast64_t level_1_product_bitsize,
                       uint_fast64_t level_1_moduli_bitsize,
                       uint_fast64_t block_size)
        : TwoPhasePargeAbstract(new GenPargeBlock(level_1_product_bitsize, level_1_moduli_bitsize, block_size), NULL)
    {
        assert(level_1_product_bitsize > level_1_moduli_bitsize * 2 && "Level 1 moduli size cannot be too large for the two-phase algorithm to be beneficial.");
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
        mpz_t input_Mi[m_level_1_moduli_count];
        mpz_t input_f[m_level_1_moduli_count];
        uint64_t input_f_expo[m_level_1_moduli_count];
        mpz_t input_work[m_level_1_moduli_count];
        mpz_t input_r[m_level_1_moduli_count];
        for (size_t i = 0; i < m_level_1_moduli_count; i++)
        {
            mpz_init(input_Mi[i]);
            mpz_init(input_r[i]);
            mpz_init(input_work[i]);
            mpz_init_set(input_f[i], m_level_1_moduli->val(i).get_mpz());
            input_f_expo[i] = (m_level_1_moduli->val(i) - 1).bitsize() - 1;
        }
        CNMA::precompute_Mi_parge_block(input_Mi, input_f, m_level_1_moduli_count);
        // recover
        for (size_t i = 0; i < out_len; i++)
        {
            Givaro::Integer &t = phase1_recovered[i];
            for (size_t f = 0; f < m_level_1_moduli_count; f++)
            {
                const Phase1_Int &in = phase2_recovered[i * m_level_1_moduli_count + f];
                mpz_set(input_r[f], in.get_mpz());
                CNMA::dc_reduce_plus(input_r[f], (m_level_1_moduli->val(f) - 1).bitsize() - 1);
            }
            CNMA::garner_parge_block(t.get_mpz(), m_level_1_moduli_count, input_r, input_f_expo, input_f, input_Mi, input_work);
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
            mpz_clear(input_Mi[i]);
            mpz_clear(input_r[i]);
            mpz_clear(input_f[i]);
            mpz_clear(input_work[i]);
        }
#if TIME_MMC
        cerr << endl;
#endif
        return phase1_recovered;
    }
};

#endif

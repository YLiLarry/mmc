#ifndef H_TWO_PHASE_PARGE_SHIFT
#define H_TWO_PHASE_PARGE_SHIFT

#include "two_phase_parge_abstract.h"
#include "gen_parge_shift.h"

#include "cnma/parge_num.h"
#include "cnma/reconstruct_parge_block.h"
#include "cnma/reconstruct_parge_shift.h"

class TwoPhasePargeShift : public TwoPhasePargeAbstract
{
    uint_fast64_t m_level_1_moduli_bitsize_coefficient;

  public:
    TwoPhasePargeShift(uint_fast64_t level_1_product_bitsize,
                       uint_fast64_t level_1_moduli_bitsize,
                       uint_fast64_t level_1_moduli_bitsize_coefficient)
        : TwoPhasePargeAbstract(new GenPargeShift(level_1_product_bitsize, level_1_moduli_bitsize, level_1_moduli_bitsize_coefficient), NULL),
          m_level_1_moduli_bitsize_coefficient(level_1_moduli_bitsize_coefficient)
    {
        assert(level_1_product_bitsize > level_1_moduli_bitsize * 2 && "Level 1 moduli size cannot be too large for the two-phase algorithm to be beneficial.");
    }


    /*
        this helper method is used by matrix_product(...)
    */
  protected:
    virtual const vector<Phase1_Int> matrix_reduce_phase_1(const vector<Givaro::Integer> &inputs) const override
    {
        size_t len_inputs = inputs.size();
        // phase 1 begins
        // p1_reduced stores multi-moduli representation of each input
        vector<Phase1_Int> p1_reduced(len_inputs * m_level_1_moduli_count);
        for (size_t i = 0; i < len_inputs; i++)
        {
            // first moduli is 2^n
            Phase1_Int &t0 = p1_reduced[i * m_level_1_moduli_count + 0];
            t0 = inputs[i];
            t0 %= m_level_1_moduli->val(0); // could be slightly better here!
            // second moduli is 2^n+3
            Phase1_Int &t1 = p1_reduced[i * m_level_1_moduli_count + 1];
            t1 = inputs[i];
            t1 %= m_level_1_moduli->val(1);
            // third moduli is a random prime
            Phase1_Int &t2 = p1_reduced[i * m_level_1_moduli_count + 2];
            t2 = inputs[i];
            t2 %= m_level_1_moduli->val(2);
            // rest moduli are 2^i+1
            for (size_t f = 3; f < m_level_1_moduli_count; f++)
            {
                Phase1_Int &t = p1_reduced[i * m_level_1_moduli_count + f];
                t = inputs[i];
                CNMA::dc_reduce_plus(t.get_mpz(), (m_level_1_moduli->val(f) - 1).bitsize() - 1);
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
            input_f_expo[i] = m_level_1_moduli->val(i).bitsize() - 1;
        }
        CNMA::precompute_Mi_parge_block(input_Mi, input_f, m_level_1_moduli_count);
        // recover
        for (size_t i = 0; i < out_len; i++)
        {
            // first moduli is 2^n
            const Phase1_Int &in0 = phase2_recovered[i * m_level_1_moduli_count + 0];
            mpz_mod(input_r[0], in0.get_mpz(), m_level_1_moduli->val(0).get_mpz());
            // second moduli is 2^n + 3
            const Phase1_Int &in1 = phase2_recovered[i * m_level_1_moduli_count + 1];
            mpz_mod(input_r[1], in1.get_mpz(), m_level_1_moduli->val(1).get_mpz());
            // third moduli is a random prime
            const Phase1_Int &in2 = phase2_recovered[i * m_level_1_moduli_count + 2];
            mpz_mod(input_r[2], in2.get_mpz(), m_level_1_moduli->val(2).get_mpz());
            // rest moduli are 2^i + 1
            for (size_t f = 3; f < m_level_1_moduli_count; f++)
            {
                const Phase1_Int &in = phase2_recovered[i * m_level_1_moduli_count + f];
                mpz_set(input_r[f], in.get_mpz());
                CNMA::dc_reduce_plus(input_r[f], m_level_1_moduli->val(f).bitsize() - 1);
            }
            Givaro::Integer &t = phase1_recovered[i];
            CNMA::garner_parge_shift_mixed(t.get_mpz(), m_level_1_moduli_count, input_r, input_f_expo, input_f, input_Mi, m_level_1_moduli_bitsize_coefficient, input_work);
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

#if !defined(H_TWO_PHASE_MARGE_ABSTRACT)
#define H_TWO_PHASE_MARGE_ABSTRACT

#include "containers.h"
#include "gen_coprime_abstract.h"
#include "sim_rns.h"
#include <iostream>
#include <gmp++/gmp++.h>
#include <fflas-ffpack/field/rns-double.h>
#include "two_phase_abstract.h"

#include "cnma/marge_num.h"
#include "cnma/reconstruct_marge.h"

// Phase 1:
// m_level_1_moduli_count is the number of co-primes moduli, each of bit length 2^B_F.
// co-primes are stored as T_F type in memory.
// Phase 2:
// N_M is the number of prime moduli, each of bit length B_M.
// primes are stored as T_M type in memory.
// The following relation must be true: N_M * B_M > 2^B_F
class TwoPhaseMargeAbstract : public TwoPhaseAbstract
{
  public:
    TwoPhaseMargeAbstract(const GenCoprimeAbstract<Givaro::Integer> *m_level_1_moduli,
                          const GenCoprimeAbstract<double> *m_level_2_moduli)
        : TwoPhaseAbstract(m_level_1_moduli, m_level_2_moduli) {}
    TwoPhaseMargeAbstract(const TwoPhaseMargeAbstract &) = delete;
    TwoPhaseMargeAbstract &operator=(const TwoPhaseMargeAbstract &) = delete;

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
        vector<Phase1_Int> p1_reduced(len_inputs * m_level_1_moduli_count);
        for (size_t i = 0; i < len_inputs; i++)
        {
            for (size_t f = 0; f < m_level_1_moduli_count; f++)
            {
                Phase1_Int &t = p1_reduced[i * m_level_1_moduli_count + f];
                t = inputs[i];
                CNMA::dc_reduce_minus(t.get_mpz(), (m_level_1_moduli->val(f) + 1).bitsize() - 1);
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
        mpz_t input_Mi[m_level_1_moduli_count]; // tmp
        mpz_t input_f[m_level_1_moduli_count];
        mpz_t input_r[m_level_1_moduli_count];
        mpz_t input_work[m_level_1_moduli_count];
        uint64_t input_f_expo[m_level_1_moduli_count];
        for (size_t i = 0; i < m_level_1_moduli_count; i++)
        {
            mpz_init(input_Mi[i]);
            mpz_init(input_r[i]);
            mpz_init(input_work[i]);
            mpz_init_set(input_f[i], m_level_1_moduli->val(i).get_mpz());
            input_f_expo[i] = (m_level_1_moduli->val(i) + 1).bitsize() - 1;
        }
        CNMA::precompute_Mi_marge(input_Mi, input_f, m_level_1_moduli_count);
        // recover
        for (size_t i = 0; i < out_len; i++)
        {
            for (size_t f = 0; f < m_level_1_moduli_count; f++)
            {
                const Phase1_Int &in = phase2_recovered[i * m_level_1_moduli_count + f];
                mpz_mod(input_r[f], in.get_mpz(), m_level_1_moduli->val(f).get_mpz());
                // mpz_set(input_r[f], in.get_mpz());
                // CNMA::dc_reduce_minus(input_r[f], (m_level_1_moduli->val(f) + 1).bitsize() - 1);
            }
            Givaro::Integer &t = phase1_recovered[i];
            CNMA::garner_marge(t.get_mpz(), m_level_1_moduli_count, input_r, input_f_expo, input_f, input_Mi, input_work);
            // CNMA::garner_simple_marge(t.get_mpz(), m_level_1_moduli_count, input_r, input_f, input_Mi);
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

#endif // H_TWO_PHASE_MARGE_ABSTRACT

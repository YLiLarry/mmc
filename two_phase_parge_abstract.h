#ifndef H_TWO_PHASE_PARGE_ABSTRACT
#define H_TWO_PHASE_PARGE_ABSTRACT

#include "two_phase_abstract.h"
#include <vector>
#include <gmp++/gmp++.h>

namespace CNMA
{
extern "C"
{
#include "cnma/parge_num.h"
}
} // namespace CNMA

class TwoPhasePargeAbstract : public TwoPhaseAbstract
{

  public:
    TwoPhasePargeAbstract(const GenCoprimeAbstract<Givaro::Integer> *m_level_1_moduli,
                          const GenCoprimeAbstract<double> *m_level_2_moduli)
        : TwoPhaseAbstract(m_level_1_moduli, m_level_2_moduli) {}
    TwoPhasePargeAbstract(const TwoPhasePargeAbstract &) = delete;
    TwoPhasePargeAbstract &operator=(const TwoPhasePargeAbstract &) = delete;

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
};

#endif
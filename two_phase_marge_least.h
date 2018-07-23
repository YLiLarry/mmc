#if !defined(H_TWO_PHASE_MARGE_LEAST)
#define H_TWO_PHASE_MARGE_LEAST

#include "containers.h"
#include "gen_marge_least.h"
#include "gen_prime.h"
#include "sim_rns.h"
#include <iostream>
#include <gmp++/gmp++.h>
#include <fflas-ffpack/field/rns-double.h>
#include "two_phase_marge_abstract.h"

namespace CNMA
{
extern "C"
{
#include "cnma/marge_num.h"
#include "cnma/reconstruct_marge.h"
}
} // namespace CNMA

// Phase 1:
// m_level_1_moduli_count is the number of co-primes moduli, each of bit length 2^B_F.
// co-primes are stored as T_F type in memory.
// Phase 2:
// N_M is the number of prime moduli, each of bit length B_M.
// primes are stored as T_M type in memory.
// The following relation must be true: N_M * B_M > 2^B_F
class TwoPhaseMargeLeast : public TwoPhaseMargeAbstract
{
  public:
    TwoPhaseMargeLeast(uint_fast64_t level_1_product_bitsize,
                  uint_fast64_t level_1_moduli_bitsize,
                  uint_fast64_t level_2_product_bitsize,
                  uint_fast64_t level_2_moduli_bitsize)
        : TwoPhaseMargeAbstract(new GenMargeLeast(level_1_product_bitsize, level_1_moduli_bitsize),
                           new GenPrimeMost<double>(level_2_product_bitsize, level_2_moduli_bitsize))
    {
    }

    ~TwoPhaseMargeLeast()
    {
        delete this->m_level_1_moduli;
        delete this->m_level_2_moduli;
    }

    TwoPhaseMargeLeast(const TwoPhaseMargeLeast &) = delete;
    TwoPhaseMargeLeast &operator=(const TwoPhaseMargeLeast &) = delete;
};

#endif // H_TWO_PHASE_MARGE_LEAST

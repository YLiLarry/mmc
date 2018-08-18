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

#include "cnma/marge_num.h"
#include "cnma/reconstruct_marge.h"


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
                       uint_fast64_t level_1_moduli_bitsize)
        : TwoPhaseMargeAbstract(new GenMargeLeast(level_1_product_bitsize, level_1_moduli_bitsize), NULL)
    {
      assert(level_1_product_bitsize > level_1_moduli_bitsize * 2 && "Level 1 moduli size cannot be too large for the two-phase algorithm to be beneficial.");
    }

    TwoPhaseMargeLeast(const TwoPhaseMargeLeast &) = delete;
    TwoPhaseMargeLeast &operator=(const TwoPhaseMargeLeast &) = delete;
};

#endif // H_TWO_PHASE_MARGE_LEAST

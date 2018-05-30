#if !defined(H_TWO_STEPS_ALGO)
#define H_TWO_STEPS_ALGO

#include <givaro/modular-integer.h>
#include <linbox/integer.h>
#include <iostream>
#include "cnma/marge_num.h"
#include "coprime_gen.h"
#include "prime_gen.h"

// Phase 1:
// N_F is the number of co-primes moduli, each of bit length 2^B_F.
// co-primes are stored as T_F type in memory.
// Phase 2:
// N_M is the number of prime moduli, each of bit length B_M.
// primes are stored as T_M type in memory.
// The following relation must be true: N_M * B_M > 2^B_F
template <typename T_F, size_t N_F, uint_fast8_t B_F, typename T_M, size_t N_M,
          uint_fast8_t B_M>
class TwoPhaseAlgo {
   private:
    MargeGenExact<T_F, N_F, B_F> _F_arr;  // co-primes in the first level
    PrimeGenExact<T_M, N_M, B_M> _M_arr;  // primes in the second level

   public:
    TwoPhaseAlgo() { assert(N_M * B_M > (1 << B_F)); }
    ~TwoPhaseAlgo() = default;
    TwoPhaseAlgo(const TwoPhaseAlgo&) = delete;
    TwoPhaseAlgo& operator=(const TwoPhaseAlgo&) = delete;

    void mult(T_F& dest, T_F& a, T_F& b) {
        assert(_F_arr.product() > a);
        assert(_F_arr.product() > b);
        // stores multimuduli representation of a in phase 1
        array<T_F, N_F> MM_A_1;
        // stores multimuduli representation of b in phase 1
        array<T_F, N_F> MM_B_1;
        // phase 1 begins
        for (size_t i = 0; i < N_F; i++) {
            MM_A_1[i] += a;  // init to a
            MM_B_1[i] += b;  // init to b
            dc_reduce_minus(MM_A_1[i], _F_arr[i]);
            dc_reduce_minus(MM_B_1[i], _F_arr[i]);
        }
        cout << MM_A_1 << endl;
        cout << MM_B_1 << endl;

        // phase 2 begins
        for (size_t i = 0; i < N_M; i++) {
        }

        // stores multimuduli representation of a in phase 2
        array<T_M, N_M> MM_A_2;
        // stores multimuduli representation of b in phase 2
        array<T_M, N_M> MM_B_2;
    }
};

#endif  // H_TWO_STEPS_ALGO

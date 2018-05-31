#if !defined(H_TWO_STEPS_ALGO)
#define H_TWO_STEPS_ALGO

#include "./sim_rns.h"
#include "cnma/marge_num.h"
#include "coprime_gen.h"
#include "prime_gen.h"
#include <givaro/modular-integer.h>
#include <iostream>
#include <linbox/integer.h>

// Phase 1:
// N_F is the number of co-primes moduli, each of bit length 2^B_E_F.
// co-primes are stored as T_F type in memory.
// Phase 2:
// N_M is the number of prime moduli, each of bit length B_M.
// primes are stored as T_M type in memory.
// The following relation must be true: N_M * B_M > 2^B_E_F
template <typename T_F, size_t N_F, uint_fast8_t B_E_F, typename T_M, size_t N_M, uint_fast64_t B_M>
class TwoPhaseAlgo {
private:
    const PrimeGenExact<T_F, N_F, B_E_F> _F_arr; // co-primes in the first level
    const PrimeGenExact<T_M, N_M, B_M> _M_arr; // primes in the second level

public:
    TwoPhaseAlgo() { assert((1 << B_E_F) < static_cast<uint_fast64_t>(_M_arr.product)); }
    ~TwoPhaseAlgo() = default;
    TwoPhaseAlgo(const TwoPhaseAlgo&) = delete;
    TwoPhaseAlgo& operator=(const TwoPhaseAlgo&) = delete;

    typedef SIM_RNS::RNS<T_F, T_F, N_F> Phase1_RNS;
    typedef SIM_RNS::RNS<T_F, T_M, N_M> Phase2_RNS;

    void mult(T_F& dest, T_F& a, T_F& b)
    {
        assert(a.bitsize() < static_cast<uint_fast64_t>(_F_arr.sum));
        assert(b.bitsize() < static_cast<uint_fast64_t>(_F_arr.sum));
        // phase 1 begins
        // MM_A_1 stores multimuduli representation of a
        typename Phase1_RNS::ReducedInt
            MM_A_1(a, _F_arr, [&](size_t i) {
                T_F* t = new T_F;
                t->expensiveCopy(a);
                dc_reduce_minus(*t, _F_arr[i]);
                return t;
            });
        // MM_B_1 stores multimuduli representation of b
        typename Phase1_RNS::ReducedInt
            MM_B_1(a, _F_arr, [&](size_t i) {
                T_F* t = new T_F;
                t->expensiveCopy(b);
                dc_reduce_minus(*t, _F_arr[i]);
                return t;
            });
        cout << MM_A_1 << endl;
        cout << MM_B_1 << endl;
        // phase 2 begins
        size_t input_size = 2 * N_F;
        vector<const T_F*> inputs(input_size);
        for (size_t i = 0; i < N_F; i++) {
            inputs[i] = &MM_A_1.residuals[i];
        }
        for (size_t i = 0; i < N_F; i++) {
            inputs[i + N_F] = &MM_B_1.residuals[i];
        }
        typename Phase2_RNS::ReducedIntVector reduced_all = Phase2_RNS::naive_reduce(inputs, _M_arr);
        cout << reduced_all << endl;
    }
};

#endif // H_TWO_STEPS_ALGO

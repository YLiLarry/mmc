#if !defined(H_TWO_STEPS_ALGO)
#define H_TWO_STEPS_ALGO

#include "cnma/marge_num.h"
#include "containers.h"
#include "coprime_gen.h"
#include "prime_gen.h"
#include "sim_rns.h"
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
    TwoPhaseAlgo() { assert((1 << B_E_F) < static_cast<uint_fast64_t>(_M_arr.sum)); }
    ~TwoPhaseAlgo() = default;
    TwoPhaseAlgo(const TwoPhaseAlgo&) = delete;
    TwoPhaseAlgo& operator=(const TwoPhaseAlgo&) = delete;

    void mult(T_F& dest, T_F& a, T_F& b)
    {
        cerr << "_F_arr.sum: " << _F_arr.sum << endl;
        cerr << "a.bitsize(): " << a.bitsize() << endl;
        assert(a.bitsize() < static_cast<uint_fast64_t>(_F_arr.sum));
        assert(b.bitsize() < static_cast<uint_fast64_t>(_F_arr.sum));
        typedef SIM_RNS::RNS<T_F, T_F, N_F> Phase1_RNS;
        cerr << "phase 1" << endl;
        // phase 1 begins
        // MM_A_1 stores multimuduli representation of a
        typename Phase1_RNS::ReducedInt
            MM_A_1(_F_arr, [&](size_t i) {
                T_F* t = new T_F;
                cerr << ".";
                t->EXPENSIVE_COPY(a);
                dc_reduce_minus(*t, _F_arr[i]);
                return t;
            });
        cerr << endl;
        // MM_B_1 stores multimuduli representation of b
        typename Phase1_RNS::ReducedInt
            MM_B_1(_F_arr, [&](size_t i) {
                T_F* t = new T_F;
                cerr << ".";
                t->EXPENSIVE_COPY(b);
                dc_reduce_minus(*t, _F_arr[i]);
                return t;
            });
        cerr << endl;
        // cout << MM_A_1 << endl;
        // cout << MM_B_1 << endl;
        cerr << "phase 2" << endl;
        typedef SIM_RNS::RNS<T_F, T_M, N_M> Phase2_RNS;
        // phase 2 begins
        vector<const T_F*> phase2_inputs(2 * N_F);
        for (size_t i = 0; i < N_F; i++) {
            phase2_inputs[i] = MM_A_1.residuals.ptr(i);
        }
        for (size_t i = 0; i < N_F; i++) {
            phase2_inputs[i + N_F] = MM_B_1.residuals.ptr(i);
        }
        typename Phase2_RNS::ReducedIntVector phase2_outputs = Phase2_RNS::sim_reduce(phase2_inputs, _M_arr);
        // A = A*B
        for (size_t i = 0; i < N_F; i++) {
            typename Phase2_RNS::ReducedInt* a = phase2_outputs[i];
            typename Phase2_RNS::ReducedInt* b = phase2_outputs[i + N_F];
            a->operator*=(*b);
        }
        phase2_outputs.erase(phase2_outputs.begin() + N_F);
        phase2_outputs.resize(N_F);
        // cout << phase2_outputs << endl;
        phase2_outputs.erase(phase2_outputs.begin());

        cerr << "phase 2 recovery" << endl;
        // phase 2 recovery begins
        vector<T_F*> recovered = Phase2_RNS::sim_recover(phase2_outputs);
        cerr << recovered << endl;
        recovered.erase(recovered.begin());
    }
};

#endif // H_TWO_STEPS_ALGO

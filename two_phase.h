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
// N_F is the number of co-primes moduli, each of bit length 2^B_F.
// co-primes are stored as T_F type in memory.
// Phase 2:
// N_M is the number of prime moduli, each of bit length B_M.
// primes are stored as T_M type in memory.
// The following relation must be true: N_M * B_M > 2^B_F
template <typename T_F, size_t N_F, uint_fast64_t B_F, typename T_M, size_t N_M, uint_fast64_t B_M>
class TwoPhaseAlgo {
private:
    const MargeGenMost<T_F, N_F, B_F> _F_arr; // co-primes in the first level
    const PrimeGenExact<T_M, N_M, B_M> _M_arr; // primes in the second level

public:
    TwoPhaseAlgo()
    {
        cerr << _F_arr << endl
             << _M_arr << endl;
        assert(B_F < _M_arr.product.bitsize() && "Your second level moduli's product is less than the input.");
    }
    ~TwoPhaseAlgo() = default;
    TwoPhaseAlgo(const TwoPhaseAlgo&) = delete;
    TwoPhaseAlgo& operator=(const TwoPhaseAlgo&) = delete;

    void mult(T_F& dest, T_F& a, T_F& b)
    {
#ifdef DEBUG_MMC
        cerr << "mult(a,b) called with inputs:" << endl
             << " - a: " << a << endl
             << " - b: " << b << endl;
#endif
        assert(static_cast<uint_least64_t>(a.bitsize()) < _F_arr.sum);
        assert(static_cast<uint_least64_t>(b.bitsize()) < _F_arr.sum);
#ifdef DEBUG_MMC
        cerr << "##### phase 1 #####" << endl;
#endif
        // phase 1 begins
        // MM_A_1 stores multi-moduli representation of a
        SIM_RNS::ReducedInt<T_F, N_F> MM_A_1(_F_arr, [&](size_t i) {
            T_F* t = new T_F;
#ifdef TIME_MMC
            cerr << ".";
#endif
            t->EXPENSIVE_COPY(a);
            dc_reduce_minus(t->get_mpz(), (_F_arr[i] + 1).bitsize() - 1);
            return t;
        });
#ifdef TIME_MMC
        cerr << endl;
#endif
        // MM_B_1 stores multi-moduli representation of b
        SIM_RNS::ReducedInt<T_F, N_F> MM_B_1(_F_arr, [&](size_t i) {
            T_F* t = new T_F;
#ifdef TIME_MMC
            cerr << ".";
#endif
            t->EXPENSIVE_COPY(b);
            dc_reduce_minus(t->get_mpz(), (_F_arr[i] + 1).bitsize() - 1);
            return t;
        });
#ifdef TIME_MMC
        cerr << endl;
#endif

#ifdef DEBUG_MMC
        cerr << "phase 1 reduced: " << endl
             << "MM_A_1 " << MM_A_1 << endl
             << "MM_B_1 " << MM_B_1 << endl;
#endif

#ifdef DEBUG_MMC
        cerr << "##### phase 2 #####" << endl;
#endif

        // phase 2 begins
        NumPtrVector<T_F> phase2_inputs(2 * N_F);
        for (size_t i = 0; i < N_F; i++) {
            phase2_inputs.ptr(i) = MM_A_1.residuals.ptr(i);
        }
        for (size_t i = 0; i < N_F; i++) {
            phase2_inputs.ptr(i + N_F) = MM_B_1.residuals.ptr(i);
        }
        NumPtrVector<SIM_RNS::ReducedInt<T_M, N_M>>* phase2_outputs = SIM_RNS::new_sim_reduce<T_F, T_M, N_M>(phase2_inputs, _M_arr);

#ifdef DEBUG_MMC
        cerr << "phase 2 reduced: " << endl
             << *phase2_outputs << endl;
#endif
        // A = A*B
        for (size_t i = 0; i < N_F; i++) {
            SIM_RNS::ReducedInt<T_M, N_M>* a = phase2_outputs->ptr(i);
            SIM_RNS::ReducedInt<T_M, N_M>* b = phase2_outputs->ptr(i + N_F);
            a->operator*=(*b);
        }
        phase2_outputs->erase(phase2_outputs->begin() + N_F);
        phase2_outputs->resize(N_F);

        // cerr << "a * b: " << phase2_outputs << endl;
#ifdef DEBUG_MMC
        cerr << "phase 2 recovery" << endl;
#endif
        // phase 2 recovery begins
        NumPtrVector<T_F>* recovered = SIM_RNS::new_sim_recover<T_F, T_M, N_M>(*phase2_outputs);
#ifdef DEBUG_MMC
        cerr << *recovered << endl;
#endif
        delete phase2_outputs;
        delete recovered;
    }
};

#endif // H_TWO_STEPS_ALGO

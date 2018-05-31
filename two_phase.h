#if !defined(H_TWO_STEPS_ALGO)
#define H_TWO_STEPS_ALGO

#include "cnma/marge_num.h"
#include "coprime_gen.h"
#include "prime_gen.h"
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
template <typename T_F, size_t N_F, uint_fast64_t B_F, typename T_M, size_t N_M,
    uint_fast64_t B_M>
class TwoPhaseAlgo {
private:
    PrimeGenExact<T_F, N_F, B_F> _F_arr; // co-primes in the first level
    PrimeGenExact<T_M, N_M, B_M> _M_arr; // primes in the second level

public:
    TwoPhaseAlgo() { assert(B_F < static_cast<uint_fast64_t>(_M_arr.sum)); }
    ~TwoPhaseAlgo() = default;
    TwoPhaseAlgo(const TwoPhaseAlgo&) = delete;
    TwoPhaseAlgo& operator=(const TwoPhaseAlgo&) = delete;

    void mult(T_F& dest, T_F& a, T_F& b)
    {
        assert(a.bitsize() < static_cast<uint_fast64_t>(_F_arr.sum));
        assert(b.bitsize() < static_cast<uint_fast64_t>(_F_arr.sum));
        // phase 1 begins
        // MM_A_1 stores multimuduli representation of a
        ConstNumPtrArray<T_F, N_F> MM_A_1([&](size_t i) {
            T_F* t = new T_F;
            t->expensiveCopy(a);
            dc_reduce_minus(*t, _F_arr[i]);
            return t;
        });
        // MM_B_1 stores multimuduli representation of b
        ConstNumPtrArray<T_F, N_F> MM_B_1([&](size_t i) {
            T_F* t = new T_F;
            t->expensiveCopy(b);
            dc_reduce_minus(*t, _F_arr[i]);
            return t;
        });
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

#endif // H_TWO_STEPS_ALGO

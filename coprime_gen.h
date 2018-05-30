#if !defined(H_COPRIME_GEN)
#define H_COPRIME_GEN

#include "givaro/modular-integer.h"
#include "prime_gen.h"
#include <cmath>

// generates N co-primes of form 2^n - 1
// where each n is exactly 2^B bits long,
// ie, all co-primes are between [2^(2^(B - 1)), 2^(2^B - 1)]
template <typename T, size_t N, uint_fast8_t B>
class MargeGenExact : public PrimeGen<T, N, B> {
public:
    MargeGenExact()
        : PrimeGen<T, N, B>([this](size_t index) {
            const uint64_t p = PrimeGen<T, N, B>::_primeItr.random_exact();
            assert(p < 3e10);
            T* t = new T{ 1 }; // init to 1
            (*t) <<= p; // = 2^p
            (*t)--; // = 2^p-1
            return t;
        })
    {
    }
    ~MargeGenExact() = default;
    MargeGenExact(MargeGenExact&) = delete;
};

// template <typename T, size_t N, uint8_t B>
// class PargeGenMost : public PrimeGen<T, N, floor(log2(B)))> {
//    public:
//     PargeGenMost() {
//         for (int i = 0; i < N; i++) {
//             T p = PrimeGen<T, N, B>::_primeItr.random();
//             PrimeGen<T, N, B>::_primes[i] = pow(2, p);
//         }
//     };
//     ~PargeGenMost() = default;
//     PargeGenMost(PargeGenMost&) = delete;
// };

#endif // H_COPRIME_GEN

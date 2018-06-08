#if !defined(H_COPRIME_GEN)
#define H_COPRIME_GEN

#include "givaro/modular-integer.h"
#include "prime_gen.h"
#include <cmath>

// generates N co-primes of form 2^n - 1
// where each n is exactly 2^B bits long,
// ie, all co-primes are between [2^(2^(B - 1)), 2^(2^B - 1)]
template <typename T, size_t N, uint_fast64_t B>
class MargeGenMost : public PrimeGen<T, N, B> {
protected:
    static LinBox::RandomPrimeIter _primeItr;

public:
    MargeGenMost()
        : PrimeGen<T, N, B>([&](size_t index) {
            const uint_fast64_t p = MargeGenMost<T, N, B>::_primeItr.random();
            T* t = new T{ 1 }; // init to 1
            t->operator<<=(p); // = 2^p
            t->operator--(); // = 2^p-1
            return t;
        })
    {
    }

    friend ostream& operator<<(ostream& out, const MargeGenMost<LInteger, N, B>& arr) {
        out << "MargeGenMost [";
        for (size_t i = 0; i < N; i++) {
            const Integer& t = arr[i];
            const uint_fast64_t p = (t+1).bitsize() - 1; // static_cast<uint_fast64_t>(std::ceil())
            out << "2^" << p << "-1=" << t;
            if (i != N - 1) {
                out << " , ";
            }
        }
        out << "]" << endl
             << " - count: " << arr.count() << endl
             << " - max bit length " << B << endl
             << " - product bit length: " << arr.product.bitsize() << endl;
        return out;
    }

    ~MargeGenMost() = default;
    MargeGenMost(MargeGenMost&) = delete;
};

template <typename T, size_t N, uint_fast64_t B>
LinBox::RandomPrimeIter MargeGenMost<T, N, B>::_primeItr{ static_cast<uint_fast64_t>(std::floor(std::log2(B))) };

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

#if !defined(H_PRIME_GEN)
#define H_PRIME_GEN

#include <linbox/randiter/random-prime.h>
#include <array>
#include <ostream>

using namespace std;

// generate N random primes each of s-bits length
// use type T
template <typename T, size_t N, size_t B>
class PrimeGen {
   protected:
    LinBox::RandomPrimeIter _primeIt;
    array<T, N> _result;

   public:
    PrimeGen() {
        LinBox::Integer s = static_cast<int>(B);
        for (int i = 0; i < N; i++) {
            _result[i] = _primeIt.random(s);
        }
    }

    friend ostream& operator<<(ostream& out, PrimeGen<int, N, B>& pg);
    T operator[](size_t i) { return _result[i]; };

    ~PrimeGen() {}
};

template <size_t N, size_t B>
ostream& operator<<(ostream& out, PrimeGen<int, N, B>& pg) {
    out << "[";
    for (size_t i = 0; i < N; i++) {
        int v = pg._result[i];
        out << v;
        if (i != N - 1) {
            out << ",";
        }
    }
    out << "]";
    return out;
}
// template<typename t, size_t n, int s>
// array<t,n> prime_gen();

#endif  // H_PRIME_GEN

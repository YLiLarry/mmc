#if !defined(H_PRIME_GEN)
#define H_PRIME_GEN

#include <array>
#include <linbox/randiter/random-prime.h>

using namespace std;

// generate N random primes each of s-bits length
// use type T
template<typename T, size_t N> 
class PrimGen<T> {
    Array<T, T>* result;
    public: 
        PrimeGen() {
            result = new Array<T>();   
        }
        ~PrimeGen() {
            delete result;
        }
}

template<typename t, size_t n, int s> 
array<t,n> prime_gen();

#endif // H_PRIME_GEN

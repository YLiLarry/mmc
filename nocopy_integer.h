#if !defined(H_NOCOPY_INTEGER)
#define H_NOCOPY_INTEGER

#include <givaro/givtimer.h>
#include <linbox/integer.h>
#include <linbox/randiter/random-prime.h>

class LInteger : public LinBox::Integer
{
  public:
    LInteger() = default;
    explicit LInteger(const double &i)
        : LinBox::Integer(i)
    {
    }
    explicit LInteger(const int &i)
        : LinBox::Integer(i)
    {
    }

    LInteger(const LInteger &) = default;
    // explicit LInteger(const LInteger &i)
    //     : LinBox::Integer(i)
    // {
    // }
    explicit LInteger(const Integer &i)
        : LinBox::Integer(i)
    {
    }
    LInteger(const LInteger &src, const uint_fast64_t &from,
             const uint_fast64_t &to)
        : LInteger{static_cast<LInteger>(src & static_cast<LInteger>(1 << ((to - from) - 1)))}
    {
    }
    // clone another LInteger
    void EXPENSIVE_COPY(const LInteger &a)
    {
        Integer::operator=(a);
    } // clone another LInteger

    void EXPENSIVE_COPY(const Integer &a)
    {
        Integer::operator=(a);
    }
    // does the same as the copy constructor, except on the (*this)
    // itself (inplace)
    void cut(const uint_fast64_t &from, const uint_fast64_t &to)
    {
        operator&=(uint_fast64_t(1 << ((to - from) - 1)));
    }

    static LInteger random_exact(const uint_fast64_t p)
    {
        LInteger a;
        a.randomize(p);
        return a;
    }

    // rewrite the current value to a random bits of length p
    void randomize(const uint_fast64_t p)
    {
#if !PSEUDO_RANDOM_MMC
        Integer::seeding();
#endif
        Integer::random_exact_2exp((*this), p);
    }

    void randomizePrime(const uint_fast64_t p)
    {
#if PSEUDO_RANDOM_MMC
        LinBox::RandomPrimeIter pit{p};
#else
        LinBox::RandomPrimeIter pit{p, (uint64_t)LinBox::BaseTimer::seed()};
#endif
        pit.random_exact(*this);
    }
};

// this is just a wrapper on LinBox::Integer (by transitivity, gmp's multi-precision integer)
// with all the copy constructors, copy assignment operators deleted
// the motivation is to represent large integer primes in memory and make sure that they are
// only referred by C++ references or pointers and the compiler does not pass the primes by copies
// when you call a function with them.
class NoCopyInteger : public LInteger
{
  public:
    NoCopyInteger() = default;
    explicit NoCopyInteger(const double &i)
        : LInteger(i)
    {
    }

    ~NoCopyInteger() = default;
    // copy constrcut the [from,to) bits from src
    // fast bit extract [from,to) (right exclusive)
    // be very careful: do not copy the entire |src| bits but only what you need (usually fewer bits)
    // assuming src is a very large NoCopyInteger (eg, src >= 2^(2^10)) this is much cheaper
    NoCopyInteger(const NoCopyInteger &src, const uint_fast64_t &from,
                  const uint_fast64_t &to)
        : LInteger{src, from, to}
    {
    }
    static NoCopyInteger *newRandom(const uint_fast64_t n)
    {
        NoCopyInteger *p = new NoCopyInteger();
        p->randomize(n);
        return p;
    }

    // NoCopyInteger& operator=(const Integer& i) { return operator=(i); }
    // NoCopyInteger& operator=(const int& i) { return operator=(i); }
    const NoCopyInteger &operator=(const uint_fast8_t &i)
    {
        Integer::operator=(i);
        return (*this);
    }

    // prevent from accidentally writing copy or assign-copy
    // do not write NoCopyInteger i = something;
    // this will deep-copy in memory and is expensive if something is big
    NoCopyInteger &operator=(const NoCopyInteger &) = delete;
    NoCopyInteger &operator=(const Integer &) = delete;
    NoCopyInteger(const Integer &) = delete;
    NoCopyInteger(const NoCopyInteger &) = delete;
};

#endif // H_NOCOPY_INTEGER

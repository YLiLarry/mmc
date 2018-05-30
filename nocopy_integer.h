#if !defined(H_NOCOPY_INTEGER)
#define H_NOCOPY_INTEGER

#include <linbox/integer.h>

using namespace LinBox;

class NoCopyInteger : public Integer {
public:
    NoCopyInteger() = default;
    explicit NoCopyInteger(const uint_fast8_t& i)
        : Integer(i)
    {
    }

    ~NoCopyInteger() = default;
    // copy the [from,to) bits from src
    // fast bit extract [from,to) (right exclusive)
    // be very careful do not copy any values that are |src| bits long
    // (ie, pass-by-reference)
    // assuming src is a very large NoCopyInteger (eg, src >= 2^(2^10))
    NoCopyInteger(const NoCopyInteger& src, const uint_fast64_t& from,
        const uint_fast64_t& to)
        : Integer{ src & uint_fast64_t(1 << ((to - from) - 1)) }
    {
    }
    // does the same as the copy constructor, except on the existing instance
    // itself (inplace)
    void cut(const uint_fast64_t& from, const uint_fast64_t& to)
    {
        operator&=(uint_fast64_t(1 << ((to - from) - 1)));
    }
    // rewrite the current value to a random bits of length 2^n
    // ie, between [2^(n-1), 2^n-1]
    void random(const uint_fast64_t p)
    {
        Integer::random_exact_2exp((*this), p);
    }
    // NoCopyInteger& operator=(const Integer& i) { return operator=(i); }
    // NoCopyInteger& operator=(const int& i) { return operator=(i); }
    const NoCopyInteger& operator=(const uint_fast8_t& i)
    {
        Integer::operator=(i);
        return (*this);
    }

    // prevent from accidentally writing copy or assign-copy
    // do not write NoCopyInteger i = something;
    // this will deep-copy in memory and is expensive if something is big
    NoCopyInteger& operator=(const NoCopyInteger&) = delete;
    NoCopyInteger& operator=(const Integer&) = delete;
    NoCopyInteger(const Integer&) = delete;
    NoCopyInteger(const NoCopyInteger&) = delete;
};

#endif // H_NOCOPY_INTEGER

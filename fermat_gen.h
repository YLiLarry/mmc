#if !defined(H_FERMAT_GEN)
#define H_FERMAT_GEN

#include <gmp++/gmp++.h>
#include "prime_gen.h"
#include <cmath>

// generates N co-primes of form 2^(2^n) + 1 where
// the largest one is at most most_bits, and
// the product of them are at most 2^product_bits bits long,
// since there are not many fermat numbers we don't generate randomly
// we just start from the largest and count down until the product is large enough
template <class T>
class FermatGenMost : public NumPtrVector<T>
{
  protected:
    Givaro::Integer _product = 1;

  public:
    const Givaro::Integer &product = _product;

  public:
    FermatGenMost(size_t most_bits, size_t product_bits) : NumPtrVector<T>(0)
    {
        _fermat_nums = new NumPtrArray<T>();
        assert(nbits > 0 && "upper bound must be a positive integer");
        for (size_t i = most_bits; product_sofar.bitsize() <= product_bits; i--)
        {
            if (i <= 0)
            {
                cerr << "Cannot generate enough fermat numbers, consider increasing bit size uppder bound." << endl;
                abort();
            }
            assert(i > 0);
            T *t = new T{1};
            t <<= i;
            t += 1;
            push_back(t); // t = 2^i + 1
        }
    }

    friend ostream &operator<<(ostream &out, const FermatGenMost<LInteger> &arr)
    {
        out << "FermatGenMost [";
        for (size_t i = 0; i < N; i++)
        {
            const Integer &t = arr[i];
            const uint_fast64_t p = (t - 1).bitsize() - 1;
            out << "2^" << p << "+1=" << t;
            if (i != N - 1)
            {
                out << " , ";
            }
        }
        out << "]" << endl
            << " - count: " << arr.count() << endl
            << " - max bit length " << B << endl
            << " - product bit length: " << arr.product.bitsize() << endl;
        return out;
    }

    ~FermatGenMost() = default;
    FermatGenMost(FermatGenMost &) = delete;
};

#endif // H_FERMAT_GEN

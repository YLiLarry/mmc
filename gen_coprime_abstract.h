#ifndef H_COPRIME_GEN_ABSTRACT
#define H_COPRIME_GEN_ABSTRACT

#include <gmp++/gmp++.h>
#include <vector>
#include <iterator>

template <class T>
class GenCoprimeAbstract : public std::vector<T>
{
  public:
    virtual const T &max() const = 0;
    virtual Givaro::Integer product() const = 0;
    virtual uint_fast64_t product_bitsize() const = 0;
    virtual uint_fast64_t max_bitsize() const = 0;
    virtual ~GenCoprimeAbstract() = default;

    inline const T &val(size_t i) const { return std::vector<T>::operator[](i); }
    inline size_t count() const { return std::vector<T>::size(); } // alias

    friend std::ostream &operator<<(std::ostream &out, const GenCoprimeAbstract<T> &arr)
    {
#if DEBUG_MMC
        size_t N = arr.count();
        out << "[";
        for (size_t i = 0; i < N; i++)
        {
            out << Givaro::Integer(arr.val(i));
            if (i != N - 1)
            {
                out << " , ";
            }
        }
        out << "]" << endl;
#endif
#if DEBUG_MMC || TIME_MMC
        cerr << " - count: " << arr.count() << endl;
#endif
#if DEBUG_MMC
        cerr << " - max: " << Givaro::Integer(arr.max()) << endl;
#endif
#if DEBUG_MMC || TIME_MMC
        cerr << " - max bit length: " << arr.max_bitsize() << endl;
        cerr << " - product bit length: " << arr.product_bitsize() << endl;
#endif
#if DEBUG_MMC
        if (arr.product_bitsize() < 1024)
        {
            out << " - product: " << arr.product() << endl;
        }
#endif
        return out;
    }
};

#endif
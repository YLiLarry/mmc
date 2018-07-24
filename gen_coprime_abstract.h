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
        out << "]" << std::endl;
#endif
#if DEBUG_MMC || TIME_MMC
        std::cerr << " - count: " << arr.count() << std::endl;
#endif
#if DEBUG_MMC
        std::cerr << " - max: " << Givaro::Integer(arr.max()) << std::endl;
#endif
#if DEBUG_MMC || TIME_MMC
        std::cerr << " - max bit length: " << arr.max_bitsize() << std::endl;
        std::cerr << " - product bit length: " << arr.product_bitsize() << std::endl;
#endif
#if DEBUG_MMC
        if (arr.product_bitsize() < 1024)
        {
            out << " - product: " << arr.product() << std::endl;
        }
#endif
        return out;
    }
};

#endif
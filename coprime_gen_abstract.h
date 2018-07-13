#ifndef H_COPRIME_GEN_ABSTRACT
#define H_COPRIME_GEN_ABSTRACT

#include <linbox/integer.h>
#include <vector>
#include <iterator>

template <class T>
class CoprimeGenAbstract : public vector<T>
{
  public:
    inline virtual const T &max() const = 0;
    inline virtual Givaro::Integer product() const = 0;
    inline virtual uint_fast64_t product_bitsize() const = 0;
    inline virtual uint_fast64_t max_bitsize() const = 0;
    virtual ~CoprimeGenAbstract() = default;

    inline const T &val(size_t i) const { return vector<T>::operator[](i); }
    inline size_t count() const { return vector<T>::size(); } // alias

    friend ostream &operator<<(ostream &out, const CoprimeGenAbstract<T> &arr)
    {
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
        out << "]" << endl
            << " - count: " << arr.count() << endl
            << " - max " << Givaro::Integer(arr.max()) << endl
            << " - max bit length " << arr.max_bitsize() << endl
            << " - product bit length: " << arr.product_bitsize() << endl;
        if (arr.product_bitsize() < 1024)
        {
            out << " - product: " << arr.product() << endl;
        }
        return out;
    }
};

#endif
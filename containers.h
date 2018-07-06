#if !defined(H_CONTAINERS)
#define H_CONTAINERS

#include "nocopy_integer.h"
#include <array>
#include <cstdint>
#include <functional>
#include <linbox/randiter/random-prime.h>
#include <ostream>
#include <fflas-ffpack/fflas/fflas.h>

using namespace std;

// a helper class that frees fflas_new allocated memory when destroyed
template <class Field>
class FFLAS_Mem
{
  public:
    const typename Field::Element_ptr data;
    FFLAS_Mem(const typename Field::Element_ptr &data)
        : data(data)
    {
    }
    ~FFLAS_Mem()
    {
        FFLAS::fflas_delete(data);
    }
};

template <class T>
class PtrVector : public vector<T *>
{
  public:
    PtrVector(size_t len)
        : vector<T *>(len)
    {
    }
    T &val(size_t i) const { return *(this->operator[](i)); }
    T *&ptr(size_t i) { return this->operator[](i); }
    T *ptr(size_t i) const { return this->operator[](i); }
    size_t length() const { return this->size(); }
    ~PtrVector()
    {
        this->erase(this->begin());
    }
};

template <class T>
class NumPtrVector : public PtrVector<T>
{
  public:
    NumPtrVector(size_t len)
        : PtrVector<T>(len)
    {
    }

    T &operator[](const size_t i) const
    {
        return this->val(i);
    }

    bool equals(const vector<T> &other) const
    {
        size_t N = this->size();
        if (N != other.size())
        {
            return false;
        }
        for (size_t i = 0; i < N; i++)
        {
            if (this->val(i) != other[i])
            {
                return false;
            }
        }
        return true;
    }
    bool equals(const NumPtrVector<T> &other) const
    {
        size_t N = this->size();
        if (N != other.size())
        {
            return false;
        }
        for (size_t i = 0; i < N; i++)
        {
            if (!this->ptr(i) || !other.ptr(i))
            {
                return false;
            }
            if (this->val(i) != other.val(i))
            {
                return false;
            }
        }
        return true;
    }

    vector<Givaro::Integer> *EXPENSIVE_NEW_INTEGER_VECTOR() const
    {
        size_t N = this->length();
        auto vec = new vector<Givaro::Integer>(N);
        for (size_t i = 0; i < N; i++)
        {
            vec->operator[](i) = this->val(i);
        }
        return vec;
    }

    friend ostream &operator<<(ostream &out, const NumPtrVector<T> &arr)
    {
        size_t N = arr.size();
        out << "NumPtrVector [";
        for (size_t i = 0; i < N; i++)
        {
            out << arr[i];
            if (i != N - 1)
            {
                out << " , ";
            }
        }
        out << "]";
        return out;
    }
};

class LIntegerPtrVector : public NumPtrVector<LInteger>
{
  public:
    LIntegerPtrVector(size_t size)
        : NumPtrVector<LInteger>(size)
    {
    }
    // LIntegerPtrVector(const NumPtrVector<Givaro::Integer> &other)
    //     : NumPtrVector<LInteger>(other.size())
    // {
    //     size_t size = other.size();
    //     for (size_t i = 0; i < size; i++)
    //     {
    //         LInteger *t = new LInteger(other.val(i));
    //         this->ptr(i) = t;
    //     }
    // }

    // bool equals(const NumPtrVector<Givaro::Integer> &other) const
    // {
    //     NumPtrVector<Givaro::Integer>
    //     LIntegerPtrVector tmp(other);
    //     return other.equals(tmp);
    // }

    static LIntegerPtrVector *new_random(size_t size, uint_fast64_t bitsize)
    {
        LIntegerPtrVector *v = new LIntegerPtrVector(size);
        for (size_t i = 0; i < size; i++)
        {
            LInteger *t = new LInteger(0);
            t->randomize(bitsize);
            v->ptr(i) = t;
        }
        return v;
    }
    static LIntegerPtrVector *new_zeros(size_t size)
    {
        LIntegerPtrVector *v = new LIntegerPtrVector(size);
        for (size_t i = 0; i < size; i++)
        {
            LInteger *t = new LInteger(0);
            v->ptr(i) = t;
        }
        return v;
    }
};

template <class T>
using PtrAllocator = std::function<T *(size_t index)>;

// used to store a list of large primes as pointers
template <class T, size_t N>
class PtrArray : public array<T *, N>
{
  public:
    PtrArray() = default;

    PtrArray(const PtrAllocator<T> &allocator)
    {
        for (size_t i = 0; i < N; i++)
        {
            this->operator[](i) = allocator(i);
        }
    }

    PtrVector<T> *EXPENSIVE_NEW_PTR_VECTOR() const
    {
        auto vec = new PtrVector<T>(N);
        for (size_t i = 0; i < N; i++)
        {
            vec->ptr(i) = this->ptr(i);
        }
        return vec;
    }

    // default C++ array<T,N> destructor frees the memory
    ~PtrArray() = default;

    PtrArray(const PtrArray &) = delete;
    const PtrArray &operator=(const PtrArray &) = delete;

    T &val(size_t i) const { return *(this->operator[](i)); }
    T *ptr(size_t i) const { return this->operator[](i); }

    size_t count() const { return N; }

    friend ostream &operator<<(ostream &out, const PtrArray<T, N> &arr)
    {
        out << "PtrArray [";
        for (size_t i = 0; i < N; i++)
        {
            out << arr[i];
            if (i != N - 1)
            {
                out << " , ";
            }
        }
        out << "]";
        return out;
    };
};

template <class T, size_t N>
class NumPtrArray : public PtrArray<T, N>
{
  public:
    NumPtrArray() = default;

    NumPtrArray(const PtrAllocator<T> &allocator)
        : PtrArray<T, N>(allocator)
    {
    }

    T &operator[](const size_t i) const
    {
        return this->val(i);
    }

    NumPtrVector<T> *EXPENSIVE_NEW_NUM_PTR_VECTOR() const
    {
        auto vec = new NumPtrVector<T>(N);
        for (size_t i = 0; i < N; i++)
        {
            vec->ptr(i) = this->ptr(i);
        }
        return vec;
    }

    vector<Givaro::Integer> *EXPENSIVE_NEW_INTEGER_VECTOR() const
    {
        auto vec = new vector<Givaro::Integer>(N);
        for (size_t i = 0; i < N; i++)
        {
            vec->operator[](i) = this->val(i);
        }
        return vec;
    }

    friend ostream &operator<<(ostream &out, const NumPtrArray<T, N> &arr)
    {
        out << "NumPtrArray [";
        for (size_t i = 0; i < N; i++)
        {
            out << arr[i];
            if (i != N - 1)
            {
                out << " , ";
            }
        }
        out << "]";
        return out;
    }

    bool equals(const NumPtrArray<T, N> &other) const
    {
        if (N != other.count())
        {
            return false;
        }
        for (size_t i = 0; i < N; i++)
        {
            if (!this->ptr(i) || !other.ptr(i))
            {
                return false;
            }
            if (this->val(i) != other.val(i))
            {
                return false;
            }
        }
        return true;
    }

    ~NumPtrArray() = default;

    NumPtrArray(const NumPtrArray &) = delete;
    const NumPtrArray &operator=(const NumPtrArray &) = delete;
};

template <class T, size_t N>
class ConstNumPtrArray : public NumPtrArray<const T, N>
{
  private:
    Givaro::Integer _product;
    Givaro::Integer _sum;

  public:
    const Givaro::Integer &product;
    const Givaro::Integer &sum;

    ConstNumPtrArray(const PtrAllocator<const T> &allocator)
        : NumPtrArray<const T, N>(allocator), _product(N > 0 ? 1 : 0), _sum(0), product(_product), sum(_sum)
    {
        for (size_t i = 0; i < N; i++)
        {
            _product *= this->operator[](i);
            _sum += this->operator[](i);
        }
    }

    const T *ptr(const size_t &i) const
    {
        return PtrArray<const T, N>::operator[](i);
    }

    const T &operator[](const size_t &i) const
    {
        return *(PtrArray<const T, N>::operator[](i));
    }

    NumPtrVector<T> *EXPENSIVE_NEW_NUM_PTR_VECTOR() const
    {
        auto vec = new NumPtrVector<T>(N);
        for (size_t i = 0; i < N; i++)
        {
            vec->ptr(i) = new T{this->val(i)};
        }
        return vec;
    }

    friend ostream &operator<<(ostream &out,
                               const ConstNumPtrArray<T, N> &arr)
    {
        out << "ConstNumPtrArray [";
        for (size_t i = 0; i < N; i++)
        {
            out << arr[i];
            if (i != N - 1)
            {
                out << " , ";
            }
        }
        out << "]";
        return out;
    }

    ~ConstNumPtrArray() = default;

    ConstNumPtrArray(const ConstNumPtrArray &) = delete;
    const ConstNumPtrArray &operator=(const ConstNumPtrArray &) = delete;
};

#endif // H_CONTAINERS

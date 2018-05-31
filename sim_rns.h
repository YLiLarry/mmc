#if !defined(H_SIM_RNS)
#define H_SIM_RNS

#include "nocopy_integer.h"
#include "prime_gen.h"
#include <fflas-ffpack/field/rns.h>
#include <linbox/integer.h>
#include <vector>

namespace SIM_RNS {

template <class T_I, class T_M, size_t N>
class RNS {

public:
    class ReducedInt {

    public:
        const T_I& original;
        const ConstNumPtrArray<T_M, N> residuals;
        const ConstNumPtrArray<T_M, N>& modulis;
        ReducedInt(const T_I& _original, const ConstNumPtrArray<T_M, N>& _modulis, const PtrAllocator<T_M>& residuals_alloc)
            : original(_original)
            , residuals(residuals_alloc)
            , modulis(_modulis)
        {
        }
        ~ReducedInt() = default;
        friend ostream& operator<<(ostream& out, const ReducedInt& r)
        {
            out << endl
                << "type ReducedInt - " << endl
                << " original:" << r.original << endl
                << " modulis:" << r.modulis << endl
                << " residuals:" << r.residuals << endl;
            return out;
        }
    };

    typedef vector<shared_ptr<ReducedInt>> ReducedIntVector;
    static ReducedIntVector naive_reduce(const vector<const NoCopyInteger*>& inputs, const ConstNumPtrArray<T_M, N>& modulis)
    {
        size_t len_inputs = inputs.size();
        ReducedIntVector vec(len_inputs);
        for (size_t input_i = 0; input_i < len_inputs; input_i++) {
            const NoCopyInteger* ith_input = inputs[input_i];
            const PtrAllocator<T_M>& residuals_alloc = [&](size_t res_index) {
                T_M* p = new T_M;
                p->expensiveCopy(*ith_input);
                const T_M& m = modulis[res_index];
                p->operator%=(m);
                return p;
            };
            vec[input_i] = shared_ptr<ReducedInt>(new ReducedInt(*ith_input, modulis, residuals_alloc));
        }
        return vec;
    }

    static void recover(T_I& output, const ReducedIntVector&, const ConstNumPtrArray<T_M, N>& modulis)
    {
        return output;
    }

    friend ostream& operator<<(ostream& out, const ReducedIntVector& arr)
    {
        size_t n = arr.size();
        out << "type vector - [";
        for (size_t i = 0; i < n; i++) {
            out << (*arr[i]);
            if (i != n - 1) {
                out << ",";
            }
        }
        out << "]";
        return out;
    }
};
}

template <class T>
ostream& operator<<(ostream& out, const vector<T*>& arr)
{
    size_t n = arr.size();
    out << "type vector - [";
    for (size_t i = 0; i < n; i++) {
        out << (*arr[i]);
        if (i != n - 1) {
            out << ",";
        }
    }
    out << "]";
    return out;
}
#endif // H_SIM_RNS
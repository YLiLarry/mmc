#if !defined(H_SIM_RNS)
#define H_SIM_RNS

#include "containers.h"
#include "nocopy_integer.h"
#include "prime_gen.h"
#include <fflas-ffpack/field/rns.h>
#include <givaro/modular-integer.h>
#include <givaro/modular-uint64.h>
#include <linbox/algorithms/cra-domain-seq.h>
#include <linbox/algorithms/cra-early-single.h>
#include <linbox/algorithms/cra-givrnsfixed.h>
#include <linbox/integer.h>
#include <memory>
#include <vector>

namespace SIM_RNS {

template <class T_I, class T_M, size_t N>
class RNS {

public:
    class ReducedInt {

    public:
        NumPtrArray<T_M, N> residuals;
        const ConstNumPtrArray<T_M, N>& modulis;
        ReducedInt(const ConstNumPtrArray<T_M, N>& _modulis, const PtrAllocator<T_M>& residuals_alloc)
            : residuals(residuals_alloc)
            , modulis(_modulis)
        {
        }
        ~ReducedInt() = default;
        friend ostream& operator<<(ostream& out, const ReducedInt& r)
        {
            out << endl
                << "type ReducedInt - " << endl
                << " modulis:" << r.modulis << endl
                << " residuals:" << r.residuals << endl;
            return out;
        }

        void operator*=(const ReducedInt& other)
        {
            for (size_t i = 0; i < N; i++) {
                assert(&modulis[i] == &other.modulis[i]);
                residuals[i] *= other.residuals[i];
                residuals[i] %= modulis[i];
            }
        }
    };

    typedef vector<ReducedInt*> ReducedIntVector;
    static ReducedIntVector naive_reduce(const vector<const NoCopyInteger*>& inputs, const ConstNumPtrArray<T_M, N>& modulis)
    {
        size_t len_inputs = inputs.size();
        ReducedIntVector vec(len_inputs);
        for (size_t input_i = 0; input_i < len_inputs; input_i++) {
            const NoCopyInteger* ith_input = inputs[input_i];
            const PtrAllocator<T_M>& residuals_alloc = [&](size_t res_index) {
                T_M* p = new T_M(*ith_input);
                const T_M& m = modulis[res_index];
                (*p) %= m;
                return p;
            };
            vec[input_i] = new ReducedInt(modulis, residuals_alloc);
        }
        return vec;
    }

    static vector<NoCopyInteger*> naive_recover(const ReducedIntVector& inputs)
    {
        size_t input_len = inputs.size();
        vector<NoCopyInteger*> outvec = vector<NoCopyInteger*>(input_len);
        if (input_len > 0) {
            vector<Integer> modulis = inputs[0]->modulis.EXPENSIVE_TO_VECTOR();
            for (size_t i = 0; i < input_len; i++) {
                ReducedInt* reduced_int = inputs[i];
                // the numbers should have the same modulis
                assert(reduced_int->modulis == inputs[0]->modulis);
                // run CRT on reduced_int to produce a T_I*
                auto prime_res_mapping = [&](auto r, auto f) {
                    return r;
                };
                auto primeiter = modulis.begin();
                typedef LinBox::EarlySingleCRA<Givaro::Modular<T_M>> CRABase;
                CRABase base{};
                LinBox::ChineseRemainderSeq<CRABase> cra{ base };
                Integer result;
                cra(result, prime_res_mapping, primeiter);
                outvec[i]->EXPENSIVE_COPY(result);
            }
        }
        return outvec;
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

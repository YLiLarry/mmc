#if !defined(H_SIM_RNS)
#define H_SIM_RNS

#include "containers.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/field/rns-double.h"
#include "fflas-ffpack/field/rns-integer.h"
#include "fflas-ffpack/utils/args-parser.h"
#include "fflas-ffpack/utils/timer.h"
#include "nocopy_integer.h"
#include "prime_gen.h"
#include <givaro/modular-integer.h>
#include <givaro/modular-uint64.h>
#include <linbox/algorithms/cra-domain-seq.h>
#include <linbox/algorithms/cra-early-single.h>
#include <linbox/algorithms/cra-givrnsfixed.h>
#include <linbox/integer.h>
#include <memory>
#include <vector>

namespace SIM_RNS {

// class T_L: type used for large integers (inputs of reduction, outputs of recovery)
// class T_M: type used for modulis and residuals
// size_t N: how many modulis
template <class T_L, class T_M, size_t N_M>
class RNS {

public:
    // an integer representation in RNS
    // contains arrays for the residuals and modulis in the same order
    class ReducedInt {

    public:
        NumPtrArray<T_M, N_M> residuals;
        const ConstNumPtrArray<T_M, N_M>& modulis;
        ReducedInt(const ConstNumPtrArray<T_M, N_M>& _modulis, const PtrAllocator<T_M>& residuals_alloc)
            : residuals(residuals_alloc)
            , modulis(_modulis)
        {
        }
        ~ReducedInt() = default;
        friend ostream& operator<<(ostream& out, const ReducedInt& r)
        {
            out << endl
                << "ReducedInt " << r.residuals << endl
                << " - from modulis:" << r.modulis << endl;
            return out;
        }

        void operator*=(const ReducedInt& other)
        {
            for (size_t i = 0; i < N_M; i++) {
                assert(&modulis[i] == &other.modulis[i]);
                residuals[i] *= other.residuals[i];
                residuals[i] %= modulis[i];
            }
        }
    };

    static vector<ReducedInt*> naive_reduce(const vector<const T_L*>& inputs, const ConstNumPtrArray<T_M, N_M>& modulis)
    {
        size_t len_inputs = inputs.size();
        vector<ReducedInt*> vec(len_inputs);
        for (size_t input_i = 0; input_i < len_inputs; input_i++) {
            const T_L* ith_input = inputs[input_i];
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
    vector<ReducedInt*> naive_reduce(const vector<const NoCopyInteger*>& inputs, const ConstNumPtrArray<NoCopyInteger, N_M>& modulis)
    {
        size_t len_inputs = inputs.size();
        vector<ReducedInt*> vec(len_inputs);
        for (size_t input_i = 0; input_i < len_inputs; input_i++) {
            const NoCopyInteger* ith_input = inputs[input_i];
            const PtrAllocator<NoCopyInteger>& residuals_alloc = [&](size_t res_index) {
                NoCopyInteger* p = new NoCopyInteger;
                p->EXPENSIVE_COPY(*ith_input);
                const NoCopyInteger& m = modulis[res_index];
                p->operator%=(m);
                return p;
            };
            vec[input_i] = new ReducedInt(modulis, residuals_alloc);
        }
        return vec;
    }
    static vector<T_L*> naive_recover(const vector<ReducedInt*>& inputs)
    {
        size_t input_len = inputs.size();
        vector<T_L*> outvec = vector<T_L*>(input_len);
        if (input_len > 0) {
            vector<Integer> modulis = inputs[0]->modulis.EXPENSIVE_TO_VECTOR();
            for (size_t i = 0; i < input_len; i++) {
                ReducedInt* reduced_int = inputs[i];
                // the numbers should have the same modulis
                assert(reduced_int->modulis == inputs[0]->modulis);
                // run CRT on reduced_int to produce a T_L*
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

    // calls the algorithm for simultaneous reduction to RNS on inputs using modulis
    static vector<ReducedInt*> sim_reduce(const vector<const T_L*>& inputs, const ConstNumPtrArray<T_M, N_M>& modulis)
    {
        size_t len_inputs = inputs.size();
        // to adopt for FFLA's interface
        // wet up input_A
        Givaro::Modular<Givaro::Integer> input_field;
        typename Givaro::Modular<Givaro::Integer>::Element_ptr input_A = FFLAS::fflas_new(input_field, len_inputs, 1);
        size_t input_max_bitsize = input_field.cardinality().bitsize();
        size_t n_16bits_chunks = (input_max_bitsize / 16) + ((input_max_bitsize % 16) ? 1 : 0);
        for (size_t r = 0; r < len_inputs; r++) {
            input_A[r] = *inputs[r];
        }
        // set up output_A
        // RNS: contains an array of primes whose product is >= P, each of primes_bits long
        FFPACK::rns_double rns(modulis.EXPENSIVE_TO_VECTOR());
        // RNSInteger is a decorator (wrapper) on FFPACK::rns_double that adds some functionalities
        FFPACK::RNSInteger<FFPACK::rns_double> zrns(rns);
        typename FFPACK::RNSInteger<FFPACK::rns_double>::Element_ptr output_A = FFLAS::fflas_new(zrns, len_inputs);
        FFLAS::finit_rns(zrns, len_inputs, 1, n_16bits_chunks, input_A, 1, output_A);
        vector<ReducedInt*> output_vec(len_inputs);
        for (size_t idx_out = 0; idx_out < len_inputs; idx_out++) {
            output_vec[idx_out] = new ReducedInt(modulis, [&](size_t idx_moduli) {
                double* res = output_A[idx_moduli * len_inputs + idx_out]._ptr;
                return new T_M{ static_cast<T_M>(std::round(*res)) };
            });
        }
        FFLAS::fflas_delete(input_A);
        FFLAS::fflas_delete(output_A);
        return output_vec;
    }

    static vector<T_L*> sim_recover(const vector<ReducedInt*>& inputs)
    {
        size_t len_inputs = inputs.size();
        vector<T_L*> output_vec(len_inputs);
        if (len_inputs == 0) {
            return output_vec;
        }
        // assuming all inputs are reduced with the same set of modulis
        FFPACK::rns_double rns(inputs[0]->modulis.EXPENSIVE_TO_VECTOR());
        // RNSInteger is a decorator (wrapper) on FFPACK::rns_double that adds some functionalities
        FFPACK::RNSInteger<FFPACK::rns_double> zrns(rns);
        // to adopt for FFLA's interface
        // create two finite fields mod the larges number representable by Givaro::Integer
        Givaro::Modular<Givaro::Integer> output_field(-1);
        auto input_A = FFLAS::fflas_new(zrns, len_inputs);
        auto output_A = FFLAS::fflas_new(output_field, len_inputs);
        FFLAS::fconvert_rns(zrns, len_inputs, 1, Givaro::Integer(0), output_A, 1, input_A);
        for (size_t idx_out = 0; idx_out < len_inputs; idx_out++) {
            Givaro::Integer& e = output_A[idx_out];
            output_vec[idx_out] = new T_L{ e };
        }
        FFLAS::fflas_delete(input_A);
        FFLAS::fflas_delete(output_A);
        return output_vec;
    }

    friend ostream& operator<<(ostream& out, const vector<ReducedInt*>& arr)
    {
        size_t n = arr.size();
        out << "[";
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
    out << "[";
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

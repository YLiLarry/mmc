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

using namespace std;

// class T_L: type used for large integers (inputs of reduction, outputs of recovery)
// class T_M: type used for moduli and residuals
// size_t N: how many moduli

// an integer representation in RNS
// contains arrays for the residuals and moduli in the same order

template <class T_M, size_t N_M>
class ReducedInt {
public:
    NumPtrArray<T_M, N_M> residuals;
    const ConstNumPtrArray<T_M, N_M>& moduli;
    ReducedInt(const ConstNumPtrArray<T_M, N_M>& _moduli, const PtrAllocator<T_M>& residuals_alloc)
        : residuals(residuals_alloc)
        , moduli(_moduli)
    {
    }
    ~ReducedInt() = default;
    friend ostream& operator<<(ostream& out, const ReducedInt& r)
    {
        out << endl
            << "ReducedInt " << r.residuals << endl
            << " - from moduli:" << r.moduli << endl;
        return out;
    }

    void operator*=(const ReducedInt& other)
    {
        for (size_t i = 0; i < N_M; i++) {
            assert(&moduli[i] == &other.moduli[i]);
            residuals[i] *= other.residuals[i];
            residuals[i] %= moduli[i];
        }
    }
};

template <class T_L, class T_M, size_t N_M>
NumPtrVector<ReducedInt<T_M, N_M>>* new_naive_reduce(const NumPtrVector<const T_L>& inputs, const ConstNumPtrArray<T_M, N_M>& moduli)
{
    size_t len_inputs = inputs.size();
    auto vec = new NumPtrVector<ReducedInt<T_M, N_M>>(len_inputs);
    for (size_t input_i = 0; input_i < len_inputs; input_i++) {
        const T_L& ith_input = *inputs[input_i];
        vec[input_i] = new ReducedInt<T_M, N_M>(moduli, [&](size_t res_index) {
            const T_M& m = moduli[res_index];
            // cerr << "ith_input %= m " << endl
            //     << ith_input << endl
            //     << m << endl;
            T_M* p = new T_M{ static_cast<T_M>(ith_input % m) };
            // cerr << *p << endl;
            return p;
        });
    }
    return vec;
}

template <class T_L, class T_M, size_t N_M>
NumPtrVector<T_L>* new_naive_recover(const NumPtrVector<ReducedInt<T_M, N_M>>& inputs)
{
    size_t input_len = inputs.size();
    auto outvec = new NumPtrVector<T_L>(input_len);
    if (input_len > 0) {
        vector<Integer> moduli = inputs[0]->moduli.EXPENSIVE_NEW_VECTOR();
        for (size_t i = 0; i < input_len; i++) {
            ReducedInt<T_M, N_M>* reduced_int = inputs[i];
            // the numbers should have the same moduli
            assert(reduced_int->moduli == inputs[0]->moduli);
            // run CRT on reduced_int to produce a T_L*
            auto prime_res_mapping = [&](auto r, auto f) {
                return r;
            };
            auto primeiter = moduli.begin();
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

template <class T_L, class T_M, size_t N_M>
NumPtrVector<T_L>* new_sim_recover(const NumPtrVector<ReducedInt<T_M, N_M>>& inputs)
{
    size_t len_inputs = inputs.size();
    auto output_vec = new NumPtrVector<T_L>(len_inputs);
    if (len_inputs == 0) {
        return output_vec;
    }
    // assuming all inputs are reduced with the same set of moduli
    vector<Integer>* input_vec = inputs[0].moduli.EXPENSIVE_NEW_INTEGER_VECTOR();
    FFPACK::rns_double rns(*input_vec);
    delete input_vec;
    // RNSInteger is a decorator (wrapper) on FFPACK::rns_double that adds some functionalities
    FFPACK::RNSInteger<FFPACK::rns_double> z_rns(rns);
    // to adopt for FFLA's interface
    // create two finite fields mod the larges number representable by Givaro::Integer
    Givaro::Modular<Givaro::Integer> output_field(-1);
    auto input_A = FFLAS::fflas_new(z_rns, len_inputs);
    auto output_A = FFLAS::fflas_new(output_field, len_inputs);
    FFLAS::fconvert_rns(z_rns, len_inputs, 1, Givaro::Integer(0), output_A, 1, input_A);
    for (size_t idx_out = 0; idx_out < len_inputs; idx_out++) {
        Givaro::Integer& e = output_A[idx_out];
        output_vec->ptr(idx_out) = new T_L{ e };
    }
    FFLAS::fflas_delete(input_A);
    FFLAS::fflas_delete(output_A);
    return output_vec;
}

// calls the algorithm for simultaneous reduction to RNS on inputs using moduli
template <class T_L, class T_M, size_t N_M>
NumPtrVector<ReducedInt<T_M, N_M>>* new_sim_reduce(const NumPtrVector<T_L>& inputs, const ConstNumPtrArray<T_M, N_M>& moduli)
{
    cerr << "sim_reduce called with" << endl
         << " - inputs: " << inputs << endl
         << " - moduli: " << moduli << endl;
    size_t len_inputs = inputs.size();
    // to adopt for FFLA's interface
    // wet up input_A
    Givaro::Modular<Givaro::Integer> input_field;
    typename Givaro::Modular<Givaro::Integer>::Element_ptr input_A = FFLAS::fflas_new(input_field, len_inputs, 1);
    size_t input_max_bitsize = input_field.cardinality().bitsize();
    size_t n_16bits_chunks = (input_max_bitsize / 16) + ((input_max_bitsize % 16) ? 1 : 0);
    for (size_t r = 0; r < len_inputs; r++) {
        input_A[r] = inputs[r];
    }
    cerr << "check: " << *input_A << endl;
    // set up output_A
    // RNS: contains an array of primes whose product is >= P, each of primes_bits long
    const vector<Givaro::Integer>* moduli_vec = moduli.EXPENSIVE_NEW_INTEGER_VECTOR();
    FFPACK::rns_double rns(*moduli_vec);
    delete moduli_vec;
    // RNSInteger is a decorator (wrapper) on FFPACK::rns_double that adds some functionalities
    FFPACK::RNSInteger<FFPACK::rns_double> z_rns(rns);
    typename FFPACK::RNSInteger<FFPACK::rns_double>::Element_ptr output_A = FFLAS::fflas_new(z_rns, len_inputs);
    FFLAS::finit_rns(z_rns, len_inputs, 1, n_16bits_chunks, input_A, 1, output_A);
    auto output_vec = new NumPtrVector<ReducedInt<T_M, N_M>>(len_inputs);
    for (size_t idx_out = 0; idx_out < len_inputs; idx_out++) {
        ReducedInt<T_M, N_M>* ptr = new ReducedInt<T_M, N_M>(moduli, [&](size_t idx_moduli) {
            double* res = output_A[idx_moduli * len_inputs + idx_out]._ptr;
            cerr << inputs[idx_out] << " mod " << moduli[idx_moduli] << " = " << *res << endl;
            return new T_M{ static_cast<T_M>(std::round(*res)) };
        });
        ReducedInt<T_M, N_M>*& _ptr = output_vec->ptr(idx_out);
        _ptr = ptr;
    }
    FFLAS::fflas_delete(input_A);
    FFLAS::fflas_delete(output_A);
    return output_vec;
}
}
#endif // H_SIM_RNS

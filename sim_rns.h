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

namespace SIM_RNS
{

using namespace std;

// class T_L: type used for large integers (inputs of reduction, outputs of
// recovery) class T_M: type used for moduli and residuals size_t N: how many
// moduli

// an integer representation in RNS
// contains arrays for the residuals and moduli in the same order

template <class T_M, size_t N_M>
class ReducedInt
{
  public:
    NumPtrArray<T_M, N_M> residuals;
    const ConstNumPtrArray<T_M, N_M> &moduli;
    ReducedInt(const ConstNumPtrArray<T_M, N_M> &_moduli,
               const PtrAllocator<T_M> &residuals_alloc)
        : residuals(residuals_alloc), moduli(_moduli) {}
    ~ReducedInt() = default;
    friend ostream &operator<<(ostream &out, const ReducedInt &r)
    {
        out << endl
            << "ReducedInt " << r.residuals << endl
            << " - from moduli: " << r.moduli << endl;
        return out;
    }

    void operator*=(const ReducedInt &other)
    {
        for (size_t i = 0; i < N_M; i++)
        {
            assert(&moduli[i] == &other.moduli[i]);
            residuals[i] *= other.residuals[i];
            residuals[i] %= moduli[i];
        }
    }
};

template <class T_L, class T_M, size_t N_M>
vector<ReducedInt<T_M, N_M>> *new_naive_reduce(
    const vector<const T_L> &inputs,
    const ConstNumPtrArray<T_M, N_M> &moduli)
{
    size_t len_inputs = inputs.size();
    auto vec = new vector<ReducedInt<T_M, N_M>>(len_inputs);
    for (size_t input_i = 0; input_i < len_inputs; input_i++)
    {
        const T_L &ith_input = *inputs[input_i];
        vec[input_i] = new ReducedInt<T_M, N_M>(moduli, [&](size_t res_index) {
            const T_M &m = moduli[res_index];
            // cerr << "ith_input %= m " << endl
            //     << ith_input << endl
            //     << m << endl;
            T_M *p = new T_M{static_cast<T_M>(ith_input % m)};
            // cerr << *p << endl;
            return p;
        });
    }
    return vec;
}

template <class T_L, class T_M, size_t N_M>
vector<T_L> *new_naive_recover(const vector<ReducedInt<T_M, N_M>> &inputs)
{
    size_t input_len = inputs.size();
    auto outvec = new vector<T_L>(input_len);
    if (input_len > 0)
    {
        vector<Integer> moduli = inputs[0]->moduli.EXPENSIVE_NEW_VECTOR();
        for (size_t i = 0; i < input_len; i++)
        {
            ReducedInt<T_M, N_M> *reduced_int = inputs[i];
            // the numbers should have the same moduli
            assert(reduced_int->moduli == inputs[0]->moduli);
            // run CRT on reduced_int to produce a T_L*
            auto prime_res_mapping = [&](auto r, auto f) { return r; };
            auto primeiter = moduli.begin();
            typedef LinBox::EarlySingleCRA<Givaro::Modular<T_M>> CRABase;
            CRABase base{};
            LinBox::ChineseRemainderSeq<CRABase> cra{base};
            Integer result;
            cra(result, prime_res_mapping, primeiter);
            outvec[i]->EXPENSIVE_COPY(result);
        }
    }
    return outvec;
}

template <class T_L, class T_M, size_t N_M>
vector<T_L> *new_sim_recover(
    const vector<ReducedInt<T_M, N_M>> &inputs)
{
#if DEBUG_MMC
    cerr << "new_sim_recover called with" << endl
         << " - inputs: " << inputs << endl;
#endif
    size_t len_inputs = inputs.size();
    auto output_vec = new vector<T_L>(len_inputs);
    if (len_inputs == 0)
    {
        return output_vec;
    }
    // assuming all inputs are reduced with the same set of moduli
    const ConstNumPtrArray<T_M, N_M> &moduli = inputs[0].moduli;
    vector<Integer> *input_vec = moduli.EXPENSIVE_NEW_INTEGER_VECTOR();
    // cerr << "check: " << *input_vec << endl;
    FFPACK::rns_double rns(*input_vec);
    delete input_vec;
    // RNSInteger is a decorator (wrapper) on FFPACK::rns_double that adds some
    // functionalities
    FFPACK::RNSInteger<FFPACK::rns_double> z_rns(rns);
    // to adopt for FFLA's interface
    // create two finite fields mod the larges number representable by
    // Givaro::Integer
    typename FFPACK::RNSInteger<FFPACK::rns_double>::Element_ptr input_A =
        FFLAS::fflas_new(z_rns, len_inputs);
    for (size_t idx_in = 0; idx_in < len_inputs; idx_in++)
    {
        for (size_t idx_mod = 0; idx_mod < N_M; idx_mod++)
        {
            *input_A[idx_mod * len_inputs + idx_in]._ptr =
                static_cast<double>(inputs[idx_in].residuals[idx_mod]);
#if TEST_MMC
            // cerr << *input_A[idx_mod * len_inputs + idx_in]._ptr << " = " <<
            // inputs[idx_in].residuals[idx_mod] << endl;
            assert(Integer(*input_A[idx_mod * len_inputs + idx_in]._ptr) =
                       Integer(inputs[idx_in].residuals[idx_mod]));
            // cerr << "input_A[" << idx_mod << " * " << len_inputs << " + " <<
            // idx_in << "] = " << *input_A[idx_mod * len_inputs + idx_in]._ptr
            // << endl;
#endif
        }
    }
    Givaro::Modular<Givaro::Integer> output_field(moduli.product);
    typename Givaro::Modular<Givaro::Integer>::Element_ptr output_A =
        FFLAS::fflas_new(output_field, len_inputs, 1);
    FFLAS::fconvert_rns(z_rns, len_inputs, 1, Givaro::Integer(0), output_A, 1,
                        input_A);
    for (size_t idx_out = 0; idx_out < len_inputs; idx_out++)
    {
        // cerr << "check: " << output_A[idx_out] << endl;
        Givaro::Integer &e = output_A[idx_out];
        output_vec->ptr(idx_out) = new T_L{e};
    }
    FFLAS::fflas_delete(input_A);
    FFLAS::fflas_delete(output_A);
    return output_vec;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <class RNS_Field, class Int_Field>
vector<typename Int_Field::Element> fflas_new_sim_recover(
    const RNS_Field &rns_field, typename RNS_Field::Element_ptr input_A,
    size_t num_integers, const Int_Field &out_field)
{
#if DEBUG_MMC
    cerr << "########## fflas_new_sim_recover ##########" << endl
         << " - moduli: " << rns_field.rns()._basis << endl;
    for (size_t i = 0; i < num_integers; i++)
    {
        rns_field.write(cerr, input_A[i]);
        cerr << endl;
    }
#endif
    vector<typename Int_Field::Element> output_vec(num_integers);
    if (num_integers == 0)
    {
        return output_vec;
    }
    typename Givaro::Modular<Givaro::Integer>::Element_ptr output_A =
        FFLAS::fflas_new(out_field, num_integers, 1);
    FFLAS::fconvert_rns(rns_field, num_integers, 1, Givaro::Integer(0),
                        output_A, 1, input_A);
    Givaro::Modular<Givaro::Integer> basis_product_int_mod_field(rns_field.rns()._M);
    for (size_t idx_out = 0; idx_out < num_integers; idx_out++)
    {
        Givaro::Integer &e = output_A[idx_out];
        output_vec[idx_out] = typename Int_Field::Element(e);
    }
    FFLAS::fflas_delete(output_A);
#if DEBUG_MMC
    cerr << "output_vec: " << output_vec << endl
         << "########## fflas_new_sim_recover ends ##########" << endl;
#endif
    return output_vec;
}

template <class Int_Field, class RNS_Field>
typename RNS_Field::Element_ptr fflas_new_sim_reduce(
    const Int_Field &int_field,
    const vector<typename Int_Field::Element> &inputs,
    const RNS_Field &rns_field)
{
#if DEBUG_MMC
    cerr << "########## fflas_new_sim_reduce ##########" << endl
         << " - inputs: " << inputs << endl
         << " - moduli: " << rns_field.rns()._basis << endl;
#endif
    size_t len_inputs = inputs.size();
    assert(len_inputs > 0);
    typename Int_Field::Element_ptr fflas_inputs =
        FFLAS::fflas_new(int_field, len_inputs);
    typename Int_Field::Element max_input = inputs[0];
    for (size_t i = 0; i < len_inputs; i++)
    {
#if TEST_MMC
        assert(inputs[i] < rns_field.rns()._M &&
               "At least of one the input integers is too large.");
#endif
        if (inputs[i] > max_input)
        {
            max_input = inputs[i];
        }
        fflas_inputs[i] = inputs[i];
    }
    size_t input_max_bitsize = max_input.bitsize();
    size_t n_16bits_chunks =
        (input_max_bitsize / 16) + ((input_max_bitsize % 16) ? 1 : 0);
    typename RNS_Field::Element_ptr output_A =
        FFLAS::fflas_new(rns_field, len_inputs);
    FFLAS::finit_rns(rns_field, len_inputs, 1, n_16bits_chunks, fflas_inputs, 1,
                     output_A);
    FFLAS::fflas_delete(fflas_inputs);
#if DEBUG_MMC
    for (size_t i = 0; i < len_inputs; i++)
    {
        rns_field.write(cerr, output_A[i]);
        cerr << endl;
    }
    cerr << "########## fflas_new_sim_reduce ##########" << endl;
#endif
    return output_A;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <class T_L, class T_M, size_t N_M>
typename FFPACK::rns_double::Element_ptr fflas_new_sim_reduce(
    const vector<T_L> &inputs,
    const ConstNumPtrArray<T_M, N_M> &moduli)
{
#if DEBUG_MMC
    cerr << "fflas_new_sim_reduce called with" << endl
         << " - inputs: " << inputs << endl
         << " - moduli: " << moduli << endl
         << "    - product: " << moduli.product << endl;
#endif
    size_t len_inputs = inputs.size();
    // to adopt for FFLA's interface
    // wet up input_A
    Givaro::Modular<Givaro::Integer> input_field;
    typename Givaro::Modular<Givaro::Integer>::Element_ptr input_A =
        FFLAS::fflas_new(input_field, len_inputs, 1);
    size_t input_max_bitsize = moduli.product.bitsize();
    size_t n_16bits_chunks =
        (input_max_bitsize / 16) + ((input_max_bitsize % 16) ? 1 : 0);
    // cerr << "input_max_bitsize: " << input_max_bitsize << endl;
    // cerr << "n_16bits_chunks: " << n_16bits_chunks << endl;
    for (size_t r = 0; r < len_inputs; r++)
    {
#if TEST_MMC
        assert(inputs[r] < moduli.product &&
               "At least of one the input integers is too large.");
#endif
        input_A[r] = inputs[r];
    }
    // cerr << "check: " << *input_A << endl;
    // set up output_A
    // RNS: contains an array of primes whose product is >= P, each of
    // primes_bits long
    const vector<Givaro::Integer> *moduli_vec =
        moduli.EXPENSIVE_NEW_INTEGER_VECTOR();
    FFPACK::rns_double rns(*moduli_vec);
    delete moduli_vec;
    // RNSInteger is a decorator (wrapper) on FFPACK::rns_double that adds some
    // functionalities
    FFPACK::RNSInteger<FFPACK::rns_double> z_rns(rns);
    typename FFPACK::RNSInteger<FFPACK::rns_double>::Element_ptr output_A =
        FFLAS::fflas_new(z_rns, len_inputs);
    FFLAS::finit_rns(z_rns, len_inputs, 1, n_16bits_chunks, input_A, 1,
                     output_A);
    FFLAS::fflas_delete(input_A);
    return output_A;
}

// calls the algorithm for simultaneous reduction to RNS on inputs using moduli
template <class T_L, class T_M, size_t N_M>
vector<ReducedInt<T_M, N_M>> *new_sim_reduce(
    const vector<T_L> &inputs, const ConstNumPtrArray<T_M, N_M> &moduli)
{
#if DEBUG_MMC
    cerr << "new_sim_reduce called with" << endl
         << " - inputs: " << inputs << endl
         << " - moduli: " << moduli << endl
         << "    - product: " << moduli.product << endl;
#endif
    size_t len_inputs = inputs.size();
    typename FFPACK::RNSInteger<FFPACK::rns_double>::Element_ptr output_A =
        fflas_new_sim_reduce<T_L, T_M, N_M>(inputs, moduli);
    auto output_vec = new vector<ReducedInt<T_M, N_M>>(len_inputs);
    for (size_t idx_out = 0; idx_out < len_inputs; idx_out++)
    {
        ReducedInt<T_M, N_M> *ptr =
            new ReducedInt<T_M, N_M>(moduli, [&](size_t idx_moduli) {
                double *res = output_A[idx_moduli * len_inputs + idx_out]._ptr;
            // cerr << inputs[idx_out] << " mod " << moduli[idx_moduli] << " = "
            // << *res << endl;
#if TEST_MMC
                assert(Integer(inputs[idx_out]) % Integer(moduli[idx_moduli]) ==
                       Integer(*res));
#endif
                return new T_M{static_cast<T_M>(std::round(*res))};
            });
        ReducedInt<T_M, N_M> *&_ptr = output_vec->ptr(idx_out);
        _ptr = ptr;
    }
    FFLAS::fflas_delete(output_A);
    return output_vec;
}

vector<Givaro::Integer> fflas_mult_integer(
    const vector<Givaro::Integer> &matrix_a,
    const vector<Givaro::Integer> &matrix_b,
    size_t dim_m, size_t dim_n, size_t dim_k)
{
    using namespace FFLAS;
    assert(dim_m && dim_n && dim_k);
    // create matrix_c to return
    typedef Givaro::ZRing<Givaro::Integer> ZRing;
    ZRing zring;
    ZRing::Element_ptr fflas_mat_a = FFLAS::fflas_new(zring, dim_m, dim_n);
    ZRing::Element_ptr fflas_mat_b = FFLAS::fflas_new(zring, dim_n, dim_k);
    ZRing::Element_ptr fflas_mat_c = FFLAS::fflas_new(zring, dim_m, dim_k);

    size_t len_a = dim_m * dim_n;
    size_t len_b = dim_n * dim_k;
    for (size_t i = 0; i < len_a; i++)
    {
        fflas_mat_a[i] = matrix_a[i];
    }

    for (size_t i = 0; i < len_b; i++)
    {
        fflas_mat_b[i] = matrix_b[i];
    }

    FFLAS::fgemm(
        zring,               // field
        FFLAS::FflasNoTrans, // transpose matrix_a?
        FFLAS::FflasNoTrans, // transpose matrix_b?
        dim_m, dim_n, dim_k,
        zring.one, // coefficient before matrix_a*matrix_b
        fflas_mat_a,
        dim_n, // row length of matrix_a
        fflas_mat_b,
        dim_k,      // row length of matrix_b
        zring.zero, // constant matrix to add
        fflas_mat_c,
        dim_k // row length of matrix_c
    );

    size_t out_len = dim_m * dim_k;
    vector<Givaro::Integer> out(out_len);
    for (size_t i = 0; i < out_len; i++)
    {
        out[i] = fflas_mat_c[i];
    }

    return out;
};

} // namespace SIM_RNS
#endif // H_SIM_RNS

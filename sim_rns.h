#if !defined(H_SIM_RNS)
#define H_SIM_RNS

#include "containers.h"
#include <fflas-ffpack/field/rns-double.h>
#include <fflas-ffpack/field/rns-integer.h>
#include <fflas-ffpack/utils/timer.h>
#include "nocopy_integer.h"
#include "gen_prime.h"
#include <gmp++/gmp++.h>
#include <memory>
#include <vector>

namespace SIM_RNS
{

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
#if TIME_MMC
    cerr << ".......... fconvert_rns .........." << endl;
#endif
    FFLAS::fconvert_rns(rns_field,
                        num_integers,       // #rows
                        1,                  // #cols
                        Givaro::Integer(1), // sign
                        output_A,
                        1, // #cols
                        input_A);
#if TIME_MMC
    cerr << ".......... fconvert_rns ends .........." << endl;
#endif
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
         << " - inputs: " << inputs << endl;
    cerr << " - moduli: " << rns_field.rns()._basis << endl;
#endif
    size_t len_inputs = inputs.size();
    assert(len_inputs > 0);
    typename Int_Field::Element_ptr fflas_inputs =
        FFLAS::fflas_new(int_field, len_inputs);
    typename Int_Field::Element max_input = inputs[0];
    for (size_t i = 0; i < len_inputs; i++)
    {
#if CHECK_MMC
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
#if TIME_MMC
    cerr << ".......... finit_rns .........." << endl;
#endif
    FFLAS::finit_rns(rns_field, len_inputs, 1, n_16bits_chunks, fflas_inputs, 1,
                     output_A);
#if TIME_MMC
    cerr << ".......... finit_rns ends .........." << endl;
#endif
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

vector<Givaro::Integer> fflas_mult_integer(
    const vector<Givaro::Integer> &matrix_a,
    const vector<Givaro::Integer> &matrix_b,
    size_t dim_m, size_t dim_n, size_t dim_k)
{
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
        dim_k, // row length of matrix_c
        FFLAS::ParSeqHelper::Sequential());

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

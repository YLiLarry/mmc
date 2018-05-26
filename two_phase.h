#if !defined(H_TWO_STEPS_ALGO)
#define H_TWO_STEPS_ALGO

namespace CNMA
{
extern "C"
{
#include "cnma/marge_num.h"
#include "cnma/reconstruct.h"
}
} // namespace CNMA
#include "containers.h"
#include "coprime_gen.h"
#include "prime_gen.h"
#include "sim_rns.h"
#include <givaro/modular-integer.h>
#include <iostream>
#include <linbox/integer.h>
#include <fflas-ffpack/fflas/fflas.h>
#include <fflas-ffpack/field/rns-double.h>

// Phase 1:
// N_F is the number of co-primes moduli, each of bit length 2^B_F.
// co-primes are stored as T_F type in memory.
// Phase 2:
// N_M is the number of prime moduli, each of bit length B_M.
// primes are stored as T_M type in memory.
// The following relation must be true: N_M * B_M > 2^B_F
template <typename T_F, size_t N_F, uint_fast64_t B_F, typename T_M, size_t N_M, uint_fast64_t B_M>
class TwoPhaseAlgo
{
  protected:
    const MargeGenMost<T_F, N_F, B_F> _F_arr;  // co-primes in the first level
    const PrimeGenExact<T_M, N_M, B_M> _M_arr; // primes in the second level
    // Phase2_RNS_Field uses _M_arr repeats N_F times as modulis
    typedef FFPACK::rns_double Phase2_RNS_Rep;
    Phase2_RNS_Rep *_phase2_rns_rep;
    Phase2_RNS_Rep *_phase2_rns_computation_rep;

    typedef FFPACK::RNSInteger<Phase2_RNS_Rep> Phase2_RNS_Field;
    typedef Phase2_RNS_Field::Element Phase2_RNS_Int;
    typedef Phase2_RNS_Field::Element_ptr Phase2_RNS_Int_Ptr;
    Phase2_RNS_Field *_phase2_rns_field;
    Phase2_RNS_Field *_phase2_rns_computation_field;

    typedef Givaro::Modular<Givaro::Integer> Phase1_Field;
    typedef Phase1_Field::Element Phase1_Int;
    typedef Phase1_Field::Element_ptr Phase1_Int_Ptr;
    Phase1_Field *_phase1_field;

    class Phase2_Matrix
    {
        // shared ptr will be deleted when no Phase2_Matrix holds the FFLAS_Mem
        // so that FFLAS::fflas_delete is called
        shared_ptr<FFLAS_Mem<Phase2_RNS_Field>> _data;

      public:
        const Phase2_RNS_Int_Ptr data;
        const size_t dim_m;
        const size_t dim_n;
        const size_t count;

        Phase2_Matrix(const Phase2_RNS_Field &f, size_t dim_m, size_t dim_n)
            : _data(make_shared<FFLAS_Mem<Phase2_RNS_Field>>(FFLAS::fflas_new(f, dim_m, dim_n))),
              data(_data->data),
              dim_m(dim_m),
              dim_n(dim_n),
              count(dim_n * dim_m)
        {
            assert(this->data._stride == this->count);
        }

        Phase2_Matrix(const Phase2_RNS_Int_Ptr &arr, size_t dim_m, size_t dim_n)
            : _data(make_shared<FFLAS_Mem<Phase2_RNS_Field>>(arr)),
              data(_data->data),
              dim_m(dim_m),
              dim_n(dim_n),
              count(dim_n * dim_m)
        {
            assert(this->data._stride == this->count);
        }

        Phase2_RNS_Int ref(size_t r, size_t c, size_t f, size_t m) const
        {
            return data[f * count * N_M + m * count + r * dim_n + c];
        }

        friend ostream &operator<<(ostream &out, const Phase2_Matrix &mat)
        {
            out << "Phase2_Matrix: " << endl;
            for (size_t r = 0; r < mat.dim_m; r++)
            {
                for (size_t f = 0; f < N_F; f++)
                {
                    if (f != 0)
                    {
                        out << " ---";
                    }
                    for (size_t m = 0; m < N_M; m++)
                    {
                        if (m != 0)
                        {
                            out << " -";
                        }
                        for (size_t c = 0; c < mat.dim_n; c++)
                        {
                            out << " " << static_cast<u_int64_t>(*mat.ref(r, c, f, m)._ptr);
                        }
                    }
                }
                out << endl;
            }
            return out;
        }
    };

  public:
    TwoPhaseAlgo()
    {
        // _M_arr repeats N_F times
        array<T_M, N_M * N_F> tmp;
        for (size_t i = 0; i < N_F; i++)
        {
            for (size_t j = 0; j < N_M; j++)
            {
                tmp[i * N_M + j] = _M_arr[j];
            }
        }
        _phase2_rns_computation_rep = new Phase2_RNS_Rep{tmp};
        _phase2_rns_computation_field = new Phase2_RNS_Field(*_phase2_rns_computation_rep);
        auto vec = _M_arr.EXPENSIVE_NEW_INTEGER_VECTOR();
        _phase2_rns_rep = new Phase2_RNS_Rep{*vec};
        _phase2_rns_field = new Phase2_RNS_Field(*_phase2_rns_rep);
        _phase1_field = new Phase1_Field(_M_arr.product);
        delete vec;
#if DEBUG_MMC
        cerr << "########## TwoPhaseAlgo constructor ##########" << endl
             << " - _F_arr: " << _F_arr << endl
             << " - _M_arr: " << _M_arr << endl
             << " - _phase2_rns_field: " << _phase2_rns_field->rns()._basis << endl
             << " - _phase2_rns_computation_field: " << _phase2_rns_computation_field->rns()._basis << endl
             << "########## TwoPhaseAlgo constructor ends ##########" << endl;

#endif
        for (size_t f = 0; f < N_F; f++)
        {
            if (_F_arr.val(f) >= _M_arr.product)
            {
                cerr << "level 1 moduli " << _F_arr[f] << " is not smaller than the product of level 2 modulis " << _M_arr.product << endl;
                assert(_F_arr.val(f) < _M_arr.product);
            }
        }
        // assert(B_F < _M_arr.product.bitsize() && "At least one of the first level moduli is greater than the product of the second level moduli.");
    }
    ~TwoPhaseAlgo()
    {
        delete _phase1_field;
        delete _phase2_rns_rep;
        delete _phase2_rns_field;
        delete _phase2_rns_computation_rep;
        delete _phase2_rns_computation_field;
    };
    TwoPhaseAlgo(const TwoPhaseAlgo &) = delete;
    TwoPhaseAlgo &operator=(const TwoPhaseAlgo &) = delete;

    ///////////////////////////////////////////////////////////////////////////////////////////

    Phase2_Matrix matrix_reduce(const NumPtrVector<LInteger> &inputs, size_t dim_m, size_t dim_n) const
    {
        size_t len_inputs = inputs.length();
        assert(len_inputs == dim_m * dim_n && "input matrix dimension is incorrect");
        assert(dim_m > 0 && "input matrix dimension is incorrect");
        assert(dim_n > 0 && "input matrix dimension is incorrect");
#if DEBUG_MMC
        cerr << "########## matrix_reduce ##########" << endl
             << inputs << endl;
#endif
#if TEST_MMC
        for (size_t i = 0; i < len_inputs; i++)
        {
            assert(inputs[i] < _F_arr.product && "Inputs must be less than the product of first level moduli.");
        }
#endif
#if DEBUG_MMC
        cerr << "..... phase 1 reduce ....." << endl;
#endif
        // phase 1 begins
        // p1_reduced stores multi-moduli representation of each input
        NumPtrVector<Phase1_Int> p1_reduced(len_inputs * N_F);
        for (size_t i = 0; i < len_inputs; i++)
        {
            for (size_t f = 0; f < N_F; f++)
            {
                Phase1_Int *t = new Phase1_Int{inputs.val(i)};
                p1_reduced.ptr(i * N_F + f) = t;
#if TIME_MMC
                cerr << ".";
#endif
                CNMA::dc_reduce_minus(t->get_mpz(), (_F_arr[f] + 1).bitsize() - 1);
            }
        }
#if TIME_MMC
        cerr << endl;
#endif

#if DEBUG_MMC
        cerr << "phase 1 reduced: " << endl
             << p1_reduced << endl;
#endif

#if DEBUG_MMC
        cerr << "..... phase 2 reduce ....." << endl;
#endif
        // phase 2 begins
        Phase2_RNS_Int_Ptr phase2_outputs = SIM_RNS::fflas_new_sim_reduce(*_phase1_field, p1_reduced, *_phase2_rns_field);

        Phase2_Matrix mat(*_phase2_rns_computation_field, dim_m, dim_n);
        for (size_t r = 0; r < dim_m; r++)
        {
            for (size_t c = 0; c < dim_n; c++)
            {
                size_t i = r * mat.dim_n + c;
                for (size_t f = 0; f < N_F; f++)
                {
                    for (size_t m = 0; m < N_M; m++)
                    {
                        mat.ref(r, c, f, m)._ptr[0] = phase2_outputs[m * (mat.count * N_F) + (i * N_F) + f]._ptr[0];
                    }
                }
            }
        }
        FFLAS::fflas_delete(phase2_outputs);
#if DEBUG_MMC
        cerr << "phase 2 reduced: " << endl
             << mat << endl
             << "########## matrix_reduce ends ##########" << endl;
#endif
        return mat;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////

    NumPtrVector<Givaro::Integer> *matrix_recover(const Phase2_Matrix &mat) const
    {
#if DEBUG_MMC
        cerr << "########## matrix_recover ##########" << endl
             << mat << endl
             << "..... phase 2 recovery ....." << endl;

        // for (size_t i = 0; i < mat.count * N_F; i++)
        // {
        //     _phase2_rns_field->write(cerr, mat.data[i]);
        // }
#endif

        Phase2_RNS_Int_Ptr phase2_inputs = FFLAS::fflas_new(*_phase2_rns_field, N_F * mat.count);
        for (size_t r = 0; r < mat.dim_m; r++)
        {
            for (size_t c = 0; c < mat.dim_n; c++)
            {
                size_t i = r * mat.dim_n + c;
                for (size_t f = 0; f < N_F; f++)
                {
                    for (size_t m = 0; m < N_M; m++)
                    {
                        // phase2_inputs[(N_M * N_F) * i + f * N_M + m]._ptr[0] = mat.ref(r, c, f, m)._ptr[0];
                        phase2_inputs[m * (mat.count * N_F) + (i * N_F) + f]._ptr[0] = mat.ref(r, c, f, m)._ptr[0];
                    }
                }
            }
        }

        // phase 2 recovery begins

        NumPtrVector<Phase1_Int> *phase2_recovered = SIM_RNS::fflas_new_sim_recover(*_phase2_rns_field, phase2_inputs, mat.count * N_F, *_phase1_field);
#if DEBUG_MMC
        cerr << *phase2_recovered << endl;
#endif
#if TEST_MMC
        // for (size_t i = 0; i < mat.count; i++)
        // {
        //     for (size_t f = 0; f < N_F; f++)
        //     {
        //         if (phase2_recovered->val(i * N_F + f) >= _F_arr[f])
        //         {
        //             cerr << "something went wrong in matrix_recover:" << endl
        //                  << phase2_recovered->val(i * N_F + f) << endl
        //                  << _F_arr[f] << endl;
        //         }
        //         assert(phase2_recovered->val(i * N_F + f) < _F_arr[f]);
        //     }
        // }
#endif

#if DEBUG_MMC
        cerr << "..... phase 2 recovery ends ....." << endl;
        cerr << "..... phase 1 recovery ....." << endl;
#endif
        // phase 1 recovery begins
        const auto phase1_recovered = new NumPtrVector<Givaro::Integer>(mat.count);

        // initialization
        mpz_t _Mi[N_F]; // tmp
        mpz_t _f[N_F];
        mpz_t _r[N_F];
        for (size_t i = 0; i < N_F; i++)
        {
            mpz_init(_Mi[i]);
            mpz_init(_r[i]);
            mpz_init_set(_f[i], _F_arr[i].get_mpz());
        }
        CNMA::precompute_Mi(_Mi, _f, N_F);
        // recover
        for (size_t i = 0; i < mat.count; i++)
        {
            Givaro::Integer *t = new Givaro::Integer(0);
            phase1_recovered->ptr(i) = t;
            for (size_t f = 0; f < N_F; f++)
            {
                Phase1_Int &in = phase2_recovered->val(i * N_F + f);
                if (in >= _F_arr.product)
                {
                    cerr << "Intermediate value overflows. "
                         << "You may need larger/more first level moduli, "
                         << "or smaller/fewer second level moduli. " << endl;
                }
                assert(in < _F_arr.product);
                mpz_mod(_r[f], in.get_mpz(), _F_arr[f].get_mpz());
            }
            CNMA::garner(t->get_mpz(), N_F, _r, _f, _Mi);
            mpz_mod(t->get_mpz(), t->get_mpz(), _F_arr.product.get_mpz());
        }

        delete phase2_recovered;
#if DEBUG_MMC
        cerr << "..... phase 1 recovery ends ....." << endl;
        cerr << "new_recover finished: " << *phase1_recovered << endl
             << "########## matrix_recover ends ##########" << endl;
#endif
        return phase1_recovered;
    }

    Phase2_Matrix phase2_mult(
        const Phase2_Matrix &matrix_a,
        const Phase2_Matrix &matrix_b) const
    {
#if DEBUG_MMC
        cerr << "########## phase2_mult ##########" << endl
             << " - matrix_a: " << endl
             << matrix_a << endl
             << " - matrix_b: " << endl
             << matrix_b << endl;
#endif
        assert(matrix_a.dim_n == matrix_b.dim_m);
        Phase2_Matrix matrix_c = phase2_matrix_fgemm(matrix_a.data, matrix_b.data, matrix_a.dim_m, matrix_a.dim_n, matrix_b.dim_m);
#if DEBUG_MMC
        cerr << " - matrix product:" << endl
             << matrix_c
             << "########## phase2_mult ends ##########" << endl;
#endif
        return matrix_c;
    }

    Phase2_Matrix phase2_matrix_fgemm(
        const Phase2_RNS_Int_Ptr &matrix_a,
        const Phase2_RNS_Int_Ptr &matrix_b,
        size_t dim_m, size_t dim_n, size_t dim_k) const
    {
        return Phase2_Matrix(fflas_new_fgemm(matrix_a, matrix_b, dim_m, dim_n, dim_k), dim_m, dim_k);
    }

    Phase2_RNS_Int_Ptr fflas_new_fgemm(
        const Phase2_RNS_Int_Ptr &matrix_a,
        const Phase2_RNS_Int_Ptr &matrix_b,
        size_t dim_m, size_t dim_n, size_t dim_k) const
    {

        assert(dim_m && dim_n && dim_k);
        // create matrix_c to return
        Phase2_RNS_Int_Ptr matrix_c = FFLAS::fflas_new(*_phase2_rns_computation_field, dim_m, dim_k);

        assert(_phase2_rns_computation_field->size() == N_M * N_F);
        assert(matrix_a._stride == dim_m * dim_n);
        assert(matrix_b._stride == dim_n * dim_k);
        assert(matrix_c._stride == dim_m * dim_k);

        // Givaro::Integer t1;
        // cerr << _phase2_rns_computation_field->convert(t1, _phase2_rns_computation_field->one) << endl;
        // _phase2_rns_computation_field->write(cerr, _phase2_rns_computation_field->one) << endl;
        // _phase2_rns_computation_field->write(cerr, _phase2_rns_computation_field->zero) << endl;
        // Givaro::Integer card;
        // cerr << _phase2_rns_computation_field->cardinality(card) << endl;

        FFLAS::fgemm(
            *_phase2_rns_computation_field, // field
            FFLAS::FflasNoTrans,            // transpose matrix_a?
            FFLAS::FflasNoTrans,            // transpose matrix_b?
            dim_m, dim_n, dim_k,
            _phase2_rns_computation_field->one, // coefficient before matrix_a*matrix_b
            matrix_a,
            dim_n, // row length of matrix_a
            matrix_b,
            dim_k,                               // row length of matrix_b
            _phase2_rns_computation_field->zero, // constant matrix to add
            matrix_c,
            dim_k // row length of matrix_c
        );

        return matrix_c;
    };
};

#endif // H_TWO_STEPS_ALGO

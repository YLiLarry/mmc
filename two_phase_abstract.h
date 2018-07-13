#ifndef H_TWO_PHASE_ABSTRACT
#define H_TWO_PHASE_ABSTRACT

namespace CNMA
{
extern "C"
{
#include "cnma/marge_num.h"
#include "cnma/reconstruct.h"
}
} // namespace CNMA

#include "containers.h"
#include "marge_gen.h"
#include "prime_gen.h"
#include "sim_rns.h"
#include <givaro/modular-integer.h>
#include <iostream>
#include <linbox/integer.h>
#include <fflas-ffpack/fflas/fflas.h>
#include <fflas-ffpack/field/rns-double.h>
#include <vector>
#include "coprime_gen_abstract.h"

class TwoPhaseAbstract
{
  protected:
    // Phase2_RNS_Field uses level_2_moduli repeats level_1_moduli_count times as modulis
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
    const CoprimeGenAbstract<Givaro::Integer> *level_1_moduli;
    const CoprimeGenAbstract<double> *level_2_moduli;
    const size_t level_1_moduli_count;
    const size_t level_2_moduli_count;

  public:
    TwoPhaseAbstract(const CoprimeGenAbstract<Givaro::Integer> *level_1_moduli,
                     const CoprimeGenAbstract<double> *level_2_moduli)
        : level_1_moduli(level_1_moduli),
          level_2_moduli(level_2_moduli),
          level_1_moduli_count(level_1_moduli->count()),
          level_2_moduli_count(level_2_moduli->count())
    {
#if DEBUG_MMC
        cerr << "########## TwoPhaseAlgo constructor ##########" << endl
             << " - level_1_moduli: " << *level_1_moduli << endl
             << " - level_2_moduli: " << *level_2_moduli << endl;
#endif
        // level_2_moduli repeats level_1_moduli_count times
        vector<Givaro::Integer> tmp(level_2_moduli_count * level_1_moduli_count);
        for (size_t i = 0; i < level_1_moduli_count; i++)
        {
            for (size_t j = 0; j < level_2_moduli_count; j++)
            {
                tmp[i * level_2_moduli_count + j] = level_2_moduli->val(j);
            }
        }
        _phase2_rns_computation_rep = new Phase2_RNS_Rep{tmp};
        _phase2_rns_computation_field = new Phase2_RNS_Field(*_phase2_rns_computation_rep);
        _phase2_rns_rep = new Phase2_RNS_Rep{*level_2_moduli};
        _phase2_rns_field = new Phase2_RNS_Field(*_phase2_rns_rep);
        _phase1_field = new Phase1_Field(level_2_moduli->product());
#if DEBUG_MMC
        cerr << " - _phase2_rns_field: " << _phase2_rns_field->rns()._basis << endl
             << " - _phase2_rns_computation_field: " << _phase2_rns_computation_field->rns()._basis << endl
             << "########## TwoPhaseAlgo constructor ends ##########" << endl;
#endif
        for (size_t f = 0; f < level_1_moduli_count; f++)
        {
            if (level_1_moduli->val(f) >= level_2_moduli->product())
            {
                cerr << "level 1 moduli " << level_1_moduli->val(f) << " is not smaller than the product of level 2 modulis " << level_2_moduli->product() << endl;
                assert(level_1_moduli->val(f) < level_2_moduli->product());
            }
        }
        // assert(B_F < level_2_moduli->product_bitsize() && "At least one of the first level moduli is greater than the product of the second level moduli.");
    };

    virtual ~TwoPhaseAbstract()
    {
        delete _phase1_field;
        delete _phase2_rns_rep;
        delete _phase2_rns_field;
        delete _phase2_rns_computation_rep;
        delete _phase2_rns_computation_field;
    };
    TwoPhaseAbstract(const TwoPhaseAbstract &) = delete;
    TwoPhaseAbstract &operator=(const TwoPhaseAbstract &) = delete;

  public:
    class Phase2_Matrix
    {
        // shared ptr will be deleted when no Phase2_Matrix holds the FFLAS_Mem
        // so that FFLAS::fflas_delete is called
        shared_ptr<FFLAS_Mem<Phase2_RNS_Field>> _data;

      public:
        size_t dim_m;
        size_t dim_n;
        size_t count;
        size_t level_1_moduli_count;
        size_t level_2_moduli_count;
        inline const Phase2_RNS_Int_Ptr &data() const { return _data->data; }

        Phase2_Matrix() = default;

        Phase2_Matrix(const TwoPhaseAbstract &f, size_t dim_m, size_t dim_n)
            : _data(make_shared<FFLAS_Mem<Phase2_RNS_Field>>(FFLAS::fflas_new(*(f._phase2_rns_computation_field), dim_m, dim_n))),
              dim_m(dim_m),
              dim_n(dim_n),
              count(dim_n * dim_m),
              level_1_moduli_count(f.level_1_moduli_count),
              level_2_moduli_count(f.level_2_moduli_count)
        {
            assert(this->data()._stride == this->count);
        }

        Phase2_Matrix(const TwoPhaseAbstract &f,
                      const Phase2_RNS_Int_Ptr &arr,
                      size_t dim_m, size_t dim_n)
            : _data(make_shared<FFLAS_Mem<Phase2_RNS_Field>>(arr)),
              dim_m(dim_m),
              dim_n(dim_n),
              count(dim_n * dim_m),
              level_1_moduli_count(f.level_1_moduli_count),
              level_2_moduli_count(f.level_2_moduli_count)
        {
            assert(this->data()._stride == this->count);
        }

        Phase2_Matrix(const Phase2_Matrix &) = default;
        Phase2_Matrix &operator=(const Phase2_Matrix &) = default;

        const Phase2_RNS_Int ref(size_t r, size_t c, size_t f, size_t m) const
        {
            return data()[f * count * level_2_moduli_count + m * count + r * dim_n + c];
        }

        friend ostream &operator<<(ostream &out, const Phase2_Matrix &mat)
        {
            out << "Phase2_Matrix: " << endl;
            for (size_t r = 0; r < mat.dim_m; r++)
            {
                for (size_t f = 0; f < mat.level_1_moduli_count; f++)
                {
                    if (f != 0)
                    {
                        out << " ---";
                    }
                    for (size_t m = 0; m < mat.level_2_moduli_count; m++)
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

  protected:
    /*
        this helper method is used by matrix_product(...)
    */
    virtual const vector<Phase1_Int> matrix_reduce_phase_1(const vector<Givaro::Integer> &inputs) const = 0;

  protected:
    /* 
        use this method to recover from a phase 1 representations to integers
    */
    virtual const vector<Givaro::Integer> matrix_recover_phase_1(const vector<Phase1_Int> &phase2_recovered) const = 0;

  public:
    /* 
        use this method to reduce a single matrix to level 2
    */
    Phase2_Matrix matrix_reduce(const vector<Givaro::Integer> &inputs, size_t dim_m, size_t dim_n) const
    {
        size_t len_inputs = inputs.size();
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
            assert(inputs[i] < level_1_moduli->product() && "Inputs must be less than the product of first level moduli.");
        }
#endif
        vector<Phase1_Int> p1_reduced = matrix_reduce_phase_1(inputs);

#if DEBUG_MMC
        cerr << "..... phase 2 reduce ....." << endl;
#endif
        // phase 2 begins
        Phase2_RNS_Int_Ptr phase2_outputs = SIM_RNS::fflas_new_sim_reduce(*_phase1_field, p1_reduced, *_phase2_rns_field);

        Phase2_Matrix mat(*this, dim_m, dim_n);
        for (size_t r = 0; r < dim_m; r++)
        {
            for (size_t c = 0; c < dim_n; c++)
            {
                size_t i = r * mat.dim_n + c;
                for (size_t f = 0; f < level_1_moduli_count; f++)
                {
                    for (size_t m = 0; m < level_2_moduli_count; m++)
                    {
                        mat.ref(r, c, f, m)._ptr[0] = phase2_outputs[m * (mat.count * level_1_moduli_count) + (i * level_1_moduli_count) + f]._ptr[0];
                    }
                }
            }
        }
        FFLAS::fflas_delete(phase2_outputs);
#if DEBUG_MMC
        cerr << "phase 2 reduced: " << endl
             << mat << endl;
#endif
#if DEBUG_MMC
        cerr << "########## matrix_reduce ends ##########" << endl;
#endif
        return mat;
    }

  public:
    /* 
        use this method to reduce multiple matrices to level 2
        shoule be faster than reducing one by one
    */
    const vector<Phase2_Matrix> matrix_reduce(const vector<Givaro::Integer> &matrices, const vector<size_t> &dimensions)
    {
        size_t len_inputs = matrices.size();
        size_t num_matrices = dimensions.size();
#if DEBUG_MMC
        cerr << "########## matrix_reduce ##########" << endl
             << "reducing " << len_inputs << " matrices" << endl;
#endif
        assert(len_inputs > 0 && num_matrices > 2 && "need at least 2 matrices");

#if TEST_MMC
        size_t len_ints = 0;
        for (size_t i = 1; i < num_matrices; i++)
        {
            size_t m = dimensions[i - 1];
            size_t n = dimensions[i];
            len_ints += m * n;
        }
        assert(len_inputs == len_ints && "supplied inputs and dimensions don't match");
#endif
        const vector<Phase1_Int> p1_reduced = matrix_reduce_phase_1(matrices);

#if DEBUG_MMC
        cerr << "..... phase 2 reduce ....." << endl;
#endif
        // phase 2 begins
        Phase2_RNS_Int_Ptr phase2_outputs = SIM_RNS::fflas_new_sim_reduce(*_phase1_field, p1_reduced, *_phase2_rns_field);

        vector<Phase2_Matrix> outputs(num_matrices);
        for (size_t o = 1; o < num_matrices; o++)
        {
            size_t dim_m = dimensions[o - 1];
            size_t dim_n = dimensions[o];
            Phase2_Matrix mat(*this, dim_m, dim_n);
            for (size_t r = 0; r < dim_m; r++)
            {
                for (size_t c = 0; c < dim_n; c++)
                {
                    size_t i = r * mat.dim_n + c;
                    for (size_t f = 0; f < level_1_moduli_count; f++)
                    {
                        for (size_t m = 0; m < level_2_moduli_count; m++)
                        {
                            mat.ref(r, c, f, m)._ptr[0] = phase2_outputs[m * (mat.count * level_1_moduli_count) + (i * level_1_moduli_count) + f]._ptr[0];
                        }
                    }
                }
            }
            outputs[o] = mat;
        }
        FFLAS::fflas_delete(phase2_outputs);
#if DEBUG_MMC
        cerr << "phase 2 reduced: " << endl
             << outputs << endl;
#endif
#if DEBUG_MMC
        cerr << "########## matrix_reduce ends ##########" << endl;
#endif
        return outputs;
    }

  public:
    /* 
        use this method to recover from a single reduced matrix to phase 1 representations
    */
    virtual vector<Phase1_Int> matrix_recover_phase_2(const Phase2_Matrix &mat) const
    {
        Phase2_RNS_Int_Ptr phase2_inputs = FFLAS::fflas_new(*_phase2_rns_field, level_1_moduli_count * mat.count);
        for (size_t r = 0; r < mat.dim_m; r++)
        {
            for (size_t c = 0; c < mat.dim_n; c++)
            {
                size_t i = r * mat.dim_n + c;
                for (size_t f = 0; f < level_1_moduli_count; f++)
                {
                    for (size_t m = 0; m < level_2_moduli_count; m++)
                    {
                        // phase2_inputs[(level_2_moduli_count * level_1_moduli_count) * i + f * level_2_moduli_count + m]._ptr[0] = mat.ref(r, c, f, m)._ptr[0];
                        phase2_inputs[m * (mat.count * level_1_moduli_count) + (i * level_1_moduli_count) + f]._ptr[0] = mat.ref(r, c, f, m)._ptr[0];
                    }
                }
            }
        }

        // phase 2 recovery begins
        return SIM_RNS::fflas_new_sim_recover(*_phase2_rns_field, phase2_inputs, mat.count * level_1_moduli_count, *_phase1_field);
    }

    /* 
        use this method to recover from a single reduced matrix 
    */

  public:
    vector<Givaro::Integer> matrix_recover(const Phase2_Matrix &mat) const
    {
#if DEBUG_MMC
        cerr << "########## matrix_recover ##########" << endl
             << mat << endl
             << "..... phase 2 recovery ....." << endl;
#endif
        const vector<Phase1_Int> phase2_recovered = matrix_recover_phase_2(mat);
#if DEBUG_MMC
        cerr << phase2_recovered << endl
             << "..... phase 2 recovery ends ....." << endl
             << "..... phase 1 recovery ....." << endl;
#endif
        const vector<Givaro::Integer> phase1_recovered = matrix_recover_phase_1(phase2_recovered);
#if DEBUG_MMC
        cerr << "..... phase 1 recovery ends ....." << endl
             << "new_recover finished: " << phase1_recovered << endl
             << "########## matrix_recover ends ##########" << endl;
#endif
        return phase1_recovered;
    }

  public:
    /* 
        use this method to multiply two reduced matrices
    */
    Phase2_Matrix phase2_mult(const Phase2_Matrix &matrix_a, const Phase2_Matrix &matrix_b) const
    {
#if DEBUG_MMC
        cerr << "########## phase2_mult ##########" << endl
             << " - matrix_a: " << endl
             << matrix_a << endl
             << " - matrix_b: " << endl
             << matrix_b << endl;
#endif
        assert(matrix_a.dim_n == matrix_b.dim_m);
        Phase2_Matrix matrix_c = phase2_matrix_fgemm(matrix_a.data(), matrix_b.data(), matrix_a.dim_m, matrix_a.dim_n, matrix_b.dim_m);
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
        return Phase2_Matrix(*this, fflas_new_fgemm(matrix_a, matrix_b, dim_m, dim_n, dim_k), dim_m, dim_k);
    }

    Phase2_RNS_Int_Ptr fflas_new_fgemm(
        const Phase2_RNS_Int_Ptr &matrix_a,
        const Phase2_RNS_Int_Ptr &matrix_b,
        size_t dim_m, size_t dim_n, size_t dim_k) const
    {

        assert(dim_m && dim_n && dim_k);
        // create matrix_c to return
        Phase2_RNS_Int_Ptr matrix_c = FFLAS::fflas_new(*_phase2_rns_computation_field, dim_m, dim_k);

        assert(_phase2_rns_computation_field->size() == level_2_moduli_count * level_1_moduli_count);
        assert(matrix_a._stride == dim_m * dim_n);
        assert(matrix_b._stride == dim_n * dim_k);
        assert(matrix_c._stride == dim_m * dim_k);

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

#endif
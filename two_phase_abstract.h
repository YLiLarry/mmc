#ifndef H_TWO_PHASE_ABSTRACT
#define H_TWO_PHASE_ABSTRACT

#include "containers.h"
#include "gen_coprime_abstract.h"
#include "sim_rns.h"
#include <iostream>
#include <memory>
#include <gmp++/gmp++.h>
#include <fflas-ffpack/field/rns-double.h>
#include <fflas-ffpack/fflas/fflas_fgemm/fgemm_classical_mp.inl>
#include <vector>
#include "gen_coprime_abstract.h"
#include <givaro/givtimer.h>

class TwoPhaseAbstract
{
  protected:
    // Phase2_RNS_Field uses m_level_2_moduli repeats m_level_1_moduli_count times as modulis
    typedef FFPACK::rns_double Phase2_RNS_Rep;
    Phase2_RNS_Rep *m_phase2_rns_rep;
    // Phase2_RNS_Rep *m_phase2_rns_computation_rep;

    typedef FFPACK::RNSInteger<Phase2_RNS_Rep> Phase2_RNS_Field;
    typedef Phase2_RNS_Field::Element Phase2_RNS_Int;
    typedef Phase2_RNS_Field::Element_ptr Phase2_RNS_Int_Ptr;
    Phase2_RNS_Field *m_phase2_rns_field;
    // Phase2_RNS_Field *m_phase2_rns_computation_field;

    typedef Givaro::Modular<Givaro::Integer> Phase1_Field;
    typedef Phase1_Field::Element Phase1_Int;
    typedef Phase1_Field::Element_ptr Phase1_Int_Ptr;
    Phase1_Field *m_phase1_field;
    const GenCoprimeAbstract<Givaro::Integer> *m_level_1_moduli;
    const GenCoprimeAbstract<double> *m_level_2_moduli;
    const size_t m_level_1_moduli_count;
    const size_t m_level_2_moduli_count;

  public:
    TwoPhaseAbstract(const GenCoprimeAbstract<Givaro::Integer> *m_level_1_moduli,
                     const GenCoprimeAbstract<double> *m_level_2_moduli)
        : m_level_1_moduli(m_level_1_moduli),
          m_level_2_moduli(m_level_2_moduli),
          m_level_1_moduli_count(m_level_1_moduli->count()),
          m_level_2_moduli_count(m_level_2_moduli->count())
    {
#if DEBUG_MMC || TIME_MMC
        cerr << "########## TwoPhaseAlgo constructor ##########" << endl;
#endif
#if DEBUG_MMC || TIME_MMC
        cerr << " - m_level_1_moduli: " << *m_level_1_moduli << endl;
        cerr << " - m_level_2_moduli: " << *m_level_2_moduli << endl;
#endif
        assert(m_level_1_moduli->product() > m_level_2_moduli->product() && "chosen RNS size must make sense.");
        // m_level_2_moduli repeats m_level_1_moduli_count times
        std::vector<Givaro::Integer> tmp(m_level_2_moduli_count * m_level_1_moduli_count);
        for (size_t i = 0; i < m_level_1_moduli_count; i++)
        {
            for (size_t j = 0; j < m_level_2_moduli_count; j++)
            {
                tmp[i * m_level_2_moduli_count + j] = m_level_2_moduli->val(j);
            }
        }
        // m_phase2_rns_computation_rep = new Phase2_RNS_Rep{tmp};
        // m_phase2_rns_computation_field = new Phase2_RNS_Field(*m_phase2_rns_computation_rep);
        m_phase2_rns_rep = new Phase2_RNS_Rep{*m_level_2_moduli};
        m_phase2_rns_field = new Phase2_RNS_Field(*m_phase2_rns_rep);
        m_phase1_field = new Phase1_Field(m_level_2_moduli->product());
#if DEBUG_MMC
        cerr << " - m_phase2_rns_field: " << m_phase2_rns_field->rns()._basis << endl;
        // cerr << " - m_phase2_rns_computation_field: " << m_phase2_rns_computation_field->rns()._basis << endl;
#endif
#if DEBUG_MMC || TIME_MMC
        cerr << "########## TwoPhaseAlgo constructor ends ##########" << endl;
#endif
#if CHECK_MMCC
        for (size_t f = 0; f < m_level_1_moduli_count; f++)
        {
            if (m_level_1_moduli->val(f) >= m_level_2_moduli->product())
            {
                cerr << "level 1 moduli " << m_level_1_moduli->val(f) << " is not smaller than the product of level 2 modulis " << m_level_2_moduli->product() << endl;
            }
            assert(m_level_1_moduli->val(f) < m_level_2_moduli->product() && "level 2 moduli product must be large enough");
        }
#endif
    };

    virtual ~TwoPhaseAbstract()
    {
        delete m_phase1_field;
        delete m_phase2_rns_rep;
        delete m_phase2_rns_field;
        // delete m_phase2_rns_computation_rep;
        // delete m_phase2_rns_computation_field;
    };
    TwoPhaseAbstract(const TwoPhaseAbstract &) = delete;
    TwoPhaseAbstract &operator=(const TwoPhaseAbstract &) = delete;

  public:
    class Phase2_Matrix
    {
        // shared ptr will be deleted when no Phase2_Matrix holds the FFLAS_Mem
        // so that FFLAS::fflas_delete is called
        shared_ptr<FFLAS_Mem<Phase2_RNS_Field>> m_data;

      public:
        size_t dim_m;
        size_t dim_n;
        size_t count;
        size_t m_level_1_moduli_count;
        size_t m_level_2_moduli_count;
        inline Phase2_RNS_Int_Ptr &data() { return m_data->data; }
        inline const Phase2_RNS_Int_Ptr &data() const { return m_data->data; }

        Phase2_Matrix() = default;

        Phase2_Matrix(const TwoPhaseAbstract &f, size_t dim_m, size_t dim_n)
            : m_data(std::make_shared<FFLAS_Mem<Phase2_RNS_Field>>(FFLAS::fflas_new(*(f.m_phase2_rns_field), f.m_level_1_moduli_count * dim_m * dim_n))),
              dim_m(dim_m),
              dim_n(dim_n),
              count(dim_n * dim_m),
              m_level_1_moduli_count(f.m_level_1_moduli_count),
              m_level_2_moduli_count(f.m_level_2_moduli_count)
        {
            m_data->data._stride = this->count;
            assert(this->data()._stride == this->count);
        }

        Phase2_Matrix(const TwoPhaseAbstract &f,
                      const Phase2_RNS_Int_Ptr &arr,
                      size_t dim_m, size_t dim_n)
            : m_data(std::make_shared<FFLAS_Mem<Phase2_RNS_Field>>(arr)),
              dim_m(dim_m),
              dim_n(dim_n),
              count(dim_n * dim_m),
              m_level_1_moduli_count(f.m_level_1_moduli_count),
              m_level_2_moduli_count(f.m_level_2_moduli_count)
        {
            m_data->data._stride = this->count;
            assert(this->data()._stride == this->count);
        }

        Phase2_Matrix(const Phase2_Matrix &) = default;
        Phase2_Matrix &operator=(const Phase2_Matrix &) = default;

        const Phase2_RNS_Int ref(size_t r, size_t c, size_t f, size_t m) const
        {
            return data()[f * count * m_level_2_moduli_count + m * count + r * dim_n + c];
        }

        friend std::ostream &operator<<(std::ostream &out, const Phase2_Matrix &mat)
        {
            out << "Phase2_Matrix: " << std::endl;
            for (size_t r = 0; r < mat.dim_m; r++)
            {
                for (size_t f = 0; f < mat.m_level_1_moduli_count; f++)
                {
                    if (f != 0)
                    {
                        out << " ---";
                    }
                    for (size_t m = 0; m < mat.m_level_2_moduli_count; m++)
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
                out << std::endl;
            }
            return out;
        }
    };

  protected:
    /*
        this helper method is used by matrix_product(...)
    */
    virtual const std::vector<Phase1_Int> matrix_reduce_phase_1(const std::vector<Givaro::Integer> &inputs) const = 0;

  protected:
    /* 
        use this method to recover from a phase 1 representations to integers
    */
    virtual const std::vector<Givaro::Integer> matrix_recover_phase_1(const std::vector<Phase1_Int> &phase2_recovered) const = 0;

  public:
    /* 
        use this method to reduce a single matrix to level 2
    */
    Phase2_Matrix matrix_reduce(const std::vector<Givaro::Integer> &inputs, size_t dim_m, size_t dim_n) const
    {
        size_t len_inputs = inputs.size();
        assert(len_inputs == dim_m * dim_n && "input matrix dimension is incorrect");
        assert(dim_m > 0 && "input matrix dimension is incorrect");
        assert(dim_n > 0 && "input matrix dimension is incorrect");
#if DEBUG_MMC || TIME_MMC
        cerr << "########## matrix_reduce ##########" << endl;
#endif
#if DEBUG_MMC
        cerr << "inputs: " << endl
             << inputs << endl;
#endif
#if CHECK_MMCC
        for (size_t i = 0; i < len_inputs; i++)
        {
            assert(inputs[i] < m_level_1_moduli->product() && "Inputs must be less than the product of first level moduli.");
        }
#endif
#if DEBUG_MMC || TIME_MMC
        cerr << "..... phase 1 reduce ....." << endl;
#endif
#if TIME_MMC
        Givaro::Timer timer;
        timer.start();
#endif
        std::vector<Phase1_Int> p1_reduced = matrix_reduce_phase_1(inputs);
#if TIME_MMC
        timer.stop();
        cerr << "Timer: " << timer << endl;
#endif
#if DEBUG_MMC
        cerr << "phase 1 reduced: " << endl
             << p1_reduced << endl;
#endif
#if DEBUG_MMC || TIME_MMC
        cerr << "..... phase 1 reduce ends ....." << endl;
        cerr << "..... phase 2 reduce ....." << endl;
#endif
#if TIME_MMC
        timer.clear();
        timer.start();
#endif
        Phase2_RNS_Int_Ptr phase2_outputs = SIM_RNS::fflas_new_sim_reduce(*m_phase1_field, p1_reduced, *m_phase2_rns_field);
#if TIME_MMC
        timer.stop();
        cerr << "Timer: " << timer << endl;
#endif
#if DEBUG_MMC || TIME_MMC
        cerr << "..... phase 2 reduce ends ....." << endl;
#endif
        Phase2_Matrix mat(*this, dim_m, dim_n);
        for (size_t r = 0; r < dim_m; r++)
        {
            for (size_t c = 0; c < dim_n; c++)
            {
                size_t i = r * mat.dim_n + c;
                for (size_t f = 0; f < m_level_1_moduli_count; f++)
                {
                    for (size_t m = 0; m < m_level_2_moduli_count; m++)
                    {
                        mat.ref(r, c, f, m)._ptr[0] = phase2_outputs[m * (mat.count * m_level_1_moduli_count) + (i * m_level_1_moduli_count) + f]._ptr[0];
                    }
                }
            }
        }
        FFLAS::fflas_delete(phase2_outputs);
#if DEBUG_MMC
        cerr << "phase 2 reduced: " << endl
             << mat << endl;
#endif
#if DEBUG_MMC || TIME_MMC
        cerr << "########## matrix_reduce ends ##########" << endl;
#endif
        return mat;
    }

  public:
    /* 
        use this method to reduce multiple matrices to level 2
        shoule be faster than reducing one by one
    */
    const std::vector<Phase2_Matrix> matrix_reduce(const std::vector<Givaro::Integer> &matrices, const std::vector<size_t> &dimensions)
    {
        size_t len_inputs = matrices.size();
        size_t num_matrices = dimensions.size();
#if DEBUG_MMC || TIME_MMC
        cerr << "########## matrix_reduce ##########" << endl;
#endif
#if DEBUG_MMC
        cerr << "reducing " << len_inputs << " matrices" << endl;
#endif
        assert(len_inputs > 0 && num_matrices > 2 && "need at least 2 matrices");

#if CHECK_MMCC
        size_t len_ints = 0;
        for (size_t i = 1; i < num_matrices; i++)
        {
            size_t m = dimensions[i - 1];
            size_t n = dimensions[i];
            len_ints += m * n;
        }
        assert(len_inputs == len_ints && "supplied inputs and dimensions don't match");
#endif

#if DEBUG_MMC || TIME_MMC
        cerr << "..... phase 1 reduce ....." << endl;
#endif
#if TIME_MMC
        Givaro::Timer timer;
        timer.start();
#endif
        const std::vector<Phase1_Int> p1_reduced = matrix_reduce_phase_1(matrices);
#if TIME_MMC
        timer.stop();
        cerr << "Timer: " << timer << endl;
#endif
#if DEBUG_MMC || TIME_MMC
        cerr << "..... phase 1 reduce ends ....." << endl;
#endif
#if DEBUG_MMC
        cerr << "phase 1 reduced: " << endl
             << p1_reduced << endl;
#endif
#if DEBUG_MMC
        cerr << "..... phase 2 reduce ....." << endl;
#endif
#if TIME_MMC
        timer.clear();
        timer.start();
#endif
        Phase2_RNS_Int_Ptr phase2_outputs = SIM_RNS::fflas_new_sim_reduce(*m_phase1_field, p1_reduced, *m_phase2_rns_field);
#if TIME_MMC
        timer.stop();
        cerr << "Timer: " << timer << endl;
#endif
#if DEBUG_MMC || TIME_MMC
        cerr << "..... phase 2 reduce ends ....." << endl;
#endif
        std::vector<Phase2_Matrix> outputs(num_matrices);
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
                    for (size_t f = 0; f < m_level_1_moduli_count; f++)
                    {
                        for (size_t m = 0; m < m_level_2_moduli_count; m++)
                        {
                            mat.ref(r, c, f, m)._ptr[0] = phase2_outputs[m * (mat.count * m_level_1_moduli_count) + (i * m_level_1_moduli_count) + f]._ptr[0];
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
#if DEBUG_MMC || TIME_MMC
        cerr << "########## matrix_reduce ends ##########" << endl;
#endif
        return outputs;
    }

  public:
    /* 
        use this method to recover from a single reduced matrix to phase 1 representations
    */
    virtual std::vector<Phase1_Int> matrix_recover_phase_2(const Phase2_Matrix &mat) const
    {
        Phase2_RNS_Int_Ptr phase2_inputs = FFLAS::fflas_new(*m_phase2_rns_field, m_level_1_moduli_count * mat.count);
        for (size_t r = 0; r < mat.dim_m; r++)
        {
            for (size_t c = 0; c < mat.dim_n; c++)
            {
                size_t i = r * mat.dim_n + c;
                for (size_t f = 0; f < m_level_1_moduli_count; f++)
                {
                    for (size_t m = 0; m < m_level_2_moduli_count; m++)
                    {
                        // phase2_inputs[(m_level_2_moduli_count * m_level_1_moduli_count) * i + f * m_level_2_moduli_count + m]._ptr[0] = mat.ref(r, c, f, m)._ptr[0];
                        phase2_inputs[m * (mat.count * m_level_1_moduli_count) + (i * m_level_1_moduli_count) + f]._ptr[0] = mat.ref(r, c, f, m)._ptr[0];
                    }
                }
            }
        }

#if DEBUG_MMC || TIME_MMC
        cerr << ".......... fflas_new_sim_recover .........." << endl;
#endif
        // phase 2 recovery begins
        auto result = SIM_RNS::fflas_new_sim_recover(*m_phase2_rns_field, phase2_inputs, mat.count * m_level_1_moduli_count, *m_phase1_field);
#if DEBUG_MMC || TIME_MMC
        cerr << ".......... fflas_new_sim_recover ends .........." << endl;
#endif
        return result;
    }

    /* 
        use this method to recover from a single reduced matrix 
    */

  public:
    std::vector<Givaro::Integer> matrix_recover(const Phase2_Matrix &mat) const
    {
#if DEBUG_MMC || TIME_MMC
        cerr << "########## matrix_recover ##########" << endl;
#endif
#if DEBUG_MMC
        cerr << mat << endl;
#endif
#if DEBUG_MMC || TIME_MMC
        cerr << "..... phase 2 recovery ....." << endl;
#endif
#if TIME_MMC
        Givaro::Timer timer;
        timer.start();
#endif
        const std::vector<Phase1_Int> phase2_recovered = matrix_recover_phase_2(mat);
#if TIME_MMC
        timer.stop();
        cerr << "Timer: " << timer << endl;
#endif
#if DEBUG_MMC
        cerr << phase2_recovered << endl;
#endif
#if DEBUG_MMC || TIME_MMC
        cerr << "..... phase 2 recovery ends ....." << endl;
        cerr << "..... phase 1 recovery ....." << endl;
#endif
#if TIME_MMC
        timer.clear();
        timer.start();
#endif
        const std::vector<Givaro::Integer> phase1_recovered = matrix_recover_phase_1(phase2_recovered);
#if TIME_MMC
        timer.stop();
        cerr << "Timer: " << timer << endl;
#endif
#if DEBUG_MMC || TIME_MMC
        cerr << "..... phase 1 recovery ends ....." << endl;
#endif
#if DEBUG_MMC
        cerr << phase1_recovered << endl;
#endif
#if DEBUG_MMC || TIME_MMC
        cerr << "########## matrix_recover ends ##########" << endl;
#endif
        return phase1_recovered;
    }

  public:
    /* 
        use this method to multiply two reduced matrices
    */
    Phase2_Matrix phase2_mult(const Phase2_Matrix &matrix_a, const Phase2_Matrix &matrix_b) const
    {
#if DEBUG_MMC || TIME_MMC
        cerr << "########## phase2_mult ##########" << endl;
#endif
#if DEBUG_MMC
        cerr << " - matrix_a: " << endl;
        cerr << matrix_a << endl;
        cerr << " - matrix_b: " << endl;
        cerr << matrix_b << endl;
#endif
        assert(matrix_a.dim_n == matrix_b.dim_m);
        Phase2_Matrix matrix_c = phase2_matrix_fgemm(matrix_a.data(), matrix_b.data(), matrix_a.dim_m, matrix_a.dim_n, matrix_b.dim_m);
#if DEBUG_MMC
        cerr << " - matrix product:" << endl;
        cerr << matrix_c;
#endif
#if DEBUG_MMC || TIME_MMC
        cerr << "########## phase2_mult ends ##########" << endl;
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
        Phase2_RNS_Int_Ptr matrix_c = FFLAS::fflas_new(*m_phase2_rns_field, m_level_1_moduli_count * dim_m * dim_k);
        matrix_c._stride = dim_m * dim_k;

        // assert(m_phase2_rns_field->size() == m_level_2_moduli_count * m_level_1_moduli_count);
        assert(matrix_a._stride == dim_m * dim_n);
        assert(matrix_b._stride == dim_n * dim_k);
        assert(matrix_c._stride == dim_m * dim_k);

#if DEBUG_MMC || TIME_MMC
        cerr << ".......... fgemm .........." << endl;
#endif
        FFLAS::MMHelper<FFPACK::RNSInteger<Phase2_RNS_Rep>,
                        FFLAS::MMHelperAlgo::Classic,
                        FFLAS::ModeCategories::DefaultTag,
                        FFLAS::ParSeqHelper::Sequential>
            tag;
        fgemm<Phase2_RNS_Rep>(
            (Phase2_RNS_Field)*m_phase2_rns_field, // field
            FFLAS::FflasNoTrans,                   // transpose matrix_a?
            FFLAS::FflasNoTrans,                   // transpose matrix_b?
            dim_m, dim_n, dim_k,
            m_phase2_rns_field->one, // coefficient before matrix_a*matrix_b
            matrix_a,
            dim_n, // row length of matrix_a
            matrix_b,
            dim_k,                    // row length of matrix_b
            m_phase2_rns_field->zero, // constant matrix to add
            matrix_c,
            dim_k, // row length of matrix_c
            tag);
#if DEBUG_MMC || TIME_MMC
        cerr << ".......... fgemm ends .........." << endl;
#endif
        return matrix_c;
    }

    // fgemm for RnsInteger sequential version
    template <typename RNS>
    inline typename FFPACK::RNSInteger<RNS>::Element_ptr
    fgemm(const FFPACK::RNSInteger<RNS> &F,
          const FFLAS::FFLAS_TRANSPOSE ta,
          const FFLAS::FFLAS_TRANSPOSE tb,
          const size_t dim_m, const size_t dim_n, const size_t dim_k,
          const typename FFPACK::RNSInteger<RNS>::Element alpha,
          typename FFPACK::RNSInteger<RNS>::ConstElement_ptr Ad, const size_t lda,
          typename FFPACK::RNSInteger<RNS>::ConstElement_ptr Bd, const size_t ldb,
          const typename FFPACK::RNSInteger<RNS>::Element beta,
          typename FFPACK::RNSInteger<RNS>::Element_ptr Cd, const size_t ldc,
          FFLAS::MMHelper<FFPACK::RNSInteger<RNS>, FFLAS::MMHelperAlgo::Classic, FFLAS::ModeCategories::DefaultTag, FFLAS::ParSeqHelper::Sequential> &H) const
    {
        // compute each fgemm componentwise
#ifdef PROFILE_FGEMM_MP
        Givaro::Timer t;
        t.start();
#endif
        for (size_t f = 0; f < m_level_1_moduli_count; f++)
        {
            for (size_t m = 0; m < F.size(); m++)
            {
                size_t i = f * F.size() + m;
                auto field = F.rns()._field_rns[m];
                FFLAS::MMHelper<typename RNS::ModField, FFLAS::MMHelperAlgo::Winograd> H2(field, H.recLevel, H.parseq);
#if CHECK_MMC
                assert(dim_m * dim_n == Ad._stride);
                assert(dim_n * dim_k == Bd._stride);
                assert(dim_m * dim_k == Cd._stride);
#endif
                FFLAS::fgemm(field, ta, tb,
                             dim_m, dim_n, dim_k,
                             alpha._ptr[m * alpha._stride],
                             Ad._ptr + i * Ad._stride, lda,
                             Bd._ptr + i * Bd._stride, ldb,
                             beta._ptr[m * beta._stride],
                             Cd._ptr + i * Cd._stride, ldc, H2);
            }
        }
#ifdef PROFILE_FGEMM_MP
        t.stop();

        std::cerr << "==========================================" << std::endl
                  << "Pointwise fgemm : " << t.realtime() << " (" << F.size() << ") moduli " << std::endl
                  << "==========================================" << std::endl;
#endif
        return Cd;
    }
};

#endif
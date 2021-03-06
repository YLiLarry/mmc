
using namespace std;

#include "nocopy_integer.h"
#include "sim_rns.h"
#include "two_phase_marge_most.h"
#include "two_phase_parge_shift.h"
#include "two_phase_parge_block.h"
#include <array>
#include <gmp++/gmp++.h>
#include <ostream>
#include <time.h>

#define TEST_MARGE_MOST 1
#define TEST_MARGE_LEAST 0
#define TEST_PARGE_BLOCK 0
#define TEST_PARGE_SHIFT 0

using namespace LinBox;
using namespace SIM_RNS;

void test()
{
    const uint_fast64_t input_bitsize = (1 << 6);

    vector<Givaro::Integer> a(4);
    a[0] = LInteger::random_exact(input_bitsize);
    a[1] = LInteger::random_exact(input_bitsize);
    a[2] = LInteger::random_exact(input_bitsize);
    a[3] = LInteger::random_exact(input_bitsize);
    vector<Givaro::Integer> b(4);
    b[0] = LInteger::random_exact(input_bitsize);
    b[1] = LInteger::random_exact(input_bitsize);
    b[2] = LInteger::random_exact(input_bitsize);
    b[3] = LInteger::random_exact(input_bitsize);
    cerr << "input a: " << a << endl;
    cerr << "input b: " << b << endl;

    clock_t _time;
    vector<Givaro::Integer> expect = SIM_RNS::fflas_mult_integer(a, b, 2, 2, 2);

#if TEST_MARGE_MOST
    cerr << "===========================================" << endl;
    cerr << "======== Testing TwoPhaseMargeMost ========" << endl;
    cerr << "===========================================" << endl;
    {
        TwoPhaseMargeMost algo(2 * input_bitsize, input_bitsize / 2);

        auto r = algo.matrix_reduce(a, 2, 2);
        vector<Givaro::Integer> a_ = algo.matrix_recover(r);
        assert(equals(a, a_));

        auto s = algo.matrix_reduce(b, 2, 2);
        vector<Givaro::Integer> s_ = algo.matrix_recover(s);
        assert(equals(b, s_));

        vector<Givaro::Integer> mm;
        mm.insert(mm.end(), a.begin(), a.end());
        mm.insert(mm.end(), b.begin(), b.end());
        assert(mm.size() == 8);
        vector<size_t> ss;
        ss.push_back(2);
        ss.push_back(2);
        ss.push_back(2);
        auto matrices = algo.matrix_reduce(mm, ss);
        cerr << matrices[0] << endl;
        cerr << r << endl;
        auto mm1 = algo.matrix_recover(matrices[0]);
        auto mm2 = algo.matrix_recover(matrices[1]);
        assert(equals(mm1, a));
        assert(equals(mm2, b));

        auto t = algo.phase2_mult(r, s);
        auto got = algo.matrix_recover(t);

        if (!equals(got, expect))
        {
            cerr << "TwoPhaseMargeMost failed" << endl
                 << " - expect: " << expect << endl
                 << " - got: " << got << endl;
            abort();
        }
        cerr << "TwoPhaseMargeMost passed!" << endl;
    }
#endif

#if TEST_PARGE_BLOCK
    cerr << "===========================================" << endl;
    cerr << "======= Testing TwoPhasePargeBlock ========" << endl;
    cerr << "===========================================" << endl;
    {
        TwoPhasePargeBlock algo_parge_block(2 * input_bitsize, input_bitsize / 2, 4);

        auto r = algo_parge_block.matrix_reduce(a, 2, 2);
        vector<Givaro::Integer> a_ = algo_parge_block.matrix_recover(r);
        assert(equals(a, a_));

        auto s = algo_parge_block.matrix_reduce(b, 2, 2);
        vector<Givaro::Integer> b_ = algo_parge_block.matrix_recover(s);
        assert(equals(b, b_));

        auto t = algo_parge_block.phase2_mult(r, s);
        auto got = algo_parge_block.matrix_recover(t);

        if (!equals(got, expect))
        {
            cerr << "TwoPhasePargeBlock failed" << endl
                 << " - expect: " << expect << endl
                 << " - got: " << got << endl;
            abort();
        }
        cerr << "TwoPhasePargeBlock passed!" << endl;
    }
#endif

#if TEST_PARGE_SHIFT
    cerr << "===========================================" << endl;
    cerr << "======= Testing TwoPhasePargeShift ========" << endl;
    cerr << "===========================================" << endl;
    {
        TwoPhasePargeShift algo_parge_shift(2 * input_bitsize, input_bitsize / 2, 1);

        auto r = algo_parge_shift.matrix_reduce(a, 2, 2);
        vector<Givaro::Integer> a_ = algo_parge_shift.matrix_recover(r);
        assert(equals(a, a_));

        auto s = algo_parge_shift.matrix_reduce(b, 2, 2);
        vector<Givaro::Integer> b_ = algo_parge_shift.matrix_recover(s);
        assert(equals(b, b_));

        auto t = algo_parge_shift.phase2_mult(r, s);
        auto got = algo_parge_shift.matrix_recover(t);

        if (!equals(got, expect))
        {
            cerr << "TwoPhasePargeShift failed" << endl
                 << " - expect: " << expect << endl
                 << " - got: " << got << endl;
            abort();
        }
        cerr << "TwoPhasePargeShift passed!" << endl;
    }
#endif

    cerr << "All tests passed!" << endl;
}

int main()
{
    for (int i = 0; i < 100; i++)
    {
        test();
    }
}

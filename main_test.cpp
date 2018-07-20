
using namespace std;

#include "nocopy_integer.h"
#include "sim_rns.h"
#include "two_phase_marge.h"
#include "two_phase_parge_shift.h"
#include "two_phase_parge_block.h"
#include <array>
#include <gmp++/gmp++.h>
#include <ostream>
#include <time.h>

using namespace LinBox;
using namespace SIM_RNS;

void test(bool test_marge, bool test_parge, bool test_fermat)
{
    const uint_fast64_t input_bitsize = (1 << 10);
    const uint_fast64_t output_bitsize = input_bitsize * 2;

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

    if (!test_marge)
    {
        goto end_test_marge;
    }

    cerr << "===========================================" << endl;
    cerr << "========== Testing TwoPhaseMarge ==========" << endl;
    cerr << "===========================================" << endl;
    {
        const uint_fast64_t phase1_moduli_bitsize = input_bitsize / 2;
        const size_t phase2_moduli_bitsize = 21;

        TwoPhaseMarge algo_marge(input_bitsize, phase1_moduli_bitsize, phase2_moduli_bitsize);

        auto r = algo_marge.matrix_reduce(a, 2, 2);
        vector<Givaro::Integer> r_ = algo_marge.matrix_recover(r);
        assert(equals(a, r_));

        auto s = algo_marge.matrix_reduce(b, 2, 2);
        vector<Givaro::Integer> s_ = algo_marge.matrix_recover(s);
        assert(equals(b, s_));

        auto t = algo_marge.phase2_mult(r, s);
        auto got = algo_marge.matrix_recover(t);

        if (!equals(got, expect))
        {
            cerr << "TwoPhaseMarge failed" << endl
                 << " - expect: " << expect << endl
                 << " - got: " << got << endl;
            abort();
        }
        cerr << "TwoPhaseMarge passed!" << endl;
    }

end_test_marge:

    if (!test_parge)
    {
        goto end_test_parge;
    }
    cerr << "===========================================" << endl;
    cerr << "========= Testing TwoPhasePargeBlock ==========" << endl;
    cerr << "===========================================" << endl;
    {
        const uint_fast64_t phase1_moduli_bitsize = input_bitsize;
        const size_t phase2_moduli_bitsize = 21;

        TwoPhasePargeBlock algo_parge_block(output_bitsize, phase1_moduli_bitsize, phase2_moduli_bitsize);

        auto r = algo_parge_block.matrix_reduce(a, 2, 2);
        vector<Givaro::Integer> r_ = algo_parge_block.matrix_recover(r);
        assert(equals(a, r_));

        auto s = algo_parge_block.matrix_reduce(b, 2, 2);
        vector<Givaro::Integer> s_ = algo_parge_block.matrix_recover(s);
        assert(equals(b, s_));

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

end_test_parge:

    if (!test_fermat)
    {
        goto end_test_fermat;
    }
    cerr << "===========================================" << endl;
    cerr << "========= Testing TwoPhasePargeShift ==========" << endl;
    cerr << "===========================================" << endl;
    {
        const uint_fast64_t phase1_moduli_bitsize = 1024;
        const size_t phase2_moduli_bitsize = 21;

        TwoPhasePargeShift algo_parge_shift(output_bitsize, phase1_moduli_bitsize, 1, phase2_moduli_bitsize);

        auto r = algo_parge_shift.matrix_reduce(a, 2, 2);
        vector<Givaro::Integer> r_ = algo_parge_shift.matrix_recover(r);
        assert(equals(a, r_));

        auto s = algo_parge_shift.matrix_reduce(b, 2, 2);
        vector<Givaro::Integer> s_ = algo_parge_shift.matrix_recover(s);
        assert(equals(b, s_));

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

end_test_fermat:

    cerr << "All tests passed!" << endl;
}

int main()
{
    for (int i = 0; i < 3; i++)
    {
        test(false, false, true);
    }
}

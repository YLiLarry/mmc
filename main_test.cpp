
using namespace std;

#include "nocopy_integer.h"
#include "sim_rns.h"
#include "two_phase_marge.h"
#include "two_phase_fermat.h"
#include "two_phase_parge.h"
#include <array>
#include <gmp++/gmp++.h>
#include <ostream>
#include <time.h>

using namespace LinBox;
using namespace SIM_RNS;

void test(bool test_marge, bool test_parge, bool test_fermat)
{
    const uint_fast64_t input_bit_length = (1 << 12);

    vector<Givaro::Integer> a(4);
    a[0] = LInteger::random_exact(input_bit_length);
    a[1] = LInteger::random_exact(input_bit_length);
    a[2] = LInteger::random_exact(input_bit_length);
    a[3] = LInteger::random_exact(input_bit_length);
    vector<Givaro::Integer> b(4);
    b[0] = LInteger::random_exact(input_bit_length);
    b[1] = LInteger::random_exact(input_bit_length);
    b[2] = LInteger::random_exact(input_bit_length);
    b[3] = LInteger::random_exact(input_bit_length);
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
        const uint_fast64_t phase1_moduli_bit_length = 512;
        const size_t phase2_moduli_bit_length = 21;

        TwoPhaseMarge algo_marge(input_bit_length, phase1_moduli_bit_length, phase2_moduli_bit_length);
        auto r = algo_marge.matrix_reduce(a, 2, 2);
        auto s = algo_marge.matrix_reduce(b, 2, 2);
        auto t = algo_marge.phase2_mult(r, s);
        auto got = algo_marge.matrix_recover(t);

        vector<Givaro::Integer> r_ = algo_marge.matrix_recover(r);
        vector<Givaro::Integer> s_ = algo_marge.matrix_recover(s);
        assert(equals(a, r_));
        assert(equals(b, s_));

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
    cerr << "========= Testing TwoPhaseParge ==========" << endl;
    cerr << "===========================================" << endl;
    {
        const uint_fast64_t phase1_moduli_bit_length = (1 << 13);
        const size_t phase2_moduli_bit_length = 21;

        TwoPhaseParge algo_parge(input_bit_length, phase1_moduli_bit_length, phase2_moduli_bit_length);
        _time = clock();
        auto r = algo_parge.matrix_reduce(a, 2, 2);
        auto s = algo_parge.matrix_reduce(b, 2, 2);
        auto t = algo_parge.phase2_mult(r, s);
        auto got = algo_parge.matrix_recover(t);

        vector<Givaro::Integer> r_ = algo_parge.matrix_recover(r);
        vector<Givaro::Integer> s_ = algo_parge.matrix_recover(s);
        assert(equals(a, r_));
        assert(equals(b, s_));

        if (!equals(got, expect))
        {
            cerr << "TwoPhaseParge failed" << endl
                 << " - expect: " << expect << endl
                 << " - got: " << got << endl;
            abort();
        }
        cerr << "TwoPhaseParge passed!" << endl;
    }

end_test_parge:

    if (!test_fermat)
    {
        goto end_test_fermat;
    }
    cerr << "===========================================" << endl;
    cerr << "========= Testing TwoPhaseFermat ==========" << endl;
    cerr << "===========================================" << endl;
    {
        const uint_fast64_t phase1_moduli_bit_length = (1 << 13);
        const size_t phase2_moduli_bit_length = 21;

        TwoPhaseFermat algo_fermat(input_bit_length, phase1_moduli_bit_length, phase2_moduli_bit_length);
        _time = clock();
        auto r = algo_fermat.matrix_reduce(a, 2, 2);
        auto s = algo_fermat.matrix_reduce(b, 2, 2);
        auto t = algo_fermat.phase2_mult(r, s);
        auto got = algo_fermat.matrix_recover(t);

        vector<Givaro::Integer> r_ = algo_fermat.matrix_recover(r);
        vector<Givaro::Integer> s_ = algo_fermat.matrix_recover(s);
        assert(equals(a, r_));
        assert(equals(b, s_));

        if (!equals(got, expect))
        {
            cerr << "TwoPhaseFermat failed" << endl
                 << " - expect: " << expect << endl
                 << " - got: " << got << endl;
            abort();
        }
        cerr << "TwoPhaseFermat passed!" << endl;
    }

end_test_fermat:

    cerr << "All tests passed!" << endl;
}

int main()
{
    for (int i = 0; i < 3; i++)
    {
        test(true, true, true);
    }
}

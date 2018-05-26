#include "coprime_gen.h"
#include "nocopy_integer.h"
#include "prime_gen.h"
#include "sim_rns.h"
#include "two_phase.h"
#include <array>
#include <linbox/integer.h>
#include <ostream>

using namespace std;
using namespace LinBox;
using namespace SIM_RNS;

void test3()
{
    const size_t phase1_moduli_count = 2;
    const uint_fast64_t phase1_moduli_bit_length = 32;

    const size_t phase2_moduli_basis_size = 2;
    const size_t phase2_moduli_bit_length = 26;

    typedef LInteger Phase1_moduli_type;
    typedef uint_fast64_t Phase2_moduli_type;

    TwoPhaseAlgo<
        Phase1_moduli_type,
        phase1_moduli_count,
        phase1_moduli_bit_length,
        Phase2_moduli_type,
        phase2_moduli_basis_size,
        phase2_moduli_bit_length>
        algo;
    LIntegerPtrVector *a = LIntegerPtrVector::new_zeros(4);
    a->val(0) = LInteger(57307);
    a->val(1) = LInteger(56520);
    a->val(2) = LInteger(35421);
    a->val(3) = LInteger(37839);
    LIntegerPtrVector *b = LIntegerPtrVector::new_zeros(4);
    b->val(0) = LInteger(58568);
    b->val(1) = LInteger(45133);
    b->val(2) = LInteger(51263);
    b->val(3) = LInteger(55291);
    cerr << "input a: " << *a << endl;
    cerr << "input b: " << *b << endl;
    auto r = algo.matrix_reduce(*a, 2, 2);
    auto s = algo.matrix_reduce(*b, 2, 2);
    auto r_ = algo.matrix_recover(r);
    auto s_ = algo.matrix_recover(s);
    assert(r_->equals(*reinterpret_cast<NumPtrVector<Givaro::Integer> *>(a)));
    assert(s_->equals(*reinterpret_cast<NumPtrVector<Givaro::Integer> *>(b)));
    delete r_;
    delete s_;
    auto t = algo.phase2_mult(r, s);
    auto p = algo.matrix_recover(t);
    vector<Givaro::Integer> expect{
        static_cast<uint64_t>(6253741136),
        static_cast<uint64_t>(5711484151),
        static_cast<uint64_t>(4014277785),
        static_cast<uint64_t>(3690812142)};
    if (!p->equals(expect))
    {
        cerr << "test failed" << endl
             << " - expect: " << expect << endl
             << " - got: " << *p << endl;
        abort();
    }
    cerr << "Passed!" << endl;
}

int main()
{
    for (int i = 0; i < 100; i++)
    {
        test3();
        // test2();
    }
}

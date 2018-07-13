#include "marge_gen.h"
#include "nocopy_integer.h"
#include "prime_gen.h"
#include "sim_rns.h"
#include "two_phase.h"
#include <array>
#include <linbox/integer.h>
#include <ostream>
#include <time.h>
#include <linbox/matrix/matrix-domain.h>
#include <linbox/algorithms/matrix-blas3/mul.h>
// #include <istringstream>

using namespace std;
using namespace LinBox;
using namespace SIM_RNS;

void test5()
{
    const uint_fast64_t input_bit_length = 100;
    const uint_fast64_t phase1_moduli_bit_length = 100;
    const size_t phase2_moduli_bit_length = 26;

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

    clock_t _time = clock();
    TwoPhaseMarge algo(input_bit_length, phase1_moduli_bit_length, phase2_moduli_bit_length);
    auto r = algo.matrix_reduce(a, 2, 2);
    auto s = algo.matrix_reduce(b, 2, 2);
    auto t = algo.phase2_mult(r, s);
    auto got = algo.matrix_recover(t);
    _time = clock() - _time;
    cerr << "It took me " << _time << " clicks (" << static_cast<float>(_time) / CLOCKS_PER_SEC << " seconds)" << endl;

    Givaro::IntegerDomain id;
    _time = clock();
    vector<Givaro::Integer> expect = SIM_RNS::fflas_mult_integer(a, b, 2, 2, 2);
    _time = clock() - _time;
    cerr << "It took me " << _time << " clicks (" << static_cast<float>(_time) / CLOCKS_PER_SEC << " seconds)" << endl;

    vector<Givaro::Integer> r_ = algo.matrix_recover(r);
    vector<Givaro::Integer> s_ = algo.matrix_recover(s);
    assert(equals(a, r_));
    assert(equals(b, s_));

    if (!equals(got, expect))
    {
        cerr << "test failed" << endl
             << " - expect: " << expect << endl
             << " - got: " << got << endl;
        abort();
    }
    cerr << "Passed!" << endl;
}
/*
void test4()
{
    const size_t phase1_moduli_count = 3;
    const uint_fast64_t phase1_moduli_bit_length = 64;

    const size_t phase2_moduli_basis_size = 5;
    const size_t phase2_moduli_bit_length = 26;

    typedef LInteger Phase1_moduli_type;
    typedef uint_fast64_t Phase2_moduli_type;

    TwoPhaseMarge<
        Phase1_moduli_type,
        phase1_moduli_count,
        phase1_moduli_bit_length,
        Phase2_moduli_type,
        phase2_moduli_basis_size,
        phase2_moduli_bit_length>
        algo;
    // stringstream ss;
    LIntegerPtrVector *a = LIntegerPtrVector::new_zeros(4);
    a->val(0) = LInteger(1);
    a->ptr(0)->operator<<=(100);
    a->val(1) = LInteger(1);
    a->ptr(1)->operator<<=(99);
    a->val(2) = LInteger(1);
    a->ptr(2)->operator<<=(98);
    a->val(3) = LInteger(1);
    a->ptr(3)->operator<<=(97);
    LIntegerPtrVector *b = LIntegerPtrVector::new_zeros(4);
    b->val(0) = LInteger(4294967295.0);
    b->val(1) = LInteger(4294967295.0);
    b->val(2) = LInteger(4294967295.0);
    b->val(3) = LInteger(4294967295.0);
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

void test3()
{
    const size_t phase1_moduli_count = 2;
    const uint_fast64_t phase1_moduli_bit_length = 32;

    const size_t phase2_moduli_basis_size = 2;
    const size_t phase2_moduli_bit_length = 26;

    typedef LInteger Phase1_moduli_type;
    typedef uint_fast64_t Phase2_moduli_type;

    TwoPhaseMarge<
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
    clock_t _time = clock();
    auto t = algo.phase2_mult(r, s);
    auto p = algo.matrix_recover(t);
    _time = clock() - _time;
    cerr << "It took me " << _time << " clicks (" << static_cast<float>(_time) / CLOCKS_PER_SEC << " seconds)" << endl;
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

    vector<Givaro::Integer> *_a = a->EXPENSIVE_NEW_INTEGER_VECTOR();
    vector<Givaro::Integer> *_b = b->EXPENSIVE_NEW_INTEGER_VECTOR();
    Givaro::IntegerDomain id;
    _time = clock();
    SIM_RNS::fflas_mult_integer(*_a, *_b, 2, 2, 2);
    _time = clock() - _time;
    cerr << "It took me " << _time << " clicks (" << static_cast<float>(_time) / CLOCKS_PER_SEC << " seconds)" << endl;
}

*/
int main()
{
    for (int i = 0; i < 4; i++)
    {
        test5();
        // test2();
    }
}

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

void test2()
{
    const size_t input_count = 100;
    const size_t input_bitsize = 50;
    const size_t moduli_count = 8;
    const size_t moduli_bitsize = 10;
    const PrimeGenExact<LInteger, input_count, input_bitsize> pg;
    const PrimeGenExact<LInteger, moduli_count, moduli_bitsize> moduli;
    const NumPtrVector<LInteger>* inputs = pg.EXPENSIVE_NEW_NUM_PTR_VECTOR();
    auto vec = new_sim_reduce(*inputs, moduli);
    auto vec2 = new_sim_recover<LInteger, LInteger, moduli_count>(*vec);
    cerr << *inputs << endl;
    cerr << *vec2 << endl;
    assert(inputs->equals(*vec2) && "test2 failed");
    delete vec;
    delete inputs;
    cerr << "test2 passed" << endl;
}

void test()
{
    const size_t phase1_moduli_count = 2;
    const uint_fast64_t phase1_moduli_bit_length = (1 << 10);

    const size_t phase2_moduli_basis_size = 100;
    const size_t phase2_moduli_bit_length = 15;

    const size_t input_bit_length = (1 << 5);

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
    Phase1_moduli_type a;
    Phase1_moduli_type b;
    Phase1_moduli_type d;
    a.randomizePrime(input_bit_length);
    b.randomizePrime(input_bit_length);

    algo.mult(d, a, b);

    // cout << a << endl;
    // cout << (a * b) << endl;
    // NoCopyInteger d;
    // algo.mult(d, a, b);

    // PrimeGenExact<NoCopyInteger, basis_size, 26> basis;
    // cout << basis << endl;
    // vector<NoCopyInteger*> inputs;
    // inputs.push_back(NoCopyInteger::newRandom(4096));
    // inputs.push_back(NoCopyInteger::newRandom(4096));
    // inputs.push_back(NoCopyInteger::newRandom(4096));
    // cout << SIM_RNS::RNS<NoCopyInteger, NoCopyInteger, NoCopyInteger, basis_size>::naive_reduce(inputs, basis) << endl;

    // inputs.erase(inputs.begin());
    // MargeGenExact<NoCopyInteger, 1, 26> a;

    // PrimeGen<NoCopyInteger, 10, 10> a([](size_t i) {
    //     cout << i << endl;
    //     return new NoCopyInteger(i);
    // });
    // cout
    // << a << endl
    // << a.product << endl;
}

int main()
{
    for (int i = 0; i < 100; i++) {
        // test();
        test2();
    }
}

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
int main()
{
    const size_t phase1_modulis_size = 2;
    const uint_fast64_t phase1_modulis_bit_length = (1 << 10);

    const size_t phase2_modulis_basis_size = 100;
    const size_t phase2_modulis_bit_length = 15;

    const size_t input_bit_length = (1 << 10);

    typedef LInteger Phase1_modulis_type;
    typedef uint_fast64_t Phase2_modulis_type;

    TwoPhaseAlgo<
        Phase1_modulis_type,
        phase1_modulis_size,
        phase1_modulis_bit_length,
        Phase2_modulis_type,
        phase2_modulis_basis_size,
        phase2_modulis_bit_length>
        algo;
    Phase1_modulis_type a;
    Phase1_modulis_type b;
    Phase1_modulis_type d;
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

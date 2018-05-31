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
    // TwoPhaseAlgo<NoCopyInteger, 11, 1024, NoCopyInteger, 1000, 26> algo;
    // NoCopyInteger a;
    // NoCopyInteger b;
    // a.random(10240);
    // b.random(10240);
    // cout << a << endl;
    // cout << b << endl;
    // NoCopyInteger d;
    // algo.mult(d, a, b);
    const size_t basis_size = 10;

    PrimeGenExact<NoCopyInteger, basis_size, 26> basis;
    cout << basis << endl;
    vector<NoCopyInteger*> inputs;
    inputs.push_back(NoCopyInteger::newRandom(4096));
    inputs.push_back(NoCopyInteger::newRandom(4096));
    inputs.push_back(NoCopyInteger::newRandom(4096));
    cout << SIM_RNS::RNS<NoCopyInteger, NoCopyInteger, NoCopyInteger, basis_size>::naive_reduce(inputs, basis) << endl;
    inputs.erase(inputs.begin());
    // MargeGenExact<NoCopyInteger, 1, 26> a;

    // PrimeGen<NoCopyInteger, 10, 10> a([](size_t i) {
    //     cout << i << endl;
    //     return new NoCopyInteger(i);
    // });
    // cout
    // << a << endl
    // << a.product << endl;
}

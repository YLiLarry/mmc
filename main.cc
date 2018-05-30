#include "coprime_gen.h"
#include "nocopy_integer.h"
#include "prime_gen.h"
#include "two_phase.h"
#include <array>
#include <linbox/integer.h>
#include <ostream>

using namespace std;
using namespace LinBox;
int main()
{
    // TwoPhaseAlgo<NoCopyInteger, 1, 14, int, 1000, 26> algo;
    // NoCopyInteger a;
    // NoCopyInteger b;
    // a.random(1000);
    // b.random(1000);
    // NoCopyInteger d;
    // algo.mult(d, a, b);
    // PrimeGenExact<NoCopyInteger, 1000, 26> a;
    MargeGenExact<NoCopyInteger, 1, 26> a;
    cout << a << endl
         << a.product << endl;
}

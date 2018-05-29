#include <array>
#include <ostream>
#include "prime_gen.h"

using namespace std;

int main() {
    PrimeGen<int, 10, 10> p{};
    // array<int,10> a = prime_gen<int,10,10>();
    cout << p << endl;
}

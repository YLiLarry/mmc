#include "marge_num.h"
#include "../nocopy_integer.h"
#include <givaro/modular-integer.h>
#include <linbox/integer.h>
#include <ostream>

using namespace LinBox;
using namespace Givaro;
using namespace std;

// void dc_reduce_minus(LInteger& a, const LInteger& moduli)
// {
//     uint_fast64_t n = (moduli+1).bitsize();
//     uint_fast64_t s = a.bitsize();
//     uint_fast64_t b;
//     uint_fast64_t bn;
//     while (s > n) {
//         // cerr << "|" << a.bitsize() << " " << n << "|";
//         b = (s - 1) / n + 1;
//         b = b >> 1; // 2b = k
//         bn = b * n;
//         LInteger t(a, bn, bn << 1); // t = a[bn, 2bn)
//         a.cut(0, bn); // a = a[0, bn)
//         a += t;
//         s = a.bitsize();
//     }
// }

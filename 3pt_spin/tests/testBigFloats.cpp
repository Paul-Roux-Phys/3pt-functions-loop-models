#include "../include/BigFloats.hpp"

using namespace mp;
using namespace mpc;

const unsigned int prec = 1024;

int main(void) {
  mpcomplex<prec> z1(1.0, 2.0);  // 1 + 2i
  mpcomplex<prec> z2(3.0, -1.0); // 3 - i
  mpcomplex<prec> z3(0.0, 0.0);
  z3 += z1;      // add
  z3 *= z2 * z2; // multiply
  cos(3.0 * z3) - sin(2.0 * z3); // cos and sin functions
  exp(2.4 * z1 * z2);            // exp function

  mpfloat<prec> r1(1.34);    // real float
  mpfloat<prec> pi = mp::pi<prec>(); // pi
  mpcomplex<prec> im(0, 1);

  z1 = r1 * z3;
  z2 = pi * r1;

  exp(im * pi / 3) * r1;

  return 0;
}

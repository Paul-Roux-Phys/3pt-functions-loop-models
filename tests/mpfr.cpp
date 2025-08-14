#include <mpfr.h>

int main () {
  mpfr_t a;
  mpfr_init_set_d(a, 1.0, MPFR_RNDN);
  mpfr_clear(a);
  return 0;
}

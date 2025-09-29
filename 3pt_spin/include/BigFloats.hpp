#pragma once
#include <iostream>
#include <mpc.h>
#include <mpfr.h>

namespace mp {

// -------------------------
// mpfloat class
// -------------------------
template <unsigned int prec> class mpfloat {
  mpfr_t val;

public:
  // Constructors / Destructor
  explicit mpfloat() {
    mpfr_init2(val, prec);
    mpfr_set_d(val, 0.0, MPFR_RNDN);
  }

  mpfloat(double x) {
    mpfr_init2(val, prec);
    mpfr_set_d(val, x, MPFR_RNDN);
  }

  mpfloat(const mpfloat &other) {
    mpfr_init2(val, prec);
    mpfr_set(val, other.val, MPFR_RNDN);
  }

  mpfloat &operator=(const mpfloat &other) {
    if (this != &other) {
      mpfr_set(val, other.val, MPFR_RNDN);
    }
    return *this;
  }

  ~mpfloat() { mpfr_clear(val); }

  // Arithmetic operators
  mpfloat &operator+=(const mpfloat &rhs) {
    mpfr_add(val, val, rhs.val, MPFR_RNDN);
    return *this;
  }

  mpfloat &operator-=(const mpfloat &rhs) {
    mpfr_sub(val, val, rhs.val, MPFR_RNDN);
    return *this;
  }

  mpfloat &operator*=(const mpfloat &rhs) {
    mpfr_mul(val, val, rhs.val, MPFR_RNDN);
    return *this;
  }

  mpfloat &operator/=(const mpfloat &rhs) {
    mpfr_div(val, val, rhs.val, MPFR_RNDN);
    return *this;
  }

  // Friend operators for mpfloat + mpfloat
  friend mpfloat operator+(mpfloat lhs, const mpfloat &rhs) {
    lhs += rhs;
    return lhs;
  }
  friend mpfloat operator-(mpfloat lhs, const mpfloat &rhs) {
    lhs -= rhs;
    return lhs;
  }
  friend mpfloat operator-(const mpfloat &x) {
    mpfloat res;
    mpfr_neg(res.get(), x.get(), MPFR_RNDN);
    return res;
  }
  friend mpfloat operator*(mpfloat lhs, const mpfloat &rhs) {
    lhs *= rhs;
    return lhs;
  }
  friend mpfloat operator/(mpfloat lhs, const mpfloat &rhs) {
    lhs /= rhs;
    return lhs;
  }

  // Friend operators for mpfloat with double
  friend mpfloat operator+(mpfloat lhs, double rhs) {
    mpfr_t tmp;
    mpfr_init2(tmp, prec);
    mpfr_set_d(tmp, rhs, MPFR_RNDN);
    mpfr_add(lhs.val, lhs.val, tmp, MPFR_RNDN);
    mpfr_clear(tmp);
    return lhs;
  }
  friend mpfloat operator+(double lhs, const mpfloat &rhs) { return rhs + lhs; }

  friend mpfloat operator-(mpfloat lhs, double rhs) {
    mpfr_t tmp;
    mpfr_init2(tmp, prec);
    mpfr_set_d(tmp, rhs, MPFR_RNDN);
    mpfr_sub(lhs.val, lhs.val, tmp, MPFR_RNDN);
    mpfr_clear(tmp);
    return lhs;
  }
  friend mpfloat operator-(double lhs, const mpfloat &rhs) {
    mpfloat res(lhs);
    res -= rhs;
    return res;
  }

  friend mpfloat operator*(mpfloat lhs, double rhs) {
    mpfr_mul_d(lhs.val, lhs.val, rhs, MPFR_RNDN);
    return lhs;
  }
  friend mpfloat operator*(double lhs, const mpfloat &rhs) { return rhs * lhs; }

  friend mpfloat operator/(mpfloat lhs, double rhs) {
    mpfr_div_d(lhs.val, lhs.val, rhs, MPFR_RNDN);
    return lhs;
  }
  friend mpfloat operator/(double lhs, const mpfloat &rhs) {
    mpfloat res(lhs);
    mpfr_div(res.val, res.val, rhs.val, MPFR_RNDN);
    return res;
  }

  // Output
  friend std::ostream &operator<<(std::ostream &os, const mpfloat &x) {
    char *s = nullptr;
    int digits = os.precision();
    mpfr_asprintf(&s, "%.*Rg", digits, x.val);
    os << s;
    mpfr_free_str(s);
    return os;
  }

  // Access
  mpfr_t &get() { return val; }
  const mpfr_t &get() const { return val; }
};

// -------------------------
// Trig / exp functions
// -------------------------
template <unsigned int prec> mpfloat<prec> sin(const mpfloat<prec> &x) {
  mpfloat<prec> res;
  mpfr_sin(res.get(), x.get(), MPFR_RNDN);
  return res;
}

template <unsigned int prec> mpfloat<prec> cos(const mpfloat<prec> &x) {
  mpfloat<prec> res;
  mpfr_cos(res.get(), x.get(), MPFR_RNDN);
  return res;
}

template <unsigned int prec> mpfloat<prec> exp(const mpfloat<prec> &x) {
  mpfloat<prec> res;
  mpfr_exp(res.get(), x.get(), MPFR_RNDN);
  return res;
}

template <unsigned int prec> mpfloat<prec> pi() {
  mpfloat<prec> res;
  mpfr_const_pi(res.get(), MPFR_RNDN);
  return res;
}

} // namespace mp

namespace mpc {

// -------------------------
// mpcomplex class
// -------------------------
template <unsigned int prec> class mpcomplex {
  mpc_t val;

public:
  // Constructors / Destructor
  explicit mpcomplex() {
    mpc_init2(val, prec);
    mpc_set_ui_ui(val, 0, 0, MPC_RNDNN);
  }

  mpcomplex(double re, double im = 0.0) {
    mpc_init2(val, prec);
    mpc_set_d_d(val, re, im, MPC_RNDNN);
  }

  mpcomplex(const mpcomplex &other) {
    mpc_init2(val, prec);
    mpc_set(val, other.val, MPC_RNDNN);
  }

  mpcomplex &operator=(const mpcomplex &other) {
    if (this != &other) {
      mpc_set(val, other.val, MPC_RNDNN);
    }
    return *this;
  }

  mpcomplex(const mp::mpfloat<prec> &re) {
    mpc_init2(val, prec);
    mpc_set_fr(val, re.get(), MPC_RNDNN); // imag part = 0
  }

  ~mpcomplex() { mpc_clear(val); }

  // Arithmetic operators
  mpcomplex &operator+=(const mpcomplex &rhs) {
    mpc_add(val, val, rhs.val, MPC_RNDNN);
    return *this;
  }

  mpcomplex &operator*=(const mpcomplex &rhs) {
    mpc_mul(val, val, rhs.val, MPC_RNDNN);
    return *this;
  }

  mpcomplex &muladd(const mpcomplex &a, const mpcomplex &b) {
    mpc_t tmp;
    mpc_init2(tmp, prec);
    mpc_mul(tmp, a.val, b.val, MPC_RNDNN);
    mpc_add(val, val, tmp, MPC_RNDNN);
    mpc_clear(tmp);
    return *this;
  }

  friend mpcomplex operator+(mpcomplex lhs, const mpcomplex &rhs) {
    lhs += rhs;
    return lhs;
  }

  friend mpcomplex operator*(mpcomplex lhs, const mpcomplex &rhs) {
    lhs *= rhs;
    return lhs;
  }

  friend std::ostream &operator<<(std::ostream &os, const mpcomplex &z) {
    char *s = mpc_get_str(10, os.precision(), z.val, MPC_RNDNN);
    os << s;
    mpc_free_str(s);
    return os;
  }

  mpc_t &get() { return val; }
  const mpc_t &get() const { return val; }

  // Scalar arithmetic
  friend mpcomplex operator*(const mpcomplex &lhs, double rhs) {
    mpcomplex res(lhs);
    mpfr_t tmp;
    mpfr_init2(tmp, prec);
    mpfr_set_d(tmp, rhs, MPFR_RNDN);
    mpc_mul_fr(res.val, lhs.val, tmp, MPC_RNDNN);
    mpfr_clear(tmp);
    return res;
  }

  friend mpcomplex operator*(double lhs, const mpcomplex &rhs) {
    return rhs * lhs;
  }

  friend mpcomplex operator+(const mpcomplex &lhs, double rhs) {
    mpcomplex res(lhs);
    mpfr_t tmp;
    mpfr_init2(tmp, prec);
    mpfr_set_d(tmp, rhs, MPFR_RNDN);
    mpc_add_fr(res.val, lhs.val, tmp, MPC_RNDNN);
    mpfr_clear(tmp);
    return res;
  }

  friend mpcomplex operator+(double lhs, const mpcomplex &rhs) {
    return rhs + lhs;
  }

  friend mpcomplex operator-(const mpcomplex &lhs, double rhs) {
    mpcomplex res(lhs);
    mpfr_t tmp;
    mpfr_init2(tmp, prec);
    mpfr_set_d(tmp, rhs, MPFR_RNDN);
    mpc_sub_fr(res.val, lhs.val, tmp, MPC_RNDNN);
    mpfr_clear(tmp);
    return res;
  }

  friend mpcomplex operator-(double lhs, const mpcomplex &rhs) {
    mpcomplex res;
    mpfr_t tmp;
    mpfr_init2(tmp, prec);
    mpfr_set_d(tmp, lhs, MPFR_RNDN);
    mpc_fr_sub(res.val, tmp, rhs.val, MPC_RNDNN);
    mpfr_clear(tmp);
    return res;
  }

  friend mpcomplex operator-(const mpcomplex &z) {
    mpcomplex res;
    mpc_neg(res.get(), z.get(), MPC_RNDNN);
    return res;
  }

  // Subtraction of two mpcomplex numbers
  friend mpcomplex operator-(const mpcomplex &lhs, const mpcomplex &rhs) {
    mpcomplex res(lhs);
    mpc_sub(res.val, lhs.val, rhs.val, MPC_RNDNN);
    return res;
  }

  friend mpcomplex operator/(const mpcomplex &lhs, const mpcomplex &rhs) {
    mpcomplex res(lhs);
    mpc_div(res.val, lhs.val, rhs.val, MPC_RNDNN);
    return res;
  }

  friend mpcomplex operator/(const mpcomplex &lhs, double rhs) {
    mpcomplex res(lhs);
    mpfr_t tmp;
    mpfr_init2(tmp, prec);
    mpfr_set_d(tmp, rhs, MPFR_RNDN);
    mpc_div_fr(res.val, lhs.val, tmp, MPC_RNDNN);
    mpfr_clear(tmp);
    return res;
  }

  friend mpcomplex operator/(double lhs, const mpcomplex &rhs) {
    mpcomplex res;
    mpfr_t tmp;
    mpfr_init2(tmp, prec);
    mpfr_set_d(tmp, lhs, MPFR_RNDN);
    mpc_fr_div(res.val, tmp, rhs.val, MPC_RNDNN);
    mpfr_clear(tmp);
    return res;
  }

  // --- Multiply by mp::mpfloat ---
  friend mpcomplex operator*(const mpcomplex &lhs,
                             const mp::mpfloat<prec> &rhs) {
    mpcomplex res(lhs);
    mpc_mul_fr(res.val, lhs.val, rhs.get(), MPC_RNDNN);
    return res;
  }

  friend mpcomplex operator*(const mp::mpfloat<prec> &lhs,
                             const mpcomplex &rhs) {
    return rhs * lhs;
  }

  // --- Add mpfloat to complex (as real part) ---
  friend mpcomplex operator+(const mpcomplex &lhs,
                             const mp::mpfloat<prec> &rhs) {
    mpcomplex res(lhs);
    mpc_add_fr(res.val, lhs.val, rhs.get(), MPC_RNDNN);
    return res;
  }

  friend mpcomplex operator+(const mp::mpfloat<prec> &lhs,
                             const mpcomplex &rhs) {
    return rhs + lhs;
  }

  // --- Subtract mpfloat ---
  friend mpcomplex operator-(const mpcomplex &lhs,
                             const mp::mpfloat<prec> &rhs) {
    mpcomplex res(lhs);
    mpc_sub_fr(res.val, lhs.val, rhs.get(), MPC_RNDNN);
    return res;
  }

  friend mpcomplex operator-(const mp::mpfloat<prec> &lhs,
                             const mpcomplex &rhs) {
    mpcomplex res;
    mpc_fr_sub(res.get(), lhs.get(), rhs.val, MPC_RNDNN);
    return res;
  }

  // --- Divide by mpfloat ---
  friend mpcomplex operator/(const mpcomplex &lhs,
                             const mp::mpfloat<prec> &rhs) {
    mpcomplex res(lhs);
    mpc_div_fr(res.val, lhs.val, rhs.get(), MPC_RNDNN);
    return res;
  }

  friend mpcomplex operator/(const mp::mpfloat<prec> &lhs,
                             const mpcomplex &rhs) {
    mpcomplex res;
    mpc_fr_div(res.val, lhs.get(), rhs.val, MPC_RNDNN);
    return res;
  }

  mp::mpfloat<prec> real() const {
    mp::mpfloat<prec> res;
    mpfr_set(res.get(), mpc_realref(val), MPFR_RNDN);
    return res;
  }

  mp::mpfloat<prec> imag() const {
    mp::mpfloat<prec> res;
    mpfr_set(res.get(), mpc_imagref(val), MPFR_RNDN);
    return res;
  }
};

// -------------------------
// Trig / exp functions
// -------------------------
template <unsigned int prec> mpcomplex<prec> cos(const mpcomplex<prec> &z) {
  mpcomplex<prec> res;
  mpc_cos(res.get(), z.get(), MPC_RNDNN);
  return res;
}

template <unsigned int prec> mpcomplex<prec> sin(const mpcomplex<prec> &z) {
  mpcomplex<prec> res;
  mpc_sin(res.get(), z.get(), MPC_RNDNN);
  return res;
}

template <unsigned int prec> mpcomplex<prec> exp(const mpcomplex<prec> &z) {
  mpcomplex<prec> res;
  mpc_exp(res.get(), z.get(), MPC_RNDNN);
  return res;
}
} // namespace mpc

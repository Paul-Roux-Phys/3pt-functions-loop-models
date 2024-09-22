#pragma once

#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;
#include <string>
#include <cstdint>
#include <array>
#include <functional>
#include <cmath>
#include<unordered_map>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/mpc.hpp>
namespace mp = boost::multiprecision;
#include "arcomp.h"


/* ________________________________________________________________________________________________________________

  Choose precision
 __________________________________________________________________________________________________________________ */
#if (!defined(STD_PREC) && !defined(MULT_PREC))
#error "You must #define either STD_PREC or MULT_PREC"
#endif

#if (defined(MULT_PREC) && !defined(DECIMAL_PLACES))
#error "In multiple precision, you must #define DECIMAL_PLACES to indicate how many decimals should be used"
#endif

#ifdef DECIMAL_PLACES
constexpr uint MULT_PREC_DECIMAL_PLACES = DECIMAL_PLACES;
#endif


/* ________________________________________________________________________________________________________________

  Choose matrix element type
 __________________________________________________________________________________________________________________ */
#if (!defined(REAL_MATRIX) && !defined(COMPLEX_MATRIX))
#error "You must #define either REAL_MATRIX or COMPLEX_MATRIX"
#endif

/* ________________________________________________________________________________________________________________

  Other parameters
 __________________________________________________________________________________________________________________ */

#ifndef SIZE
#error "You must #define SIZE to the required size of keys representing the system's states"
#endif

constexpr std::size_t KEY_SIZE = SIZE;

#ifndef HASH_SIZE
#error "You must #define HASH_SIZE to the required hash table size." \
       "It is recommended to set it to around 0.7*(expected vector size)" \
       "for optimal performance."
#endif

constexpr std::size_t HASH_TABLE_SIZE = HASH_SIZE;
/* ________________________________________________________________________________________________________________

  Create type aliases
 __________________________________________________________________________________________________________________ */
namespace transfermatrices {

using Key = std::array<int8_t, KEY_SIZE>;

#ifdef STD_PREC
  using Real = double;
#else
  using Real = mp::number<mp::mpfr_float_backend<MULT_PREC_DECIMAL_PLACES>>;
#endif
  using Complex = arcomplex<Real>;

#ifdef COMPLEX_MATRIX
  using Weight = Complex;
#else
  using Weight = Real;
#endif

template <typename T>
class Vector;

template <typename W, typename BV, typename Vec>
class VectorPair;

#include "../transfer_matrices/keys.hpp"
#include "../transfer_matrices/vectors.hpp"
#include "../transfer_matrices/r_matrices.hpp"
#include "../transfer_matrices/vector_pair.hpp"
#include "../transfer_matrices/diagonalise.hpp"

} // namespace
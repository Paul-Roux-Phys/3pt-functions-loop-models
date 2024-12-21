#ifndef DLARFG_H
#define DLARFG_H

#include "cdgemv.h"
#include <cln/real.h>
#include <cln/output.h>
#include <cln/real_io.h>

/* CLN version */

template<typename T>
void dlarf(const std::string& side, int m, int n, T* v, int incv, T tau, T* c, int ldc, T* work, int digits)
{
  /*
*  Purpose
*  =======
*
*  DLARF applies a real elementary reflector H to a real m by n matrix
*  C, from either the left or the right. H is represented in the form
*
*        H = I - tau * v * v'
*
*  where tau is a real scalar and v is a real vector.
*
*  If tau = 0, then H is taken to be the unit matrix.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': form  H * C
*          = 'R': form  C * H
*
*  M       (input) INTEGER
*          The number of rows of the matrix C.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C.
*
*  V       (input) DOUBLE PRECISION array, dimension
*                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
*                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
*          The vector v in the representation of H. V is not used if
*          TAU = 0.
*
*  INCV    (input) INTEGER
*          The increment between elements of v. INCV <> 0.
*
*  TAU     (input) DOUBLE PRECISION
*          The value tau in the representation of H.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the m by n matrix C.
*          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
*          or C * H if SIDE = 'R'.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension
*                         (N) if SIDE = 'L'
*                      or (M) if SIDE = 'R'
*
*  =====================================================================
*
  */
  cln::float_format_t precision=cln::float_format(digits);
  const T one = cln::cl_float(1,precision);
  const T zero = cln::cl_float(0,precision);
  
  if(side=="L")
  {
    /* Form  H * C */
    if(tau!=zero)
    {
      /* w := C' * v */
      dgemv<T>("T",m,n,one,c,ldc,v,incv,zero,work,1,digits);
      
      /* C := C - v * w' */
      dger<T>(m,n,-tau,v,incv,work,1,c,ldc,digits);
    }
  }
  else
  {
    /* Form  C * H */
    if(tau!=zero)
    {
      /* w := C * v */
      dgemv<T>("N",m,n,one,c,ldc,v,incv,zero,work,1,digits);
      
      /* C := C - w * v' */
      dger<T>(m,n,-tau,work,1,v,incv,c,ldc,digits);
    }
  }
  return;
  
  /* End of DLARF */
}

template<typename T>
void dlarf(const std::string& side, int m, int n, T* v, int incv, T tau, T* c, int ldc, T* work)
{
  /*
*  Purpose
*  =======
*
*  DLARF applies a real elementary reflector H to a real m by n matrix
*  C, from either the left or the right. H is represented in the form
*
*        H = I - tau * v * v'
*
*  where tau is a real scalar and v is a real vector.
*
*  If tau = 0, then H is taken to be the unit matrix.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': form  H * C
*          = 'R': form  C * H
*
*  M       (input) INTEGER
*          The number of rows of the matrix C.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C.
*
*  V       (input) DOUBLE PRECISION array, dimension
*                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
*                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
*          The vector v in the representation of H. V is not used if
*          TAU = 0.
*
*  INCV    (input) INTEGER
*          The increment between elements of v. INCV <> 0.
*
*  TAU     (input) DOUBLE PRECISION
*          The value tau in the representation of H.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the m by n matrix C.
*          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
*          or C * H if SIDE = 'R'.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension
*                         (N) if SIDE = 'L'
*                      or (M) if SIDE = 'R'
*
*  =====================================================================
*
  */
  T one=1.0;
  T zero=0.0;
  if(side=="L")
  {
    /* Form  H * C */
    if(tau!=zero)
    {
      /* w := C' * v */
      dgemv<T>("T",m,n,one,c,ldc,v,incv,zero,work,1);
      
      /* C := C - v * w' */
      dger<T>(m,n,-tau,v,incv,work,1,c,ldc);
    }
  }
  else
  {
    /* Form  C * H */
    if(tau!=zero)
    {
      /* w := C * v */
      dgemv<T>("N",m,n,one,c,ldc,v,incv,zero,work,1);
      
      /* C := C - w * v' */
      dger<T>(m,n,-tau,work,1,v,incv,c,ldc);
    }
  }
  return;
  
  /* End of DLARF */
}

#endif

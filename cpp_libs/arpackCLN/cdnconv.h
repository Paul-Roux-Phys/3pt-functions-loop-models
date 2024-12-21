#ifndef DNCONV_H
#define DNCONV_H

/*
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dnconv
c
c\Description: 
c  Convergence testing for the nonsymmetric Arnoldi eigenvalue routine.
c
c\Usage:
c  call dnconv
c     ( N, RITZR, RITZI, BOUNDS, TOL, NCONV )
c
c\Arguments
c  N       Integer.  (INPUT)
c          Number of Ritz values to check for convergence.
c
c  RITZR,  Double precision arrays of length N.  (INPUT)
c  RITZI   Real and imaginary parts of the Ritz values to be checked
c          for convergence.

c  BOUNDS  Double precision array of length N.  (INPUT)
c          Ritz estimates for the Ritz values in RITZR and RITZI.
c
c  TOL     Double precision scalar.  (INPUT)
c          Desired backward error for a Ritz value to be considered
c          "converged".
c
c  NCONV   Integer scalar.  (OUTPUT)
c          Number of "converged" Ritz values.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     second  ARPACK utility routine for timing.
c     dlamch  LAPACK routine that determines machine constants.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University 
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics 
c     Rice University           
c     Houston, Texas    
c
c\Revision history:
c     xx/xx/92: Version ' 2.1'
c
c\SCCS Information: @(#) 
c FILE: nconv.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
c
c\Remarks
c     1. xxxx
c
c\EndLib
c
c-----------------------------------------------------------------------
c
*/

#include "cdebug.h"
#include "cstat.h"
#include "fortranfuncs.h"
#include <math.h>
#include <cln/real.h>
#include <cln/output.h>
#include <cln/real_io.h>

/* CLN version */

template<typename T>
void dnconv(int n, T* ritzr, T* ritzi, T* bounds, T tol, int& nconv, int digits)
{
  int i;
  T temp,eps23;
  
  cln::float_format_t precision = cln::float_format(digits);
  const cln::cl_F two=cln::cl_float(2,precision);
  const cln::cl_F three=cln::cl_float(3,precision);
  
  /*
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c 
c     %-------------------------------------------------------------%
c     | Convergence test: unlike in the symmetric code, I am not    |
c     | using things like refined error bounds and gap condition    |
c     | because I don't know the exact equivalent concept.          |
c     |                                                             |
c     | Instead the i-th Ritz value is considered "converged" when: |
c     |                                                             |
c     |     bounds(i) .le. ( TOL * | ritz | )                       |
c     |                                                             |
c     | for some appropriate choice of norm.                        |
c     %-------------------------------------------------------------%
c
  */
  
  second(timing.t0);
  
  /*
c
c     %---------------------------------%
c     | Get machine dependent constant. |
c     %---------------------------------%
c
  */
  
  eps23=dlamch<T>("E",digits);
//  std::cout << "eps = " << eps23 << "\n";
  eps23=cln::cl_float(cln::exp(two/three*cln::ln(eps23)),precision); //  eps23=pow(eps23,2.0/3.0);
//  std::cout << "eps23 = " << eps23 << "\n";
//  std::cout << "tol = " << tol << "\n";

  nconv=0;
  for(i=1;i<=n;i++)
  {
    temp=cln::max(eps23,dlapy2<T>(ritzr[i-1],ritzi[i-1]));
    if(bounds[i-1]<=tol*temp) nconv+=1;
  }
  
  second(timing.t1);
  timing.tnconv+=timing.t1-timing.t0;
  
  return;
  /*
c
c     %---------------%
c     | End of dnconv |
c     %---------------%
c
  */
}

template<typename T>
void dnconv(int n, T* ritzr, T* ritzi, T* bounds, T tol, int& nconv)
{
  int i;
  T temp,eps23;
  
  /*
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c 
c     %-------------------------------------------------------------%
c     | Convergence test: unlike in the symmetric code, I am not    |
c     | using things like refined error bounds and gap condition    |
c     | because I don't know the exact equivalent concept.          |
c     |                                                             |
c     | Instead the i-th Ritz value is considered "converged" when: |
c     |                                                             |
c     |     bounds(i) .le. ( TOL * | ritz | )                       |
c     |                                                             |
c     | for some appropriate choice of norm.                        |
c     %-------------------------------------------------------------%
c
  */
  
  second(timing.t0);
  
  /*
c
c     %---------------------------------%
c     | Get machine dependent constant. |
c     %---------------------------------%
c
  */
  
  eps23=dlamch<T>("E");
  eps23=pow(eps23,(2.0/3.0));
  
  nconv=0;
  for(i=1;i<=n;i++)
  {
    temp=std::max(eps23,dlapy2<T>(ritzr[i-1],ritzi[i-1]));
    if(bounds[i-1]<=tol*temp) nconv+=1;
  }
  
  second(timing.t1);
  timing.tnconv+=timing.t1-timing.t0;
  
  return;
  /*
c
c     %---------------%
c     | End of dnconv |
c     %---------------%
c
  */
}


#endif
  
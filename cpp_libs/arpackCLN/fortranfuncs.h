#ifndef FORTRAN_FUNCS_H
#define FORTRAN_FUNCS_H

/* Fortran dcopy */

#include <iostream>
#include <math.h>
#include <string>
#include <algorithm>
#include <iomanip>
#include <time.h>
#include <stdlib.h>
#include <cln/real.h>
#include <cln/output.h>
#include <cln/real_io.h>
#include <cln/float_io.h>
//#include <random>

template<typename T>
void dlamc2(int&, int&, bool, T&, int&, T&, int&, T&);
template<typename T>
void dlamc2(int&, int&, bool, T&, int&, T&, int&, T&, int);
template<typename T>
T dlamc3(T,T);
template<typename T>
void dlamc4(int&, T, int);
template<typename T>
void dlamc4(int&, T, int, int);
template<typename T>
void dlamc5(int, int, int, bool, int, T&);
template<typename T>
void dlamc5(int, int, int, bool, int, T&, int);
template<typename T>
void dlassq(int, T*, int, T&, T&);
template<typename T>
void dlassq(int, T*, int, T&, T&, int);
template<typename T>
void dger(int, int, T, T*, int, T*, int, T*, int);

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template <> int sgn <cln::cl_F> (cln::cl_F val) {
    cln::float_format_t precision = cln::float_format(cln::float_digits(val) );
    /* Parameters */
    const cln::cl_F zero=cln::cl_float(0,precision);
    return ((zero < val) - (val < zero));
}

template <typename T> T sign(T A, T B) {
    int sgnB=sgn<T>(B);
    if(sgnB==0)
      return sgn<T>(A)*A;
    else
      return sgn<T>(A)*sgn<T>(B)*A;
}

template <> cln::cl_F sign(cln::cl_F A, cln::cl_F B) {
    int sgnB=sgn<cln::cl_F>(B);
    cln::float_format_t precision = cln::float_format(cln::float_digits(A) );
    cln::cl_F sgnAfloat=cln::cl_float(sgn<cln::cl_F>(A),precision);
    cln::cl_F sgnBfloat;
    if(sgnB==0)
      return sgnAfloat*A;
    else
    {
      sgnBfloat=cln::cl_float(sgn<cln::cl_F>(B),precision);
      return sgnAfloat*sgnBfloat*A;
    }
}

template<typename T>
void copy(int n, T* dx, int incx, T* dy, int incy)
{
  int i,ix,iy,m,mp1;
  
  for(i=0;i<n;i++)
    dy[i]=dx[i];
  
 /* if(n<0) return;
  if((incx==1)&&(incy==1))
  {
    m=n%7;
    if(m!=0)
    {
      for(i=0; i<m; i++)
        dy[i]=dx[i];
      if(n<7) return;
    }

    mp1=m;
    for(i=mp1;i<n;i=i+7)
    {
      dy[i]=dx[i];
      dy[i+1]=dx[i+1];
      dy[i+2]=dx[i+2];
      dy[i+3]=dx[i+3];
      dy[i+4]=dx[i+4];
      dy[i+5]=dx[i+5];
      dy[i+6]=dx[i+6];
    }    
  }
  else
  {
    ix=1;
    iy=1;
    if(incx<0)ix = (-n+1)*incx+1;
    if(incy<0)iy = (-n+1)*incy+1;
    for(i=0;i<n;i++)
    {
      dy[iy]=dx[ix];
      ix+=incx;
      iy+=incy;
    }    
  }*/
  
  return;
}

template<typename T>
T dlamch(const std::string cmach, int digits)
{
  
  /* CLN version*/
  
  cln::float_format_t precision = cln::float_format(digits);
    /* Parameters */
  const cln::cl_F one=cln::cl_float(1,precision);
  const cln::cl_F zero=cln::cl_float(0,precision);
  
  /* Local Scalars */
  static bool first=true;
  bool lrnd;
  int beta,imax,imin,it;
  cln::cl_F rmach,small;
  static cln::cl_F eps,sfmin,base,t,rnd,emin,rmin,emax,rmax,prec;
  
  const std::string valuetype=cmach;
  
//  std::cout << "digits = " << digits << "\n";
  /* This section put in for debugging purposes */
  
  
  
  /* Executable Statements */
  
/*  if(first)
  {
    first=false;
    dlamc2<cln::cl_F>(beta,it,lrnd,eps,imin,rmin,imax,rmax,digits);
    base=cl_float(beta,precision);
    cln::print_float(std::cout,cln::default_print_flags,base);
    std::cout << "\n";
    t=cl_float(it,precision);
    if(lrnd)
    {
      rnd=one;
//      eps=pow(base,1-it)/2;
      eps=cln::cl_float(expt(base,1-it)/2,precision);
    }
    else
    {
      rnd=zero;
      eps=cln::cl_float(expt(base,1-it),precision);
    }
    prec=eps*base;
    emin=cl_float(imin,precision);
    emax=cl_float(imax,precision);
    sfmin=rmin;
    small=one/rmax;
    if(small>=sfmin)
    {
      /*
*
*           Use SMALL plus a bit, to avoid the possibility of rounding
*           causing overflow when computing  1/sfmin.
*
*     */
      
/*      sfmin=small*(one+eps);
    }
  }*/

  if(first)
  {
    eps=cln::float_epsilon(precision);
//    beta=cln::float_radix(precision);
    base=cln::cl_float(2,precision);
//    std::cout << "base = ";
//    print_float(std::cout,cln::default_print_flags,base);
//    std::cout << "\n";
    prec=eps*base;
    rmin=cln::least_positive_float(precision);
    sfmin=rmin;
    rmax=cln::most_positive_float(precision);
 //   std::cout << "rmax = " << rmax << "\n";
 //   std::cout << "rmin = " << rmin << "\n";
 //   std::cout << "one = " << one << "\n";
    small=one/rmax;
    if(small>=sfmin)
    {
      /*
*
*           Use SMALL plus a bit, to avoid the possibility of rounding
*           causing overflow when computing  1/sfmin.
*
*     */
      
      sfmin=small*(one+eps);
    }
  }    
  
  if((valuetype=="E")||(valuetype=="e"))
    rmach=eps;
  else if((valuetype=="S")||(valuetype=="s"))
    rmach=sfmin;
  else if((valuetype=="B")||(valuetype=="b"))
    rmach=base;
  else if((valuetype=="P")||(valuetype=="p"))
    rmach=prec;
//  else if((valuetype=="N")||(valuetype=="n"))
//    rmach=t;
//  else if((valuetype=="R")||(valuetype=="r"))
//    rmach=rnd;
//  else if((valuetype=="M")||(valuetype=="m"))
//    rmach=emin;
  else if((valuetype=="U")||(valuetype=="u"))
    rmach=rmin;
//  else if((valuetype=="L")||(valuetype=="l"))
//    rmach=emax;
  else if((valuetype=="O")||(valuetype=="o"))
    rmach=rmax;
  
  return rmach;
  
  /* End of CLN DLAMCH */
}

template<typename T>
T dlamch(const std::string cmach)
{
  /*
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          CMACH
*     ..
*
*  Purpose
*  =======
*
*  DLAMCH determines double precision machine parameters.
*
*  Arguments
*  =========
*
*  CMACH   (input) CHARACTER*1
*          Specifies the value to be returned by DLAMCH:
*          = 'E' or 'e',   DLAMCH := eps
*          = 'S' or 's ,   DLAMCH := sfmin
*          = 'B' or 'b',   DLAMCH := base
*          = 'P' or 'p',   DLAMCH := eps*base
*          = 'N' or 'n',   DLAMCH := t
*          = 'R' or 'r',   DLAMCH := rnd
*          = 'M' or 'm',   DLAMCH := emin
*          = 'U' or 'u',   DLAMCH := rmin
*          = 'L' or 'l',   DLAMCH := emax
*          = 'O' or 'o',   DLAMCH := rmax
*
*          where
*
*          eps   = relative machine precision
*          sfmin = safe minimum, such that 1/sfmin does not overflow
*          base  = base of the machine
*          prec  = eps*base
*          t     = number of (base) digits in the mantissa
*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*          emin  = minimum exponent before (gradual) underflow
*          rmin  = underflow threshold - base**(emin-1)
*          emax  = largest exponent before overflow
*          rmax  = overflow threshold  - (base**emax)*(1-eps)
*
* =====================================================================
*
  */
  
  /* Parameters */
  const T one=1.0000000000000000000000;
  const T zero=0.000000000000000000000;
  
  /* Local Scalars */
  static bool first=true;
  bool lrnd;
  int beta,imax,imin,it;
  T rmach,small;
  static T eps,sfmin,base,t,rnd,emin,rmin,emax,rmax,prec;
  
  const std::string valuetype=cmach;
  
  
  /* This section put in for debugging purposes */
  
  
  
  /* Executable Statements */
  
  if(first)
  {
    first=false;
    dlamc2<T>(beta,it,lrnd,eps,imin,rmin,imax,rmax);
    base=beta;
    t=it;
    if(lrnd)
    {
      rnd=1.0;
      eps=pow(base,1-it)/2;
    }
    else
    {
      rnd=0.0;
      eps=pow(base,1-it);
    }
    prec=eps*base;
    emin=imin;
    emax=imax;
    sfmin=rmin;
    small=1.0/rmax;
    if(small>=sfmin)
    {
      /*
*
*           Use SMALL plus a bit, to avoid the possibility of rounding
*           causing overflow when computing  1/sfmin.
*
*     */
      
      sfmin=small*(one+eps);
    }
  }
  
  if((valuetype=="E")||(valuetype=="e"))
    rmach=eps;
  else if((valuetype=="S")||(valuetype=="s"))
    rmach=sfmin;
  else if((valuetype=="B")||(valuetype=="b"))
    rmach=base;
  else if((valuetype=="P")||(valuetype=="p"))
    rmach=prec;
  else if((valuetype=="N")||(valuetype=="n"))
    rmach=t;
  else if((valuetype=="R")||(valuetype=="r"))
    rmach=rnd;
  else if((valuetype=="M")||(valuetype=="m"))
    rmach=emin;
  else if((valuetype=="U")||(valuetype=="u"))
    rmach=rmin;
  else if((valuetype=="L")||(valuetype=="l"))
    rmach=emax;
  else if((valuetype=="O")||(valuetype=="o"))
    rmach=rmax;
  
  return rmach;
  
  /* End of DLAMCH */
}

template<typename T>
void dlamc1(int& beta, int& t, bool& rnd, bool& ieee1)
{
  /*
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*  Purpose
*  =======
*
*  DLAMC1 determines the machine parameters given by BETA, T, RND, and
*  IEEE1.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  IEEE1   (output) LOGICAL
*          Specifies whether rounding appears to be done in the IEEE
*          'round to nearest' style.
*
*  Further Details
*  ===============
*
*  The routine is based on the routine  ENVRON  by Malcolm and
*  incorporates suggestions by Gentleman and Marovich. See
*
*     Malcolm M. A. (1972) Algorithms to reveal properties of
*        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
*
*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
*        that reveal properties of floating point arithmetic units.
*        Comms. of the ACM, 17, 276-277.
*
* =====================================================================
*
  */
  
  /* Local Scalars */
  static bool first=true;
  static bool lieee1,lrnd;
  static int lbeta,lt;
  T a,b,c,f,one,qtr,savec,t1,t2;
  
  /* Executable Statements */
  
  if(first)
  {
    first=false;
    one=1.0;
    /*
*
*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
*        IEEE1, T and RND.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are  stored and not held in registers,  or
*        are not affected by optimizers.
*
*        Compute  a = 2.0**m  with the  smallest positive integer m such
*        that
*
*           fl( a + 1.0 ) = a.
*
    */
    a=1.0;
    c=1.0;
    
label10:
    if(c==one)
    {
      a=2*a;
      c=dlamc3<T>(a,one);
      c=dlamc3<T>(c,-a);
      goto label10;
    }
    /*
*
*        Now compute  b = 2.0**m  with the smallest positive integer m
*        such that
*
*           fl( a + b ) .gt. a.
*
    */
    
    b=1.0;
    c=dlamc3(a,b);
    
label20:
    if(c==a)
    {
      b=2*b;
      c=dlamc3(a,b);
      goto label20;
    }
    
    /*
*
*        Now compute the base.  a and c  are neighbouring floating point
*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
*        their difference is beta. Adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).
*
    */
    qtr=one/4;
    savec=c;
    c=dlamc3(c,-a);
    lbeta=c+qtr;
    /*
*
*        Now determine whether rounding or chopping occurs,  by adding a
*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
*
    */
    b=lbeta;
    f=dlamc3(b/2,-b/100);
    c=dlamc3(f,a);
    if(c==a)
      lrnd=true;
    else
      lrnd=false;
    f=dlamc3<T>(b/2,b/100);
    c=dlamc3<T>(f,a);
    if(lrnd&&(c==a))
      lrnd=false;
    /*
*
*        Try and decide whether rounding is done in the  IEEE  'round to
*        nearest' style. B/2 is half a unit in the last place of the two
*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
*        A, but adding B/2 to SAVEC should change SAVEC.
*
    */
    t1=dlamc3<T>(b/2,a);
    t2=dlamc3<T>(b/2,savec);
    lieee1=(t1==a)&&(t2>savec)&&lrnd;
    /*
*
*        Now find  the  mantissa, t.  It should  be the  integer part of
*        log to the base beta of a,  however it is safer to determine  t
*        by powering.  So we find t as the smallest positive integer for
*        which
*
*           fl( beta**t + 1.0 ) = 1.0.
*
    */
    lt=0;
    a=1;
    c=1;

label30:
    if(c==one)
    {
      lt+=1;
      a=a*lbeta;
      c=dlamc3<T>(a,one);
      c=dlamc3<T>(c,-a);
      goto label30;
    }    
  }
  
  beta=lbeta;
  t=lt;
  rnd=lrnd;
  ieee1=lieee1;
  
  return;
  /* End of DLAMC1 */
}

/* CLN version */

template<typename T>
void dlamc1(int& beta, int& t, bool& rnd, bool& ieee1, int digits)
{
  /*
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*  Purpose
*  =======
*
*  DLAMC1 determines the machine parameters given by BETA, T, RND, and
*  IEEE1.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  IEEE1   (output) LOGICAL
*          Specifies whether rounding appears to be done in the IEEE
*          'round to nearest' style.
*
*  Further Details
*  ===============
*
*  The routine is based on the routine  ENVRON  by Malcolm and
*  incorporates suggestions by Gentleman and Marovich. See
*
*     Malcolm M. A. (1972) Algorithms to reveal properties of
*        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
*
*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
*        that reveal properties of floating point arithmetic units.
*        Comms. of the ACM, 17, 276-277.
*
* =====================================================================
*
  */
  
  /* Local Scalars */
  static bool first=true;
  static bool lieee1,lrnd;
  static int lbeta,lt;
  T a,b,c,f,one,qtr,savec,t1,t2;
  T two;
  cln::float_format_t precision=cln::float_format(digits);
  
  /* Executable Statements */
  
  if(first)
  {
    first=false;
    one = cln::cl_float(1,precision);
    two = cln::cl_float(2,precision);
    /*
*
*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
*        IEEE1, T and RND.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are  stored and not held in registers,  or
*        are not affected by optimizers.
*
*        Compute  a = 2.0**m  with the  smallest positive integer m such
*        that
*
*           fl( a + 1.0 ) = a.
*
    */
    a=one;
    c=one;
    
label10:
    if(c==one)
    {
      a=two*a;
      c=dlamc3<T>(a,one);
      c=dlamc3<T>(c,-a);
      goto label10;
    }
    /*
*
*        Now compute  b = 2.0**m  with the smallest positive integer m
*        such that
*
*           fl( a + b ) .gt. a.
*
    */
    
    b=one;
    c=dlamc3(a,b);
    
label20:
    if(c==a)
    {
      b=two*b;
      c=dlamc3(a,b);
      goto label20;
    }
    
    /*
*
*        Now compute the base.  a and c  are neighbouring floating point
*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
*        their difference is beta. Adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).
*
    */
    qtr=one/4;
    savec=c;
    c=dlamc3(c,-a);
    lbeta=int(cln::double_approx(c+qtr));
    /*
*
*        Now determine whether rounding or chopping occurs,  by adding a
*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
*
    */
    b=cln::cl_float(lbeta,precision);
    f=dlamc3(b/two,-b/100);
    c=dlamc3(f,a);
    if(c==a)
      lrnd=true;
    else
      lrnd=false;
    f=dlamc3<T>(b/two,b/100);
    c=dlamc3<T>(f,a);
    if(lrnd&&(c==a))
      lrnd=false;
    /*
*
*        Try and decide whether rounding is done in the  IEEE  'round to
*        nearest' style. B/2 is half a unit in the last place of the two
*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
*        A, but adding B/2 to SAVEC should change SAVEC.
*
    */
    t1=dlamc3<T>(b/two,a);
    t2=dlamc3<T>(b/two,savec);
    lieee1=(t1==a)&&(t2>savec)&&lrnd;
    /*
*
*        Now find  the  mantissa, t.  It should  be the  integer part of
*        log to the base beta of a,  however it is safer to determine  t
*        by powering.  So we find t as the smallest positive integer for
*        which
*
*           fl( beta**t + 1.0 ) = 1.0.
*
    */
    lt=0;
    a=one;
    c=one;

label30:
    if(c==one)
    {
      lt+=1;
      a=a*cln::cl_float(lbeta,precision);
      c=dlamc3<T>(a,one);
      c=dlamc3<T>(c,-a);
      goto label30;
    }    
  }
  
  beta=lbeta;
  t=lt;
  rnd=lrnd;
  ieee1=lieee1;
  
  return;
  /* End of CLN DLAMC1 */
}

template<typename T>
void dlamc2(int& beta, int& t, bool rnd, T& eps, int& emin, T& rmin, int& emax, T& rmax)
{
  /*
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
  */
  
  /*
*
*  Purpose
*  =======
*
*  DLAMC2 determines the machine parameters specified in its argument
*  list.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  EPS     (output) DOUBLE PRECISION
*          The smallest positive number such that
*
*             fl( 1.0 - EPS ) .LT. 1.0,
*
*          where fl denotes the computed value.
*
*  EMIN    (output) INTEGER
*          The minimum exponent before (gradual) underflow occurs.
*
*  RMIN    (output) DOUBLE PRECISION
*          The smallest normalized number for the machine, given by
*          BASE**( EMIN - 1 ), where  BASE  is the floating point value
*          of BETA.
*
*  EMAX    (output) INTEGER
*          The maximum exponent before overflow occurs.
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest positive number for the machine, given by
*          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
*          value of BETA.
*
*  Further Details
*  ===============
*
*  The computation of  EPS  is based on a routine PARANOIA by
*  W. Kahan of the University of California at Berkeley.
*
* =====================================================================
*
  */
  
  /* Local Scalars */
  static bool first=true;
  static bool iwarn=false;
  bool ieee,lieee1,lrnd;
  int gnmin,gpmin,i,ngnmin,ngpmin;
  static int lemin,lemax,lt,lbeta;
  T a,b,c,half,one,rbase,sixth,small,third,two,zero;
  static T lrmax,lrmin,leps;
  
  /* Executable Statements */
  
  if(first)
  {
    first=false;
    zero=0;
    one=1;
    two=2;
  
  /*
*
*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
*        BETA, T, RND, EPS, EMIN and RMIN.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are stored  and not held in registers,  or
*        are not affected by optimizers.
*
*        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
*
  */
    dlamc1<T>(lbeta,lt,lrnd,lieee1);
    
    /* Start to find EPS. */
  
    b=lbeta;
    a=pow(b,-lt);
    leps=a;
  
    /* Try some tricks to see whether or not this is the correct  EPS. */
  
    b=two/3.0;
    half=one/2.0;
    sixth=dlamc3<T>(b,-half);
    third=dlamc3<T>(sixth,sixth);
    b=dlamc3<T>(third,-half);
    b=dlamc3<T>(b,sixth);
    b=abs(b);
    if(b<leps)
      b=leps;
    
    leps=1.0;
  
label10:
    if((leps>b)&&(b>zero))
    {
      leps=b;
      c=dlamc3<T>(half*leps,pow(two,5)*leps*leps);
      c=dlamc3<T>(half,-c);
      b=dlamc3<T>(half,c);
      c=dlamc3<T>(half,-b);
      b=dlamc3<T>(half,c);
      goto label10;
    }
  
    if(a<leps)
      leps=a;
  
  //  std::cout << "eps = " << leps << "\n";
    /*
*
*        Computation of EPS complete.
*
*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
*        Keep dividing  A by BETA until (gradual) underflow occurs. This
*        is detected when we cannot recover the previous A.
*
    */
  
    rbase=one/lbeta;
    small=one;
    for(i=1;i<=3;i++)
      small=dlamc3<T>(small*rbase,zero);
    a=dlamc3<T>(one,small);
    dlamc4<T>(ngpmin,one,lbeta);
    dlamc4<T>(ngnmin,-one,lbeta);
    dlamc4<T>(gpmin,a,lbeta);
    dlamc4<T>(gnmin,-a,lbeta);
    ieee=false;
  
    if((ngpmin==ngnmin)&&(gpmin==gnmin))
    {
      if(ngpmin==gpmin)
        lemin=ngpmin;
        /*
*            ( Non twos-complement machines, no gradual underflow;
*              e.g.,  VAX )
        */
      else if((gpmin-ngpmin)==3)
      {
        lemin=ngpmin-1+lt;
        ieee=true;
        /*
*            ( Non twos-complement machines, with gradual underflow;
*              e.g., IEEE standard followers )
*       */
      }
      else
      {
        lemin=std::min(ngpmin,gpmin);
        /* ( A guess; no known machine ) */
        iwarn=true;
      }
    }
    else if((ngpmin==gpmin)&&(ngnmin==gnmin))
    {
      if(abs(ngpmin-ngnmin)==1)
      {
        lemin=std::max(ngpmin,ngnmin);
        /* ( Twos-complement machines, no gradual underflow;
               e.g., CYBER 205 ) */
      }
      else
      {
        lemin=std::min(ngpmin,ngnmin);
        /* ( A guess; no known machine ) */
        iwarn=true;
      }
    }
    else if((abs(ngpmin-ngnmin)==1)&&(gpmin==gnmin))
    {
      if((gpmin-std::min(ngpmin,ngnmin))==3)
      {
        lemin=std::max(ngpmin,ngnmin)-1+lt;
        /*
*            ( Twos-complement machines with gradual underflow;
*              no known machine )
        */
      }
      else
      {
        lemin=std::min(ngpmin,ngnmin);
        /* ( A guess; no known machine ) */
        iwarn=true;
      }
    }
    else
    {
      lemin=std::min(std::min(ngpmin,ngnmin),std::min(gpmin,gnmin));
      /* ( A guess; no known machine ) */
      iwarn=true;
    }
    /*
***
* Comment out this if block if EMIN is ok
    */
    if(iwarn)
    {
      first=true;
      std::cout << "WARNING. The value EMIN may be incorrect:-EMIN = " << lemin << "\n";
      std::cout << "If, after inspection, the value EMIN looks acceptable please comment out \n";
      std::cout << "the IF block as marked within the code of routine DLAMC2, otherwise supply EMIN explicitly. \n";
    }
  
    /*
***
*
*        Assume IEEE arithmetic if we found denormalised  numbers above,
*        or if arithmetic seems to round in the  IEEE style,  determined
*        in routine DLAMC1. A true IEEE machine should have both  things
*        true; however, faulty machines may have one or the other.
*
    */
    ieee=ieee||lieee1;
    /*
*
*        Compute  RMIN by successive division by  BETA. We could compute
*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
*        this computation.
*
    */
    lrmin=1.0;
    for(i=1;i<=(1-lemin);i++)
      lrmin=dlamc3<T>(lrmin*rbase,zero);
  
    /* Finally, call DLAMC5 to compute EMAX and RMAX. */
    dlamc5<T>(lbeta,lt,lemin,ieee,lemax,lrmax);
  }
  
  beta=lbeta;
  t=lt;
  rnd=lrnd;
  eps=leps;
  emin=lemin;
  rmin=lrmin;
  emax=lemax;
  rmax=lrmax;

 // std::cout << "got to the end of dlamc2 \n";
  return;
  
  /* End of DLAMC2 */
}

  /* CLN version */

template<typename T>
void dlamc2(int& beta, int& t, bool rnd, T& eps, int& emin, T& rmin, int& emax, T& rmax, int digits)
{
  
  /*
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
  */
  
  /*
*
*  Purpose
*  =======
*
*  DLAMC2 determines the machine parameters specified in its argument
*  list.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  EPS     (output) DOUBLE PRECISION
*          The smallest positive number such that
*
*             fl( 1.0 - EPS ) .LT. 1.0,
*
*          where fl denotes the computed value.
*
*  EMIN    (output) INTEGER
*          The minimum exponent before (gradual) underflow occurs.
*
*  RMIN    (output) DOUBLE PRECISION
*          The smallest normalized number for the machine, given by
*          BASE**( EMIN - 1 ), where  BASE  is the floating point value
*          of BETA.
*
*  EMAX    (output) INTEGER
*          The maximum exponent before overflow occurs.
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest positive number for the machine, given by
*          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
*          value of BETA.
*
*  Further Details
*  ===============
*
*  The computation of  EPS  is based on a routine PARANOIA by
*  W. Kahan of the University of California at Berkeley.
*
* =====================================================================
*
  */
  
  /* Local Scalars */
  static bool first=true;
  static bool iwarn=false;
  bool ieee,lieee1,lrnd;
  int gnmin,gpmin,i,ngnmin,ngpmin;
  static int lemin,lemax,lt,lbeta;
  T a,b,c,half,one,rbase,sixth,small,third,two,zero,three;
  static T lrmax,lrmin,leps;
  
  cln::float_format_t precision = cln::float_format(digits);
    /* Parameters */
//  const cln::cl_F one=cl_float(1,precision);
 // const cln::cl_F zero=cl_float(0,precision);
  /* Executable Statements */
  
  if(first)
  {
    first=false;
    zero=cl_float(0,precision);
    one=cl_float(1,precision);
    two=cl_float(2,precision);
    three=cl_float(3,precision);
  
  /*
*
*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
*        BETA, T, RND, EPS, EMIN and RMIN.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are stored  and not held in registers,  or
*        are not affected by optimizers.
*
*        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
*
  */
    dlamc1<T>(lbeta,lt,lrnd,lieee1,digits);
    
    /* Start to find EPS. */
  
    b=cln::cl_float(lbeta,precision);
    //a=pow(b,-lt);
    a=cln::cl_float(expt(b,-lt),precision);
    leps=a;
  
    /* Try some tricks to see whether or not this is the correct  EPS. */
  
    //b=cln::cl_float(2.0/3.0,precision);
    b=two/three;
    half=one/two;
    sixth=dlamc3<T>(b,-half);
    third=dlamc3<T>(sixth,sixth);
    b=dlamc3<T>(third,-half);
    b=dlamc3<T>(b,sixth);
    b=cln::abs(b);
    if(b<leps)
      b=leps;
    
    leps=one;
  
label10:
    if((leps>b)&&(b>zero))
    {
      leps=b;
      //c=dlamc3<T>(half*leps,pow(two,5)*leps*leps);
      c=dlamc3<T>(half*leps,cl_float(expt(two,5))*leps*leps);
      c=dlamc3<T>(half,-c);
      b=dlamc3<T>(half,c);
      c=dlamc3<T>(half,-b);
      b=dlamc3<T>(half,c);
      goto label10;
    }
  
    if(a<leps)
      leps=a;
  
  //  std::cout << "eps = " << leps << "\n";
    /*
*
*        Computation of EPS complete.
*
*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
*        Keep dividing  A by BETA until (gradual) underflow occurs. This
*        is detected when we cannot recover the previous A.
*
    */
  
    rbase=one/lbeta;
    small=one;
    for(i=1;i<=3;i++)
      small=dlamc3<T>(small*rbase,zero);
    a=dlamc3<T>(one,small);
    dlamc4<T>(ngpmin,one,lbeta,digits);
    dlamc4<T>(ngnmin,-one,lbeta,digits);
    dlamc4<T>(gpmin,a,lbeta,digits);
    dlamc4<T>(gnmin,-a,lbeta,digits);
    ieee=false;
  
    if((ngpmin==ngnmin)&&(gpmin==gnmin))
    {
      if(ngpmin==gpmin)
        lemin=ngpmin;
        /*
*            ( Non twos-complement machines, no gradual underflow;
*              e.g.,  VAX )
        */
      else if((gpmin-ngpmin)==3)
      {
        lemin=ngpmin-1+lt;
        ieee=true;
        /*
*            ( Non twos-complement machines, with gradual underflow;
*              e.g., IEEE standard followers )
*       */
      }
      else
      {
        lemin=std::min(ngpmin,gpmin);
        /* ( A guess; no known machine ) */
        iwarn=true;
      }
    }
    else if((ngpmin==gpmin)&&(ngnmin==gnmin))
    {
      if(abs(ngpmin-ngnmin)==1)
      {
        lemin=std::max(ngpmin,ngnmin);
        /* ( Twos-complement machines, no gradual underflow;
               e.g., CYBER 205 ) */
      }
      else
      {
        lemin=std::min(ngpmin,ngnmin);
        /* ( A guess; no known machine ) */
        iwarn=true;
      }
    }
    else if((abs(ngpmin-ngnmin)==1)&&(gpmin==gnmin))
    {
      if((gpmin-std::min(ngpmin,ngnmin))==3)
      {
        lemin=std::max(ngpmin,ngnmin)-1+lt;
        /*
*            ( Twos-complement machines with gradual underflow;
*              no known machine )
        */
      }
      else
      {
        lemin=std::min(ngpmin,ngnmin);
        /* ( A guess; no known machine ) */
        iwarn=true;
      }
    }
    else
    {
      lemin=std::min(std::min(ngpmin,ngnmin),std::min(gpmin,gnmin));
      /* ( A guess; no known machine ) */
      iwarn=true;
    }
    /*
***
* Comment out this if block if EMIN is ok
    */
    if(iwarn)
    {
      first=true;
      std::cout << "WARNING. The value EMIN may be incorrect:-EMIN = " << lemin << "\n";
      std::cout << "If, after inspection, the value EMIN looks acceptable please comment out \n";
      std::cout << "the IF block as marked within the code of routine DLAMC2, otherwise supply EMIN explicitly. \n";
    }
  
    /*
***
*
*        Assume IEEE arithmetic if we found denormalised  numbers above,
*        or if arithmetic seems to round in the  IEEE style,  determined
*        in routine DLAMC1. A true IEEE machine should have both  things
*        true; however, faulty machines may have one or the other.
*
    */
    ieee=ieee||lieee1;
    /*
*
*        Compute  RMIN by successive division by  BETA. We could compute
*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
*        this computation.
*
    */
    lrmin=one;
    for(i=1;i<=(1-lemin);i++)
      lrmin=dlamc3<T>(lrmin*rbase,zero);
  
    /* Finally, call DLAMC5 to compute EMAX and RMAX. */
    dlamc5<T>(lbeta,lt,lemin,ieee,lemax,lrmax,digits);
  }
  
  beta=lbeta;
  t=lt;
  rnd=lrnd;
  eps=leps;
  emin=lemin;
  rmin=lrmin;
  emax=lemax;
  rmax=lrmax;

 // std::cout << "got to the end of dlamc2 \n";
  return;
  
  /* End of CLN DLAMC2 */
}
    
template<typename T>
T dlamc3(T a,T b)
{
  return a+b;
}

template<typename T>
void dlamc4(int& emin, T start, int base)
{
  /*
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*
*  Purpose
*  =======
*
*  DLAMC4 is a service routine for DLAMC2.
*
*  Arguments
*  =========
*
*  EMIN    (output) EMIN
*          The minimum exponent before (gradual) underflow, computed by
*          setting A = START and dividing by BASE until the previous A
*          can not be recovered.
*
*  START   (input) DOUBLE PRECISION
*          The starting point for determining EMIN.
*
*  BASE    (input) INTEGER
*          The base of the machine.
*
* =====================================================================
*
  */
    
  /* Local Scalars */
  int i;
  T a,b1,b2,c1,c2,d1,d2,one,rbase,zero;
  
  /* Executable Statements */
  
  a=start;
  one=1.0;
  rbase=one/base;
  zero=0;
  emin=1;
  b1=dlamc3<T>(a*rbase,zero);
  c1=a;
  c2=a;
  d1=a;
  d2=a;
  
label10:
  if((c1==a)&&(c2==a)&&(d1==a)&&(d2==a))
  {
    emin-=1;
    a=b1;
    b1=dlamc3<T>(a/base,zero);
    c1=dlamc3<T>(b1*base,zero);
    d1=zero;
    for(i=1;i<=base;i++)
      d1+=b1;
    b2=dlamc3<T>(a*rbase,zero);
    c2=dlamc3<T>(b2/rbase,zero);
    d2=zero;
    for(i=1;i<=base;i++)
      d2+=b2;
    goto label10;
  }
  
  return;
  
  /* End of DLAMC4 */
}

/* CLN version */

template<typename T>
void dlamc4(int& emin, T start, int base, int digits)
{
  /*
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*
*  Purpose
*  =======
*
*  DLAMC4 is a service routine for DLAMC2.
*
*  Arguments
*  =========
*
*  EMIN    (output) EMIN
*          The minimum exponent before (gradual) underflow, computed by
*          setting A = START and dividing by BASE until the previous A
*          can not be recovered.
*
*  START   (input) DOUBLE PRECISION
*          The starting point for determining EMIN.
*
*  BASE    (input) INTEGER
*          The base of the machine.
*
* =====================================================================
*
  */
    
  /* Local Scalars */
  int i;
  T a,b1,b2,c1,c2,d1,d2,one,rbase,zero;
  cln::float_format_t precision=cln::float_format(digits);
  
  /* Executable Statements */
  
  a=start;
  zero=cln::cl_float(0,precision);
  one=cln::cl_float(1,precision);
  //one="1.121591834751928357";
  rbase=one/base;
  cln::print_float(std::cout,cln::default_print_flags,rbase);
  std::cout << "\n";
  emin=1;
  b1=dlamc3<T>(a*rbase,zero);
  c1=a;
  c2=a;
  d1=a;
  d2=a;
  cln::print_float(std::cout,cln::default_print_flags,c1);
  std::cout << "\n";
  cln::print_float(std::cout,cln::default_print_flags,c2);
  std::cout << "\n";
  cln::print_float(std::cout,cln::default_print_flags,d1);
  std::cout << "\n";
  cln::print_float(std::cout,cln::default_print_flags,d2);
  std::cout << "\n";
  
label10:
  if((c1==a)&&(c2==a)&&(d1==a)&&(d2==a))
  {
 //   print_float(std::cout,cln::default_print_flags,a/base);
 //   std::cout << "\n";
    emin-=1;
    a=b1;
    b1=dlamc3<T>(a/base,zero);
    c1=dlamc3<T>(cln::cl_float(b1*base,precision),zero);
    d1=zero;
    for(i=1;i<=base;i++)
      d1+=b1;
    b2=dlamc3<T>(cln::cl_float(a*rbase,precision),zero);
    c2=dlamc3<T>(b2/rbase,zero);
    d2=zero;
    for(i=1;i<=base;i++)
      d2+=b2;
    goto label10;
  }
  
  return;
  
  /* End of CLN DLAMC4 */
}

template<typename T>
void dlamc5(int beta, int p, int emin, bool ieee, int emax, T& rmax)
{
  /*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*  Purpose
*  =======
*
*  DLAMC5 attempts to compute RMAX, the largest machine floating-point
*  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
*  approximately to a power of 2.  It will fail on machines where this
*  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
*  EMAX = 28718).  It will also fail if the value supplied for EMIN is
*  too large (i.e. too close to zero), probably with overflow.
*
*  Arguments
*  =========
*
*  BETA    (input) INTEGER
*          The base of floating-point arithmetic.
*
*  P       (input) INTEGER
*          The number of base BETA digits in the mantissa of a
*          floating-point value.
*
*  EMIN    (input) INTEGER
*          The minimum exponent before (gradual) underflow.
*
*  IEEE    (input) LOGICAL
*          A logical flag specifying whether or not the arithmetic
*          system is thought to comply with the IEEE standard.
*
*  EMAX    (output) INTEGER
*          The largest exponent before overflow
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest machine floating-point number.
*
* =====================================================================
*
  */
  
  const T zero=0.0;
  const T one=1.0;
  
  /* Local Scalars */
  int exbits,expsum,i,lexp,nbits,ctry,uexp; // note: ctry = try in original
  T oldy,recbas,y,z;
  
  /*
*     .. Executable Statements ..
*
*     First compute LEXP and UEXP, two powers of 2 that bound
*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
*     approximately to the bound that is closest to abs(EMIN).
*     (EMAX is the exponent of the required number RMAX).
*
  */
  
  lexp=1;
  exbits=1;
label10:
  ctry=lexp*2;
 // std::cout << "ctry = " << ctry << "\n";
  if(ctry<=(-emin))
  {
 //   std::cout << "ctry = " << ctry << "\n";
    lexp=ctry;
    exbits+=1;
    goto label10;
  }
  if(lexp==(-emin))
    uexp=lexp;
  else
  {
    uexp=ctry;
    exbits+=1;
  }
  
//  std::cout << "got to here \n";
  
  /*
*
*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
*     than or equal to EMIN. EXBITS is the number of bits needed to
*     store the exponent.
*
  */
  
  if((uexp+emin)>(-lexp-emin))
    expsum=2*lexp;
  else
    expsum=2*uexp;
  
  /*
*
*     EXPSUM is the exponent range, approximately equal to
*     EMAX - EMIN + 1 .
*
  */
  
  emax=expsum+emin-1;
  nbits=1+exbits+p;
  
  /*
*
*     NBITS is the total number of bits needed to store a
*     floating-point number.
*
  */
  
  if(((nbits%2)==1)&&(beta==2))
  {
    /*
*
*        Either there are an odd number of bits used to store a
*        floating-point number, which is unlikely, or some bits are
*        not used in the representation of numbers, which is possible,
*        (e.g. Cray machines) or the mantissa has an implicit bit,
*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
*        most likely. We have to assume the last alternative.
*        If this is true, then we need to reduce EMAX by one because
*        there must be some way of representing zero in an implicit-bit
*        system. On machines like Cray, we are reducing EMAX by one
*        unnecessarily.
*
    */
    
    emax-=1;
  }
  
  if(ieee)
  {
    /*
*
*        Assume we are on an IEEE machine which reserves one exponent
*        for infinity and NaN.
*
    */
    emax-=1;
  }
  
  /*
*
*     Now create RMAX, the largest machine number, which should
*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
*
*     First compute 1.0 - BETA**(-P), being careful that the
*     result is less than 1.0 .
*
  */
  recbas=one/beta;
  z=beta-one;
  y=zero;
  for(i=1;i<=p;i++)
  {
    z=z*recbas;
    if(y<one)
      oldy=y;
    y=dlamc3<T>(y,z);
  }
  if(y>=one)
    y=oldy;
  
  /*
*
*     Now multiply by BETA**EMAX to get RMAX.
*
  */
  for(i=1;i<=emax;i++)
  {
 //   std::cout << "i=" << i << "\n";
    y=dlamc3<T>(y*beta,zero);
  }
  
  rmax=y;
  
  /* End of dlamc5 */
  
  return;
}

/* CLN version */

template<typename T>
void dlamc5(int beta, int p, int emin, bool ieee, int emax, T& rmax, int digits)
{
  /*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*  Purpose
*  =======
*
*  DLAMC5 attempts to compute RMAX, the largest machine floating-point
*  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
*  approximately to a power of 2.  It will fail on machines where this
*  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
*  EMAX = 28718).  It will also fail if the value supplied for EMIN is
*  too large (i.e. too close to zero), probably with overflow.
*
*  Arguments
*  =========
*
*  BETA    (input) INTEGER
*          The base of floating-point arithmetic.
*
*  P       (input) INTEGER
*          The number of base BETA digits in the mantissa of a
*          floating-point value.
*
*  EMIN    (input) INTEGER
*          The minimum exponent before (gradual) underflow.
*
*  IEEE    (input) LOGICAL
*          A logical flag specifying whether or not the arithmetic
*          system is thought to comply with the IEEE standard.
*
*  EMAX    (output) INTEGER
*          The largest exponent before overflow
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest machine floating-point number.
*
* =====================================================================
*
  */
  
  cln::float_format_t precision=cln::float_format(digits);
  const T zero = cln::cl_float(0,precision);
  const T one = cln::cl_float(1,precision);
  
  /* Local Scalars */
  int exbits,expsum,i,lexp,nbits,ctry,uexp; // note: ctry = try in original
  T oldy,recbas,y,z;
  
  /*
*     .. Executable Statements ..
*
*     First compute LEXP and UEXP, two powers of 2 that bound
*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
*     approximately to the bound that is closest to abs(EMIN).
*     (EMAX is the exponent of the required number RMAX).
*
  */
  
  lexp=1;
  exbits=1;
label10:
  ctry=lexp*2;
 // std::cout << "ctry = " << ctry << "\n";
  if(ctry<=(-emin))
  {
 //   std::cout << "ctry = " << ctry << "\n";
    lexp=ctry;
    exbits+=1;
    goto label10;
  }
  if(lexp==(-emin))
    uexp=lexp;
  else
  {
    uexp=ctry;
    exbits+=1;
  }
  
//  std::cout << "got to here \n";
  
  /*
*
*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
*     than or equal to EMIN. EXBITS is the number of bits needed to
*     store the exponent.
*
  */
  
  if((uexp+emin)>(-lexp-emin))
    expsum=2*lexp;
  else
    expsum=2*uexp;
  
  /*
*
*     EXPSUM is the exponent range, approximately equal to
*     EMAX - EMIN + 1 .
*
  */
  
  emax=expsum+emin-1;
  nbits=1+exbits+p;
  
  /*
*
*     NBITS is the total number of bits needed to store a
*     floating-point number.
*
  */
  
  if(((nbits%2)==1)&&(beta==2))
  {
    /*
*
*        Either there are an odd number of bits used to store a
*        floating-point number, which is unlikely, or some bits are
*        not used in the representation of numbers, which is possible,
*        (e.g. Cray machines) or the mantissa has an implicit bit,
*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
*        most likely. We have to assume the last alternative.
*        If this is true, then we need to reduce EMAX by one because
*        there must be some way of representing zero in an implicit-bit
*        system. On machines like Cray, we are reducing EMAX by one
*        unnecessarily.
*
    */
    
    emax-=1;
  }
  
  if(ieee)
  {
    /*
*
*        Assume we are on an IEEE machine which reserves one exponent
*        for infinity and NaN.
*
    */
    emax-=1;
  }
  
  /*
*
*     Now create RMAX, the largest machine number, which should
*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
*
*     First compute 1.0 - BETA**(-P), being careful that the
*     result is less than 1.0 .
*
  */
  recbas=one/beta;
  z=beta-one;
  y=zero;
  for(i=1;i<=p;i++)
  {
    z=z*recbas;
    if(y<one)
      oldy=y;
    y=dlamc3<T>(y,z);
  }
  if(y>=one)
    y=oldy;
  
  /*
*
*     Now multiply by BETA**EMAX to get RMAX.
*
  */
  for(i=1;i<=emax;i++)
  {
 //   std::cout << "i=" << i << "\n";
    y=dlamc3<T>(cln::cl_float(y*beta,precision),zero);
  }
  
  rmax=y;
  
  /* End of CLN dlamc5 */
  
  return;
}
    
    
template<typename T>
T ddot(int n,T* dx, int incx, T* dy, int incy)
{
  /* Comment from Fortran code:
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
 */
  
  T result;
  int i,ix,iy,m,mp1;
  
  result=0.0;
  
  if(n<=0) return result;
  if((incx==1)&&(incy==1)) goto label20;
  
  /*
c
c        code for unequal increments or equal increments
c          not equal to 1
c
  */
  
  ix=1;
  iy=1;
  if(incx<0)ix=(-n+1)*incx+1;
  if(incy<0)iy=(-n+1)*incy+1;
  for(i=1;i<=n;i++)
  {
    result+=dx[ix-1]*dy[iy-1];
    ix+=incx;
    iy+=incy;
  }
  
  return result;
  
  /*
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
  */
  
label20:
  m=n%5;
  if(m==0) goto label40;
  for(i=1;i<=m;i++)
    result+=dx[i-1]*dy[i-1];
  if(n<5) return result;
label40:
  mp1=m+1;
  for(i=mp1;i<=n;i=i+5)
    result+=dx[i-1]*dy[i-1]+dx[i]*dy[i]+dx[i+1]*dy[i+1]+dx[i+2]*dy[i+2]+dx[i+3]*dy[i+3];
  
  return result;
}
 
/* CLN version */ 

template<typename T>
T ddot(int n,T* dx, int incx, T* dy, int incy, int digits)
{
  /* Comment from Fortran code:
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
 */
  cln::float_format_t precision=cln::float_format(digits);
  T result;
  int i,ix,iy,m,mp1;
  
  result=cln::cl_float(0,precision);
  
  if(n<=0) return result;
  if((incx==1)&&(incy==1)) goto label20;
  
  /*
c
c        code for unequal increments or equal increments
c          not equal to 1
c
  */
  
  ix=1;
  iy=1;
  if(incx<0)ix=(-n+1)*incx+1;
  if(incy<0)iy=(-n+1)*incy+1;
  for(i=1;i<=n;i++)
  {
    result+=dx[ix-1]*dy[iy-1];
    ix+=incx;
    iy+=incy;
  }
  
  return result;
  
  /*
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
  */
  
label20:
  m=n%5;
  if(m==0) goto label40;
  for(i=1;i<=m;i++)
    result+=dx[i-1]*dy[i-1];
  if(n<5) return result;
label40:
  mp1=m+1;
  for(i=mp1;i<=n;i=i+5)
    result+=dx[i-1]*dy[i-1]+dx[i]*dy[i]+dx[i+1]*dy[i+1]+dx[i+2]*dy[i+2]+dx[i+3]*dy[i+3];
  
  return result;
}
  
template<typename T>
void dvout(int lout, int n, T* sx, int idigit, const std::string& ifmt)
{
  int i,j;
  int length;
  int upper;
  
  std::cout.precision(10);
  length = ifmt.length();
  std::cout << ifmt << "\n";
  for(i=0;i<length;i++)
    std::cout << "-";
  std::cout << "\n";
  j=0;
  while(j<n)
  {
    upper=std::min(n-1,j+9);
    std::cout << j << " - " << upper << ": ";
    for(i=j;i<=upper;i++)
    {
      std::cout << std::setw(20) << sx[i];
    }
    std::cout << "\n";
    j=j+10;
  }
  
  std::cout << "\n";
  return;
}

template<>
void dvout<cln::cl_F>(int lout, int n, cln::cl_F* sx, int idigit, const std::string& ifmt)
{
  int i,j;
  int length;
  int upper;
  
  cln::float_format_t precision = cln::float_format(cln::float_digits(sx[0]) );
  
  length = ifmt.length();
  std::cout << ifmt << "\n";
  for(i=0;i<length;i++)
    std::cout << "-";
  std::cout << "\n";
  j=0;
  while(j<n)
  {
    upper=std::min(n-1,j+9);
    std::cout << j << " - " << upper << ": ";
    for(i=j;i<=upper;i++)
    {
      std::cout << std::setw(15) << sx[i];
    }
    std::cout << "\n";
    j=j+10;
  }
  
  std::cout << "\n";
  return;
}

template<typename T>
void ivout(int lout, int n, int* ix, int idigit, const std::string ifmt)
{
  int i,j;
  int length;
  int upper;
      
  length = ifmt.length();
  std::cout << ifmt << "\n";
  for(i=0;i<length;i++)
    std::cout << "-";
  std::cout << "\n";
  j=0;
  while(j<n)
  {
    upper=std::min(n-1,j+9);
    std::cout << j << " - " << upper << ": ";
    for(i=j;i<=upper;i++)
    {
      std::cout << std::setw(10) << ix[i];
    }
    std::cout << "\n";
    j=j+10;
  }
  
  std::cout << "\n";
  return;
}

template<typename T>
void dmout(int lout, int m, int n, T* A, int lda, int idigit, const std::string& ifmt)
{
  int i,j;
  int length;
  
  length = ifmt.length();
  std::cout << ifmt << "\n";
  for(i=0;i<length;i++)
    std::cout << "-";
  std::cout << "\n";
  
  std::cout << "          ";
  for(j=1;j<=m;j++)
    std::cout << std::setw(15) << "Col " << j;
  
  std::cout << "\n";
  
  for(i=1;i<=n;i++)
  {
    std::cout << "Row " << i << ":";
    for(j=1;j<=m;j++)
    {
      std::cout << std::setw(15) << A[i-1+lda*(j-1)];
 //     cln::print_float(std::cout,cln::default_print_flags,A[i-1+lda*(j-1)]);
 //     std::cout << "\n";
    }
    std::cout << "\n";
  }
  
  std::cout << "\n";
  
  return;
}

template<typename T>
T dlapy2(T x, T y)
{
  /* adapted from LAPACK:
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
* 
*  Purpose
*  =======
*
*  DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
*  overflow.
*
*  Arguments
*  =========
*
*  X       (input) DOUBLE PRECISION
*  Y       (input) DOUBLE PRECISION
*          X and Y specify the values x and y.
*
*  =====================================================================
*
  */
  
  const T zero=0.0;
  const T one=1.0;
  T w,xabs,yabs,z;
  T result;
  
  xabs=fabs(x);
  yabs=fabs(y);
  w=std::max(xabs,yabs);
  z=std::min(xabs,yabs);
  if(z==zero)
    result=w;
  else
    result=w*sqrt(one+(z/w)*(z/w));
  
  //return sqrt(w*w+z*z);
  return result;
}

template<>
cln::cl_F dlapy2<cln::cl_F>(cln::cl_F x, cln::cl_F y)
{
  /* adapted from LAPACK:
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
* 
*  Purpose
*  =======
*
*  DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
*  overflow.
*
*  Arguments
*  =========
*
*  X       (input) DOUBLE PRECISION
*  Y       (input) DOUBLE PRECISION
*          X and Y specify the values x and y.
*
*  =====================================================================
*
  */
  
  cln::float_format_t precision = cln::float_format(cln::float_digits(x) );

  const cln::cl_F zero=cln::cl_float(0,precision);
  const cln::cl_F one=cln::cl_float(1,precision);
  cln::cl_F w,xabs,yabs,z;
  cln::cl_F result;
  
  xabs=cln::abs(x);
  yabs=cln::abs(y);
  w=std::max(xabs,yabs);
  z=std::min(xabs,yabs);
  if(z==zero)
    result=w;
  else
    result=w*cln::sqrt(one+(z/w)*(z/w));
  
  //return sqrt(w*w+z*z);
  return result;
}

void second(float& t)
{
  time_t t1;
  
  t1=time(NULL);
  
  t=(float) t1;
  return;
}

template<typename T>
T dnrm2(int n, T* dx, int incx)
{
  T sum;
  int i;
  
  T zero=0.0;
 // T cutlo = 0.1415686533102923e-145;
//  T cuthi = 0.1340780792994260e+155;
  
  if(n<=0) return 0;
  
  sum=zero;
  for(i=0;i<n;i+=incx)
    sum+=dx[i]*dx[i];
  
  return sqrt(sum);
}
 /* CLN version */
template<>
cln::cl_F dnrm2(int n, cln::cl_F* dx, int incx)
{
  cln::cl_F sum;
  int i;
//  cln::float_format_t precision=cln::float_format(digits);
  cln::float_format_t precision = cln::float_format(cln::float_digits(dx[0]) );
  
  cln::cl_F zero = cln::cl_float(0,precision);
  
//  std::cout << "dnrm2 digits = " << cln::float_digits(dx[0]) << "\n";
//  if(cln::float_digits(dx[0])==53)
//    std::cout << "too few digits \n";
  
  if(n<=0) return zero;
  
  sum=zero;
  for(i=0;i<n;i+=incx)
    sum+=dx[i]*dx[i];
  
  return cln::sqrt(sum);
}

/* Populate array x with uniform random numbers */
/* in the range [-1,1]                          */
template<typename T>
void dlarnv(int* iseed, int n, T* x)
{
  int i;
  double randomnum;
  //std::default_random_engine generator(iseed[0]);
 // std::uniform_real_distribution<double> distribution(-1.0,1.0);
  
  for(i=0;i<n;i++)
  {
  //  randomnum=distribution(generator);
    x[i]=(T) 2.0*((T) rand())/((T) RAND_MAX)-1.0;
   // std::cout << "x[" << i << "] = " << x[i] << "\n";
  }
}

/* CLN version */
template<typename T>
void dlarnv(int* iseed, int n, T* x, int digits)
{
  int i;
  double randomnum;
  //std::default_random_engine generator(iseed[0]);
 // std::uniform_real_distribution<double> distribution(-1.0,1.0);
  cln::float_format_t precision=cln::float_format(digits);
  
  for(i=0;i<n;i++)
  {
  //  randomnum=distribution(generator);
    x[i]=cln::cl_float((double) 2.0*((double) rand())/((double) RAND_MAX)-1.0,precision);
   // std::cout << "x[" << i << "] = " << x[i] << "\n";
  }
}

template<typename T>
void dlacpy(const std::string& uplo, int m, int n, T* A, int lda, T* B, int ldb)
{
  /*
   *
*  Purpose
*  =======
*
*  DLACPY copies all or part of a two-dimensional matrix A to another
*  matrix B.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies the part of the matrix A to be copied to B.
*          = 'U':      Upper triangular part
*          = 'L':      Lower triangular part
*          Otherwise:  All of the matrix A
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The m by n matrix A.  If UPLO = 'U', only the upper triangle
*          or trapezoid is accessed; if UPLO = 'L', only the lower
*          triangle or trapezoid is accessed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  B       (output) DOUBLE PRECISION array, dimension (LDB,N)
*          On exit, B = A in the locations specified by UPLO.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,M).
*
*  =====================================================================
* */
  
  int i,j;
  
  if(uplo=="U")
  {
    for(j=1;j<=n;j++)
      for(i=1;i<=std::min(j,m);i++)
        B[i-1+ldb*(j-1)]=A[i-1+lda*(j-1)];
  }
  else if(uplo=="L")
  {
    for(j=1;j<=n;j++)
      for(i=j;i<=m;i++)
        B[i-1+ldb*(j-1)]=A[i-1+lda*(j-1)];
  }
  else
  {
    for(j=1;j<=n;j++)
      for(i=1;i<=m;i++)
        B[i-1+ldb*(j-1)]=A[i-1+lda*(j-1)];
  }
  
  return;
  
  /*
*
*     End of DLACPY
*
  */
}

template<typename T>
T dlanhs(const std::string& norm, int n, T* A, int lda, T* work)
{
  /*
*
*  Purpose
*  =======
*
*  DLANHS  returns the value of the one norm,  or the Frobenius norm, or
*  the  infinity norm,  or the  element of  largest absolute value  of a
*  Hessenberg matrix A.
*
*  Description
*  ===========
*
*  DLANHS returns the value
*
*     DLANHS = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*              (
*              ( norm1(A),         NORM = '1', 'O' or 'o'
*              (
*              ( normI(A),         NORM = 'I' or 'i'
*              (
*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*
*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normF  denotes the  Frobenius norm of a matrix (square root of sum of
*  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies the value to be returned in DLANHS as described
*          above.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.  When N = 0, DLANHS is
*          set to zero.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The n by n upper Hessenberg matrix A; the part of A below the
*          first sub-diagonal is not referenced.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(N,1).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),
*          where LWORK >= N when NORM = 'I'; otherwise, WORK is not
*          referenced.
*
* =====================================================================
*
*   
  */
  int i,j;
  T scale, sum, value;
  T zero=0.0;
  T one=1.0;
  
  if(n==0)
    value=zero;
  else if(norm=="M")
  {
    /* Find max(abs(A(i,j))). */
    value=zero;
    for(j=1;j<=n;j++)
      for(i=1;i<=std::min(n,j+1);i++)
        value=std::max(value,(T) fabs(A[i-1+lda*(j-1)]) );
  }
  else if( (norm=="O")||(norm=="1") )
  {
    /* Find norm1(A). */
    value=zero;
    for(j=1;j<=n;j++)
    {
      sum=zero;
      for(i=1;i<=std::min(n,j+1);i++)
        sum+=fabs(A[i-1+lda*(j-1)]);
      value=std::max(value,sum);
    }
  }
  else if(norm=="I")
  {
    /* Find normI(A). */
    for(i=1;i<=n;i++)
      work[i-1]=zero;
    for(j=1;j<=n;j++)
    {
      for(i=1;i<=std::min(n,j+1);i++)
        work[i-1]+=fabs(A[i-1+lda*(j-1)]);
    }
    value=zero;
    for(i=1;i<=n;i++)
      value=std::max(value,work[i-1]);
  }
  else if( (norm=="F")||(norm=="E") )
  {
    /* Find normF(A). */
    scale=zero;
    sum=one;
    for(j=1;j<=n;j++)
      dlassq<T>(std::min(n,j+1),&A[lda*(j-1)],1,scale,sum);
    value=scale*sqrt(sum);
  }
  
  return value;
}

/* CLN version */

template<typename T>
T dlanhs(const std::string& norm, int n, T* A, int lda, T* work, int digits)
{
  /*
*
*  Purpose
*  =======
*
*  DLANHS  returns the value of the one norm,  or the Frobenius norm, or
*  the  infinity norm,  or the  element of  largest absolute value  of a
*  Hessenberg matrix A.
*
*  Description
*  ===========
*
*  DLANHS returns the value
*
*     DLANHS = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*              (
*              ( norm1(A),         NORM = '1', 'O' or 'o'
*              (
*              ( normI(A),         NORM = 'I' or 'i'
*              (
*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*
*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normF  denotes the  Frobenius norm of a matrix (square root of sum of
*  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies the value to be returned in DLANHS as described
*          above.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.  When N = 0, DLANHS is
*          set to zero.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The n by n upper Hessenberg matrix A; the part of A below the
*          first sub-diagonal is not referenced.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(N,1).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),
*          where LWORK >= N when NORM = 'I'; otherwise, WORK is not
*          referenced.
*
* =====================================================================
*
*   
  */
  int i,j;
  T scale, sum, value;
  cln::float_format_t precision=cln::float_format(digits);
  const T zero = cln::cl_float(0,precision);
  const T one = cln::cl_float(1,precision);
  
  if(n==0)
    value=zero;
  else if(norm=="M")
  {
    /* Find max(abs(A(i,j))). */
    value=zero;
    for(j=1;j<=n;j++)
      for(i=1;i<=std::min(n,j+1);i++)
        value=cln::max(value,cln::abs(A[i-1+lda*(j-1)]) );
  }
  else if( (norm=="O")||(norm=="1") )
  {
    /* Find norm1(A). */
    value=zero;
    for(j=1;j<=n;j++)
    {
      sum=zero;
      for(i=1;i<=std::min(n,j+1);i++)
        sum+=cln::abs(A[i-1+lda*(j-1)]);
      value=cln::max(value,sum);
    }
  }
  else if(norm=="I")
  {
    /* Find normI(A). */
    for(i=1;i<=n;i++)
      work[i-1]=zero;
    for(j=1;j<=n;j++)
    {
      for(i=1;i<=std::min(n,j+1);i++)
        work[i-1]+=cln::abs(A[i-1+lda*(j-1)]);
    }
    value=zero;
    for(i=1;i<=n;i++)
      value=cln::max(value,work[i-1]);
  }
  else if( (norm=="F")||(norm=="E") )
  {
    /* Find normF(A). */
    scale=zero;
    sum=one;
    for(j=1;j<=n;j++)
      dlassq<T>(std::min(n,j+1),&A[lda*(j-1)],1,scale,sum,digits);
    value=scale*sqrt(sum);
  }
  
  return value;
}

template<typename T>
void dscal(int n, T da, T* dx, int incx)
{
  /* original comment:
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified to correct problem with negative increment, 8/21/90.
c
  */
  int i,ix,m,mp1;
  
  if(n<=0)
    return;
  if(incx==1) 
    goto label20;
  /*
c
c        code for increment not equal to 1
c
  */
  ix=1;
  if(incx<0)
    ix=(-n+1)*incx+1;
  for(i=1;i<=n;i++)
  {
    dx[ix-1]=da*dx[ix-1];
    ix+=incx;
  }
  
  return;
  
  /*
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
  */
label20:
  m=n%5;
  if(m==0)
    goto label40;
  
  for(i=1;i<=m;i++)
    dx[i-1]=da*dx[i-1];
  if(n<5)
    return;
label40:
  
  mp1=m+1;
  for(i=mp1;i<=n;i+=5)
  {
    dx[i-1]=da*dx[i-1];
    dx[i]=da*dx[i];
    dx[i+1]=da*dx[i+1];
    dx[i+2]=da*dx[i+2];
    dx[i+3]=da*dx[i+3];
  }
  
  return;
}
          
template<typename T>
void dlarfg(int n, T& alpha, T* x, int incx, T& tau)
{
  /*
*
*  Purpose
*  =======
*
*  DLARFG generates a real elementary reflector H of order n, such
*  that
*
*        H * ( alpha ) = ( beta ),   H' * H = I.
*            (   x   )   (   0  )
*
*  where alpha and beta are scalars, and x is an (n-1)-element real
*  vector. H is represented in the form
*
*        H = I - tau * ( 1 ) * ( 1 v' ) ,
*                      ( v )
*
*  where tau is a real scalar and v is a real (n-1)-element
*  vector.
*
*  If the elements of x are all zero, then tau = 0 and H is taken to be
*  the unit matrix.
*
*  Otherwise  1 <= tau <= 2.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the elementary reflector.
*
*  ALPHA   (input/output) DOUBLE PRECISION
*          On entry, the value alpha.
*          On exit, it is overwritten with the value beta.
*
*  X       (input/output) DOUBLE PRECISION array, dimension
*                         (1+(N-2)*abs(INCX))
*          On entry, the vector x.
*          On exit, it is overwritten with the vector v.
*
*  INCX    (input) INTEGER
*          The increment between elements of X. INCX > 0.
*
*  TAU     (output) DOUBLE PRECISION
*          The value tau.
*
*  =====================================================================
  */
  int j,knt;
  T beta,rsafmn,safmin,xnorm;
  T zero=0.0;
  T one=1.0;
  
  /* .. Executable Statements .. */
  if(n<=1)
  {
    tau=zero;
    return;
  }
  
  xnorm=dnrm2<T>(n-1,x,incx);
 // std::cout << "xnorm = " << xnorm << "\n";
  
  if(xnorm==zero)
    tau=zero;
  else
  {
    /* general case */
    beta=-sign<T>(dlapy2<T>(alpha,xnorm),alpha);
 //   std::cout << "beta = " << beta << "\n";
    safmin=dlamch<T>("S")/dlamch<T>("E");
    if(fabs(beta)<safmin)
    {
      /* XNORM, BETA may be inaccurate; scale X and recompute them */
      rsafmn=one/safmin;
      knt=0;
label10:
      knt+=1;
      dscal<T>(n-1,rsafmn,x,incx);
      beta=beta*rsafmn;
      alpha=alpha*rsafmn;
      if(fabs(beta)<safmin)
        goto label10;
      
      /* New BETA is at most 1, at least SAFMIN */
      
      xnorm=dnrm2<T>(n-1,x,incx);
      beta=-sign<T>(dlapy2<T>(alpha,xnorm),alpha);
      tau=(beta-alpha)/beta;
      dscal<T>(n-1,one/(alpha-beta),x,incx);
      
      /* If ALPHA is subnormal, it may lose relative accuracy */
      
      alpha=beta;
      for(j=1;j<=knt;j++)
        alpha=alpha*safmin;
    }
    else
    {
      tau=(beta-alpha)/beta;
      dscal<T>(n-1,one/(alpha-beta),x,incx);
      alpha=beta;
    }
  }
  
  return;
  
  /* End of DLARFG */
  
}

/* CLN version */

template<typename T>
void dlarfg(int n, T& alpha, T* x, int incx, T& tau, int digits)
{
  /*
*
*  Purpose
*  =======
*
*  DLARFG generates a real elementary reflector H of order n, such
*  that
*
*        H * ( alpha ) = ( beta ),   H' * H = I.
*            (   x   )   (   0  )
*
*  where alpha and beta are scalars, and x is an (n-1)-element real
*  vector. H is represented in the form
*
*        H = I - tau * ( 1 ) * ( 1 v' ) ,
*                      ( v )
*
*  where tau is a real scalar and v is a real (n-1)-element
*  vector.
*
*  If the elements of x are all zero, then tau = 0 and H is taken to be
*  the unit matrix.
*
*  Otherwise  1 <= tau <= 2.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the elementary reflector.
*
*  ALPHA   (input/output) DOUBLE PRECISION
*          On entry, the value alpha.
*          On exit, it is overwritten with the value beta.
*
*  X       (input/output) DOUBLE PRECISION array, dimension
*                         (1+(N-2)*abs(INCX))
*          On entry, the vector x.
*          On exit, it is overwritten with the vector v.
*
*  INCX    (input) INTEGER
*          The increment between elements of X. INCX > 0.
*
*  TAU     (output) DOUBLE PRECISION
*          The value tau.
*
*  =====================================================================
  */
  int j,knt;
  T beta,rsafmn,safmin,xnorm;
  cln::float_format_t precision=cln::float_format(digits);
  const T zero = cln::cl_float(0,precision);
  const T one = cln::cl_float(1,precision);
  
  /* .. Executable Statements .. */
  if(n<=1)
  {
    tau=zero;
    return;
  }
  
  xnorm=dnrm2<T>(n-1,x,incx);
 // std::cout << "xnorm = " << xnorm << "\n";
  
  if(xnorm==zero)
    tau=zero;
  else
  {
    /* general case */
    beta=-sign<T>(dlapy2<T>(alpha,xnorm),alpha);
 //   std::cout << "beta = " << beta << "\n";
    safmin=dlamch<T>("S",digits)/dlamch<T>("E",digits);
    if(cln::abs(beta)<safmin)
    {
      /* XNORM, BETA may be inaccurate; scale X and recompute them */
      rsafmn=one/safmin;
      knt=0;
label10:
      knt+=1;
      dscal<T>(n-1,rsafmn,x,incx);
      beta=beta*rsafmn;
      alpha=alpha*rsafmn;
      if(cln::abs(beta)<safmin)
        goto label10;
      
      /* New BETA is at most 1, at least SAFMIN */
      
      xnorm=dnrm2<T>(n-1,x,incx);
      beta=-sign<T>(dlapy2<T>(alpha,xnorm),alpha);
      tau=(beta-alpha)/beta;
      dscal<T>(n-1,one/(alpha-beta),x,incx);
      
      /* If ALPHA is subnormal, it may lose relative accuracy */
      
      alpha=beta;
      for(j=1;j<=knt;j++)
        alpha=alpha*safmin;
    }
    else
    {
      tau=(beta-alpha)/beta;
      dscal<T>(n-1,one/(alpha-beta),x,incx);
      alpha=beta;
    }
  }
  
  return;
  
  /* End of DLARFG */
  
}
    
template<typename T>
void dlartg(T f, T g, T& cs, T& sn, T& r)
{
  /*
*  Purpose
*  =======
*
*  DLARTG generate a plane rotation so that
*
*     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
*     [ -SN  CS  ]     [ G ]     [ 0 ]
*
*  This is a slower, more accurate version of the BLAS1 routine DROTG,
*  with the following other differences:
*     F and G are unchanged on return.
*     If G=0, then CS=1 and SN=0.
*     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
*        floating point operations (saves work in DBDSQR when
*        there are zeros on the diagonal).
*
*  If F exceeds G in magnitude, CS will be positive.
*
*  Arguments
*  =========
*
*  F       (input) DOUBLE PRECISION
*          The first component of vector to be rotated.
*
*  G       (input) DOUBLE PRECISION
*          The second component of vector to be rotated.
*
*  CS      (output) DOUBLE PRECISION
*          The cosine of the rotation.
*
*  SN      (output) DOUBLE PRECISION
*          The sine of the rotation.
*
*  R       (output) DOUBLE PRECISION
*          The nonzero component of the rotated vector.
*
*  =====================================================================
*
  */
  static bool first=true;
  int count,i;
  T eps,f1,g1,scale;
  static T safmx2,safmin,safmn2;
  T zero=0.0;
  T one=1.0;
  T two=2.0;
  
  /* .. Executable Statements .. */
  if(first)
  {
    first=false;
    safmin=dlamch<T>("S");
    eps=dlamch<T>("E");
    safmn2=pow(dlamch<T>("B"),int(log(safmin/eps)/log(dlamch<T>("B"))/two) );
    safmx2=one/safmn2;
  }
  if(g==zero)
  {
    cs=one;
    sn=zero;
    r=f;
  }
  else if(f==zero)
  {
    cs=zero;
    sn=one;
    r=g;
  }
  else
  {
    f1=f;
    g1=g;
    scale=std::max(abs(f1),abs(g1));
    if(scale>=safmx2)
    {
      count=0;
label10:
      count+=1;
      f1=f1*safmn2;
      g1=g1*safmn2;
      scale=std::max(fabs(f1),fabs(g1));
      if(scale>=safmx2)
        goto label10;
      r=sqrt(f1*f1+g1*g1);
      cs=f1/r;
      sn=g1/r;
      for(i=1;i<=count;i++)
        r=r*safmx2;
    }
    else if(scale<=safmn2)
    {
      count=0;
label30:
      count+=1;
      f1=f1*safmx2;
      g1=g1*safmx2;
      scale=std::max(fabs(f1),fabs(g1));
      if(scale<=safmn2)
        goto label30;
      r=sqrt(f1*f1+g1*g1);
      cs=f1/r;
      sn=g1/r;
      for(i=1;i<=count;i++)
        r=r*safmn2;
    }
    else
    {
      r=sqrt(f1*f1+g1*g1);
      cs=f1/r;
      sn=g1/r;
    }
    if((fabs(f)>fabs(g)) && (cs<zero) )
    {
      cs=-cs;
      sn=-sn;
      r=-r;
    }
  }
  
  /* End of DLARTG */
  
  return;
}

/* CLN dlartg */

template<>
void dlartg<cln::cl_F>(cln::cl_F f, cln::cl_F g, cln::cl_F& cs, cln::cl_F& sn, cln::cl_F& r)
{
  /*
*  Purpose
*  =======
*
*  DLARTG generate a plane rotation so that
*
*     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
*     [ -SN  CS  ]     [ G ]     [ 0 ]
*
*  This is a slower, more accurate version of the BLAS1 routine DROTG,
*  with the following other differences:
*     F and G are unchanged on return.
*     If G=0, then CS=1 and SN=0.
*     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
*        floating point operations (saves work in DBDSQR when
*        there are zeros on the diagonal).
*
*  If F exceeds G in magnitude, CS will be positive.
*
*  Arguments
*  =========
*
*  F       (input) DOUBLE PRECISION
*          The first component of vector to be rotated.
*
*  G       (input) DOUBLE PRECISION
*          The second component of vector to be rotated.
*
*  CS      (output) DOUBLE PRECISION
*          The cosine of the rotation.
*
*  SN      (output) DOUBLE PRECISION
*          The sine of the rotation.
*
*  R       (output) DOUBLE PRECISION
*          The nonzero component of the rotated vector.
*
*  =====================================================================
*
  */
  static bool first=true;
  int count,i;
  cln::cl_F eps,f1,g1,scale;
  static cln::cl_F safmx2,safmin,safmn2;
  
  int digits = cln::float_digits(f);
  cln::float_format_t precision=cln::float_format(digits );
  const cln::cl_F zero = cln::cl_float(0,precision);
  const cln::cl_F one = cln::cl_float(1,precision);
  const cln::cl_F two = cln::cl_float(2,precision);
  cln::cl_I expnt;
  
  /* .. Executable Statements .. */
  if(first)
  {
    first=false;
    safmin=dlamch<cln::cl_F>("S",digits);
    eps=dlamch<cln::cl_F>("E",digits);
    expnt = cln::floor1(cln::ln(safmin/eps)/cln::ln(dlamch<cln::cl_F>("B",digits))/two);
    safmn2= cln::cl_float( cln::expt(dlamch<cln::cl_F>("B",digits),expnt),precision );
   // safmn2=pow(dlamch<cln::cl_F>("B",digits),int(log(safmin/eps)/log(dlamch<cln::cl_F>("B",digits))/two) );
    safmx2=one/safmn2;
  }
  if(g==zero)
  {
    cs=one;
    sn=zero;
    r=f;
  }
  else if(f==zero)
  {
    cs=zero;
    sn=one;
    r=g;
  }
  else
  {
    f1=f;
    g1=g;
    scale=std::max(cln::abs(f1),cln::abs(g1));
    if(scale>=safmx2)
    {
      count=0;
label10:
      count+=1;
      f1=f1*safmn2;
      g1=g1*safmn2;
      scale=std::max(cln::abs(f1),cln::abs(g1));
      if(scale>=safmx2)
        goto label10;
      r=sqrt(f1*f1+g1*g1);
      cs=f1/r;
      sn=g1/r;
      for(i=1;i<=count;i++)
        r=r*safmx2;
    }
    else if(scale<=safmn2)
    {
      count=0;
label30:
      count+=1;
      f1=f1*safmx2;
      g1=g1*safmx2;
      scale=std::max(cln::abs(f1),cln::abs(g1));
      if(scale<=safmn2)
        goto label30;
      r=sqrt(f1*f1+g1*g1);
      cs=f1/r;
      sn=g1/r;
      for(i=1;i<=count;i++)
        r=r*safmn2;
    }
    else
    {
      r=sqrt(f1*f1+g1*g1);
      cs=f1/r;
      sn=g1/r;
    }
    if((cln::abs(f)>cln::abs(g)) && (cs<zero) )
    {
      cs=-cs;
      sn=-sn;
      r=-r;
    }
  }
  
  /* End of DLARTG */
  
  return;
}
  
template<typename T>
void dlaset(const std::string& uplo, int m, int n, T alpha, T beta, T* A, int lda)
{
  /*
*
*  Purpose
*  =======
*
*  DLASET initializes an m-by-n matrix A to BETA on the diagonal and
*  ALPHA on the offdiagonals.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies the part of the matrix A to be set.
*          = 'U':      Upper triangular part is set; the strictly lower
*                      triangular part of A is not changed.
*          = 'L':      Lower triangular part is set; the strictly upper
*                      triangular part of A is not changed.
*          Otherwise:  All of the matrix A is set.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  ALPHA   (input) DOUBLE PRECISION
*          The constant to which the offdiagonal elements are to be set.
*
*  BETA    (input) DOUBLE PRECISION
*          The constant to which the diagonal elements are to be set.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On exit, the leading m-by-n submatrix of A is set as follows:
*
*          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
*          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
*          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,
*
*          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
* =====================================================================
*
  */
  int i,j;
  
  /* .. Executable Statements .. */
  if(uplo=="U")
  {
    /* Set the strictly upper triangular or trapezoidal part of the */
    /* array to ALPHA. */
    for(j=2;j<=n;j++)
      for(i=1;i<=std::min(j-1,m);i++)
        A[i-1+lda*(j-1)]=alpha;
  }
  else if(uplo=="L")
  {
    /* Set the strictly lower triangular or trapezoidal part of the */
    /* array to ALPHA. */    
    for(j=1;j<=std::min(m,n);j++)
      for(i=j+1;i<=m;i++)
        A[i-1+lda*(j-1)]=alpha;
  }
  else
  {
    /* Set the leading m-by-n submatrix to ALPHA. */
    for(j=1;j<=n;j++)
      for(i=1;i<=m;i++)
        A[i-1+lda*(j-1)]=alpha;
  }
  
  /* Set the first min(M,N) diagonal elements to BETA. */
  for(i=1;i<=std::min(m,n);i++)
    A[i-1+lda*(i-1)]=beta;
  
  return;
  
  /* End of DLASET */
}

template<typename T>
void dlassq(int n, T* x, int incx, T& scale, T& sumsq)
{
  /*
*
*  Purpose
*  =======
*
*  DLASSQ  returns the values  scl  and  smsq  such that
*
*     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
*
*  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
*  assumed to be non-negative and  scl  returns the value
*
*     scl = max( scale, abs( x( i ) ) ).
*
*  scale and sumsq must be supplied in SCALE and SUMSQ and
*  scl and smsq are overwritten on SCALE and SUMSQ respectively.
*
*  The routine makes only one pass through the vector x.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of elements to be used from the vector X.
*
*  X       (input) DOUBLE PRECISION
*          The vector for which a scaled sum of squares is computed.
*             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
*
*  INCX    (input) INTEGER
*          The increment between successive values of the vector X.
*          INCX > 0.
*
*  SCALE   (input/output) DOUBLE PRECISION
*          On entry, the value  scale  in the equation above.
*          On exit, SCALE is overwritten with  scl , the scaling factor
*          for the sum of squares.
*
*  SUMSQ   (input/output) DOUBLE PRECISION
*          On entry, the value  sumsq  in the equation above.
*          On exit, SUMSQ is overwritten with  smsq , the basic sum of
*          squares from which  scl  has been factored out.
*
* =====================================================================
  */
  int ix;
  T absxi;
  
  const T zero=0.0;
  
  /* .. Executable Statements .. */
  
  if(n>0)
  {
    for(ix=1;ix<=(1+(n-1)*incx);ix+=incx)
    {
      if(x[ix-1]!=zero)
      {
        absxi=fabs(x[ix-1]);
        if(scale<absxi)
        {
          sumsq=1.0+sumsq*(scale/absxi)*(scale/absxi);
          scale=absxi;
        }
        else
          sumsq=sumsq+(absxi/scale)*(absxi/scale);
      }
    }
  }
  return;
  
  /* End of DLASSQ */
}

/* CLN version */

template<typename T>
void dlassq(int n, T* x, int incx, T& scale, T& sumsq, int digits)
{
  /*
*
*  Purpose
*  =======
*
*  DLASSQ  returns the values  scl  and  smsq  such that
*
*     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
*
*  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
*  assumed to be non-negative and  scl  returns the value
*
*     scl = max( scale, abs( x( i ) ) ).
*
*  scale and sumsq must be supplied in SCALE and SUMSQ and
*  scl and smsq are overwritten on SCALE and SUMSQ respectively.
*
*  The routine makes only one pass through the vector x.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of elements to be used from the vector X.
*
*  X       (input) DOUBLE PRECISION
*          The vector for which a scaled sum of squares is computed.
*             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
*
*  INCX    (input) INTEGER
*          The increment between successive values of the vector X.
*          INCX > 0.
*
*  SCALE   (input/output) DOUBLE PRECISION
*          On entry, the value  scale  in the equation above.
*          On exit, SCALE is overwritten with  scl , the scaling factor
*          for the sum of squares.
*
*  SUMSQ   (input/output) DOUBLE PRECISION
*          On entry, the value  sumsq  in the equation above.
*          On exit, SUMSQ is overwritten with  smsq , the basic sum of
*          squares from which  scl  has been factored out.
*
* =====================================================================
  */
  int ix;
  T absxi;
  
  cln::float_format_t precision=cln::float_format(digits);
  const T zero = cln::cl_float(0,precision);
  const T one = cln::cl_float(1,precision);
  
  /* .. Executable Statements .. */
  
  if(n>0)
  {
    for(ix=1;ix<=(1+(n-1)*incx);ix+=incx)
    {
      if(x[ix-1]!=zero)
      {
        absxi=cln::abs(x[ix-1]);
        if(scale<absxi)
        {
          sumsq=one+sumsq*(scale/absxi)*(scale/absxi);
          scale=absxi;
        }
        else
          sumsq=sumsq+(absxi/scale)*(absxi/scale);
      }
    }
  }
  return;
  
  /* End of DLASSQ */
}

template<typename T>
void axpy(int n, T da, T* dx, int incx, T* dy, int incy)
{
  /* original comment:
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
  */
  int i,ix,iy,m,mp1;
  
  //std::cout << "calling axpy \n";
  
  if(n<=0)
    return;
//  if(da==0.0)
//    return;
  if((incx==1)&&(incy==1))
    goto label20;
  
  /*  code for unequal increments or equal increments */
  /*   not equal to 1 */
  ix=1;
  iy=1;
  if(incx<0)
    ix=(-n+1)*incx+1;
  if(incy<0)
    iy=(-n+1)*incy+1;
  for(i=1;i<=n;i++)
  {
    dy[iy-1]=dy[iy-1]+da*dx[ix-1];
    ix+=incx;
    iy+=incy;
  }
  return;
  
/*
c        code for both increments equal to 1
c
c
c        clean-up loop
*/
label20:
  m=n%4;
  if(m==0) 
    goto label40;
  for(i=1;i<=m;i++)
    dy[i-1]=dy[i-1]+da*dx[i-1];
  if(n<4)
    return;
  
label40:
  mp1=m+1;
  for(i=mp1;i<=n;i=i+4)
  {
    dy[i-1]=dy[i-1]+da*dx[i-1];
    dy[i]=dy[i]+da*dx[i];
    dy[i+1]=dy[i+1]+da*dx[i+1];
    dy[i+2]=dy[i+2]+da*dx[i+2];
  }
  
  return;
}

template<typename T>
void dlascl(const std::string& type, int kl, int ku, T cfrom, T cto, int m, int n, T* A, int lda, int& info)
{
  /* Multiplies a matrix A by the constant cto/cfrom */
  T mul;
  int i, j;
  
  mul=cto/cfrom;
  
  for(j=1;j<=n;j++)
    for(i=1;i<=m;i++)
      A[i-1+lda*(j-1)]=mul*A[i-1+lda*(j-1)];
    
  return;
}

template<typename T>
void dlabad(T small, T large)
{
  /*
*
*  Purpose
*  =======
*
*  DLABAD takes as input the values computed by SLAMCH for underflow and
*  overflow, and returns the square root of each of these values if the
*  log of LARGE is sufficiently large.  This subroutine is intended to
*  identify machines with a large exponent range, such as the Crays, and
*  redefine the underflow and overflow limits to be the square roots of
*  the values computed by DLAMCH.  This subroutine is needed because
*  DLAMCH does not compensate for poor arithmetic in the upper half of
*  the exponent range, as is found on a Cray.
*
*  Arguments
*  =========
*
*  SMALL   (input/output) DOUBLE PRECISION
*          On entry, the underflow threshold as computed by DLAMCH.
*          On exit, if LOG10(LARGE) is sufficiently large, the square
*          root of SMALL, otherwise unchanged.
*
*  LARGE   (input/output) DOUBLE PRECISION
*          On entry, the overflow threshold as computed by DLAMCH.
*          On exit, if LOG10(LARGE) is sufficiently large, the square
*          root of LARGE, otherwise unchanged.
*
*  =====================================================================
  */
  
  /*
*     If it looks like we're on a Cray, take the square root of
*     SMALL and LARGE to avoid overflow and underflow problems.
  */
  if(log10(large)>2000.0)
  {
    small=sqrt(small);
    large=sqrt(large);
  }
  
  return;
  
  /* End of DLABAD */
}

template<typename T>
int idamax(int n, T* dx, int incx)
{
  /* original comment:
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified to correct problem with negative increment, 8/21/90. */
  
  T dmax;
  int i,ix;
  
  int result=0;
  if(n<1)
    return result;
  result=1;
  if(n==1)
    return result;
  if(incx==1)
    goto label20;
  
  /* code for increment not equal to 1 */
  
  ix=1;
  if(incx<0)
    ix=(-n+1)*incx+1;
  dmax=fabs(dx[ix-1]);
  ix+=incx;
  for(i=2;i<=n;i++)
  {
    if(fabs(dx[ix-1])<=dmax)
      goto label5;
    result=i;
    dmax=fabs(dx[ix-1]);
label5:
    ix+=incx;
  }
  return result;
  
  /* code for increment equal to 1 */
  
label20:
  dmax=fabs(dx[0]);
  for(i=2;i<=n;i++)
  {
    if(fabs(dx[i-1])<=dmax)
      continue;
    result=i;
    dmax=fabs(dx[i-1]);
label30:;
  }
  return result;
}

/* CLN version */

template<>
int idamax <cln::cl_F> (int n, cln::cl_F* dx, int incx)
{
  /* original comment:
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified to correct problem with negative increment, 8/21/90. */
  
  cln::cl_F dmax;
  int i,ix;
  
  int result=0;
  if(n<1)
    return result;
  result=1;
  if(n==1)
    return result;
  if(incx==1)
    goto label20;
  
  /* code for increment not equal to 1 */
  
  ix=1;
  if(incx<0)
    ix=(-n+1)*incx+1;
  dmax=cln::abs(dx[ix-1]);
  ix+=incx;
  for(i=2;i<=n;i++)
  {
    if(cln::abs(dx[ix-1])<=dmax)
      goto label5;
    result=i;
    dmax=cln::abs(dx[ix-1]);
label5:
    ix+=incx;
  }
  return result;
  
  /* code for increment equal to 1 */
  
label20:
  dmax=cln::abs(dx[0]);
  for(i=2;i<=n;i++)
  {
    if(cln::abs(dx[i-1])<=dmax)
      continue;
    result=i;
    dmax=cln::abs(dx[i-1]);
label30:;
  }
  return result;
}

template<typename T>
void drot(int n, T* dx,int incx, T* dy, int incy, T c, T s)
{
  /* original comment:
c
c     applies a plane rotation.
c     jack dongarra, linpack, 3/11/78.
c
  */
  T dtemp;
  int i,ix,iy;
  
  if(n<=0)
    return;
  if((incx==1)&&(incy==1))
    goto label20;
  
  /*
c
c       code for unequal increments or equal increments not equal
c         to 1
c
  */
  ix=1;
  iy=1;
  if(incx<0)
    ix=(-n+1)*incx+1;
  if(incy<0)
    iy=(-n+1)*incy+1;
  for(i=1;i<=n;i++)
  {
    dtemp=c*dx[ix-1]+s*dy[iy-1];
    dy[iy-1]=c*dy[iy-1]-s*dx[ix-1];
    dx[ix-1]=dtemp;
    ix+=incx;
    iy+=incy;
  }
  return;
  
label20:
  /*
c
c       code for both increments equal to 1
c
  */
  for(i=1;i<=n;i++)
  {
    dtemp=c*dx[i-1]+s*dy[i-1];
    dy[i-1]=c*dy[i-1]-s*dx[i-1];
    dx[i-1]=dtemp;
  }
  return;
}

template<typename T>
void dladiv(T a, T b, T c, T d, T& p, T& q)
{
  /*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
* *
*  Purpose
*  =======
*
*  DLADIV performs complex division in  real arithmetic
*
*                        a + i*b
*             p + i*q = ---------
*                        c + i*d
*
*  The algorithm is due to Robert L. Smith and can be found
*  in D. Knuth, The art of Computer Programming, Vol.2, p.195
*
*  Arguments
*  =========
*
*  A       (input) DOUBLE PRECISION
*  B       (input) DOUBLE PRECISION
*  C       (input) DOUBLE PRECISION
*  D       (input) DOUBLE PRECISION
*          The scalars a, b, c, and d in the above expression.
*
*  P       (output) DOUBLE PRECISION
*  Q       (output) DOUBLE PRECISION
*          The scalars p and q in the above expression.
*
*  =====================================================================
*/
  
  /* .. Local Scalars .. */
  T e, f;
  
  /* .. Executable Statements .. */
  
  if(fabs(d)<fabs(c))
  {
    e=d/c;
    f=c+d*e;
    p=(a+b*e)/f;
    q=(b-a*e)/f;
  }
  else
  {
    e=c/d;
    f=d+c*e;
    p=(b+a*e)/f;
    q=(-a+b*e)/f;
  }
  
  return;
  
  /* End of DLADIV */
}

/* CLN version */

template<>
void dladiv<cln::cl_F>(cln::cl_F a, cln::cl_F b, cln::cl_F c, cln::cl_F d, cln::cl_F& p, cln::cl_F& q)
{
  /*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
* *
*  Purpose
*  =======
*
*  DLADIV performs complex division in  real arithmetic
*
*                        a + i*b
*             p + i*q = ---------
*                        c + i*d
*
*  The algorithm is due to Robert L. Smith and can be found
*  in D. Knuth, The art of Computer Programming, Vol.2, p.195
*
*  Arguments
*  =========
*
*  A       (input) DOUBLE PRECISION
*  B       (input) DOUBLE PRECISION
*  C       (input) DOUBLE PRECISION
*  D       (input) DOUBLE PRECISION
*          The scalars a, b, c, and d in the above expression.
*
*  P       (output) DOUBLE PRECISION
*  Q       (output) DOUBLE PRECISION
*          The scalars p and q in the above expression.
*
*  =====================================================================
*/
  
  /* .. Local Scalars .. */
  cln::cl_F e, f;
  
  /* .. Executable Statements .. */
  
  if(cln::abs(d)<cln::abs(c))
  {
    e=d/c;
    f=c+d*e;
    p=(a+b*e)/f;
    q=(b-a*e)/f;
  }
  else
  {
    e=c/d;
    f=d+c*e;
    p=(b+a*e)/f;
    q=(-a+b*e)/f;
  }
  
  return;
  
  /* End of DLADIV */
}

template<typename T>
void dger(int m, int n, T alpha, T* x, int incx, T* y, int incy, T* A, int lda)
{
  /* Purpose
*  =======
*
*  DGER   performs the rank 1 operation
*
*     A := alpha*x*y' + A,
*
*  where alpha is a scalar, x is an m element vector, y is an n element
*  vector and A is an m by n matrix.
*
*  Parameters
*  ==========
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( m - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the m
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients. On exit, A is
*           overwritten by the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*/
  T zero=0.0;
  /* .. Local Scalars .. */
  T temp;
  int i,info,ix,j,jy,kx;
  
  /* .. Executable Statements ..
*
*     Test the input parameters. */
  info=0;
  if(m<0)
    info=1;
  else if(n<0)
    info=2;
  else if(incx==0)
    info=5;
  else if(incy==0)
    info=7;
  else if(lda<std::max(1,m))
    info=9;
  
  if(info!=0)
  {
    std::cout << "Error in dger \n";
    exit(1);
  }
  
  /* Start the operations. In this version the elements of A are */
  /* accessed sequentially with one pass through A. */
  
  if(incy>0)
    jy=1;
  else
    jy=1-(n-1)*incy;
  
  if(incx==1)
  {
    for(j=1;j<=n;j++)
    {
      if(y[jy-1]!=zero)
      {
        temp=alpha*y[jy-1];
        for(i=1;i<=m;i++)
          A[i-1+lda*(j-1)]=A[i-1+lda*(j-1)]+x[i-1]*temp;
      }
      jy+=incy;
    }
  }
  else
  {
    if(incx>0)
      kx=1;
    else
      kx=1-(m-1)*incx;
    for(j=1;j<=n;j++)
    {
      if(y[jy-1]!=zero)
      {
        temp=alpha*y[jy-1];
        ix=kx;
        for(i=1;i<=m;i++)
        {
          A[i-1+lda*(j-1)]=A[i-1+lda*(j-1)]+x[ix-1]*temp;
          ix+=incx;
        }
      }
      jy+=incy;
    }
  }
  
  return;
  
  /* End of DGER  . */
}

/* CLN version */

template<typename T>
void dger(int m, int n, T alpha, T* x, int incx, T* y, int incy, T* A, int lda, int digits)
{
  /* Purpose
*  =======
*
*  DGER   performs the rank 1 operation
*
*     A := alpha*x*y' + A,
*
*  where alpha is a scalar, x is an m element vector, y is an n element
*  vector and A is an m by n matrix.
*
*  Parameters
*  ==========
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( m - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the m
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients. On exit, A is
*           overwritten by the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*/
  cln::float_format_t precision=cln::float_format(digits);
  const T zero = cln::cl_float(0,precision);

  /* .. Local Scalars .. */
  T temp;
  int i,info,ix,j,jy,kx;
  
  /* .. Executable Statements ..
*
*     Test the input parameters. */
  info=0;
  if(m<0)
    info=1;
  else if(n<0)
    info=2;
  else if(incx==0)
    info=5;
  else if(incy==0)
    info=7;
  else if(lda<std::max(1,m))
    info=9;
  
  if(info!=0)
  {
    std::cout << "Error in dger \n";
    exit(1);
  }
  
  /* Start the operations. In this version the elements of A are */
  /* accessed sequentially with one pass through A. */
  
  if(incy>0)
    jy=1;
  else
    jy=1-(n-1)*incy;
  
  if(incx==1)
  {
    for(j=1;j<=n;j++)
    {
      if(y[jy-1]!=zero)
      {
        temp=alpha*y[jy-1];
        for(i=1;i<=m;i++)
          A[i-1+lda*(j-1)]=A[i-1+lda*(j-1)]+x[i-1]*temp;
      }
      jy+=incy;
    }
  }
  else
  {
    if(incx>0)
      kx=1;
    else
      kx=1-(m-1)*incx;
    for(j=1;j<=n;j++)
    {
      if(y[jy-1]!=zero)
      {
        temp=alpha*y[jy-1];
        ix=kx;
        for(i=1;i<=m;i++)
        {
          A[i-1+lda*(j-1)]=A[i-1+lda*(j-1)]+x[ix-1]*temp;
          ix+=incx;
        }
      }
      jy+=incy;
    }
  }
  
  return;
  
  /* End of DGER  . */
}

int mindigits(cln::cl_F* A, int lda, int digits)
{
  int i,cur,min;
  min=1000*digits;
  
  for(i=0;i<lda;i++)
  {
    cur=cln::float_digits(A[i]);
    if(cln::float_digits(A[i])<min)
      min=cur;
  }
  
  return min;
}

#endif

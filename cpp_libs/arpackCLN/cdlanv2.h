/* 
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994 
*
*
*  Purpose
*  =======
*
*  DLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric
*  matrix in standard form:
*
*       [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]
*       [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ]
*
*  where either
*  1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or
*  2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex
*  conjugate eigenvalues.
*
*  Arguments
*  =========
*
*  A       (input/output) DOUBLE PRECISION
*  B       (input/output) DOUBLE PRECISION
*  C       (input/output) DOUBLE PRECISION
*  D       (input/output) DOUBLE PRECISION
*          On entry, the elements of the input matrix.
*          On exit, they are overwritten by the elements of the
*          standardised Schur form.
*
*  RT1R    (output) DOUBLE PRECISION
*  RT1I    (output) DOUBLE PRECISION
*  RT2R    (output) DOUBLE PRECISION
*  RT2I    (output) DOUBLE PRECISION
*          The real and imaginary parts of the eigenvalues. If the
*          eigenvalues are both real, abs(RT1R) >= abs(RT2R); if the
*          eigenvalues are a complex conjugate pair, RT1I > 0.
*
*  CS      (output) DOUBLE PRECISION
*  SN      (output) DOUBLE PRECISION
*          Parameters of the rotation matrix.
*
*  =====================================================================
*
*/

#include "fortranfuncs.h"
#include <math.h>
#include <cln/real.h>
#include <cln/output.h>
#include <cln/real_io.h>

/* CLN version */

template<typename T>
void dlanv2(T& a, T& b, T& c, T& d, T& rt1r, T& rt1i, T& rt2r, T& rt2i, T& cs, T& sn, int digits)
{
  
  cln::float_format_t precision = cln::float_format(digits);
    /* Parameters */
  const cln::cl_F zero=cln::cl_float(0,precision);
  const cln::cl_F one=cln::cl_float(1,precision);
  const cln::cl_F half=cln::cl_float(0.5,precision);
  const cln::cl_F multpl=cln::cl_float(0.4,precision);
  
  T aa,bb,bcmax,bcmis,cc,cs1,dd,eps,p,sab,sac,scale;
  T sigma,sn1,tau,temp,z;
  
  
  /* .. Executable Statements ..
*
*     Initialize CS and SN */
  
  eps=dlamch<T>("P",digits); // precision
  
  if(c==zero)
  {
    cs=one;
    sn=zero;
    goto label10;
  }
  else if(b==zero)
  {
    /* Swap rows and columns */
    cs=zero;
    sn=one;
    temp=d;
    d=a;
    a=temp;
    b=-c;
    c=zero;
    goto label10;
  }
  else if( ((a-d)==zero)&&(sign<T>(one,b)!=sign<T>(one,c)) )
  {
    cs=one;
    sn=zero;
    goto label10;
  }
  else
  {
    /* Make diagonal elements equal */
    
    temp=a-d;
    p=half*temp;
    bcmax=cln::max(cln::abs(b),cln::abs(c));
    bcmis=cln::min(cln::abs(b),cln::abs(c))*sign(one,b)*sign(one,c);
    scale=cln::max(cln::abs(p),bcmax);
    z=(p/scale)*p+(bcmax/scale)*bcmis;
    
    /*
*        If Z is of the order of the machine accuracy, postpone the
*        decision on the nature of eigenvalues
    */
    if(z>=multpl*eps)
    {
      /* Real eigenvalues. Compute A and D. */
      
      z=p+sign(sqrt(scale)*sqrt(z),p);
      a=d+z;
      d=d-(bcmax/z)*bcmis;
      
      /* Compute B and the rotation matrix */
      
      tau=dlapy2<T>(c,z);
      cs=z/tau;
      sn=c/tau;
      b=b-c;
      c=zero;
    }
    else
    {
    /* Complex eigenvalues, or real (almost) equal eigenvalues.
       Make diagonal elements equal. */
      
    sigma=b+c;
    tau=dlapy2<T>(sigma,temp);
    cs=sqrt( half*(one+cln::abs(sigma)/tau) );
    sn=-(p/(tau*cs))*sign(one,sigma);
    
    /*
*
*        Compute [ AA  BB ] = [ A  B ] [ CS1 -SN1 ]
*                [ CC  DD ]   [ C  D ] [ SN1  CS1 ]
*
    */
    aa=a*cs+b*sn;
    bb=-a*sn+b*cs;
    cc=c*cs+d*sn;
    dd=-c*sn+d*cs;
    /*
*
*        Compute [ A  B ] = [ CS1  SN1 ] [ AA  BB ]
*                [ C  D ]   [-SN1  CS1 ] [ CC  DD ]
*
    */
    a=aa*cs+cc*sn;
    b=bb*cs+dd*sn;
    c=-aa*sn+cc*cs;
    d=-bb*sn+dd*cs;
    /*
*
*        Accumulate transformation
*
    */
//    temp=cs*cs1-sn*sn1;
//    sn=cs*sn1+sn*cs1;
 //   cs=temp;
    
    temp=half*(a+d);
    a=temp;
    d=temp;
    
    if(c!=zero)
    {
      if(b!=zero)
      {
        if(sign(one,b)==sign(one,c))
        {
          /* Real eigenvalues: reduce to upper triangular form */
          sab=sqrt(cln::abs(b));
          sac=sqrt(cln::abs(c));
          p=sign(sab*sac,c);
    //      std::cout << "sac*sab = " << sab*sac << ", c = " << c << ", p=" << p << "\n";
          tau=one/sqrt(cln::abs(b+c));
          a=temp+p;
          d=temp-p;
          b=b-c;
          c=zero;
          cs1=sab*tau;
          sn1=sac*tau;
          temp=cs*cs1-sn*sn1;
          sn=cs*sn1+sn*cs1;
          cs=temp;
        }
      }
      else
      {
        b=-c;
        c=zero;
        temp=cs;
        cs=-sn;
        sn=temp;
      }
    }
  }
  
  }
label10:

  /* Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I). */
          
  rt1r=a;
  rt2r=d;
  if(c==zero)
  {
    rt1i=zero;
    rt2i=zero;
  }
  else
  {
    rt1i=sqrt(cln::abs(b))*sqrt(cln::abs(c));
    rt2i=-rt1i;
  }
  return;
  
  /* End of CLN DLANV2 */
}

template<typename T>
void dlanv2(T& a, T& b, T& c, T& d, T& rt1r, T& rt1i, T& rt2r, T& rt2i, T& cs, T& sn)
{
  T zero=0.0;
  T one=1.0;
  T half=0.5;
  T multpl=4.0;
  
  T aa,bb,bcmax,bcmis,cc,cs1,dd,eps,p,sab,sac,scale;
  T sigma,sn1,tau,temp,z;
  
  
  /* .. Executable Statements ..
*
*     Initialize CS and SN */
//  std::cout << "(a,b,c,d)=(" << a << "," << b << "," << c << "," << d << ")\n";
  
  eps=dlamch<T>("P"); // precision
  
  if(c==zero)
  {
    cs=one;
    sn=zero;
    goto label10;
  }
  else if(b==zero)
  {
    /* Swap rows and columns */
    cs=zero;
    sn=one;
    temp=d;
    d=a;
    a=temp;
    b=-c;
    c=zero;
    goto label10;
  }
  else if( ((a-d)==zero)&&(sign<T>(one,b)!=sign<T>(one,c)) )
  {
    cs=one;
    sn=zero;
    goto label10;
  }
  else
  {
    /* Make diagonal elements equal */
    
    temp=a-d;
    p=half*temp;
    bcmax=std::max(T(fabs(b)),T(fabs(c)));
    bcmis=std::min(T(fabs(b)),T(fabs(c)))*sign<T>(one,b)*sign<T>(one,c);
    scale=std::max(T(fabs(p)),bcmax);
    z=(p/scale)*p+(bcmax/scale)*bcmis;
    
    /*
*        If Z is of the order of the machine accuracy, postpone the
*        decision on the nature of eigenvalues
    */
    
//    std::cout << "z = " << z << "\n";
    if(z>=multpl*eps)
    {
      /* Real eigenvalues. Compute A and D. */
//      std::cout << "sqrt(scale)*sqrt(z)=" << sqrt(scale)*sqrt(z) << "\n";
//      std::cout << "p=" << p << "\n";
//      std::cout << "bcmax=" << bcmax << "\n";
//      std::cout << "bcmis=" << bcmis << "\n";
//      std::cout << "scale=" << scale << "\n";
//      std::cout << "sign<T>(sqrt(scale)*sqrt(z),p)=" << sign<T>(sqrt(scale)*sqrt(z),p) << "\n";
      z=p+sign<T>(sqrt(scale)*sqrt(z),p);
//      std::cout << "z before add= " << z << "\n";
      a=d+z;
      d=d-(bcmax/z)*bcmis;
      
      /* Compute B and the rotation matrix */
      
      tau=dlapy2<T>(c,z);
      cs=z/tau;
      sn=c/tau;
      b=b-c;
      c=zero;
//      std::cout << "(a,b,c,d)=(" << a << "," << b << "," << c << "," << d << ")\n";
    }
    else
    {
    /* Complex eigenvalues, or real (almost) equal eigenvalues.
       Make diagonal elements equal. */
      
      sigma=b+c;
      tau=dlapy2<T>(sigma,temp);
      cs=sqrt( half*(one+fabs(sigma)/tau) );
      sn=-(p/(tau*cs))*sign(one,sigma);
    
    /*
*
*        Compute [ AA  BB ] = [ A  B ] [ CS1 -SN1 ]
*                [ CC  DD ]   [ C  D ] [ SN1  CS1 ]
*
    */
      aa=a*cs+b*sn;
      bb=-a*sn+b*cs;
      cc=c*cs+d*sn;
      dd=-c*sn+d*cs;
    /*
*
*        Compute [ A  B ] = [ CS1  SN1 ] [ AA  BB ]
*                [ C  D ]   [-SN1  CS1 ] [ CC  DD ]
*
    */
      a=aa*cs+cc*sn;
      b=bb*cs+dd*sn;
      c=-aa*sn+cc*cs;
      d=-bb*sn+dd*cs;
    /*
*
*        Accumulate transformation
*
    */
//    temp=cs*cs1-sn*sn1;
//    sn=cs*sn1+sn*cs1;
 //   cs=temp;
    
      temp=half*(a+d);
      a=temp;
      d=temp;
    
      if(c!=zero)
      {
        if(b!=zero)
        {
          if(sign(one,b)==sign(one,c))
          {
            /* Real eigenvalues: reduce to upper triangular form */
            sab=sqrt(fabs(b));
            sac=sqrt(fabs(c));
            p=sign(sab*sac,c);
      //      std::cout << "sac*sab = " << sab*sac << ", c = " << c << ", p=" << p << "\n";
            tau=one/sqrt(fabs(b+c));
            a=temp+p;
            d=temp-p;
            b=b-c;
            c=zero;
            cs1=sab*tau;
            sn1=sac*tau;
            temp=cs*cs1-sn*sn1;
            sn=cs*sn1+sn*cs1;
            cs=temp;
          }
        }
        else
        {
          b=-c;
          c=zero;
          temp=cs;
          cs=-sn;
          sn=temp;
        }
      }
    }
  
  }
label10:

  /* Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I). */
          
  rt1r=a;
  rt2r=d;
  if(c==zero)
  {
    rt1i=zero;
    rt2i=zero;
  }
  else
  {
    rt1i=sqrt(fabs(b))*sqrt(fabs(c));
    rt2i=-rt1i;
  }
  return;
  
  /* End of DLANV2 */
}
#ifndef C_DGEMV
#define C_DGEMV

/*
*  Purpose
*  =======
*
*  DGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
*
*           Unchanged on exit.
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
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
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
*     Converted to C++ by Chris Scullard, LLNL, July 2014
*/

#include "cdebug.h"
#include "cstat.h"
#include "fortranfuncs.h"
#include <cln/real.h>
#include <cln/output.h>
#include <cln/real_io.h>

/* CLN version */

template<typename T>
void dgemv(const std::string& trans, int m, int n, T alpha, T* A, int lda, T* x, int incx, T beta, T* y, int incy, int digits)
{
  cln::float_format_t precision=cln::float_format(digits);
  const T one = cln::cl_float(1,precision);
  const T zero = cln::cl_float(0,precision);
  
  T temp;
  int i,info,ix,iy,j,jx,jy,kx,ky,lenx,leny;
  
  /*
*     .. Executable Statements ..
*
*     Test the input parameters.
*
  */
  //std::cout << "trans = " << trans << "\n";
  info=0;
  if( (trans!="N")&&(trans!="T")&&(trans!="C") )
    info = 1;
  else if(m<0)
    info=2;
  else if(n<0)
    info=3;
  else if(lda<std::max(1,m) )
    info=6;
  else if(incx==0)
    info=8;
  else if(incy==0)
    info=11;
  
  if(info!=0)
  {
    std::cout << "Error in dgemv, info = " << info << "\n";
    exit(1);
    return;
  }
  
  /*
*
*     Quick return if possible.
*
  */
  
  if((m==0)||(n==0)||( (alpha==zero)&&(beta==one) ))
    return;
  
  /*
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
  */
  
  if(trans=="N")
  {
    lenx=n;
    leny=m;
  }
  else
  {
    lenx=m;
    leny=n;
  }
  
  if(incx>0)
    kx=1;
  else
    kx=1-(lenx-1)*incx;
  
  if(incy>0)
    ky=1;
  else
    ky=1-(leny-1)*incy;
  
  /*
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
  */
  //dvout<T>(debug.logfil,n,y,debug.ndigit,"_dgemv: y before y=beta*y");
  
  if(beta!=one)
  {
    if(incy==1)
    {
      if(beta==zero)
      {
        for(i=1;i<=leny;i++)
          y[i-1]=zero;
      }
      else
      {
        for(i=1;i<=leny;i++)
          y[i-1]=beta*y[i-1];
      }
    }
    else
    {
      iy=ky;
      if(beta==zero)
      {
        for(i=1;i<=leny;i++)
        {
          y[iy-1]=zero;
          iy+=incy;
        }
      }
      else
      {
        for(i=1;i<=leny;i++)
        {
          y[iy-1]=beta*y[iy-1];
          iy+=incy;
        }
      }
    }
  }
  if(alpha==zero)
    return;
  
  //dvout<T>(debug.logfil,n,y,debug.ndigit,"_dgemv: y after y=beta*y");
  
  if(trans=="N")
  {
    /*
*
*        Form  y := alpha*A*x + y.
*
    */
    jx=kx;
    if(incy==1)
    {
      for(j=1;j<=n;j++)
      {
        if(x[jx-1]!=zero)
    //    if(cln::abs(x[jx-1])>cln::float_epsilon(precision))
        {
          temp=alpha*x[jx-1];
     //     std::cout << alpha << "\n";
          for(i=1;i<=m;i++)
          {
    //        std::cout << "y[i-1]=" << y[i-1] << ", temp = " << temp << ", A=" << A[i-1+lda*(j-1)] << "\n";
     //       std::cout << "temp*A=" << temp*A[i-1+lda*(j-1)] << "\n";
            y[i-1]+=temp*A[i-1+lda*(j-1)];
          }
        }
        jx+=incx;
      }
    }
    else
    {
      for(j=1;j<=n;j++)
      {
        if(x[jx-1]!=zero)
        {
          temp=alpha*x[jx-1];
          iy=ky;
          for(i=1;i<=m;i++)
          {
            y[iy-1]+=temp*A[i-1+lda*(j-1)];
            iy+=incy;
          }
        }
        jx+=incx;
      }
    }
  }
  else
  {
    /*
*
*        Form  y := alpha*A'*x + y.
*
    */
    jy=ky;
    if(incx==1)
    {
      
      for(j=1;j<=n;j++)
      {
    //    ivout<T>(debug.logfil,1,&j,debug.ndigit,"_dgemv: j");
    //    dvout<T>(debug.logfil,leny,y,debug.ndigit,"_dgemv: y at this j");
    //    dvout<T>(debug.logfil,lenx,x,debug.ndigit,"_dgemv: x at this j");
    //    dmout<T>(debug.logfil,n,n,A,lda,debug.ndigit,"_dgemv: A at this j");
        temp=zero;
        for(i=1;i<=m;i++)
          temp+=A[i-1+lda*(j-1)]*x[i-1];
        y[jy-1]+=alpha*temp;
        jy+=incy;
      }
    }
    else
    {
      for(j=1;j<=n;j++)
      {
        temp=zero;
        ix=kx;
        for(i=1;i<=m;i++)
        {
          temp+=A[i-1+lda*(j-1)]*x[ix-1];
          ix+=incx;
        }
        y[jy-1]+=alpha*temp;
        jy+=incy;
      }
    }
  }
  
  return;
  /*
*
*     End of DGEMV .
*
  */
}

template<typename T>
void dgemv(const std::string& trans, int m, int n, T alpha, T* A, int lda, T* x, int incx, T beta, T* y, int incy)
{
  T one=1.0;
  T zero=0.0;
  
  T temp;
  int i,info,ix,iy,j,jx,jy,kx,ky,lenx,leny;
  
 // std::cout << "calling wrong dgemv \n";
 // exit(EXIT_FAILURE);
  
  /*
*     .. Executable Statements ..
*
*     Test the input parameters.
*
  */
//  std::cout << "beta = " << beta << "\n";
//  std::cout << "alpha = " << alpha << "\n";
  info=0;
  if( (trans!="N")&&(trans!="T")&&(trans!="C") )
    info = 1;
  else if(m<0)
    info=2;
  else if(n<0)
    info=3;
  else if(lda<std::max(1,m) )
    info=6;
  else if(incx==0)
    info=8;
  else if(incy==0)
    info=11;
  
  if(info!=0)
  {
    std::cout << "Error in dgemv, info = " << info << "\n";
    exit(1);
    return;
  }
  
  /*
*
*     Quick return if possible.
*
  */
  
  if((m==0)||(n==0)||( (alpha==zero)&&(beta==one) ))
    return;
  
  /*
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
  */
  
  if(trans=="N")
  {
    lenx=n;
    leny=m;
  }
  else
  {
    lenx=m;
    leny=n;
  }
  
  if(incx>0)
    kx=1;
  else
    kx=1-(lenx-1)*incx;
  
  if(incy>0)
    ky=1;
  else
    ky=1-(leny-1)*incy;
  
  /*
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
  */
  //dvout<T>(debug.logfil,n,y,debug.ndigit,"_dgemv: y before y=beta*y");
  
  if(incy==1)
    for(i=1;i<=leny;i++)
      y[i-1]=beta*y[i-1];
  else
  {
    iy=ky;
    for(i=1;i<=leny;i++)
    {
      y[iy-1]=beta*y[iy-1];
      iy+=incy;
    }
  }
  
  //dvout<T>(debug.logfil,n,y,debug.ndigit,"_dgemv: y after y=beta*y");
  
  if(trans=="N")
  {
    /*
*
*        Form  y := alpha*A*x + y.
*
    */
    jx=kx;
    if(incy==1)
    {
      for(j=1;j<=n;j++)
      {
        if(x[jx-1]!=zero)
        {
          temp=alpha*x[jx-1];
          for(i=1;i<=m;i++)
            y[i-1]+=temp*A[i-1+lda*(j-1)];
        }
        jx+=incx;
      }
    }
    else
    {
      for(j=1;j<=n;j++)
      {
        if(x[jx-1]!=zero)
        {
          temp=alpha*x[jx-1];
          iy=ky;
          for(i=1;i<=m;i++)
          {
            y[iy-1]+=temp*A[i-1+lda*(j-1)];
            iy+=incy;
          }
        }
        jx+=incx;
      }
    }
  }
  else
  {
    /*
*
*        Form  y := alpha*A'*x + y.
*
    */
    jy=ky;
    if(incx==1)
    {
      
      for(j=1;j<=n;j++)
      {
    //    ivout<T>(debug.logfil,1,&j,debug.ndigit,"_dgemv: j");
    //    dvout<T>(debug.logfil,leny,y,debug.ndigit,"_dgemv: y at this j");
    //    dvout<T>(debug.logfil,lenx,x,debug.ndigit,"_dgemv: x at this j");
    //    dmout<T>(debug.logfil,n,n,A,lda,debug.ndigit,"_dgemv: A at this j");
        temp=zero;
        for(i=1;i<=m;i++)
          temp+=A[i-1+lda*(j-1)]*x[i-1];
        y[jy-1]+=alpha*temp;
        jy+=incy;
      }
    }
    else
    {
      for(j=1;j<=n;j++)
      {
        temp=zero;
        ix=kx;
        for(i=1;i<=m;i++)
        {
          temp+=A[i-1+lda*(j-1)]*x[ix-1];
          ix+=incx;
        }
        y[jy-1]+=alpha*temp;
        jy+=incy;
      }
    }
  }
  
  return;
  /*
*
*     End of DGEMV .
*
  */
}

#endif
    
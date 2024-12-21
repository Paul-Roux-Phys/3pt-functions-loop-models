/*
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dneigh
c
c\Description:
c  Compute the eigenvalues of the current upper Hessenberg matrix
c  and the corresponding Ritz estimates given the current residual norm.
c
c\Usage:
c  call dneigh
c     ( RNORM, N, H, LDH, RITZR, RITZI, BOUNDS, Q, LDQ, WORKL, IERR )
c
c\Arguments
c  RNORM   Double precision scalar.  (INPUT)
c          Residual norm corresponding to the current upper Hessenberg 
c          matrix H.
c
c  N       Integer.  (INPUT)
c          Size of the matrix H.
c
c  H       Double precision N by N array.  (INPUT)
c          H contains the current upper Hessenberg matrix.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.
c
c  RITZR,  Double precision arrays of length N.  (OUTPUT)
c  RITZI   On output, RITZR(1:N) (resp. RITZI(1:N)) contains the real 
c          (respectively imaginary) parts of the eigenvalues of H.
c
c  BOUNDS  Double precision array of length N.  (OUTPUT)
c          On output, BOUNDS contains the Ritz estimates associated with
c          the eigenvalues RITZR and RITZI.  This is equal to RNORM 
c          times the last components of the eigenvectors corresponding 
c          to the eigenvalues in RITZR and RITZI.
c
c  Q       Double precision N by N array.  (WORKSPACE)
c          Workspace needed to store the eigenvectors of H.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  WORKL   Double precision work array of length N**2 + 3*N.  (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  This is needed to keep the full Schur form
c          of H and also in the calculation of the eigenvectors of H.
c
c  IERR    Integer.  (OUTPUT)
c          Error exit flag from dlaqrb or dtrevc.
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
c     dlaqrb  ARPACK routine to compute the real Schur form of an
c             upper Hessenberg matrix and last row of the Schur vectors.
c     second  ARPACK utility routine for timing.
c     dmout   ARPACK utility routine that prints matrices
c     dvout   ARPACK utility routine that prints vectors.
c     dlacpy  LAPACK matrix copy routine.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dtrevc  LAPACK routine to compute the eigenvectors of a matrix
c             in upper quasi-triangular form
c     dgemv   Level 2 BLAS routine for matrix vector multiplication.
c     dcopy   Level 1 BLAS that copies one vector to another .
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c     dscal   Level 1 BLAS that scales a vector.
c     
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
c FILE: neigh.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
c
c\Remarks
c     None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
*/
#include "cdebug.h"
#include "cstat.h"
#include "fortranfuncs.h"
#include "cdgemv.h"
#include "cdtrevc.h"
#include "cdlaqrb.h"
#include <cln/real.h>
#include <cln/output.h>
#include <cln/real_io.h>

/* CLN version */

template<typename T>
void dneigh(T& rnorm, int& n, T* h, int& ldh, T* ritzr, T* ritzi, T* bounds, T* q, int ldq, T* workl, int& ierr, int digits)
{
  /*
c
c     %------------%
c     | Parameters |
c     %------------%
c
  */
  
  cln::float_format_t precision=cln::float_format(digits);
  const T one = cln::cl_float(1,precision);
  const T zero = cln::cl_float(0,precision);
  
  /*
c 
c     %------------------------%
c     | Local Scalars & Arrays |
c     %------------------------%
c
  */
  
  bool select[1];
  int i,iconj,msglvl;
  T temp;
  T vl[1];
  
  /*
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%
c
  */
  
  //dmout<T>(debug.logfil,n,n,q,ldq,debug.ndigit,"_neigh: q upon starting dneigh");
  //dvout<T>(debug.logfil,n*n+3*n,workl,debug.ndigit,"_dneigh: workl upon starting dneigh");
  second(timing.t0);
  msglvl=debug.mneigh;
  
  if(msglvl>2)
  {
    dmout<T>(debug.logfil,n,n,h,ldh,debug.ndigit,"_neigh: Entering upper Hessenberg matrix H");
  }
  
  /*
c 
c     %-----------------------------------------------------------%
c     | 1. Compute the eigenvalues, the last components of the    |
c     |    corresponding Schur vectors and the full Schur form T  |
c     |    of the current upper Hessenberg matrix H.              |
c     | dlaqrb returns the full Schur form of H in WORKL(1:N**2)  |
c     | and the last components of the Schur vectors in BOUNDS.   |
c     %-----------------------------------------------------------%
c
  */
  
  dlacpy<T>("All",n,n,h,ldh,workl,n);
//  dmout<T>(debug.logfil,n,n,workl,ldh,debug.ndigit,"_neigh: matrix copied into workl");
//  dvout<T>(debug.logfil,n,bounds,debug.ndigit,"_dneigh: bounds before call to dlaqrb");
  dlaqrb<T>(true,n,1,n,workl,n,ritzr,ritzi,bounds,ierr,digits);
//  dvout<T>(debug.logfil,n,bounds,debug.ndigit,"_dneigh: bounds after call to dlaqrb");
//  dmout<T>(debug.logfil,n,n,workl,ldh,debug.ndigit,"_neigh: workl after call to dlaqrb");
  if(ierr!=0) return;
  
  if(msglvl>1)
    dvout<T>(debug.logfil,n,bounds,debug.ndigit,"_neigh: last row of the Schur matrix for H");
  
  /*
c
c     %-----------------------------------------------------------%
c     | 2. Compute the eigenvectors of the full Schur form T and  |
c     |    apply the last components of the Schur vectors to get  |
c     |    the last components of the corresponding eigenvectors. |
c     | Remember that if the i-th and (i+1)-st eigenvalues are    |
c     | complex conjugate pairs, then the real & imaginary part   |
c     | of the eigenvector components are split across adjacent   |
c     | columns of Q.                                             |
c     %-----------------------------------------------------------%
c
  */
//  dvout<T>(debug.logfil,n*n+3*n,workl,debug.ndigit,"_dneigh: workl before call to dtrevc");
//  dmout<T>(debug.logfil,n,n,q,ldq,debug.ndigit,"_neigh: q before calling dtrevc");
  dtrevc<T>("R","A",select,n,workl,n,vl,n,q,ldq,n,n,&workl[n*n],ierr,digits);
//  dmout<T>(debug.logfil,n,n,q,ldq,debug.ndigit,"_neigh: q after calling dtrevc");
//  dvout<T>(debug.logfil,n*n+3*n,workl,debug.ndigit,"_dneigh: workl after call to dtrevc");
  
  if(ierr!=0) return;
  
  /*
c
c     %------------------------------------------------%
c     | Scale the returning eigenvectors so that their |
c     | euclidean norms are all one. LAPACK subroutine |
c     | dtrevc returns each eigenvector normalized so  |
c     | that the element of largest magnitude has      |
c     | magnitude 1; here the magnitude of a complex   |
c     | number (x,y) is taken to be |x| + |y|.         |
c     %------------------------------------------------%
c
  */
  
  iconj=0;
  for(i=1;i<=n;i++)
    if(cln::abs(ritzi[i-1])<=zero)
    {
      /*
c
c           %----------------------%
c           | Real eigenvalue case |
c           %----------------------%
c 
      */
      temp=dnrm2<T>(n,&q[0+ldq*(i-1)],1);
      dscal<T>(n,one/temp,&q[0+ldq*(i-1)],1);
    }
    else
    {
      /*
c
c           %-------------------------------------------%
c           | Complex conjugate pair case. Note that    |
c           | since the real and imaginary part of      |
c           | the eigenvector are stored in consecutive |
c           | columns, we further normalize by the      |
c           | square root of two.                       |
c           %-------------------------------------------%
c
      */
      
      if(iconj==0)
      {
        temp=dlapy2<T>(dnrm2<T>(n,&q[0+ldq*(i-1)],1),dnrm2<T>(n,&q[0+ldq*i],1));
        dscal<T>(n,one/temp,&q[0+ldq*(i-1)],1);
        dscal<T>(n,one/temp,&q[0+ldq*i],1);
        iconj=1;
      }
      else
        iconj=0;
    }
    
    //dvout<T>(debug.logfil,n,bounds,debug.ndigit,"_dneigh: bounds before call to dgemv");
   // dvout<T>(debug.logfil,n*n+3*n,workl,debug.ndigit,"_dneigh: workl before call to dgemv");
   // dmout<T>(debug.logfil,n,n,q,ldq,debug.ndigit,"_neigh: q matrix before calling dgemv");
    dgemv<T>("T",n,n,one,q,ldq,bounds,1,zero,workl,1,digits);
   // dvout<T>(debug.logfil,n*n+3*n,workl,debug.ndigit,"_dneigh: workl after call to dgemv");
    //dvout<T>(debug.logfil,n,bounds,debug.ndigit,"_dneigh: bounds after call to dgemv");
    
    if(msglvl>1)
      dvout<T>(debug.logfil,n,workl,debug.ndigit,"_neigh: Last row of the eigenvector matrix for H");
    
    /*
c
c     %----------------------------%
c     | Compute the Ritz estimates |
c     %----------------------------%
c
    */
    
    iconj=0;
    for(i=1;i<=n;i++)
    {
     // std::cout << "i=" << i << "\n";
    //  dvout<T>(debug.logfil,n,bounds,debug.ndigit,"_dneigh: bounds at this i");
      if(cln::abs(ritzi[i-1])<=zero)
      {
        /*
c
c           %----------------------%
c           | Real eigenvalue case |
c           %----------------------%
c 
        */
   //     std::cout << "workl[" << i-1 << "]=" << workl[i-1] << "\n";
        bounds[i-1]=rnorm*cln::abs(workl[i-1]);
      }
      else
      {
        /*
c
c           %-------------------------------------------%
c           | Complex conjugate pair case. Note that    |
c           | since the real and imaginary part of      |
c           | the eigenvector are stored in consecutive |
c           | columns, we need to take the magnitude    |
c           | of the last components of the two vectors |
c           %-------------------------------------------%
c
        */
        if(iconj==0)
        {
          bounds[i-1]=rnorm*dlapy2<T>(workl[i-1],workl[i]);
          bounds[i]=bounds[i-1];
          iconj=1;
        }
        else
          iconj=0;
      }
    }
      
      if(msglvl>2)
      {
        dvout<T>(debug.logfil,n,ritzr,debug.ndigit,"_neigh: Real part of the eigenvalues of H");
        dvout<T>(debug.logfil,n,ritzi,debug.ndigit,"_neigh: Imaginary part of the eigenvalues of H");
        dvout<T>(debug.logfil,n,bounds,debug.ndigit,"_neigh: Ritz estimates for the eigenvalues of H");
      }
      
      second(timing.t1);
      timing.tneigh+=timing.t1-timing.t0;
      
      return;
      /*
c
c     %---------------%
c     | End of dneigh |
c     %---------------%
c
      */
}

template<typename T>
void dneigh(T& rnorm, int& n, T* h, int& ldh, T* ritzr, T* ritzi, T* bounds, T* q, int ldq, T* workl, int& ierr)
{
  /*
c
c     %------------%
c     | Parameters |
c     %------------%
c
  */
  
  const T one=1.0;
  const T zero=0.0;
  
//  std::cout << "Called wrong dneigh \n";
//  exit(1);
  
  /*
c 
c     %------------------------%
c     | Local Scalars & Arrays |
c     %------------------------%
c
  */
  
  bool select[1];
  int i,iconj,msglvl;
  T temp;
  T vl[1];
  
  /*
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%
c
  */
  
  //dmout<T>(debug.logfil,n,n,q,ldq,debug.ndigit,"_neigh: q upon starting dneigh");
  //dvout<T>(debug.logfil,n*n+3*n,workl,debug.ndigit,"_dneigh: workl upon starting dneigh");
  second(timing.t0);
  msglvl=debug.mneigh;
  
  if(msglvl>2)
  {
    dmout<T>(debug.logfil,n,n,h,ldh,debug.ndigit,"_neigh: Entering upper Hessenberg matrix H");
  }
  
  /*
c 
c     %-----------------------------------------------------------%
c     | 1. Compute the eigenvalues, the last components of the    |
c     |    corresponding Schur vectors and the full Schur form T  |
c     |    of the current upper Hessenberg matrix H.              |
c     | dlaqrb returns the full Schur form of H in WORKL(1:N**2)  |
c     | and the last components of the Schur vectors in BOUNDS.   |
c     %-----------------------------------------------------------%
c
  */
  
  dlacpy<T>("All",n,n,h,ldh,workl,n);
//  dmout<T>(debug.logfil,n,n,workl,ldh,debug.ndigit,"_neigh: matrix copied into workl");
//  dvout<T>(debug.logfil,n,bounds,debug.ndigit,"_dneigh: bounds before call to dlaqrb");
  dlaqrb<T>(true,n,1,n,workl,n,ritzr,ritzi,bounds,ierr);
//  dvout<T>(debug.logfil,n,bounds,debug.ndigit,"_dneigh: bounds after call to dlaqrb");
//  dmout<T>(debug.logfil,n,n,workl,ldh,debug.ndigit,"_neigh: workl after call to dlaqrb");
  if(ierr!=0) return;
  
  if(msglvl>1)
    dvout<T>(debug.logfil,n,bounds,debug.ndigit,"_neigh: last row of the Schur matrix for H");
  
  /*
c
c     %-----------------------------------------------------------%
c     | 2. Compute the eigenvectors of the full Schur form T and  |
c     |    apply the last components of the Schur vectors to get  |
c     |    the last components of the corresponding eigenvectors. |
c     | Remember that if the i-th and (i+1)-st eigenvalues are    |
c     | complex conjugate pairs, then the real & imaginary part   |
c     | of the eigenvector components are split across adjacent   |
c     | columns of Q.                                             |
c     %-----------------------------------------------------------%
c
  */
//  dvout<T>(debug.logfil,n*n+3*n,workl,debug.ndigit,"_dneigh: workl before call to dtrevc");
//  dmout<T>(debug.logfil,n,n,q,ldq,debug.ndigit,"_neigh: q before calling dtrevc");
  dtrevc<T>("R","A",select,n,workl,n,vl,n,q,ldq,n,n,&workl[n*n],ierr);
//  dmout<T>(debug.logfil,n,n,q,ldq,debug.ndigit,"_neigh: q after calling dtrevc");
//  dvout<T>(debug.logfil,n*n+3*n,workl,debug.ndigit,"_dneigh: workl after call to dtrevc");
  
  if(ierr!=0) return;
  
  /*
c
c     %------------------------------------------------%
c     | Scale the returning eigenvectors so that their |
c     | euclidean norms are all one. LAPACK subroutine |
c     | dtrevc returns each eigenvector normalized so  |
c     | that the element of largest magnitude has      |
c     | magnitude 1; here the magnitude of a complex   |
c     | number (x,y) is taken to be |x| + |y|.         |
c     %------------------------------------------------%
c
  */
  
  iconj=0;
  for(i=1;i<=n;i++)
    if(fabs(ritzi[i-1])<=zero)
    {
      /*
c
c           %----------------------%
c           | Real eigenvalue case |
c           %----------------------%
c 
      */
      temp=dnrm2<T>(n,&q[0+ldq*(i-1)],1);
      dscal<T>(n,one/temp,&q[0+ldq*(i-1)],1);
    }
    else
    {
      /*
c
c           %-------------------------------------------%
c           | Complex conjugate pair case. Note that    |
c           | since the real and imaginary part of      |
c           | the eigenvector are stored in consecutive |
c           | columns, we further normalize by the      |
c           | square root of two.                       |
c           %-------------------------------------------%
c
      */
      
      if(iconj==0)
      {
	temp=dlapy2<T>(dnrm2<T>(n,&q[0+ldq*(i-1)],1),dnrm2<T>(n,&q[0+ldq*i],1));
	dscal<T>(n,one/temp,&q[0+ldq*(i-1)],1);
	dscal<T>(n,one/temp,&q[0+ldq*i],1);
	iconj=1;
      }
      else
	iconj=0;
    }
    
    //dvout<T>(debug.logfil,n,bounds,debug.ndigit,"_dneigh: bounds before call to dgemv");
   // dvout<T>(debug.logfil,n*n+3*n,workl,debug.ndigit,"_dneigh: workl before call to dgemv");
   // dmout<T>(debug.logfil,n,n,q,ldq,debug.ndigit,"_neigh: q matrix before calling dgemv");
    dgemv<T>("T",n,n,one,q,ldq,bounds,1,zero,workl,1);
   // dvout<T>(debug.logfil,n*n+3*n,workl,debug.ndigit,"_dneigh: workl after call to dgemv");
    //dvout<T>(debug.logfil,n,bounds,debug.ndigit,"_dneigh: bounds after call to dgemv");
    
    if(msglvl>1)
      dvout<T>(debug.logfil,n,workl,debug.ndigit,"_neigh: Last row of the eigenvector matrix for H");
    
    /*
c
c     %----------------------------%
c     | Compute the Ritz estimates |
c     %----------------------------%
c
    */
    
    iconj=0;
    for(i=1;i<=n;i++)
    {
     // std::cout << "i=" << i << "\n";
    //  dvout<T>(debug.logfil,n,bounds,debug.ndigit,"_dneigh: bounds at this i");
      if(fabs(ritzi[i-1])<=zero)
      {
	/*
c
c           %----------------------%
c           | Real eigenvalue case |
c           %----------------------%
c 
        */
   //     std::cout << "workl[" << i-1 << "]=" << workl[i-1] << "\n";
	bounds[i-1]=rnorm*fabs(workl[i-1]);
      }
      else
      {
	/*
c
c           %-------------------------------------------%
c           | Complex conjugate pair case. Note that    |
c           | since the real and imaginary part of      |
c           | the eigenvector are stored in consecutive |
c           | columns, we need to take the magnitude    |
c           | of the last components of the two vectors |
c           %-------------------------------------------%
c
        */
	if(iconj==0)
	{
	  bounds[i-1]=rnorm*dlapy2<T>(workl[i-1],workl[i]);
	  bounds[i]=bounds[i-1];
	  iconj=1;
	}
	else
	  iconj=0;
      }
    }
      
      if(msglvl>2)
      {
	dvout<T>(debug.logfil,n,ritzr,debug.ndigit,"_neigh: Real part of the eigenvalues of H");
	dvout<T>(debug.logfil,n,ritzi,debug.ndigit,"_neigh: Imaginary part of the eigenvalues of H");
	dvout<T>(debug.logfil,n,bounds,debug.ndigit,"_neigh: Ritz estimates for the eigenvalues of H");
      }
      
      second(timing.t1);
      timing.tneigh+=timing.t1-timing.t0;
      
      return;
      /*
c
c     %---------------%
c     | End of dneigh |
c     %---------------%
c
      */
}
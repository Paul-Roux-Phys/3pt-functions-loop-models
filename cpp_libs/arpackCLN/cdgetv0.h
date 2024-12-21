/* Conversion of dgetv0 to C++ by Chris Scullard May 18, 2014 */


/* 
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dgetv0
c
c\Description: 
c  Generate a random initial residual vector for the Arnoldi process.
c  Force the residual vector to be in the range of the operator OP.  
c
c\Usage:
c  call dgetv0
c     ( IDO, BMAT, ITRY, INITV, N, J, V, LDV, RESID, RNORM, 
c       IPNTR, WORKD, IERR )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.  IDO must be zero on the first
c          call to dgetv0.
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    This is for the initialization phase to force the
c                    starting vector into the range of OP.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c          IDO = 99: done
c          -------------------------------------------------------------
c
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B in the (generalized)
c          eigenvalue problem A*x = lambda*B*x.
c          B = 'I' -> standard eigenvalue problem A*x = lambda*x
c          B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
c
c  ITRY    Integer.  (INPUT)
c          ITRY counts the number of times that dgetv0 is called.  
c          It should be set to 1 on the initial call to dgetv0.
c
c  INITV   Logical variable.  (INPUT)
c          .TRUE.  => the initial residual vector is given in RESID.
c          .FALSE. => generate a random initial residual vector.
c
c  N       Integer.  (INPUT)
c          Dimension of the problem.
c
c  J       Integer.  (INPUT)
c          Index of the residual vector to be generated, with respect to
c          the Arnoldi process.  J > 1 in case of a "restart".
c
c  V       Double precision N by J array.  (INPUT)
c          The first J-1 columns of V contain the current Arnoldi basis
c          if this is a "restart".
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling 
c          program.
c
c  RESID   Double precision array of length N.  (INPUT/OUTPUT)
c          Initial residual vector to be generated.  If RESID is 
c          provided, force RESID into the range of the operator OP.
c
c  RNORM   Double precision scalar.  (OUTPUT)
c          B-norm of the generated residual.
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c
c  WORKD   Double precision work array of length 2*N.  (REVERSE COMMUNICATION).
c          On exit, WORK(1:N) = B*RESID to be used in SSAITR.
c
c  IERR    Integer.  (OUTPUT)
c          =  0: Normal exit.
c          = -1: Cannot generate a nontrivial restarted residual vector
c                in the range of the operator OP.
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
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c
c\Routines called:
c     second  ARPACK utility routine for timing.
c     dvout   ARPACK utility routine for vector output.
c     dlarnv  LAPACK routine for generating a random vector.
c     dgemv   Level 2 BLAS routine for matrix vector multiplication.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     ddot    Level 1 BLAS that computes the scalar product of two vectors. 
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c
*/

#include "cdebug.h"
#include "cstat.h"
#include "fortranfuncs.h"
#include "cdgemv.h"
#include <cln/real.h>
#include <cln/output.h>
#include <cln/real_io.h>

/* CLN version */

template<typename T>
void dgetv0(ARint& ido, char bmat, int itry, bool initv, ARint n, int j, T* v, int ldv, T* resid, T& rnorm, int* ipntr, T* workd, int& ierr, int digits)
{
  //debugstruct debug;
  /* Parameters */
  
  cln::float_format_t precision=cln::float_format(digits);
  const T one = cln::cl_float(1,precision);
  const T zero = cln::cl_float(0,precision);
  
  /* Local Scalars & Arrays */
  
  static bool first, inits, orth;
  int idist,jj;
  static int msglvl, iter;
  static int iseed[4];
  static T rnorm0;
  
  /* Initialize the seed of the LAPACK */
  /* random number generator           */
  
  if(inits)
  {
    iseed[0]=1;
    iseed[1]=3;
    iseed[2]=5;
    iseed[3]=7;
    inits=false;
  }
  
  if(ido==0)
  {
    
    /* Initialize timing statistics  */
    /* & message level for debugging */
    
    // call second(t0)
    msglvl=debug.mgetv0;
    
    ierr=0;
    iter=0;
    first=false;
    orth=false;
    
    /*
c
c        %-----------------------------------------------------%
c        | Possibly generate a random starting vector in RESID |
c        | Use a LAPACK random number generator used by the    |
c        | matrix generation routines.                         |
c        |    idist = 1: uniform (0,1)  distribution;          |
c        |    idist = 2: uniform (-1,1) distribution;          |
c        |    idist = 3: normal  (0,1)  distribution;          |
c        %-----------------------------------------------------%
c
    */
    
    if(!initv)
    {
      //idist=2; // leftover from old fortran call. Uniform random number wanted.
      dlarnv<T>(iseed,n,resid,digits);
   //   for(jj=1;jj<=n;jj++)
   //     resid[jj-1]=one;  /* !!! Take this out. Set all components to 1 for debugging !!! */
   //   rnorm=dnrm2<T>(n,resid,1); // this too!!!
    }
    
 //   dvout<T>(debug.logfil, n, resid, debug.ndigit,"_getv0: resid");
    
    /*
c 
c        %----------------------------------------------------------%
c        | Force the starting vector into the range of OP to handle |
c        | the generalized problem when B is possibly (singular).   |
c        %----------------------------------------------------------%
c
    */
    //std::cout << "bmat=" << bmat << "\n";
    // call second(t2)
    if(bmat=='G')
    {
      timing.nopx+=1;
      ipntr[1]=1; // fixed indicees
      ipntr[2]=n+1;
      copy<T>(n,resid,1,workd,1);
      ido=-1;
      return;
    }
  }
  
  /*
c 
c     %-----------------------------------------%
c     | Back from computing OP*(initial-vector) |
c     %-----------------------------------------%
c
  */
  //std::cout << "first = " << first << "\n";
  if(first) goto label20;
  
  /*
c
c     %-----------------------------------------------%
c     | Back from computing B*(orthogonalized-vector) |
c     %-----------------------------------------------%
c
  */
  //std::cout << "orth = " << orth << "\n";
  if(orth) goto label40;
  //std::cout << "after orth \n";
  
  if(bmat=='G')
  {
    // call second(t3);
    timing.tmvopx=timing.tmvopx+(timing.t3-timing.t2);
  }
  
  /*
c 
c     %------------------------------------------------------%
c     | Starting vector is now in the range of OP; r = OP*r; |
c     | Compute B-norm of starting vector.                   |
c     %------------------------------------------------------%
c
  */
  
  // call second(t2);
  first=true;
  if(bmat=='G')
  {
    timing.nbx+=1;
    copy<T>(n,&workd[n],1,resid,1);
    ipntr[1]=n+1; // fixed indices
    ipntr[2]=1;
    ido=2;
    return;
  }
  else if(bmat=='I')
  {
    //std::cout << "doing the copy \n";
    copy<T>(n,resid,1,workd,1);
  }
  //std::cout << "back from copy \n";
label20:
  if(bmat=='G')
  {
    //call second(t3);
    timing.tmvbx+=timing.t3-timing.t2;
  }
  
  first=false;
  if(bmat=='G')
  {
    rnorm0=ddot<T>(n,resid,1,workd,1,digits);
    rnorm0=sqrt(cln::abs(rnorm0));
  }
  else if (bmat=='I')
    rnorm0=dnrm2<T>(n,resid,1);
  
  //std::cout << "computed norm \n";
  
  rnorm=rnorm0;
  
  /*
c
c     %---------------------------------------------%
c     | Exit if this is the very first Arnoldi step |
c     %---------------------------------------------%
c
  */
  
  if(j==1) goto label50;
  
  /*
c 
c     %----------------------------------------------------------------
c     | Otherwise need to B-orthogonalize the starting vector against |
c     | the current Arnoldi basis using Gram-Schmidt with iter. ref.  |
c     | This is the case where an invariant subspace is encountered   |
c     | in the middle of the Arnoldi factorization.                   |
c     |                                                               |
c     |       s = V^{T}*B*r;   r = r - V*s;                           |
c     |                                                               |
c     | Stopping criteria used for iter. ref. is discussed in         |
c     | Parlett's book, page 107 and in Gragg & Reichel TOMS paper.   |
c     %---------------------------------------------------------------%
c
  */
  
  orth=true;
  
label30:
  
  dgemv<T>("T",n,j-1,one,v,ldv,workd,1,zero,&workd[n],1,digits);
  dgemv<T>("N",n,j-1,-one,v,ldv,&workd[n],1,one,resid,1,digits);
  
  /*
c 
c     %----------------------------------------------------------%
c     | Compute the B-norm of the orthogonalized starting vector |
c     %----------------------------------------------------------%
c
  */
  
  // call second(t2);
  if(bmat=='G')
  {
    timing.nbx+=1;
    copy<T>(n,resid,1,&workd[n],1);
    ipntr[1]=n+1; // fixed indices
    ipntr[2]=1;
    ido=2;
    return;
  }
  else if(bmat=='I')
    copy<T>(n,resid,1,workd,1);
  
label40:
  
  if(bmat=='G')
  {
    // call second(t3);
    timing.tmvbx+=timing.t3-timing.t2;
  }
  
  
  if(bmat=='G')
  {
    rnorm=ddot<T>(n,resid,1,workd,1,digits);
    rnorm=sqrt(cln::abs(rnorm));
  }
  else if(bmat=='I')
    rnorm=dnrm2<T>(n,resid,1);
  
  //std::cout << "rnorm = " << rnorm << "\n"; 
  
  /*
c
c     %--------------------------------------%
c     | Check for further orthogonalization. |
c     %--------------------------------------%
c
  */
  
  if(msglvl>2)
  {
    dvout<T>(debug.logfil,1,&rnorm0,debug.ndigit,"_getv0: re-orthonalization ; rnorm0 is");
    dvout<T>(debug.logfil,1,&rnorm,debug.ndigit,"_getv0: re-orthonalization ; rnorm0 is");
  }
  
  if(rnorm>(0.717*rnorm0)) goto label50;
  
  iter+=1;
  if(iter<=1)
  {
    /*
c
c        %-----------------------------------%
c        | Perform iterative refinement step |
c        %-----------------------------------%
c
    */
    
    rnorm0=rnorm;
    goto label30;
  }
  else
  {
    /*
c
c        %------------------------------------%
c        | Iterative refinement step "failed" |
c        %------------------------------------%
c
    */
    for(jj=1;jj<=n;jj++)
      resid[jj-1]=zero;
    rnorm=zero;
    ierr=-1;
  }
  
label50:
  
  if(msglvl>0)
    dvout<T>(debug.logfil, 1, &rnorm, debug.ndigit,"_getv0: B-norm of initial / restarted starting vector");
  if(msglvl>2)
    dvout<T>(debug.logfil, n, resid, debug.ndigit,"_getv0: initial / restarted starting vector");
  
  ido=99;
  
  // call second(t1);
 // timing.tgetv0+=timing.t1-timing.t0;
  
  return;
  
  /*
c
c     %---------------%
c     | End of dgetv0 |
c     %---------------%
c
  */
}

template<typename T>
void dgetv0(ARint& ido, char bmat, int itry, bool initv, ARint n, int j, T* v, int ldv, T* resid, T& rnorm, int* ipntr, T* workd, int& ierr)
{
  //debugstruct debug;
  /* Parameters */
  
  const T one=1.0000000000000;
  const T zero=0.0000000000000;
  
  /* Local Scalars & Arrays */
  
  static bool first, inits, orth;
  int idist,jj;
  static int msglvl, iter;
  static int iseed[4];
  static T rnorm0;
  
  /* Initialize the seed of the LAPACK */
  /* random number generator           */
  
  if(inits)
  {
    iseed[0]=1;
    iseed[1]=3;
    iseed[2]=5;
    iseed[3]=7;
    inits=false;
  }
  
  if(ido==0)
  {
    
    /* Initialize timing statistics  */
    /* & message level for debugging */
    
    // call second(t0)
    msglvl=debug.mgetv0;
    
    ierr=0;
    iter=0;
    first=false;
    orth=false;
    
    /*
c
c        %-----------------------------------------------------%
c        | Possibly generate a random starting vector in RESID |
c        | Use a LAPACK random number generator used by the    |
c        | matrix generation routines.                         |
c        |    idist = 1: uniform (0,1)  distribution;          |
c        |    idist = 2: uniform (-1,1) distribution;          |
c        |    idist = 3: normal  (0,1)  distribution;          |
c        %-----------------------------------------------------%
c
    */
    
    if(!initv)
    {
      //idist=2; // leftover from old fortran call. Uniform random number wanted.
      dlarnv<T>(iseed,n,resid);
      for(jj=1;jj<=n;jj++)
        resid[jj-1]=0.5;  /* !!! Take this out. Set all components to 1 for debugging !!! */
      rnorm=dnrm2<T>(n,resid,1); // this too!!!
    }
    
 //   dvout<T>(debug.logfil, n, resid, debug.ndigit,"_getv0: resid");
    
    /*
c 
c        %----------------------------------------------------------%
c        | Force the starting vector into the range of OP to handle |
c        | the generalized problem when B is possibly (singular).   |
c        %----------------------------------------------------------%
c
    */
    //std::cout << "bmat=" << bmat << "\n";
    // call second(t2)
    if(bmat=='G')
    {
      timing.nopx+=1;
      ipntr[1]=1; // fixed indicees
      ipntr[2]=n+1;
      copy<T>(n,resid,1,workd,1);
      ido=-1;
      return;
    }
  }
  
  /*
c 
c     %-----------------------------------------%
c     | Back from computing OP*(initial-vector) |
c     %-----------------------------------------%
c
  */
  //std::cout << "first = " << first << "\n";
  if(first) goto label20;
  
  /*
c
c     %-----------------------------------------------%
c     | Back from computing B*(orthogonalized-vector) |
c     %-----------------------------------------------%
c
  */
  //std::cout << "orth = " << orth << "\n";
  if(orth) goto label40;
  //std::cout << "after orth \n";
  
  if(bmat=='G')
  {
    // call second(t3);
    timing.tmvopx=timing.tmvopx+(timing.t3-timing.t2);
  }
  
  /*
c 
c     %------------------------------------------------------%
c     | Starting vector is now in the range of OP; r = OP*r; |
c     | Compute B-norm of starting vector.                   |
c     %------------------------------------------------------%
c
  */
  
  // call second(t2);
  first=true;
  if(bmat=='G')
  {
    timing.nbx+=1;
    copy(n,&workd[n],1,resid,1);
    ipntr[1]=n+1; // fixed indices
    ipntr[2]=1;
    ido=2;
    return;
  }
  else if(bmat=='I')
  {
    //std::cout << "doing the copy \n";
    copy<T>(n,resid,1,workd,1);
  }
  //std::cout << "back from copy \n";
label20:
  if(bmat=='G')
  {
    //call second(t3);
    timing.tmvbx+=timing.t3-timing.t2;
  }
  
  first=false;
  if(bmat=='G')
  {
    rnorm0=ddot<T>(n,resid,1,workd,1);
    rnorm0=sqrt(fabs(rnorm0));
  }
  else if (bmat=='I')
    rnorm0=dnrm2<T>(n,resid,1);
  
  //std::cout << "computed norm \n";
  
  rnorm=rnorm0;
  
  /*
c
c     %---------------------------------------------%
c     | Exit if this is the very first Arnoldi step |
c     %---------------------------------------------%
c
  */
  
  if(j==1) goto label50;
  
  /*
c 
c     %----------------------------------------------------------------
c     | Otherwise need to B-orthogonalize the starting vector against |
c     | the current Arnoldi basis using Gram-Schmidt with iter. ref.  |
c     | This is the case where an invariant subspace is encountered   |
c     | in the middle of the Arnoldi factorization.                   |
c     |                                                               |
c     |       s = V^{T}*B*r;   r = r - V*s;                           |
c     |                                                               |
c     | Stopping criteria used for iter. ref. is discussed in         |
c     | Parlett's book, page 107 and in Gragg & Reichel TOMS paper.   |
c     %---------------------------------------------------------------%
c
  */
  
  orth=true;
  
label30:
  
  dgemv<T>("T",n,j-1,one,v,ldv,workd,1,zero,&workd[n],1);
  dgemv<T>("N",n,j-1,-one,v,ldv,&workd[n],1,one,resid,1);
  
  /*
c 
c     %----------------------------------------------------------%
c     | Compute the B-norm of the orthogonalized starting vector |
c     %----------------------------------------------------------%
c
  */
  
  // call second(t2);
  if(bmat=='G')
  {
    timing.nbx+=1;
    copy<T>(n,resid,1,&workd[n],1);
    ipntr[1]=n+1; // fixed indices
    ipntr[2]=1;
    ido=2;
    return;
  }
  else if(bmat=='I')
    copy<T>(n,resid,1,workd,1);
  
label40:
  
  if(bmat=='G')
  {
    // call second(t3);
    timing.tmvbx+=timing.t3-timing.t2;
  }
  
  
  if(bmat=='G')
  {
    rnorm=ddot<T>(n,resid,1,workd,1);
    rnorm=sqrt(abs(rnorm));
  }
  else if(bmat=='I')
    rnorm=dnrm2<T>(n,resid,1);
  
  //std::cout << "rnorm = " << rnorm << "\n"; 
  
  /*
c
c     %--------------------------------------%
c     | Check for further orthogonalization. |
c     %--------------------------------------%
c
  */
  
  if(msglvl>2)
  {
    dvout<T>(debug.logfil,1,&rnorm0,debug.ndigit,"_getv0: re-orthonalization ; rnorm0 is");
    dvout<T>(debug.logfil,1,&rnorm,debug.ndigit,"_getv0: re-orthonalization ; rnorm0 is");
  }
  
  if(rnorm>(0.717*rnorm0)) goto label50;
  
  iter+=1;
  if(iter<=1)
  {
    /*
c
c        %-----------------------------------%
c        | Perform iterative refinement step |
c        %-----------------------------------%
c
    */
    
    rnorm0=rnorm;
    goto label30;
  }
  else
  {
    /*
c
c        %------------------------------------%
c        | Iterative refinement step "failed" |
c        %------------------------------------%
c
    */
    for(jj=1;jj<=n;jj++)
      resid[jj-1]=zero;
    rnorm=zero;
    ierr=-1;
  }
  
label50:
  
  if(msglvl>0)
    dvout<T>(debug.logfil, 1, &rnorm, debug.ndigit,"_getv0: B-norm of initial / restarted starting vector");
  if(msglvl>2)
    dvout<T>(debug.logfil, n, resid, debug.ndigit,"_getv0: initial / restarted starting vector");
  
  ido=99;
  
  // call second(t1);
 // timing.tgetv0+=timing.t1-timing.t0;
  
  return;
  
  /*
c
c     %---------------%
c     | End of dgetv0 |
c     %---------------%
c
  */
}
    

/*
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dnaitr
c
c\Description: 
c  Reverse communication interface for applying NP additional steps to 
c  a K step nonsymmetric Arnoldi factorization.
c
c  Input:  OP*V_{k}  -  V_{k}*H = r_{k}*e_{k}^T
c
c          with (V_{k}^T)*B*V_{k} = I, (V_{k}^T)*B*r_{k} = 0.
c
c  Output: OP*V_{k+p}  -  V_{k+p}*H = r_{k+p}*e_{k+p}^T
c
c          with (V_{k+p}^T)*B*V_{k+p} = I, (V_{k+p}^T)*B*r_{k+p} = 0.
c
c  where OP and B are as in dnaupd.  The B-norm of r_{k+p} is also
c  computed and returned.
c
c\Usage:
c  call dnaitr
c     ( IDO, BMAT, N, K, NP, NB, RESID, RNORM, V, LDV, H, LDH, 
c       IPNTR, WORKD, INFO )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y.
c                    This is for the restart phase to force the new
c                    starting vector into the range of OP.
c          IDO =  1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y,
c                    IPNTR(3) is the pointer into WORK for B * X.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y.
c          IDO = 99: done
c          -------------------------------------------------------------
c          When the routine is used in the "shift-and-invert" mode, the
c          vector B * Q is already available and do not need to be
c          recompute in forming OP * Q.
c
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B that defines the
c          semi-inner product for the operator OP.  See dnaupd.
c          B = 'I' -> standard eigenvalue problem A*x = lambda*x
c          B = 'G' -> generalized eigenvalue problem A*x = lambda*M**x
c
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  K       Integer.  (INPUT)
c          Current size of V and H.
c
c  NP      Integer.  (INPUT)
c          Number of additional Arnoldi steps to take.
c
c  NB      Integer.  (INPUT)
c          Blocksize to be used in the recurrence.          
c          Only work for NB = 1 right now.  The goal is to have a 
c          program that implement both the block and non-block method.
c
c  RESID   Double precision array of length N.  (INPUT/OUTPUT)
c          On INPUT:  RESID contains the residual vector r_{k}.
c          On OUTPUT: RESID contains the residual vector r_{k+p}.
c
c  RNORM   Double precision scalar.  (INPUT/OUTPUT)
c          B-norm of the starting residual on input.
c          B-norm of the updated residual r_{k+p} on output.
c
c  V       Double precision N by K+NP array.  (INPUT/OUTPUT)
c          On INPUT:  V contains the Arnoldi vectors in the first K 
c          columns.
c          On OUTPUT: V contains the new NP Arnoldi vectors in the next
c          NP columns.  The first K columns are unchanged.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling 
c          program.
c
c  H       Double precision (K+NP) by (K+NP) array.  (INPUT/OUTPUT)
c          H is used to store the generated upper Hessenberg matrix.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling 
c          program.
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c          Pointer to mark the starting locations in the WORK for 
c          vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X.
c          IPNTR(2): pointer to the current result vector Y.
c          IPNTR(3): pointer to the vector B * X when used in the 
c                    shift-and-invert mode.  X is the current operand.
c          -------------------------------------------------------------
c          
c  WORKD   Double precision work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The calling program should not 
c          use WORKD as temporary workspace during the iteration !!!!!!
c          On input, WORKD(1:N) = B*RESID and is used to save some 
c          computation at the first step.
c
c  INFO    Integer.  (OUTPUT)
c          = 0: Normal exit.
c          > 0: Size of the spanning invariant subspace of OP found.
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
c     dgetv0  ARPACK routine to generate the initial vector.
c     ivout   ARPACK utility routine that prints integers.
c     second  ARPACK utility routine for timing.
c     dmout   ARPACK utility routine that prints matrices
c     dvout   ARPACK utility routine that prints vectors.
c     dlabad  LAPACK routine that computes machine constants.
c     dlamch  LAPACK routine that determines machine constants.
c     dlascl  LAPACK routine for careful scaling of a matrix.
c     dlanhs  LAPACK routine that computes various norms of a matrix.
c     dgemv   Level 2 BLAS routine for matrix vector multiplication.
c     daxpy   Level 1 BLAS that computes a vector triad.
c     dscal   Level 1 BLAS that scales a vector.
c     dcopy   Level 1 BLAS that copies one vector to another .
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
c\Revision history:
c     xx/xx/92: Version ' 2.4'
c
c\SCCS Information: @(#) 
c FILE: naitr.F   SID: 2.4   DATE OF SID: 8/27/96   RELEASE: 2
c
c\Remarks
c  The algorithm implemented is:
c  
c  restart = .false.
c  Given V_{k} = [v_{1}, ..., v_{k}], r_{k}; 
c  r_{k} contains the initial residual vector even for k = 0;
c  Also assume that rnorm = || B*r_{k} || and B*r_{k} are already 
c  computed by the calling program.
c
c  betaj = rnorm ; p_{k+1} = B*r_{k} ;
c  For  j = k+1, ..., k+np  Do
c     1) if ( betaj < tol ) stop or restart depending on j.
c        ( At present tol is zero )
c        if ( restart ) generate a new starting vector.
c     2) v_{j} = r(j-1)/betaj;  V_{j} = [V_{j-1}, v_{j}];  
c        p_{j} = p_{j}/betaj
c     3) r_{j} = OP*v_{j} where OP is defined as in dnaupd
c        For shift-invert mode p_{j} = B*v_{j} is already available.
c        wnorm = || OP*v_{j} ||
c     4) Compute the j-th step residual vector.
c        w_{j} =  V_{j}^T * B * OP * v_{j}
c        r_{j} =  OP*v_{j} - V_{j} * w_{j}
c        H(:,j) = w_{j};
c        H(j,j-1) = rnorm
c        rnorm = || r_(j) ||
c        If (rnorm > 0.717*wnorm) accept step and go back to 1)
c     5) Re-orthogonalization step:
c        s = V_{j}'*B*r_{j}
c        r_{j} = r_{j} - V_{j}*s;  rnorm1 = || r_{j} ||
c        alphaj = alphaj + s_{j};   
c     6) Iterative refinement step:
c        If (rnorm1 > 0.717*rnorm) then
c           rnorm = rnorm1
c           accept step and go back to 1)
c        Else
c           rnorm = rnorm1
c           If this is the first time in step 6), go to 5)
c           Else r_{j} lies in the span of V_{j} numerically.
c              Set r_{j} = 0 and rnorm = 0; go to 1)
c        EndIf 
c  End Do
c
c\EndLib
c
c-----------------------------------------------------------------------
c
*/

#include "cdebug.h"
#include "cstat.h"
#include "fortranfuncs.h"
#include <algorithm>
#include "cdgemv.h"
#include <cln/real.h>
#include <cln/output.h>
#include <cln/real_io.h>

/* CLN version */

template <typename T>
void dnaitr(int& ido, char bmat, int n, int k, int np, int nb, T* resid, T& rnorm, T* v, int ldv, T* h, int ldh, int* ipntr, T* workd, int& info, int digits)
{
  /*
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
  */
  static bool first = true;
  static bool orth1, orth2, rstart, step3, step4;
  int i, infol, jj;
  static int ierr, ipj, irj, ivj, iter, itry, j, msglvl;
  static T betaj, ovfl, rnorm1, smlnum, ulp, unfl, wnorm;
  T temp1, tst1;
  cln::float_format_t precision=cln::float_format(digits);
  const T one = cln::cl_float(1,precision);
  const T zero = cln::cl_float(0,precision);
  const T n0717 = cln::cl_float(0.717,precision); // CLN representation of 0.717
  
  /*
c
c     %-----------------------%
c     | Local Array Arguments | 
c     %-----------------------%
c
  */
  T xtemp[2];
  
  /*
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
  */
  
 // std::cout << "dnaitr: step3 = " << step3 << "\n";
 // dvout<T>(debug.logfil, n, workd, debug.ndigit,"_dnaitr: workd");
  
  if(first)
  {
    
    /*
c
c        %-----------------------------------------%
c        | Set machine-dependent constants for the |
c        | the splitting and deflation criterion.  |
c        | If norm(H) <= sqrt(OVFL),               |
c        | overflow should not occur.              |
c        | REFERENCE: LAPACK subroutine dlahqr     |
c        %-----------------------------------------%
c
    */
    
    unfl = dlamch<T>("S",digits); // safe minimum
 //   print_float(std::cout,cln::default_print_flags,unfl);
 //   std::cout << "\n";
    ovfl = one/unfl;
 //   dlabad<T>(unfl,ovfl);
    ulp = dlamch<T>("P",digits); // precision
 //   print_float(std::cout,cln::default_print_flags,ulp);
 //   std::cout << "\n";
    smlnum = unfl*(cln::cl_float(n,precision)/ulp);
 /*   dvout<T>(debug.logfil, 1, &smlnum, debug.ndigit,"_dnaitr: smlnum");
    dvout<T>(debug.logfil, 1, &unfl, debug.ndigit,"_dnaitr: unfl");
    dvout<T>(debug.logfil, 1, &ovfl, debug.ndigit,"_dnaitr: ovfl");
    dvout<T>(debug.logfil, 1, &ulp, debug.ndigit,"_dnaitr: ulp");*/
    first=false;
  }
  
  if(ido==0)
  {
    /*
c 
c        %-------------------------------%
c        | Initialize timing statistics  |
c        | & message level for debugging |
c        %-------------------------------%
c
    */
    
    second(timing.t0);
    msglvl = debug.mnaitr;
    
    /*
c 
c        %------------------------------%
c        | Initial call to this routine |
c        %------------------------------%
c
    */
    
    info=0;
    step3=false;
    step4=false;
    rstart=false;
    orth1=false;
    orth2=false;
    j=k+1;
    ipj=1;
    irj=ipj+n;
    ivj=irj+n;
  }
  
  /*
c 
c     %-------------------------------------------------%
c     | When in reverse communication mode one of:      |
c     | STEP3, STEP4, ORTH1, ORTH2, RSTART              |
c     | will be .true. when ....                        |
c     | STEP3: return from computing OP*v_{j}.          |
c     | STEP4: return from computing B-norm of OP*v_{j} |
c     | ORTH1: return from computing B-norm of r_{j+1}  |
c     | ORTH2: return from computing B-norm of          |
c     |        correction to the residual vector.       |
c     | RSTART: return from OP computations needed by   |
c     |         dgetv0.                                 |
c     %-------------------------------------------------%
c
  */
  
  if(step3) goto label50;
  if(step4) goto label60;
  if(orth1) goto label70;
  if(orth2) goto label90;
  if(rstart) goto label30;
  /*
c
c     %-----------------------------%
c     | Else this is the first step |
c     %-----------------------------%
c
c     %--------------------------------------------------------------%
c     |                                                              |
c     |        A R N O L D I     I T E R A T I O N     L O O P       |
c     |                                                              |
c     | Note:  B*r_{j-1} is already in WORKD(1:N)=WORKD(IPJ:IPJ+N-1) |
c     %--------------------------------------------------------------%
  */
 
label1000:

  if(msglvl>1)
  {
    ivout<T>(debug.logfil, 1, &j, debug.ndigit,"_naitr: generating Arnoldi vector number");
    dvout<T>(debug.logfil, 1, &rnorm, debug.ndigit,"_naitr: B-norm of the current residual is");
  }
  /*
c 
c        %---------------------------------------------------%
c        | STEP 1: Check if the B norm of j-th residual      |
c        | vector is zero. Equivalent to determing whether   |
c        | an exact j-step Arnoldi factorization is present. |
c        %---------------------------------------------------%
c
  */
  betaj=rnorm;
  if(rnorm>zero) goto label40;
  /*
c
c           %---------------------------------------------------%
c           | Invariant subspace found, generate a new starting |
c           | vector which is orthogonal to the current Arnoldi |
c           | basis and continue the iteration.                 |
c           %---------------------------------------------------%
c
  */
  
  if(msglvl>0)
    ivout<T>(debug.logfil, 1, &j, debug.ndigit,"_naitr: ****** RESTART AT STEP ******");
    
  /*
c 
c           %---------------------------------------------%
c           | ITRY is the loop variable that controls the |
c           | maximum amount of times that a restart is   |
c           | attempted. NRSTRT is used by stat.h         |
c           %---------------------------------------------%
c 
  */
  
  betaj=zero;
  timing.nrstrt+=1;
  itry=1;
  
label20:

  rstart=true;
  ido=0;
  
label30:

  /*
c
c           %--------------------------------------%
c           | If in reverse communication mode and |
c           | RSTART = .true. flow returns here.   |
c           %--------------------------------------%
c
  */
//  dvout(debug.logfil,n,resid,debug.ndigit,"_naitr: resid before getv0");
  dgetv0<T>(ido,bmat,itry,false,n,j,v,ldv,resid,rnorm,ipntr,workd,ierr,digits);
//  dvout(debug.logfil,n,resid,debug.ndigit,"_naitr: resid after getv0");
  
  if(ido!=99) return;
  if(ierr<0)
  {
    itry+=1;
    if(itry<3) goto label20;
    
    /*
c
c              %------------------------------------------------%
c              | Give up after several restart attempts.        |
c              | Set INFO to the size of the invariant subspace |
c              | which spans OP and exit.                       |
c              %------------------------------------------------%
c
    */
    
    info=j-1;
    // call second(t1);
    timing.tnaitr+=timing.t1-timing.t0;
    ido=99;
    return;
  }
  
label40:

  /*
c
c        %---------------------------------------------------------%
c        | STEP 2:  v_{j} = r_{j-1}/rnorm and p_{j} = p_{j}/rnorm  |
c        | Note that p_{j} = B*r_{j-1}. In order to avoid overflow |
c        | when reciprocating a small RNORM, test against lower    |
c        | machine bound.                                          |
c        %---------------------------------------------------------%
c
  */
  
  copy<T>(n,resid,1,&v[0+ldv*(j-1)],1);
  if(rnorm>=unfl)
  {
    temp1=one/rnorm;
 //   dvout<T>(debug.logfil, n, &v[0+ldv*(j-1)], debug.ndigit,"_dnaitr: v before dscal");
    dscal<T>(n,temp1,&v[0+ldv*(j-1)],1);
 //   dvout<T>(debug.logfil, n, &v[0+ldv*(j-1)], debug.ndigit,"_dnaitr: v after dscal");
 //   dvout<T>(debug.logfil, n, &workd[ipj-1], debug.ndigit,"_dnaitr: workd before dscal");
    dscal<T>(n,temp1,&workd[ipj-1],1);
 //   dvout<T>(debug.logfil, n, &workd[ipj-1], debug.ndigit,"_dnaitr: workd after dscal");
  }
  else
  {
    /*
c
c            %-----------------------------------------%
c            | To scale both v_{j} and p_{j} carefully |
c            | use LAPACK routine SLASCL               |
c            %-----------------------------------------%
c
    */
    
    dlascl<T>("General",i,i,rnorm,one,n,1,&v[0+ldv*(j-1)],n,infol);
    dlascl<T>("General",i,i,rnorm,one,n,1,&workd[ipj-1],n,infol);
    
  }
  
  /*
c
c        %------------------------------------------------------%
c        | STEP 3:  r_{j} = OP*v_{j}; Note that p_{j} = B*v_{j} |
c        | Note that this is not quite yet r_{j}. See STEP 4    |
c        %------------------------------------------------------%
c
  */
  
  step3=true;
  timing.nopx+=1;
  second(timing.t2);
  copy<T>(n,&v[0+ldv*(j-1)],1,&workd[ivj-1],1);
  ipntr[1]=ivj; // fixed indices
  ipntr[2]=irj;
  ipntr[3]=ipj;
  ido=1;
  
  /*
c 
c        %-----------------------------------%
c        | Exit in order to compute OP*v_{j} |
c        %-----------------------------------%
c 
  */
  
  return;
  
label50:
  
  /*
c 
c        %----------------------------------%
c        | Back from reverse communication; |
c        | WORKD(IRJ:IRJ+N-1) := OP*v_{j}   |
c        | if step3 = .true.                |
c        %----------------------------------%
c
  */
  
  second(timing.t3);
  timing.tmvopx+=timing.t3-timing.t2;
  
  step3=false;
  /*
c
c        %------------------------------------------%
c        | Put another copy of OP*v_{j} into RESID. |
c        %------------------------------------------%
c
  */
  
  copy<T>(n,&workd[irj-1],1,resid,1);
  
  /*
c 
c        %---------------------------------------%
c        | STEP 4:  Finish extending the Arnoldi |
c        |          factorization to length j.   |
c        %---------------------------------------%
c
  */
  
  second(timing.t2);
  if(bmat=='G')
  {
    timing.nbx+=1;
    step4=true;
    ipntr[1]=irj; // fixed indices
    ipntr[2]=ipj;
    ido=2;
    
    /*
c 
c           %-------------------------------------%
c           | Exit in order to compute B*OP*v_{j} |
c           %-------------------------------------%
c 
    */
    
    return;
  }
  else if(bmat=='I')
    copy<T>(n,resid,1,&workd[ipj-1],1);
  
label60:

  /*
c 
c        %----------------------------------%
c        | Back from reverse communication; |
c        | WORKD(IPJ:IPJ+N-1) := B*OP*v_{j} |
c        | if step4 = .true.                |
c        %----------------------------------%
c
  */
  
  if(bmat=='G')
  {
    second(timing.t3);
    timing.tmvbx+=timing.t3-timing.t2;
  }
  
  step4=false;
  
  /*
c
c        %-------------------------------------%
c        | The following is needed for STEP 5. |
c        | Compute the B-norm of OP*v_{j}.     |
c        %-------------------------------------%
c
  */
  
  if(bmat=='G')
  {
    wnorm=ddot<T>(n,resid,1,&workd[ipj-1],1,digits);
    wnorm=sqrt(abs(wnorm));
  }
  else if(bmat=='I')
    wnorm=dnrm2<T>(n,resid,1);
  
  /*
c
c        %-----------------------------------------%
c        | Compute the j-th residual corresponding |
c        | to the j step factorization.            |
c        | Use Classical Gram Schmidt and compute: |
c        | w_{j} <-  V_{j}^T * B * OP * v_{j}      |
c        | r_{j} <-  OP*v_{j} - V_{j} * w_{j}      |
c        %-----------------------------------------%
c
c
c        %------------------------------------------%
c        | Compute the j Fourier coefficients w_{j} |
c        | WORKD(IPJ:IPJ+N-1) contains B*OP*v_{j}.  |
c        %------------------------------------------%
c 
  */
 // dvout<T>(debug.logfil, n, &workd[ipj-1], debug.ndigit,"_dnaitr: workd before dgemv");
 // dmout<T>(debug.logfil,k+np,k+np,h,ldh,debug.ndigit,"_naitr: H before dgemv");
  dgemv<T>("T",n,j,one,v,ldv,&workd[ipj-1],1,zero,&h[0+ldh*(j-1)],1,digits);
 // dmout<T>(debug.logfil,k+np,k+np,h,ldh,debug.ndigit,"_naitr: H after dgemv");
 // dvout<T>(debug.logfil, n, &workd[ipj-1], debug.ndigit,"_dnaitr: workd after first dgemv");
  
  /*
c
c        %--------------------------------------%
c        | Orthogonalize r_{j} against V_{j}.   |
c        | RESID contains OP*v_{j}. See STEP 3. | 
c        %--------------------------------------%
c
  */
  
  dgemv<T>("N",n,j,-one,v,ldv,&h[0+ldh*(j-1)],1,one,resid,1,digits);
 // dvout<T>(debug.logfil, n, &workd[ipj-1], debug.ndigit,"_dnaitr: workd after second dgemv");
  
  if(j>1) h[j-1+ldh*(j-2)]=betaj;
  
  second(timing.t4);
  
  orth1=true;
  
  second(timing.t2);
  if(bmat=='G')
  {
    timing.nbx+=1;
    copy<T>(n,resid,1,&workd[irj-1],1);
    ipntr[1]=irj; // fixed indices
    ipntr[2]=ipj;
    ido=2;
    
    /*
c 
c           %----------------------------------%
c           | Exit in order to compute B*r_{j} |
c           %----------------------------------%
c 
    */
    
    return;
  }
  else if(bmat=='I')
    copy<T>(n,resid,1,&workd[ipj-1],1);
  
label70:

  /*
c 
c        %---------------------------------------------------%
c        | Back from reverse communication if ORTH1 = .true. |
c        | WORKD(IPJ:IPJ+N-1) := B*r_{j}.                    |
c        %---------------------------------------------------%
c
  */
  
  if(bmat=='G')
  {
    second(timing.t3);
    timing.tmvbx+=timing.t3-timing.t2;
  }
  
  orth1=false;
  
  /*
c
c        %------------------------------%
c        | Compute the B-norm of r_{j}. |
c        %------------------------------%
c
  */
  
  if(bmat=='G')
  {
    rnorm=ddot<T>(n,resid,1,&workd[ipj-1],1,digits);
    rnorm=sqrt(abs(rnorm));
  }
  else if(bmat=='I')
    rnorm = dnrm2<T>(n,resid,1);
  
  /*
c 
c        %-----------------------------------------------------------%
c        | STEP 5: Re-orthogonalization / Iterative refinement phase |
c        | Maximum NITER_ITREF tries.                                |
c        |                                                           |
c        |          s      = V_{j}^T * B * r_{j}                     |
c        |          r_{j}  = r_{j} - V_{j}*s                         |
c        |          alphaj = alphaj + s_{j}                          |
c        |                                                           |
c        | The stopping criteria used for iterative refinement is    |
c        | discussed in Parlett's book SEP, page 107 and in Gragg &  |
c        | Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.         |
c        | Determine if we need to correct the residual. The goal is |
c        | to enforce ||v(:,1:j)^T * r_{j}|| .le. eps * || r_{j} ||  |
c        | The following test determines whether the sine of the     |
c        | angle between  OP*x and the computed residual is less     |
c        | than or equal to 0.717.                                   |
c        %-----------------------------------------------------------%
c
  */
  
  if(rnorm>n0717*wnorm) goto label100;
  iter=0;
  timing.nrorth+=1;
  
  /*
c 
c        %---------------------------------------------------%
c        | Enter the Iterative refinement phase. If further  |
c        | refinement is necessary, loop back here. The loop |
c        | variable is ITER. Perform a step of Classical     |
c        | Gram-Schmidt using all the Arnoldi vectors V_{j}  |
c        %---------------------------------------------------%
c 
  */
  
label80:

  if(msglvl>2)
  {
    xtemp[0]=wnorm;
    xtemp[1]=rnorm;
    dvout<T>(debug.logfil,2,xtemp,debug.ndigit,"_naitr: re-orthonalization; wnorm and rnorm are");
    dvout<T>(debug.logfil,j,&h[0+ldh*(j-1)],debug.ndigit,"_naitr: j-th column of H");
  }
  
  /*
c
c        %----------------------------------------------------%
c        | Compute V_{j}^T * B * r_{j}.                       |
c        | WORKD(IRJ:IRJ+J-1) = v(:,1:J)'*WORKD(IPJ:IPJ+N-1). |
c        %----------------------------------------------------%
c
  */
  
  dgemv<T>("T",n,j,one,v,ldv,&workd[ipj-1],1,zero,&workd[irj-1],1,digits);
  
  /*
c
c        %---------------------------------------------%
c        | Compute the correction to the residual:     |
c        | r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1). |
c        | The correction to H is v(:,1:J)*H(1:J,1:J)  |
c        | + v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j.         |
c        %---------------------------------------------%
c
  */
  
  dgemv<T>("N",n,j,-one,v,ldv,&workd[irj-1],1,one,resid,1,digits);
  axpy<T>(j,one,&workd[irj-1],1,&h[0+ldh*(j-1)],1);
  
  orth2=true;
  second(timing.t2);
  if(bmat=='G')
  {
    timing.nbx+=1;
    copy<T>(n,resid,1,&workd[irj-1],1);
    ipntr[1]=irj; //fixed indices
    ipntr[2]=ipj;
    ido=2;
    
    /*
c 
c           %-----------------------------------%
c           | Exit in order to compute B*r_{j}. |
c           | r_{j} is the corrected residual.  |
c           %-----------------------------------%
c 
    */
    
    return;
  }
  else if(bmat=='I')
    copy<T>(n,resid,1,&workd[ipj-1],1);

label90:

  /*
c
c        %---------------------------------------------------%
c        | Back from reverse communication if ORTH2 = .true. |
c        %---------------------------------------------------%
c
  */
  
  if(bmat=='G')
  {
    second(timing.t3);
    timing.tmvbx+=timing.t3-timing.t2;
  }
  
  /*
c
c        %-----------------------------------------------------%
c        | Compute the B-norm of the corrected residual r_{j}. |
c        %-----------------------------------------------------%
c 
  */
  
  if(bmat=='G')
  {
    rnorm1=ddot<T>(n,resid,1,&workd[ipj-1],1,digits);
    rnorm1=cln::sqrt(cln::abs(rnorm1));
  }
  else if(bmat=='I')
    rnorm1=dnrm2<T>(n,resid,1);
  
  if((msglvl>0)&&(iter>0))
  {
    ivout<T>(debug.logfil,1,&j,debug.ndigit,"_naitr: Iterative refinement for Arnoldi residual");
    if(msglvl>2)
    {
      xtemp[0]=rnorm;
      xtemp[1]=rnorm1;
      dvout(debug.logfil,2,xtemp,debug.ndigit,"_naitr: iterative refinement ; rnorm and rnorm1 are");
    }
  }
  
  /*
c
c        %-----------------------------------------%
c        | Determine if we need to perform another |
c        | step of re-orthogonalization.           |
c        %-----------------------------------------%
c
  */
  
  if(rnorm1>(n0717*rnorm) )
  {
    /*
c
c           %---------------------------------------%
c           | No need for further refinement.       |
c           | The cosine of the angle between the   |
c           | corrected residual vector and the old |
c           | residual vector is greater than 0.717 |
c           | In other words the corrected residual |
c           | and the old residual vector share an  |
c           | angle of less than arcCOS(0.717)      |
c           %---------------------------------------%
c
    */
    
    rnorm=rnorm1;
  }
  else
  {
    /*
c
c           %-------------------------------------------%
c           | Another step of iterative refinement step |
c           | is required. NITREF is used by stat.h     |
c           %-------------------------------------------%
c
   */
    
    timing.nitref+=1;
    rnorm=rnorm1;
    iter+=1;
    if(iter<=1) goto label80;
    
    /*
c
c           %-------------------------------------------------%
c           | Otherwise RESID is numerically in the span of V |
c           %-------------------------------------------------%
c
    */
    for(jj=1;jj<=n;jj++)
      resid[jj-1]=zero;
    rnorm=zero;
  }
  
  /*
c 
c        %----------------------------------------------%
c        | Branch here directly if iterative refinement |
c        | wasn't necessary or after at most NITER_REF  |
c        | steps of iterative refinement.               |
c        %----------------------------------------------%
c 
  */
  
label100:

  rstart=false;
  orth2=false;
  
  second(timing.t5);
  timing.titref+=timing.t5-timing.t4;
  /*
c 
c        %------------------------------------%
c        | STEP 6: Update  j = j+1;  Continue |
c        %------------------------------------%
c
  */
  j+=1;
  if(j>(k+np))
  {
    second(timing.t1);
    timing.tnaitr+=timing.t1-timing.t0;
    ido=99;
    for(i=std::max(1,k);i<=(k+np-1);i++)
    {
      /*
c     
c              %--------------------------------------------%
c              | Check for splitting and deflation.         |
c              | Use a standard test as in the QR algorithm |
c              | REFERENCE: LAPACK subroutine dlahqr        |
c              %--------------------------------------------%
c  
      */
      tst1=cln::abs(h[i-1+ldh*(i-1)])+cln::abs(h[i+ldh*i]);
      if(tst1==zero)
	tst1=dlanhs<T>("1",k+np,h,ldh,&workd[n],digits);
      if(cln::abs(h[i+ldh*(i-1)])<=std::max(ulp*tst1,smlnum))
	h[i+ldh*(i-1)]=zero;
    }
    
    if(msglvl>2)
    {
      dmout<T>(debug.logfil,k+np,k+np,h,ldh,debug.ndigit,"_naitr: Final upper Hessenberg matrix H of order K+NP");
    }
    
    return;
  }
  /*
c
c        %--------------------------------------------------------%
c        | Loop back to extend the factorization by another step. |
c        %--------------------------------------------------------%
c
  */ 
  goto label1000;
  /*
c 
c     %---------------------------------------------------------------%
c     |                                                               |
c     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
c     |                                                               |
c     %---------------------------------------------------------------%
c
  */
  return;
  /*
c
c     %---------------%
c     | End of dnaitr |
c     %---------------%
c
  */
}

template <typename T>
void dnaitr(int& ido, char bmat, int n, int k, int np, int nb, T* resid, T& rnorm, T* v, int ldv, T* h, int ldh, int* ipntr, T* workd, int& info)
{
  /*
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
  */
  static bool first = true;
  static bool orth1, orth2, rstart, step3, step4;
  int i, infol, jj;
  static int ierr, ipj, irj, ivj, iter, itry, j, msglvl;
  static T betaj, ovfl, rnorm1, smlnum, ulp, unfl, wnorm;
  T temp1, tst1;
  const T one=1.0;
  const T zero=0.0;
  
  /*
c
c     %-----------------------%
c     | Local Array Arguments | 
c     %-----------------------%
c
  */
  T xtemp[2];
  
  /*
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
  */
  
 // std::cout << "dnaitr: step3 = " << step3 << "\n";
 // dvout<T>(debug.logfil, n, workd, debug.ndigit,"_dnaitr: workd");
//  dvout<T>(debug.logfil,n,resid,debug.ndigit,"_naitr: resid at the start");
//  dvout<T>(debug.logfil,3*n,workd,debug.ndigit,"_naitr: workd at the start");
  
  if(first)
  {
    
    /*
c
c        %-----------------------------------------%
c        | Set machine-dependent constants for the |
c        | the splitting and deflation criterion.  |
c        | If norm(H) <= sqrt(OVFL),               |
c        | overflow should not occur.              |
c        | REFERENCE: LAPACK subroutine dlahqr     |
c        %-----------------------------------------%
c
    */
    
    unfl = dlamch<T>("S"); // safe minimum
    ovfl = one/unfl;
    dlabad<T>(unfl,ovfl);
    ulp = dlamch<T>("P"); // precision
    smlnum = unfl*(n/ulp);
 /*   dvout<T>(debug.logfil, 1, &smlnum, debug.ndigit,"_dnaitr: smlnum");
    dvout<T>(debug.logfil, 1, &unfl, debug.ndigit,"_dnaitr: unfl");
    dvout<T>(debug.logfil, 1, &ovfl, debug.ndigit,"_dnaitr: ovfl");
    dvout<T>(debug.logfil, 1, &ulp, debug.ndigit,"_dnaitr: ulp");*/
    first=false;
  }
  
  if(ido==0)
  {
    /*
c 
c        %-------------------------------%
c        | Initialize timing statistics  |
c        | & message level for debugging |
c        %-------------------------------%
c
    */
    
    second(timing.t0);
    msglvl = debug.mnaitr;
    
    /*
c 
c        %------------------------------%
c        | Initial call to this routine |
c        %------------------------------%
c
    */
    
    info=0;
    step3=false;
    step4=false;
    rstart=false;
    orth1=false;
    orth2=false;
    j=k+1;
    ipj=1;
    irj=ipj+n;
    ivj=irj+n;
  }
  
  /*
c 
c     %-------------------------------------------------%
c     | When in reverse communication mode one of:      |
c     | STEP3, STEP4, ORTH1, ORTH2, RSTART              |
c     | will be .true. when ....                        |
c     | STEP3: return from computing OP*v_{j}.          |
c     | STEP4: return from computing B-norm of OP*v_{j} |
c     | ORTH1: return from computing B-norm of r_{j+1}  |
c     | ORTH2: return from computing B-norm of          |
c     |        correction to the residual vector.       |
c     | RSTART: return from OP computations needed by   |
c     |         dgetv0.                                 |
c     %-------------------------------------------------%
c
  */
  
  if(step3) goto label50;
  if(step4) goto label60;
  if(orth1) goto label70;
  if(orth2) goto label90;
  if(rstart) goto label30;
  /*
c
c     %-----------------------------%
c     | Else this is the first step |
c     %-----------------------------%
c
c     %--------------------------------------------------------------%
c     |                                                              |
c     |        A R N O L D I     I T E R A T I O N     L O O P       |
c     |                                                              |
c     | Note:  B*r_{j-1} is already in WORKD(1:N)=WORKD(IPJ:IPJ+N-1) |
c     %--------------------------------------------------------------%
  */
 
label1000:

  if(msglvl>1)
  {
    ivout<T>(debug.logfil, 1, &j, debug.ndigit,"_naitr: generating Arnoldi vector number");
    dvout<T>(debug.logfil, 1, &rnorm, debug.ndigit,"_naitr: B-norm of the current residual is");
  }
  /*
c 
c        %---------------------------------------------------%
c        | STEP 1: Check if the B norm of j-th residual      |
c        | vector is zero. Equivalent to determing whether   |
c        | an exact j-step Arnoldi factorization is present. |
c        %---------------------------------------------------%
c
  */
  betaj=rnorm;
  if(rnorm>zero) goto label40;
  /*
c
c           %---------------------------------------------------%
c           | Invariant subspace found, generate a new starting |
c           | vector which is orthogonal to the current Arnoldi |
c           | basis and continue the iteration.                 |
c           %---------------------------------------------------%
c
  */
  
  if(msglvl>0)
    ivout<T>(debug.logfil, 1, &j, debug.ndigit,"_naitr: ****** RESTART AT STEP ******");
    
  /*
c 
c           %---------------------------------------------%
c           | ITRY is the loop variable that controls the |
c           | maximum amount of times that a restart is   |
c           | attempted. NRSTRT is used by stat.h         |
c           %---------------------------------------------%
c 
  */
  
  betaj=zero;
  timing.nrstrt+=1;
  itry=1;
  
label20:

  rstart=true;
  ido=0;
  
label30:

  /*
c
c           %--------------------------------------%
c           | If in reverse communication mode and |
c           | RSTART = .true. flow returns here.   |
c           %--------------------------------------%
c
  */
//  dvout(debug.logfil,n,resid,debug.ndigit,"_naitr: resid before getv0");
  dgetv0<T>(ido,bmat,itry,false,n,j,v,ldv,resid,rnorm,ipntr,workd,ierr);
//  dvout(debug.logfil,n,resid,debug.ndigit,"_naitr: resid after getv0");
  
  if(ido!=99) return;
  if(ierr<0)
  {
    itry+=1;
    if(itry<3) goto label20;
    
    /*
c
c              %------------------------------------------------%
c              | Give up after several restart attempts.        |
c              | Set INFO to the size of the invariant subspace |
c              | which spans OP and exit.                       |
c              %------------------------------------------------%
c
    */
    
    info=j-1;
    // call second(t1);
    timing.tnaitr+=timing.t1-timing.t0;
    ido=99;
    return;
  }
  
label40:

  /*
c
c        %---------------------------------------------------------%
c        | STEP 2:  v_{j} = r_{j-1}/rnorm and p_{j} = p_{j}/rnorm  |
c        | Note that p_{j} = B*r_{j-1}. In order to avoid overflow |
c        | when reciprocating a small RNORM, test against lower    |
c        | machine bound.                                          |
c        %---------------------------------------------------------%
c
  */
  
  copy<T>(n,resid,1,&v[0+ldv*(j-1)],1);
  if(rnorm>=unfl)
  {
    temp1=one/rnorm;
 //   dvout<T>(debug.logfil, n, &v[0+ldv*(j-1)], debug.ndigit,"_dnaitr: v before dscal");
    dscal<T>(n,temp1,&v[0+ldv*(j-1)],1);
 //   dvout<T>(debug.logfil, n, &v[0+ldv*(j-1)], debug.ndigit,"_dnaitr: v after dscal");
 //   dvout<T>(debug.logfil, n, &workd[ipj-1], debug.ndigit,"_dnaitr: workd before dscal");
    dscal<T>(n,temp1,&workd[ipj-1],1);
 //   dvout<T>(debug.logfil, n, &workd[ipj-1], debug.ndigit,"_dnaitr: workd after dscal");
  }
  else
  {
    /*
c
c            %-----------------------------------------%
c            | To scale both v_{j} and p_{j} carefully |
c            | use LAPACK routine SLASCL               |
c            %-----------------------------------------%
c
    */
    
    dlascl<T>("General",i,i,rnorm,one,n,1,&v[0+ldv*(j-1)],n,infol);
    dlascl<T>("General",i,i,rnorm,one,n,1,&workd[ipj-1],n,infol);
    
  }
  
  /*
c
c        %------------------------------------------------------%
c        | STEP 3:  r_{j} = OP*v_{j}; Note that p_{j} = B*v_{j} |
c        | Note that this is not quite yet r_{j}. See STEP 4    |
c        %------------------------------------------------------%
c
  */
  
  step3=true;
  timing.nopx+=1;
  second(timing.t2);
  copy<T>(n,&v[0+ldv*(j-1)],1,&workd[ivj-1],1);
  ipntr[1]=ivj; // fixed indices
  ipntr[2]=irj;
  ipntr[3]=ipj;
  ido=1;
  
  /*
c 
c        %-----------------------------------%
c        | Exit in order to compute OP*v_{j} |
c        %-----------------------------------%
c 
  */
  
  return;
  
label50:
  
  /*
c 
c        %----------------------------------%
c        | Back from reverse communication; |
c        | WORKD(IRJ:IRJ+N-1) := OP*v_{j}   |
c        | if step3 = .true.                |
c        %----------------------------------%
c
  */
  
  second(timing.t3);
  timing.tmvopx+=timing.t3-timing.t2;
  
  step3=false;
  /*
c
c        %------------------------------------------%
c        | Put another copy of OP*v_{j} into RESID. |
c        %------------------------------------------%
c
  */
  
  copy<T>(n,&workd[irj-1],1,resid,1);
  
  /*
c 
c        %---------------------------------------%
c        | STEP 4:  Finish extending the Arnoldi |
c        |          factorization to length j.   |
c        %---------------------------------------%
c
  */
  
  second(timing.t2);
  if(bmat=='G')
  {
    timing.nbx+=1;
    step4=true;
    ipntr[1]=irj; // fixed indices
    ipntr[2]=ipj;
    ido=2;
    
    /*
c 
c           %-------------------------------------%
c           | Exit in order to compute B*OP*v_{j} |
c           %-------------------------------------%
c 
    */
    
    return;
  }
  else if(bmat=='I')
    copy<T>(n,resid,1,&workd[ipj-1],1);
  
label60:

  /*
c 
c        %----------------------------------%
c        | Back from reverse communication; |
c        | WORKD(IPJ:IPJ+N-1) := B*OP*v_{j} |
c        | if step4 = .true.                |
c        %----------------------------------%
c
  */
  
  if(bmat=='G')
  {
    second(timing.t3);
    timing.tmvbx+=timing.t3-timing.t2;
  }
  
  step4=false;
  
  /*
c
c        %-------------------------------------%
c        | The following is needed for STEP 5. |
c        | Compute the B-norm of OP*v_{j}.     |
c        %-------------------------------------%
c
  */
  
  if(bmat=='G')
  {
    wnorm=ddot<T>(n,resid,1,&workd[ipj-1],1);
    wnorm=sqrt(abs(wnorm));
  }
  else if(bmat=='I')
    wnorm=dnrm2<T>(n,resid,1);
//  dvout<T>(debug.logfil, 1, &wnorm, debug.ndigit,"_dnaitr: wnorm before dgemv");
  
  /*
c
c        %-----------------------------------------%
c        | Compute the j-th residual corresponding |
c        | to the j step factorization.            |
c        | Use Classical Gram Schmidt and compute: |
c        | w_{j} <-  V_{j}^T * B * OP * v_{j}      |
c        | r_{j} <-  OP*v_{j} - V_{j} * w_{j}      |
c        %-----------------------------------------%
c
c
c        %------------------------------------------%
c        | Compute the j Fourier coefficients w_{j} |
c        | WORKD(IPJ:IPJ+N-1) contains B*OP*v_{j}.  |
c        %------------------------------------------%
c 
  */
 // dvout<T>(debug.logfil, n, &workd[ipj-1], debug.ndigit,"_dnaitr: workd before dgemv");
 // dmout<T>(debug.logfil,k+np,k+np,h,ldh,debug.ndigit,"_naitr: H before dgemv");
  dgemv<T>("T",n,j,one,v,ldv,&workd[ipj-1],1,zero,&h[0+ldh*(j-1)],1);
 // dmout<T>(debug.logfil,k+np,k+np,h,ldh,debug.ndigit,"_naitr: H after dgemv");
 // dvout<T>(debug.logfil, n, &workd[ipj-1], debug.ndigit,"_dnaitr: workd after first dgemv");
  
  /*
c
c        %--------------------------------------%
c        | Orthogonalize r_{j} against V_{j}.   |
c        | RESID contains OP*v_{j}. See STEP 3. | 
c        %--------------------------------------%
c
  */
//  dmout<T>(debug.logfil,k+np,k+np,h,ldh,debug.ndigit,"_naitr: h before second dgemv");
//  dmout<T>(debug.logfil,n,k+np,v,ldv,debug.ndigit,"_naitr: v before second dgemv");
 // dvout<T>(debug.logfil,n,resid,debug.ndigit,"_naitr: resid before second dgemv");
 // dvout<T>(debug.logfil,n,&h[0+ldh*(j-1)],debug.ndigit,"_naitr: x vector going into second dgemv");
  dgemv<T>("N",n,j,-one,v,ldv,&h[0+ldh*(j-1)],1,one,resid,1);
//  dmout<T>(debug.logfil,n,k+np,v,ldv,debug.ndigit,"_naitr: v after second dgemv");
//  dvout<T>(debug.logfil,n,resid,debug.ndigit,"_naitr: resid after second dgemv");
 // dvout<T>(debug.logfil, n, &workd[ipj-1], debug.ndigit,"_dnaitr: workd after second dgemv");
  
  if(j>1) h[j-1+ldh*(j-2)]=betaj;
  
  second(timing.t4);
  
  orth1=true;
  
  second(timing.t2);
  if(bmat=='G')
  {
    timing.nbx+=1;
    copy<T>(n,resid,1,&workd[irj-1],1);
    ipntr[1]=irj; // fixed indices
    ipntr[2]=ipj;
    ido=2;
    
    /*
c 
c           %----------------------------------%
c           | Exit in order to compute B*r_{j} |
c           %----------------------------------%
c 
    */
    
    return;
  }
  else if(bmat=='I')
    copy<T>(n,resid,1,&workd[ipj-1],1);
//  dvout<T>(debug.logfil,n,resid,debug.ndigit,"_naitr: resid after copy");
  
label70:

  /*
c 
c        %---------------------------------------------------%
c        | Back from reverse communication if ORTH1 = .true. |
c        | WORKD(IPJ:IPJ+N-1) := B*r_{j}.                    |
c        %---------------------------------------------------%
c
  */
  
  if(bmat=='G')
  {
    second(timing.t3);
    timing.tmvbx+=timing.t3-timing.t2;
  }
  
  orth1=false;
  
  /*
c
c        %------------------------------%
c        | Compute the B-norm of r_{j}. |
c        %------------------------------%
c
  */
  
  if(bmat=='G')
  {
    rnorm=ddot<T>(n,resid,1,&workd[ipj-1],1);
    rnorm=sqrt(fabs(rnorm));
  }
  else if(bmat=='I')
    rnorm = dnrm2<T>(n,resid,1);
  
  /*
c 
c        %-----------------------------------------------------------%
c        | STEP 5: Re-orthogonalization / Iterative refinement phase |
c        | Maximum NITER_ITREF tries.                                |
c        |                                                           |
c        |          s      = V_{j}^T * B * r_{j}                     |
c        |          r_{j}  = r_{j} - V_{j}*s                         |
c        |          alphaj = alphaj + s_{j}                          |
c        |                                                           |
c        | The stopping criteria used for iterative refinement is    |
c        | discussed in Parlett's book SEP, page 107 and in Gragg &  |
c        | Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.         |
c        | Determine if we need to correct the residual. The goal is |
c        | to enforce ||v(:,1:j)^T * r_{j}|| .le. eps * || r_{j} ||  |
c        | The following test determines whether the sine of the     |
c        | angle between  OP*x and the computed residual is less     |
c        | than or equal to 0.717.                                   |
c        %-----------------------------------------------------------%
c
  */
  
  if(rnorm>T(0.717)*wnorm) goto label100;
  iter=0;
  timing.nrorth+=1;
  
  /*
c 
c        %---------------------------------------------------%
c        | Enter the Iterative refinement phase. If further  |
c        | refinement is necessary, loop back here. The loop |
c        | variable is ITER. Perform a step of Classical     |
c        | Gram-Schmidt using all the Arnoldi vectors V_{j}  |
c        %---------------------------------------------------%
c 
  */
  
label80:

  if(msglvl>2)
  {
    xtemp[0]=wnorm;
    xtemp[1]=rnorm;
    dvout<T>(debug.logfil,2,xtemp,debug.ndigit,"_naitr: re-orthonalization; wnorm and rnorm are");
    dvout<T>(debug.logfil,j,&h[0+ldh*(j-1)],debug.ndigit,"_naitr: j-th column of H");
  }
  
  /*
c
c        %----------------------------------------------------%
c        | Compute V_{j}^T * B * r_{j}.                       |
c        | WORKD(IRJ:IRJ+J-1) = v(:,1:J)'*WORKD(IPJ:IPJ+N-1). |
c        %----------------------------------------------------%
c
  */
  
  dgemv<T>("T",n,j,one,v,ldv,&workd[ipj-1],1,zero,&workd[irj-1],1);
  
  /*
c
c        %---------------------------------------------%
c        | Compute the correction to the residual:     |
c        | r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1). |
c        | The correction to H is v(:,1:J)*H(1:J,1:J)  |
c        | + v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j.         |
c        %---------------------------------------------%
c
  */
  
  dgemv<T>("N",n,j,-one,v,ldv,&workd[irj-1],1,one,resid,1);
  axpy<T>(j,one,&workd[irj-1],1,&h[0+ldh*(j-1)],1);
  
  orth2=true;
  second(timing.t2);
  if(bmat=='G')
  {
    timing.nbx+=1;
    copy<T>(n,resid,1,&workd[irj-1],1);
    ipntr[1]=irj; //fixed indices
    ipntr[2]=ipj;
    ido=2;
    
    /*
c 
c           %-----------------------------------%
c           | Exit in order to compute B*r_{j}. |
c           | r_{j} is the corrected residual.  |
c           %-----------------------------------%
c 
    */
    
    return;
  }
  else if(bmat=='I')
    copy<T>(n,resid,1,&workd[ipj-1],1);

label90:

  /*
c
c        %---------------------------------------------------%
c        | Back from reverse communication if ORTH2 = .true. |
c        %---------------------------------------------------%
c
  */
  
  if(bmat=='G')
  {
    second(timing.t3);
    timing.tmvbx+=timing.t3-timing.t2;
  }
  
  /*
c
c        %-----------------------------------------------------%
c        | Compute the B-norm of the corrected residual r_{j}. |
c        %-----------------------------------------------------%
c 
  */
  
  if(bmat=='G')
  {
    rnorm1=ddot<T>(n,resid,1,&workd[ipj-1],1);
    rnorm1=sqrt(abs(rnorm1));
  }
  else if(bmat=='I')
    rnorm1=dnrm2<T>(n,resid,1);
  
  if((msglvl>0)&&(iter>0))
  {
    ivout<T>(debug.logfil,1,&j,debug.ndigit,"_naitr: Iterative refinement for Arnoldi residual");
    if(msglvl>2)
    {
      xtemp[0]=rnorm;
      xtemp[1]=rnorm1;
      dvout(debug.logfil,2,xtemp,debug.ndigit,"_naitr: iterative refinement ; rnorm and rnorm1 are");
    }
  }
  
  /*
c
c        %-----------------------------------------%
c        | Determine if we need to perform another |
c        | step of re-orthogonalization.           |
c        %-----------------------------------------%
c
  */
  
  if(rnorm1>(T(0.717)*rnorm) )
  {
    /*
c
c           %---------------------------------------%
c           | No need for further refinement.       |
c           | The cosine of the angle between the   |
c           | corrected residual vector and the old |
c           | residual vector is greater than 0.717 |
c           | In other words the corrected residual |
c           | and the old residual vector share an  |
c           | angle of less than arcCOS(0.717)      |
c           %---------------------------------------%
c
    */
    
    rnorm=rnorm1;
  }
  else
  {
    /*
c
c           %-------------------------------------------%
c           | Another step of iterative refinement step |
c           | is required. NITREF is used by stat.h     |
c           %-------------------------------------------%
c
   */
    
    timing.nitref+=1;
    rnorm=rnorm1;
    iter+=1;
    if(iter<=1) goto label80;
    
    /*
c
c           %-------------------------------------------------%
c           | Otherwise RESID is numerically in the span of V |
c           %-------------------------------------------------%
c
    */
    for(jj=1;jj<=n;jj++)
      resid[jj-1]=zero;
    rnorm=zero;
  }
  
  /*
c 
c        %----------------------------------------------%
c        | Branch here directly if iterative refinement |
c        | wasn't necessary or after at most NITER_REF  |
c        | steps of iterative refinement.               |
c        %----------------------------------------------%
c 
  */
  
label100:

  rstart=false;
  orth2=false;
  
  second(timing.t5);
  timing.titref+=timing.t5-timing.t4;
  /*
c 
c        %------------------------------------%
c        | STEP 6: Update  j = j+1;  Continue |
c        %------------------------------------%
c
  */
  j+=1;
  if(j>(k+np))
  {
    second(timing.t1);
    timing.tnaitr+=timing.t1-timing.t0;
    ido=99;
    for(i=std::max(1,k);i<=(k+np-1);i++)
    {
      /*
c     
c              %--------------------------------------------%
c              | Check for splitting and deflation.         |
c              | Use a standard test as in the QR algorithm |
c              | REFERENCE: LAPACK subroutine dlahqr        |
c              %--------------------------------------------%
c  
      */
      tst1=fabs(h[i-1+ldh*(i-1)])+fabs(h[i+ldh*i]);
      if(tst1==zero)
        tst1=dlanhs<T>("1",k+np,h,ldh,&workd[n]);
      if(fabs(h[i+ldh*(i-1)])<=std::max(ulp*tst1,smlnum))
        h[i+ldh*(i-1)]=zero;
    }
    
    if(msglvl>2)
    {
      dmout<T>(debug.logfil,k+np,k+np,h,ldh,debug.ndigit,"_naitr: Final upper Hessenberg matrix H of order K+NP");
    }
    
    return;
  }
  /*
c
c        %--------------------------------------------------------%
c        | Loop back to extend the factorization by another step. |
c        %--------------------------------------------------------%
c
  */ 
  goto label1000;
  /*
c 
c     %---------------------------------------------------------------%
c     |                                                               |
c     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
c     |                                                               |
c     %---------------------------------------------------------------%
c
  */
  return;
  /*
c
c     %---------------%
c     | End of dnaitr |
c     %---------------%
c
  */
}
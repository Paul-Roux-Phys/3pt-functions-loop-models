/*
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dnapps
c
c\Description:
c  Given the Arnoldi factorization
c
c     A*V_{k} - V_{k}*H_{k} =gr_{k+p}*e_{k+p}^T,
c
c  apply NP implicit shifts resulting in
c
c     A*(V_{k}*Q) - (V_{k}*Q)*(Q^T* H_{k}*Q) = r_{k+p}*e_{k+p}^T * Q
c
c  where Q is an orthogonal matrix which is the product of rotations
c  and reflections resulting from the NP bulge chage sweeps.
c  The updated Arnoldi factorization becomes:
c
c     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.
c
c\Usage:
c  call dnapps
c     ( N, KEV, NP, SHIFTR, SHIFTI, V, LDV, H, LDH, RESID, Q, LDQ, 
c       WORKL, WORKD )
c
c\Arguments
c  N       Integer.  (INPUT)
c          Problem size, i.e. size of matrix A.
c
c  KEV     Integer.  (INPUT/OUTPUT)
c          KEV+NP is the size of the input matrix H.
c          KEV is the size of the updated matrix HNEW.  KEV is only 
c          updated on ouput when fewer than NP shifts are applied in
c          order to keep the conjugate pair together.
c
c  NP      Integer.  (INPUT)
c          Number of implicit shifts to be applied.
c
c  SHIFTR, Double precision array of length NP.  (INPUT)
c  SHIFTI  Real and imaginary part of the shifts to be applied.
c          Upon, entry to dnapps, the shifts must be sorted so that the 
c          conjugate pairs are in consecutive locations.
c
c  V       Double precision N by (KEV+NP) array.  (INPUT/OUTPUT)
c          On INPUT, V contains the current KEV+NP Arnoldi vectors.
c          On OUTPUT, V contains the updated KEV Arnoldi vectors
c          in the first KEV columns of V.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  H       Double precision (KEV+NP) by (KEV+NP) array.  (INPUT/OUTPUT)
c          On INPUT, H contains the current KEV+NP by KEV+NP upper 
c          Hessenber matrix of the Arnoldi factorization.
c          On OUTPUT, H contains the updated KEV by KEV upper Hessenberg
c          matrix in the KEV leading submatrix.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.
c
c  RESID   Double precision array of length N.  (INPUT/OUTPUT)
c          On INPUT, RESID contains the the residual vector r_{k+p}.
c          On OUTPUT, RESID is the update residual vector rnew_{k} 
c          in the first KEV locations.
c
c  Q       Double precision KEV+NP by KEV+NP work array.  (WORKSPACE)
c          Work array used to accumulate the rotations and reflections
c          during the bulge chase sweep.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  WORKL   Double precision work array of length (KEV+NP).  (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.
c
c  WORKD   Double precision work array of length 2*N.  (WORKSPACE)
c          Distributed array used in the application of the accumulated
c          orthogonal matrix Q.
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
c
c\Routines called:
c     ivout   ARPACK utility routine that prints integers.
c     second  ARPACK utility routine for timing.
c     dmout   ARPACK utility routine that prints matrices.
c     dvout   ARPACK utility routine that prints vectors.
c     dlabad  LAPACK routine that computes machine constants.
c     dlacpy  LAPACK matrix copy routine.
c     dlamch  LAPACK routine that determines machine constants. 
c     dlanhs  LAPACK routine that computes various norms of a matrix.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dlarf   LAPACK routine that applies Householder reflection to
c             a matrix.
c     dlarfg  LAPACK Householder reflection construction routine.
c     dlartg  LAPACK Givens rotation construction routine.
c     dlaset  LAPACK matrix initialization routine.
c     dgemv   Level 2 BLAS routine for matrix vector multiplication.
c     daxpy   Level 1 BLAS that computes a vector triad.
c     dcopy   Level 1 BLAS that copies one vector to another .
c     dscal   Level 1 BLAS that scales a vector.
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
c FILE: napps.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
c
c\Remarks
c  1. In this version, each shift is applied to all the sublocks of
c     the Hessenberg matrix H and not just to the submatrix that it
c     comes from. Deflation as in LAPACK routine dlahqr (QR algorithm
c     for upper Hessenberg matrices ) is used.
c     The subdiagonals of H are enforced to be non-negative.
c
c\EndLib
c
c-----------------------------------------------------------------------
c
*/

#include "cdebug.h"
#include "cstat.h"
#include "fortranfuncs.h"
#include "cdlarf.h"
#include <cln/real.h>
#include <cln/output.h>
#include <cln/real_io.h>

/* CLN version */

template <typename T>
void dnapps(int n, int& kev, int np, T* shiftr, T* shifti, T* v, int ldv, T* h, int ldh, T* resid, T* q, int ldq, T* workl, T* workd, int digits)
{
  /*
c
c     %------------%
c     | Parameters |
c     %------------%
c
  */
  cln::float_format_t precision=cln::float_format(digits);
  const T zero = cln::cl_float(0,precision);
  const T one = cln::cl_float(1,precision);
  const T two = cln::cl_float(2,precision);
  //std::cout << "entering dnapps \n";
  /*
c
c     %------------------------%
c     | Local Scalars & Arrays |
c     %------------------------%
c
  */
  
  static bool first=true;
  static T ovfl, smlnum, ulp, unfl;
  bool cconj;
  int i,iend,ir,istart,j,jj,kplusp,msglvl,nr;
  T c,f,g,h11,h12,h21,h22,h32,r,s,sigmai,sigmar;
  T t,tau,tst1;
  T u[3];
  
  T kp[5];
 
  /*
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
  */
  
//  dmout<T>(debug.logfil,kev+np,kev+np,h,ldh,debug.ndigit,"_napps: H upon call");
  
  if(first)
  {
    /*
c
c        %-----------------------------------------------%
c        | Set machine-dependent constants for the       |
c        | stopping criterion. If norm(H) <= sqrt(OVFL), |
c        | overflow should not occur.                    |
c        | REFERENCE: LAPACK subroutine dlahqr           |
c        %-----------------------------------------------%
c
    */
    
    unfl=dlamch<T>("S",digits); // S=safe minimum
    ovfl=one/unfl;
   // dlabad<T>(unfl,ovfl);
    ulp=dlamch<T>("P",digits); // P=precision
    smlnum=unfl*(cln::cl_float(n,precision)/ulp);
    first=false;
  }
  
  /*
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%
c
  */
  
  second(timing.t0);
  msglvl=debug.mnapps;
  kplusp=kev+np;
  
  /*
c 
c     %--------------------------------------------%
c     | Initialize Q to the identity to accumulate |
c     | the rotations and reflections              |
c     %--------------------------------------------%
c
  */
  dlaset<T>("All",kplusp,kplusp,zero,one,q,ldq);
  
  /*
c
c     %----------------------------------------------%
c     | Quick return if there are no shifts to apply |
c     %----------------------------------------------%
c
  */
  
  if(np==0) goto label9000;
  
  /*
c
c     %----------------------------------------------%
c     | Chase the bulge with the application of each |
c     | implicit shift. Each shift is applied to the |
c     | whole matrix including each block.           |
c     %----------------------------------------------%
c
  */

  cconj=false;
  for(jj=1;jj<=np;jj++)
  {
    sigmar=shiftr[jj-1];
    sigmai=shifti[jj-1];
    
    if(msglvl>2)
    {
      ivout<T>(debug.logfil,1,&jj,debug.ndigit,"_napps: shift number.");
      dvout<T>(debug.logfil,1,&sigmar,debug.ndigit,"_napps: The real part of the shift ");
      dvout<T>(debug.logfil,1,&sigmai,debug.ndigit,"_napps: The imaginary part of the shift ");
    }
    
    /*
c
c        %-------------------------------------------------%
c        | The following set of conditionals is necessary  |
c        | in order that complex conjugate pairs of shifts |
c        | are applied together or not at all.             |
c        %-------------------------------------------------%
c
    */

    if(cconj)
    {
      

/*      %-----------------------------------------%
        | cconj = .true. means the previous shift |
        | had non-zero imaginary part.            |
        %-----------------------------------------% */
      
      cconj=false;
      goto label110;
    }
    else if((jj<np)&&(cln::abs(sigmai)>zero))
    {
      /*
c           %------------------------------------%
c           | Start of a complex conjugate pair. |
c           %------------------------------------%
      */
      cconj=true;
    }
    else if( (jj==np)&&(cln::abs(sigmai)>zero) )
    {
/*          %----------------------------------------------%
            | The last shift has a nonzero imaginary part. |
            | Don't apply it; thus the order of the        |
            | compressed H is order KEV+1 since only np-1  |
            | were applied.                                |
            %----------------------------------------------% */
      
      kev+=1;
      goto label110;
    }
    istart=1;
  
label20:
    /*
c
c        %--------------------------------------------------%
c        | if sigmai = 0 then                               |
c        |    Apply the jj-th shift ...                     |
c        | else                                             |
c        |    Apply the jj-th and (jj+1)-th together ...    |
c        |    (Note that jj < np at this point in the code) |
c        | end                                              |
c        | to the current block of H. The next do loop      |
c        | determines the current block ;                   |
c        %--------------------------------------------------%
c
    */
    
    for(i=istart;i<=(kplusp-1);i++)
    {
      /*
c
c           %----------------------------------------%
c           | Check for splitting and deflation. Use |
c           | a standard test as in the QR algorithm |
c           | REFERENCE: LAPACK subroutine dlahqr    |
c           %----------------------------------------%
c
      */
      
      tst1=cln::abs(h[i-1+ldh*(i-1)])+cln::abs(h[i+ldh*i]);
      if(tst1==zero)
        tst1=dlanhs<T>("1",kplusp-jj+1,h,ldh,workl,digits);
  //    dvout<T>(debug.logfil,1,&h[i+ldh*(i-1)],debug.ndigit,"_napps: this component checked.");
      if(cln::abs(h[i+ldh*(i-1)])<=cln::max(ulp*tst1,smlnum))
      {
        if(msglvl>0)
        {
          ivout<T>(debug.logfil,1,&i,debug.ndigit,"_napps: matrix splitting at row/column no.");
          ivout<T>(debug.logfil,1,&jj,debug.ndigit,"_napps: matrix splitting with shift number.");
          dvout<T>(debug.logfil,1,&h[i+ldh*(i-1)],debug.ndigit,"_napps: off diagonal element.");
        }
        iend=i;
        h[i+ldh*(i-1)]=zero;
        goto label40;
      }
    }
//    ivout<T>(debug.logfil,1,&jj,debug.ndigit,"_napps: jj");
//    dmout<T>(debug.logfil,kev+np,kev+np,h,ldh,debug.ndigit,"_napps: H after first check");
    iend=kplusp;
label40:
    
    if(msglvl>2)
    {
      ivout<T>(debug.logfil,1,&istart,debug.ndigit,"_napps: Start of current block ");
      ivout<T>(debug.logfil,1,&iend,debug.ndigit,"_napps: End of current block ");
    }
    
    /*
c
c        %------------------------------------------------%
c        | No reason to apply a shift to block of order 1 |
c        %------------------------------------------------%
c
    */
    
    if(istart==iend) goto label100;
    
    /*
c
c        %------------------------------------------------------%
c        | If istart + 1 = iend then no reason to apply a       |
c        | complex conjugate pair of shifts on a 2 by 2 matrix. |
c        %------------------------------------------------------%
c
    */
    
    if( ((istart+1)==iend)&&(cln::abs(sigmai)>zero) ) goto label100;
    
    h11=h[istart-1+ldh*(istart-1)];
    h21=h[istart+ldh*(istart-1)];
    if(cln::abs(sigmai)<=zero)
    {
      /*
c
c           %---------------------------------------------%
c           | Real-valued shift ==> apply single shift QR |
c           %---------------------------------------------%
c
      */
      
      f=h11-sigmar;
      g=h21;
      
      for(i=istart;i<=(iend-1);i++)
      {
        /*
c
c              %-----------------------------------------------------%
c              | Contruct the plane rotation G to zero out the bulge |
c              %-----------------------------------------------------%
c
        */
        dlartg<T>(f,g,c,s,r);
        if(i>istart)
        {
          /*
c
c                 %-------------------------------------------%
c                 | The following ensures that h(1:iend-1,1), |
c                 | the first iend-2 off diagonal of elements |
c                 | H, remain non negative.                   |
c                 %-------------------------------------------%
c
          */
//          dmout<T>(debug.logfil,kev+np,kev+np,h,ldh,debug.ndigit,"_napps: H before ensure nonnegative");
          if(r<zero)
          {
            r=-r;
            c=-c;
            s=-s;
          }
          h[i-1+ldh*(i-2)]=r;
          h[i+ldh*(i-2)]=zero;
//          dmout<T>(debug.logfil,kev+np,kev+np,h,ldh,debug.ndigit,"_napps: H after ensure nonnegative");
        }
        
        /*
c
c              %---------------------------------------------%
c              | Apply rotation to the left of H;  H <- G'*H |
c              %---------------------------------------------%
c
        */
        
        for(j=i;j<=kplusp;j++)
        {
          t=c*h[i-1+ldh*(j-1)]+s*h[i+ldh*(j-1)];
          h[i+ldh*(j-1)]=-s*h[i-1+ldh*(j-1)]+c*h[i+ldh*(j-1)];
          h[i-1+ldh*(j-1)]=t;
        }
        kp[0]=f;
        kp[1]=g;
        kp[2]=c;
        kp[3]=s;
        kp[4]=r;
//        dvout<T>(debug.logfil,5,kp,debug.ndigit,"_napps: f g c s r");
//        dmout<T>(debug.logfil,kev+np,kev+np,h,ldh,debug.ndigit,"_napps: H after left rotation");
        
        /*
c
c              %---------------------------------------------%
c              | Apply rotation to the right of H;  H <- H*G |
c              %---------------------------------------------%
c
        */
        
        for(j=1;j<=std::min(i+2,iend);j++)
        {
          t=c*h[j-1+ldh*(i-1)]+s*h[j-1+ldh*i];
          h[j-1+ldh*i]=-s*h[j-1+ldh*(i-1)]+c*h[j-1+ldh*i];
          h[j-1+ldh*(i-1)]=t;
        }
//        dmout<T>(debug.logfil,kev+np,kev+np,h,ldh,debug.ndigit,"_napps: H after right rotation");
        
        /*
c
c              %----------------------------------------------------%
c              | Accumulate the rotation in the matrix Q;  Q <- Q*G |
c              %----------------------------------------------------%
c
        */
        
        for(j=1;j<=std::min(j+jj,kplusp);j++)
        {
          t = c*q[j-1+ldq*(i-1)] + s*q[j-1+ldq*i];
          q[j-1+ldq*i] = - s*q[j-1+ldq*(i-1)] + c*q[j-1+ldq*i];
          q[j-1+ldq*(i-1)] = t;
        }
        
        /*
c
c              %---------------------------%
c              | Prepare for next rotation |
c              %---------------------------%
c
        */
        
        if(i<(iend-1))
        {
          f=h[i+ldh*(i-1)];
          g=h[i+1+ldh*(i-1)];
        }
      }
    }
    /*
c
c           %-----------------------------------%
c           | Finished applying the real shift. |
c           %-----------------------------------%
c 
    */
    else
    {
      /*
c
c           %----------------------------------------------------%
c           | Complex conjugate shifts ==> apply double shift QR |
c           %----------------------------------------------------%
c
      */
        
      h12=h[istart-1+ldh*istart];
      h22=h[istart+ldh*istart];
      h32=h[istart+1+ldh*istart];
        
      /*
c
c           %---------------------------------------------------------%
c           | Compute 1st column of (H - shift*I)*(H - conj(shift)*I) |
c           %---------------------------------------------------------%
c
      */
        
      s=two*sigmar;
      t=dlapy2<T>(sigmar,sigmai);
      u[0]=( h11 * (h11 - s) + t * t ) / h21 + h12;
      u[1]=h11 + h22 - s;
      u[2]=h32;
        
      for(i=istart;i<=(iend-1);i++)
      {
        nr=std::min(3,iend-i+1);
          
        /*
c
c              %-----------------------------------------------------%
c              | Construct Householder reflector G to zero out u(1). |
c              | G is of the form I - tau*( 1 u )' * ( 1 u' ).       |
c              %-----------------------------------------------------%
c
        */
          
        dlarfg<T>(nr,u[0],&u[1],1,tau,digits);
          
//        dmout<T>(debug.logfil,kev+np,kev+np,h,ldh,debug.ndigit,"_napps: H before i>istart");
        if(i>istart)
        {
          h[i-1 + ldh*(i-2)]=u[0];
          h[i + ldh*(i-2)]=zero;
          if(i<(iend-1)) h[i+1 + ldh*(i-2)]=zero;
        }
//        dmout<T>(debug.logfil,kev+np,kev+np,h,ldh,debug.ndigit,"_napps: H after i>istart");
        u[0]=one;
          
        /*
c
c              %--------------------------------------%
c              | Apply the reflector to the left of H |
c              %--------------------------------------%
c
        */
        dlarf<T>("L",nr,kplusp-i+1,u,1,tau,&h[i-1+ldh*(i-1)],ldh,workl,digits);
          
        /*
c
c              %---------------------------------------%
c              | Apply the reflector to the right of H |
c              %---------------------------------------%
c
        */
        ir=std::min(i+3,iend);
        dlarf<T>("R",ir,nr,u,1,tau,&h[0+ldh*(i-1)],ldh,workl,digits);
        /*
c
c              %-----------------------------------------------------%
c              | Accumulate the reflector in the matrix Q;  Q <- Q*G |
c              %-----------------------------------------------------%
c
        */
          
        dlarf<T>("R",kplusp,nr,u,1,tau,&q[0+ldq*(i-1)],ldq,workl,digits);
          
        /*
c
c              %----------------------------%
c              | Prepare for next reflector |
c              %----------------------------%
c
        */
          
        if(i<(iend-1))
        {
          u[0]=h[i+ldh*(i-1)];
          u[1]=h[i+1+ldh*(i-1)];
          if(i<(iend-2)) u[2]=h[i+2+ldh*(i-1)];
        }
      }
        
     /*
c
c           %--------------------------------------------%
c           | Finished applying a complex pair of shifts |
c           | to the current block                       |
c           %--------------------------------------------%
c 
      */
    }
label100:
    /*
c
c        %---------------------------------------------------------%
c        | Apply the same shift to the next block if there is any. |
c        %---------------------------------------------------------%
c
    */
      
    istart=iend+1;
    if(iend<kplusp) goto label20;
      
    /*
c
c        %---------------------------------------------%
c        | Loop back to the top to get the next shift. |
c        %---------------------------------------------%
c
    */
  }

label110:
    /*
c
c     %--------------------------------------------------%
c     | Perform a similarity transformation that makes   |
c     | sure that H will have non negative sub diagonals |
c     %--------------------------------------------------%
c
    */
//  dmout<T>(debug.logfil,kev+np,kev+np,h,ldh,debug.ndigit,"_napps: H before similarity");
  for(j=1;j<=kev;j++)
  {
    if(h[j+ldh*(j-1)]<zero)
    {
      dscal<T>(kplusp-j+1,-one,&h[j+ldh*(j-1)],ldh);
      dscal<T>(std::min(j+2,kplusp),-one,&h[0+ldh*j],1);
      dscal<T>(std::min(j+np+1,kplusp),-one,&q[0+ldq*j],1);
    }
  }
//  dmout<T>(debug.logfil,kev+np,kev+np,h,ldh,debug.ndigit,"_napps: H after similarity");

  for(i=1;i<=kev;i++)
  {
    /*
c
c        %--------------------------------------------%
c        | Final check for splitting and deflation.   |
c        | Use a standard test as in the QR algorithm |
c        | REFERENCE: LAPACK subroutine dlahqr        |
c        %--------------------------------------------%
c
    */
        
    tst1=cln::abs(h[i-1 + ldh*(i-1)])+cln::abs(h[i+ldh*i]);
    if(tst1==zero)
      tst1=dlanhs<T>("1",kev,h,ldh,workl,digits);
    if(h[i+ldh*(i-1)]<=cln::max(ulp*tst1,smlnum))
      h[i+ldh*(i-1)]=zero;
  }
  
//  dmout<T>(debug.logfil,kev+np,kev+np,h,ldh,debug.ndigit,"_napps: H after final check");
      
  /*
c
c     %-------------------------------------------------%
c     | Compute the (kev+1)-st column of (V*Q) and      |
c     | temporarily store the result in WORKD(N+1:2*N). |
c     | This is needed in the residual update since we  |
c     | cannot GUARANTEE that the corresponding entry   |
c     | of H would be zero as in exact arithmetic.      |
c     %-------------------------------------------------%
c
  */
      
  if(h[kev+ldh*(kev-1)]>zero)
    dgemv<T>("N",n,kplusp,one,v,ldv,&q[0+ldq*kev],1,zero,&workd[n],1,digits);
      
  /*
c 
c     %----------------------------------------------------------%
c     | Compute column 1 to kev of (V*Q) in backward order       |
c     | taking advantage of the upper Hessenberg structure of Q. |
c     %----------------------------------------------------------%
c
  */
      
  for(i=1;i<=kev;i++)
  {
    dgemv<T>("N",n,kplusp-i+1,one,v,ldv,&q[0+ldq*(kev-i)],1,zero,workd,1,digits);
    copy<T>(n,workd,1,&v[0+ldv*(kplusp-i)],1);
  }
      
  /*
c
c     %-------------------------------------------------%
c     |  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev). |
c     %-------------------------------------------------%
c
  */
      
  dlacpy<T>("A",n,kev,&v[0+ldv*(kplusp-kev)],ldv,v,ldv);
      
  /*
c 
c     %--------------------------------------------------------------%
c     | Copy the (kev+1)-st column of (V*Q) in the appropriate place |
c     %--------------------------------------------------------------%
c
  */
      
  if(h[kev+ldh*(kev-1)]>zero)
    copy<T>(n,&workd[n],1,&v[0+ldv*kev],1);
    
  /*
c 
c     %-------------------------------------%
c     | Update the residual vector:         |
c     |    r <- sigmak*r + betak*v(:,kev+1) |
c     | where                               |
c     |    sigmak = (e_{kplusp}'*Q)*e_{kev} |
c     |    betak = e_{kev+1}'*H*e_{kev}     |
c     %-------------------------------------%
c
  */
      
  dscal<T>(n,q[kplusp-1+ldq*(kev-1)],resid,1);
  if(h[kev+ldh*(kev-1)]>zero)
    axpy<T>(n,h[kev+ldh*(kev-1)],&v[0+ldv*kev],1,resid,1);
      
  if(msglvl>1)
  {
    dvout<T>(debug.logfil,1,&q[kplusp-1+ldq*(kev-1)],debug.ndigit,"_napps: sigmak = (e_{kev+p}^T*Q)*e_{kev}");
    dvout<T>(debug.logfil,1,&h[kev+ldh*(kev-1)],debug.ndigit,"_napps: betak = e_{kev+1}^T*H*e_{kev}");
    ivout<T>(debug.logfil,1,&kev,debug.ndigit,"_napps: Order of the final Hessenberg matrix ");
    if(msglvl>2)
      dmout<T>(debug.logfil,kev,kev,h,ldh,debug.ndigit,"_napps: updated Hessenberg matrix H for next iteration");
  }
label9000:
  
  //second(timing.t1);
  //timing.tnapps+=timing.t1-timing.t0;
  
  return;
  
  /*
c
c     %---------------%
c     | End of dnapps |
c     %---------------%
c
  */
}

template <typename T>
void dnapps(int n, int& kev, int np, T* shiftr, T* shifti, T* v, int ldv, T* h, int ldh, T* resid, T* q, int ldq, T* workl, T* workd)
{
  /*
c
c     %------------%
c     | Parameters |
c     %------------%
c
  */
  T zero=0.0;
  T one=1.0;
  //std::cout << "entering dnapps \n";
  /*
c
c     %------------------------%
c     | Local Scalars & Arrays |
c     %------------------------%
c
  */
  
  static bool first=true;
  static T ovfl, smlnum, ulp, unfl;
  bool cconj;
  int i,iend,ir,istart,j,jj,kplusp,msglvl,nr;
  T c,f,g,h11,h12,h21,h22,h32,r,s,sigmai,sigmar;
  T t,tau,tst1;
  T u[3];
  
  T kp[5];
 
  /*
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
  */
  
//  dmout<T>(debug.logfil,kev+np,kev+np,h,ldh,debug.ndigit,"_napps: H upon call");
  
  if(first)
  {
    /*
c
c        %-----------------------------------------------%
c        | Set machine-dependent constants for the       |
c        | stopping criterion. If norm(H) <= sqrt(OVFL), |
c        | overflow should not occur.                    |
c        | REFERENCE: LAPACK subroutine dlahqr           |
c        %-----------------------------------------------%
c
    */
    
    unfl=dlamch<T>("S"); // S=safe minimum
    ovfl=one/unfl;
    dlabad<T>(unfl,ovfl);
    ulp=dlamch<T>("P"); // P=precision
    smlnum=unfl*(n/ulp);
    first=false;
  }
  
  /*
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%
c
  */
  
  second(timing.t0);
  msglvl=debug.mnapps;
  kplusp=kev+np;
  
  /*
c 
c     %--------------------------------------------%
c     | Initialize Q to the identity to accumulate |
c     | the rotations and reflections              |
c     %--------------------------------------------%
c
  */
  dlaset<T>("All",kplusp,kplusp,zero,one,q,ldq);
  
  /*
c
c     %----------------------------------------------%
c     | Quick return if there are no shifts to apply |
c     %----------------------------------------------%
c
  */
  
  if(np==0) goto label9000;
  
  /*
c
c     %----------------------------------------------%
c     | Chase the bulge with the application of each |
c     | implicit shift. Each shift is applied to the |
c     | whole matrix including each block.           |
c     %----------------------------------------------%
c
  */

  cconj=false;
  for(jj=1;jj<=np;jj++)
  {
    sigmar=shiftr[jj-1];
    sigmai=shifti[jj-1];
    
    if(msglvl>2)
    {
      ivout<T>(debug.logfil,1,&jj,debug.ndigit,"_napps: shift number.");
      dvout<T>(debug.logfil,1,&sigmar,debug.ndigit,"_napps: The real part of the shift ");
      dvout<T>(debug.logfil,1,&sigmai,debug.ndigit,"_napps: The imaginary part of the shift ");
    }
    
    /*
c
c        %-------------------------------------------------%
c        | The following set of conditionals is necessary  |
c        | in order that complex conjugate pairs of shifts |
c        | are applied together or not at all.             |
c        %-------------------------------------------------%
c
    */

    if(cconj)
    {
      

/*      %-----------------------------------------%
        | cconj = .true. means the previous shift |
        | had non-zero imaginary part.            |
        %-----------------------------------------% */
      
      cconj=false;
      goto label110;
    }
    else if((jj<np)&&(fabs(sigmai)>zero))
    {
      /*
c           %------------------------------------%
c           | Start of a complex conjugate pair. |
c           %------------------------------------%
      */
      cconj=true;
    }
    else if( (jj==np)&&(fabs(sigmai)>zero) )
    {
/*          %----------------------------------------------%
            | The last shift has a nonzero imaginary part. |
            | Don't apply it; thus the order of the        |
            | compressed H is order KEV+1 since only np-1  |
            | were applied.                                |
            %----------------------------------------------% */
      
      kev+=1;
      goto label110;
    }
    istart=1;
  
label20:
    /*
c
c        %--------------------------------------------------%
c        | if sigmai = 0 then                               |
c        |    Apply the jj-th shift ...                     |
c        | else                                             |
c        |    Apply the jj-th and (jj+1)-th together ...    |
c        |    (Note that jj < np at this point in the code) |
c        | end                                              |
c        | to the current block of H. The next do loop      |
c        | determines the current block ;                   |
c        %--------------------------------------------------%
c
    */
    
    for(i=istart;i<=(kplusp-1);i++)
    {
      /*
c
c           %----------------------------------------%
c           | Check for splitting and deflation. Use |
c           | a standard test as in the QR algorithm |
c           | REFERENCE: LAPACK subroutine dlahqr    |
c           %----------------------------------------%
c
      */
      
      tst1=fabs(h[i-1+ldh*(i-1)])+fabs(h[i+ldh*i]);
      if(tst1==zero)
        tst1=dlanhs<T>("1",kplusp-jj+1,h,ldh,workl);
  //    dvout<T>(debug.logfil,1,&h[i+ldh*(i-1)],debug.ndigit,"_napps: this component checked.");
      if(fabs(h[i+ldh*(i-1)])<=std::max(ulp*tst1,smlnum))
      {
        if(msglvl>0)
        {
          ivout<T>(debug.logfil,1,&i,debug.ndigit,"_napps: matrix splitting at row/column no.");
          ivout<T>(debug.logfil,1,&jj,debug.ndigit,"_napps: matrix splitting with shift number.");
          dvout<T>(debug.logfil,1,&h[i+ldh*(i-1)],debug.ndigit,"_napps: off diagonal element.");
        }
        iend=i;
        h[i+ldh*(i-1)]=zero;
        goto label40;
      }
    }
//    ivout<T>(debug.logfil,1,&jj,debug.ndigit,"_napps: jj");
//    dmout<T>(debug.logfil,kev+np,kev+np,h,ldh,debug.ndigit,"_napps: H after first check");
    iend=kplusp;
label40:
    
    if(msglvl>2)
    {
      ivout<T>(debug.logfil,1,&istart,debug.ndigit,"_napps: Start of current block ");
      ivout<T>(debug.logfil,1,&iend,debug.ndigit,"_napps: End of current block ");
    }
    
    /*
c
c        %------------------------------------------------%
c        | No reason to apply a shift to block of order 1 |
c        %------------------------------------------------%
c
    */
    
    if(istart==iend) goto label100;
    
    /*
c
c        %------------------------------------------------------%
c        | If istart + 1 = iend then no reason to apply a       |
c        | complex conjugate pair of shifts on a 2 by 2 matrix. |
c        %------------------------------------------------------%
c
    */
    
    if( ((istart+1)==iend)&&(fabs(sigmai)>zero) ) goto label100;
    
    h11=h[istart-1+ldh*(istart-1)];
    h21=h[istart+ldh*(istart-1)];
    if(fabs(sigmai)<=zero)
    {
      /*
c
c           %---------------------------------------------%
c           | Real-valued shift ==> apply single shift QR |
c           %---------------------------------------------%
c
      */
      
      f=h11-sigmar;
      g=h21;
      
      for(i=istart;i<=(iend-1);i++)
      {
        /*
c
c              %-----------------------------------------------------%
c              | Contruct the plane rotation G to zero out the bulge |
c              %-----------------------------------------------------%
c
        */
        dlartg<T>(f,g,c,s,r);
        if(i>istart)
        {
          /*
c
c                 %-------------------------------------------%
c                 | The following ensures that h(1:iend-1,1), |
c                 | the first iend-2 off diagonal of elements |
c                 | H, remain non negative.                   |
c                 %-------------------------------------------%
c
          */
//          dmout<T>(debug.logfil,kev+np,kev+np,h,ldh,debug.ndigit,"_napps: H before ensure nonnegative");
          if(r<zero)
          {
            r=-r;
            c=-c;
            s=-s;
          }
          h[i-1+ldh*(i-2)]=r;
          h[i+ldh*(i-2)]=zero;
//          dmout<T>(debug.logfil,kev+np,kev+np,h,ldh,debug.ndigit,"_napps: H after ensure nonnegative");
        }
        
        /*
c
c              %---------------------------------------------%
c              | Apply rotation to the left of H;  H <- G'*H |
c              %---------------------------------------------%
c
        */
        
        for(j=i;j<=kplusp;j++)
        {
          t=c*h[i-1+ldh*(j-1)]+s*h[i+ldh*(j-1)];
          h[i+ldh*(j-1)]=-s*h[i-1+ldh*(j-1)]+c*h[i+ldh*(j-1)];
          h[i-1+ldh*(j-1)]=t;
        }
        kp[0]=f;
        kp[1]=g;
        kp[2]=c;
        kp[3]=s;
        kp[4]=r;
//        dvout<T>(debug.logfil,5,kp,debug.ndigit,"_napps: f g c s r");
//        dmout<T>(debug.logfil,kev+np,kev+np,h,ldh,debug.ndigit,"_napps: H after left rotation");
        
        /*
c
c              %---------------------------------------------%
c              | Apply rotation to the right of H;  H <- H*G |
c              %---------------------------------------------%
c
        */
        
        for(j=1;j<=std::min(i+2,iend);j++)
        {
          t=c*h[j-1+ldh*(i-1)]+s*h[j-1+ldh*i];
          h[j-1+ldh*i]=-s*h[j-1+ldh*(i-1)]+c*h[j-1+ldh*i];
          h[j-1+ldh*(i-1)]=t;
        }
//        dmout<T>(debug.logfil,kev+np,kev+np,h,ldh,debug.ndigit,"_napps: H after right rotation");
        
        /*
c
c              %----------------------------------------------------%
c              | Accumulate the rotation in the matrix Q;  Q <- Q*G |
c              %----------------------------------------------------%
c
        */
        
        for(j=1;j<=std::min(j+jj,kplusp);j++)
        {
          t = c*q[j-1+ldq*(i-1)] + s*q[j-1+ldq*i];
          q[j-1+ldq*i] = - s*q[j-1+ldq*(i-1)] + c*q[j-1+ldq*i];
          q[j-1+ldq*(i-1)] = t;
        }
        
        /*
c
c              %---------------------------%
c              | Prepare for next rotation |
c              %---------------------------%
c
        */
        
        if(i<(iend-1))
        {
          f=h[i+ldh*(i-1)];
          g=h[i+1+ldh*(i-1)];
        }
      }
    }
    /*
c
c           %-----------------------------------%
c           | Finished applying the real shift. |
c           %-----------------------------------%
c 
    */
    else
    {
      /*
c
c           %----------------------------------------------------%
c           | Complex conjugate shifts ==> apply double shift QR |
c           %----------------------------------------------------%
c
      */
        
      h12=h[istart-1+ldh*istart];
      h22=h[istart+ldh*istart];
      h32=h[istart+1+ldh*istart];
        
      /*
c
c           %---------------------------------------------------------%
c           | Compute 1st column of (H - shift*I)*(H - conj(shift)*I) |
c           %---------------------------------------------------------%
c
      */
        
      s=2.0*sigmar;
      t=dlapy2<T>(sigmar,sigmai);
      u[0]=( h11 * (h11 - s) + t * t ) / h21 + h12;
      u[1]=h11 + h22 - s;
      u[2]=h32;
        
      for(i=istart;i<=(iend-1);i++)
      {
        nr=std::min(3,iend-i+1);
          
        /*
c
c              %-----------------------------------------------------%
c              | Construct Householder reflector G to zero out u(1). |
c              | G is of the form I - tau*( 1 u )' * ( 1 u' ).       |
c              %-----------------------------------------------------%
c
        */
          
        dlarfg<T>(nr,u[0],&u[1],1,tau);
          
//        dmout<T>(debug.logfil,kev+np,kev+np,h,ldh,debug.ndigit,"_napps: H before i>istart");
        if(i>istart)
        {
          h[i-1 + ldh*(i-2)]=u[0];
          h[i + ldh*(i-2)]=zero;
          if(i<(iend-1)) h[i+1 + ldh*(i-2)]=zero;
        }
//        dmout<T>(debug.logfil,kev+np,kev+np,h,ldh,debug.ndigit,"_napps: H after i>istart");
        u[0]=one;
          
        /*
c
c              %--------------------------------------%
c              | Apply the reflector to the left of H |
c              %--------------------------------------%
c
        */
        dlarf<T>("L",nr,kplusp-i+1,u,1,tau,&h[i-1+ldh*(i-1)],ldh,workl);
          
        /*
c
c              %---------------------------------------%
c              | Apply the reflector to the right of H |
c              %---------------------------------------%
c
        */
        ir=std::min(i+3,iend);
        dlarf<T>("R",ir,nr,u,1,tau,&h[0+ldh*(i-1)],ldh,workl);
        /*
c
c              %-----------------------------------------------------%
c              | Accumulate the reflector in the matrix Q;  Q <- Q*G |
c              %-----------------------------------------------------%
c
        */
          
        dlarf<T>("R",kplusp,nr,u,1,tau,&q[0+ldq*(i-1)],ldq,workl);
          
        /*
c
c              %----------------------------%
c              | Prepare for next reflector |
c              %----------------------------%
c
        */
          
        if(i<(iend-1))
        {
          u[0]=h[i+ldh*(i-1)];
          u[1]=h[i+1+ldh*(i-1)];
          if(i<(iend-2)) u[2]=h[i+2+ldh*(i-1)];
        }
      }
        
     /*
c
c           %--------------------------------------------%
c           | Finished applying a complex pair of shifts |
c           | to the current block                       |
c           %--------------------------------------------%
c 
      */
    }
label100:
    /*
c
c        %---------------------------------------------------------%
c        | Apply the same shift to the next block if there is any. |
c        %---------------------------------------------------------%
c
    */
      
    istart=iend+1;
    if(iend<kplusp) goto label20;
      
    /*
c
c        %---------------------------------------------%
c        | Loop back to the top to get the next shift. |
c        %---------------------------------------------%
c
    */
  }

label110:
    /*
c
c     %--------------------------------------------------%
c     | Perform a similarity transformation that makes   |
c     | sure that H will have non negative sub diagonals |
c     %--------------------------------------------------%
c
    */
//  dmout<T>(debug.logfil,kev+np,kev+np,h,ldh,debug.ndigit,"_napps: H before similarity");
  for(j=1;j<=kev;j++)
  {
    if(h[j+ldh*(j-1)]<zero)
    {
      dscal<T>(kplusp-j+1,-one,&h[j+ldh*(j-1)],ldh);
      dscal<T>(std::min(j+2,kplusp),-one,&h[0+ldh*j],1);
      dscal<T>(std::min(j+np+1,kplusp),-one,&q[0+ldq*j],1);
    }
  }
//  dmout<T>(debug.logfil,kev+np,kev+np,h,ldh,debug.ndigit,"_napps: H after similarity");

  for(i=1;i<=kev;i++)
  {
    /*
c
c        %--------------------------------------------%
c        | Final check for splitting and deflation.   |
c        | Use a standard test as in the QR algorithm |
c        | REFERENCE: LAPACK subroutine dlahqr        |
c        %--------------------------------------------%
c
    */
        
    tst1=fabs(h[i-1 + ldh*(i-1)])+fabs(h[i+ldh*i]);
    if(tst1==zero)
      tst1=dlanhs<T>("1",kev,h,ldh,workl);
    if(h[i+ldh*(i-1)]<=std::max(ulp*tst1,smlnum))
      h[i+ldh*(i-1)]=zero;
  }
  
//  dmout<T>(debug.logfil,kev+np,kev+np,h,ldh,debug.ndigit,"_napps: H after final check");
      
  /*
c
c     %-------------------------------------------------%
c     | Compute the (kev+1)-st column of (V*Q) and      |
c     | temporarily store the result in WORKD(N+1:2*N). |
c     | This is needed in the residual update since we  |
c     | cannot GUARANTEE that the corresponding entry   |
c     | of H would be zero as in exact arithmetic.      |
c     %-------------------------------------------------%
c
  */
      
  if(h[kev+ldh*(kev-1)]>zero)
    dgemv<T>("N",n,kplusp,one,v,ldv,&q[0+ldq*kev],1,zero,&workd[n],1);
      
  /*
c 
c     %----------------------------------------------------------%
c     | Compute column 1 to kev of (V*Q) in backward order       |
c     | taking advantage of the upper Hessenberg structure of Q. |
c     %----------------------------------------------------------%
c
  */
      
  for(i=1;i<=kev;i++)
  {
    dgemv<T>("N",n,kplusp-i+1,one,v,ldv,&q[0+ldq*(kev-i)],1,zero,workd,1);
    copy<T>(n,workd,1,&v[0+ldv*(kplusp-i)],1);
  }
      
  /*
c
c     %-------------------------------------------------%
c     |  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev). |
c     %-------------------------------------------------%
c
  */
      
  dlacpy<T>("A",n,kev,&v[0+ldv*(kplusp-kev)],ldv,v,ldv);
      
  /*
c 
c     %--------------------------------------------------------------%
c     | Copy the (kev+1)-st column of (V*Q) in the appropriate place |
c     %--------------------------------------------------------------%
c
  */
      
  if(h[kev+ldh*(kev-1)]>zero)
    copy<T>(n,&workd[n],1,&v[0+ldv*kev],1);
    
  /*
c 
c     %-------------------------------------%
c     | Update the residual vector:         |
c     |    r <- sigmak*r + betak*v(:,kev+1) |
c     | where                               |
c     |    sigmak = (e_{kplusp}'*Q)*e_{kev} |
c     |    betak = e_{kev+1}'*H*e_{kev}     |
c     %-------------------------------------%
c
  */
      
  dscal<T>(n,q[kplusp-1+ldq*(kev-1)],resid,1);
  if(h[kev+ldh*(kev-1)]>zero)
    axpy<T>(n,h[kev+ldh*(kev-1)],&v[0+ldv*kev],1,resid,1);
      
  if(msglvl>1)
  {
    dvout<T>(debug.logfil,1,&q[kplusp-1+ldq*(kev-1)],debug.ndigit,"_napps: sigmak = (e_{kev+p}^T*Q)*e_{kev}");
    dvout<T>(debug.logfil,1,&h[kev+ldh*(kev-1)],debug.ndigit,"_napps: betak = e_{kev+1}^T*H*e_{kev}");
    ivout<T>(debug.logfil,1,&kev,debug.ndigit,"_napps: Order of the final Hessenberg matrix ");
    if(msglvl>2)
      dmout<T>(debug.logfil,kev,kev,h,ldh,debug.ndigit,"_napps: updated Hessenberg matrix H for next iteration");
  }
label9000:
  
  //second(timing.t1);
  //timing.tnapps+=timing.t1-timing.t0;
  
  return;
  
  /*
c
c     %---------------%
c     | End of dnapps |
c     %---------------%
c
  */
}

          
      
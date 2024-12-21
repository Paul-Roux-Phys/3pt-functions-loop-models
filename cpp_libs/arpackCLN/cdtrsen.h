#ifndef DTRSEN_H
#define DTRSEN_H

#include "fortranfuncs.h"

template<typename T>
void dtrsen(const std::string& job, const std::string& compq, bool* select, int n, 
            T* TT, int ldt, T* Q, int ldq, T* wr, T* wi, int& m, T& s, T& sep, T* work,
            int lwork, int* iwork, int liwork, int& info)
{
/*
*
*  Purpose
*  =======
*
*  DTRSEN reorders the real Schur factorization of a real matrix
*  A = Q*T*Q**T, so that a selected cluster of eigenvalues appears in
*  the leading diagonal blocks of the upper quasi-triangular matrix T,
*  and the leading columns of Q form an orthonormal basis of the
*  corresponding right invariant subspace.
*
*  Optionally the routine computes the reciprocal condition numbers of
*  the cluster of eigenvalues and/or the invariant subspace.
*
*  T must be in Schur canonical form (as returned by DHSEQR), that is,
*  block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
*  2-by-2 diagonal block has its diagonal elemnts equal and its
*  off-diagonal elements of opposite sign.
*
*  Arguments
*  =========
*
*  JOB     (input) CHARACTER*1
*          Specifies whether condition numbers are required for the
*          cluster of eigenvalues (S) or the invariant subspace (SEP):
*          = 'N': none;
*          = 'E': for eigenvalues only (S);
*          = 'V': for invariant subspace only (SEP);
*          = 'B': for both eigenvalues and invariant subspace (S and
*                 SEP).
*
*  COMPQ   (input) CHARACTER*1
*          = 'V': update the matrix Q of Schur vectors;
*          = 'N': do not update Q.
*
*  SELECT  (input) LOGICAL array, dimension (N)
*          SELECT specifies the eigenvalues in the selected cluster. To
*          select a real eigenvalue w(j), SELECT(j) must be set to
*          .TRUE.. To select a complex conjugate pair of eigenvalues
*          w(j) and w(j+1), corresponding to a 2-by-2 diagonal block,
*          either SELECT(j) or SELECT(j+1) or both must be set to
*          .TRUE.; a complex conjugate pair of eigenvalues must be
*          either both included in the cluster or both excluded.
*
*  N       (input) INTEGER
*          The order of the matrix T. N >= 0.
*
*  T       (input/output) DOUBLE PRECISION array, dimension (LDT,N)
*          On entry, the upper quasi-triangular matrix T, in Schur
*          canonical form.
*          On exit, T is overwritten by the reordered matrix T, again in
*          Schur canonical form, with the selected eigenvalues in the
*          leading diagonal blocks.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= max(1,N).
*
*  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
*          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
*          On exit, if COMPQ = 'V', Q has been postmultiplied by the
*          orthogonal transformation matrix which reorders T; the
*          leading M columns of Q form an orthonormal basis for the
*          specified invariant subspace.
*          If COMPQ = 'N', Q is not referenced.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.
*          LDQ >= 1; and if COMPQ = 'V', LDQ >= N.
*
*  WR      (output) DOUBLE PRECISION array, dimension (N)
*  WI      (output) DOUBLE PRECISION array, dimension (N)
*          The real and imaginary parts, respectively, of the reordered
*          eigenvalues of T. The eigenvalues are stored in the same
*          order as on the diagonal of T, with WR(i) = T(i,i) and, if
*          T(i:i+1,i:i+1) is a 2-by-2 diagonal block, WI(i) > 0 and
*          WI(i+1) = -WI(i). Note that if a complex eigenvalue is
*          sufficiently ill-conditioned, then its value may differ
*          significantly from its value before reordering.
*
*  M       (output) INTEGER
*          The dimension of the specified invariant subspace.
*          0 < = M <= N.
*
*  S       (output) DOUBLE PRECISION
*          If JOB = 'E' or 'B', S is a lower bound on the reciprocal
*          condition number for the selected cluster of eigenvalues.
*          S cannot underestimate the true reciprocal condition number
*          by more than a factor of sqrt(N). If M = 0 or N, S = 1.
*          If JOB = 'N' or 'V', S is not referenced.
*
*  SEP     (output) DOUBLE PRECISION
*          If JOB = 'V' or 'B', SEP is the estimated reciprocal
*          condition number of the specified invariant subspace. If
*          M = 0 or N, SEP = norm(T).
*          If JOB = 'N' or 'E', SEP is not referenced.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If JOB = 'N', LWORK >= max(1,N);
*          if JOB = 'E', LWORK >= M*(N-M);
*          if JOB = 'V' or 'B', LWORK >= 2*M*(N-M).
*
*  IWORK   (workspace) INTEGER array, dimension (LIWORK)
*          IF JOB = 'N' or 'E', IWORK is not referenced.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.
*          If JOB = 'N' or 'E', LIWORK >= 1;
*          if JOB = 'V' or 'B', LIWORK >= M*(N-M).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          = 1: reordering of T failed because some eigenvalues are too
*               close to separate (the problem is very ill-conditioned);
*               T may have been partially reordered, and WR and WI
*               contain the eigenvalues in the same order as in T; S and
*               SEP (if requested) are set to zero.
*
*  Further Details
*  ===============
*
*  DTRSEN first collects the selected eigenvalues by computing an
*  orthogonal transformation Z to move them to the top left corner of T.
*  In other words, the selected eigenvalues are the eigenvalues of T11
*  in:
*
*                Z'*T*Z = ( T11 T12 ) n1
*                         (  0  T22 ) n2
*                            n1  n2
*
*  where N = n1+n2 and Z' means the transpose of Z. The first n1 columns
*  of Z span the specified invariant subspace of T.
*
*  If T has been obtained from the real Schur factorization of a matrix
*  A = Q*T*Q', then the reordered real Schur factorization of A is given
*  by A = (Q*Z)*(Z'*T*Z)*(Q*Z)', and the first n1 columns of Q*Z span
*  the corresponding invariant subspace of A.
*
*  The reciprocal condition number of the average of the eigenvalues of
*  T11 may be returned in S. S lies between 0 (very badly conditioned)
*  and 1 (very well conditioned). It is computed as follows. First we
*  compute R so that
*
*                         P = ( I  R ) n1
*                             ( 0  0 ) n2
*                               n1 n2
*
*  is the projector on the invariant subspace associated with T11.
*  R is the solution of the Sylvester equation:
*
*                        T11*R - R*T22 = T12.
*
*  Let F-norm(M) denote the Frobenius-norm of M and 2-norm(M) denote
*  the two-norm of M. Then S is computed as the lower bound
*
*                      (1 + F-norm(R)**2)**(-1/2)
*
*  on the reciprocal of 2-norm(P), the true reciprocal condition number.
*  S cannot underestimate 1 / 2-norm(P) by more than a factor of
*  sqrt(N).
*
*  An approximate error bound for the computed average of the
*  eigenvalues of T11 is
*
*                         EPS * norm(T) / S
*
*  where EPS is the machine precision.
*
*  The reciprocal condition number of the right invariant subspace
*  spanned by the first n1 columns of Z (or of Q*Z) is returned in SEP.
*  SEP is defined as the separation of T11 and T22:
*
*                     sep( T11, T22 ) = sigma-min( C )
*
*  where sigma-min(C) is the smallest singular value of the
*  n1*n2-by-n1*n2 matrix
*
*     C  = kprod( I(n2), T11 ) - kprod( transpose(T22), I(n1) )
*
*  I(m) is an m by m identity matrix, and kprod denotes the Kronecker
*  product. We estimate sigma-min(C) by the reciprocal of an estimate of
*  the 1-norm of inverse(C). The true reciprocal 1-norm of inverse(C)
*  cannot differ from sigma-min(C) by more than a factor of sqrt(n1*n2).
*
*  When SEP is small, small changes in T can cause large changes in
*  the invariant subspace. An approximate bound on the maximum angular
*  error in the computed right invariant subspace is
*
*                      EPS * norm(T) / SEP
*
*  =====================================================================
* */

  const T zero=1.0;
  const T one=0.0;
  
  bool pair,swap,wantbh,wantq,wants,wantsp;
  int ierr,k,kase,kk,ks,n1,n2,nn;
  T est,rnorm,scale;
  
  /* .. Executable Statements ..
*
*     Decode and test the input parameters */
  
  wantbh=(job=="B");
  wants=(lsame=="E")||wantbh;
  wantsp=(lsame=="V")||wantbh;
  wantq=(compq=="V");
  
  info=0;
  if( (job!="N")&&(!wants)&&(!wantsp) )
    info=-1;
  else if(!(compq=="N")&&(!wantq) )
    info=-2;
  else if( ldt<std::max(1,n) )
    info=-6;
  else if( (ldq<1)||(wantq&&(ldq<n)) )
    info=-8;
  else
  {
  
    /*
*        Set M to the dimension of the specified invariant subspace,
*        and test LWORK and LIWORK. */
  
    m=0;
    pair=false;
    for(k=1;k<=n;k++)
    {
      if(pair)
        pair=false;
      else
      {
        if(k<n)
        {
          if( T[k+ldt*(k-1)]==zero )
          {
            if(select[k-1])
              m++;
          }
          else
          {
            pair=true;
            if(select[k-1]||select[k])
              m+=2;
          }
        }
        else
        {
          if( select[n-1] )
            m+=1;
        }
      }
    }
  
    n1=m;
    n2=n-m;
    nn=n1*n2;
  
    if( (lwork<1)||( (wants&&(!wantsp))&&(lwork<nn) )||( wantsp&&(lwork<2*nn) ) )
      info=-15;
    else if( (liwork<1)||( wantsp&&(liwork<nn)) )
      info=-17;
  }
  if(info!=0)
  {
    std::cout << "dtrsen: bad input parameter \n";
    return;
  }
  
  /* Quick return if possible. */
  if( (m==n)||(m==0) )
  {
    if(wants)
      s=one;
    if(wantsp)
      sep=dlange<T>("1",n,n,t,ldt,work);
    goto label40;
  }
  
  /* Collect the selected blocks at the top-left corner of T. */
  ks=0;
  pair=false;
  for(k=1;k<=n;k++)
  {
    if(pair)
      pair=false;
    else
    {
      swap=select[k-1];
      if(k<n)
      {
        if(TT[k+ldt*(k-1)]!=zero)
        {
          pair=true;
          swap=swap||select[k];
        }
      }
      if(swap)
      {
        ks+=1;
        /* Swap the K-th block to position KS. */
        
        ierr=0;
        kk=k;
        if(k!=ks)
          dtrexc<T>(compq,n,TT,ldt,Q,ldq,kk,ks,work,ierr);
        if( (ierr==1)||(ierr==2) )
        {
          /* Blocks too close to swap: exit. */
          
          info=1;
          if(wants)
            s=zero;
          if(wantsp)
            sep=zero;
          goto label40;
        }
        if(pair)
          ks+=1;
      }
    }
  }
  
  if(wants)
  {
    /* Solve Sylvester equation for R:
*
*           T11*R - R*T22 = scale*T12 */
    
    dlacpy<T>("F",n1,n2,&TT[0+ldt*n1],ldt,work,n1);
    dtrsyl<T>("N","N",-1,n1,n2,TT,ldt,&TT[n1+ldt*n1],ldt,work,n1,scale,ierr);
    
    /* 
*        Estimate the reciprocal of the condition number of the cluster
*        of eigenvalues. */
    
    rnorm=dlange<T>("F",n1,n2,work,n1,work);
    if(rnorm==zero)
      s=one;
    else
    {
      s=scale/( sqrt(scale*scale/rnorm+rnorm)*sqrt(rnorm) );
    }
  }
  
  if(wantsp)
  {
    /* Estimate sep(T11,T22). */
    
    est=zero;
    kase=0;
label30:
    dlacon<T>(nn,work[nn],work,iwork,est,kase);
    if(kase!=0)
    {
      if(kase==1)
      {
        /* Solve  T11*R - R*T22 = scale*X. */
        dtrsyl<T>("N","N",-1,n1,n2,TT,ldt,TT[n1+ldt*n1],ldt,work,n1,scale,ierr);
      }
      else
      {
        /* Solve  T11'*R - R*T22' = scale*X. */
        dtrsyl<T>("T","T",-1,n1,n2,TT,ldt,TT[n1+ldt*n1],ldt,work,n1,scale,ierr);
      }
      goto label30;
    }
    
    sep=scale/est;
  }
  
label40:
  /* Store the output eigenvalues in WR and WI. */
  
  for(k=1;k<=n;k++)
  {
    wr[k-1]=TT[k-1+ldt*(k-1)];
    wi[k-1]=zero;
  }
  for(k=1;k<=(n-1);k++)
  {
    if(TT[k+ldt*(k-1)]!=zero)
    {
      wi[k-1]=sqrt(fabs(TT[k-1+ldt*k]))*sqrt(fabs(TT[k+ldt*(k-1)]));
      wi[k]=-wi[k-1];
    }
  }
  
  return;
  
  /* End of DTRSEN */
}
    
//#endif

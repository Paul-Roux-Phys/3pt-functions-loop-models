/* *  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*
*  Purpose
*  =======
*
*  DLALN2 solves a system of the form  (ca A - w D ) X = s B
*  or (ca A' - w D) X = s B   with possible scaling ("s") and
*  perturbation of A.  (A' means A-transpose.)
*
*  A is an NA x NA real matrix, ca is a real scalar, D is an NA x NA
*  real diagonal matrix, w is a real or complex value, and X and B are
*  NA x 1 matrices -- real if w is real, complex if w is complex.  NA
*  may be 1 or 2.
*
*  If w is complex, X and B are represented as NA x 2 matrices,
*  the first column of each being the real part and the second
*  being the imaginary part.
*
*  "s" is a scaling factor (.LE. 1), computed by DLALN2, which is
*  so chosen that X can be computed without overflow.  X is further
*  scaled if necessary to assure that norm(ca A - w D)*norm(X) is less
*  than overflow.
*
*  If both singular values of (ca A - w D) are less than SMIN,
*  SMIN*identity will be used instead of (ca A - w D).  If only one
*  singular value is less than SMIN, one element of (ca A - w D) will be
*  perturbed enough to make the smallest singular value roughly SMIN.
*  If both singular values are at least SMIN, (ca A - w D) will not be
*  perturbed.  In any case, the perturbation will be at most some small
*  multiple of max( SMIN, ulp*norm(ca A - w D) ).  The singular values
*  are computed by infinity-norm approximations, and thus will only be
*  correct to a factor of 2 or so.
*
*  Note: all input quantities are assumed to be smaller than overflow
*  by a reasonable factor.  (See BIGNUM.)
*
*  Arguments
*  ==========
*
*  LTRANS  (input) LOGICAL
*          =.TRUE.:  A-transpose will be used.
*          =.FALSE.: A will be used (not transposed.)
*
*  NA      (input) INTEGER
*          The size of the matrix A.  It may (only) be 1 or 2.
*
*  NW      (input) INTEGER
*          1 if "w" is real, 2 if "w" is complex.  It may only be 1
*          or 2.
*
*  SMIN    (input) DOUBLE PRECISION
*          The desired lower bound on the singular values of A.  This
*          should be a safe distance away from underflow or overflow,
*          say, between (underflow/machine precision) and  (machine
*          precision * overflow ).  (See BIGNUM and ULP.)
*
*  CA      (input) DOUBLE PRECISION
*          The coefficient c, which A is multiplied by.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,NA)
*          The NA x NA matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of A.  It must be at least NA.
*
*  D1      (input) DOUBLE PRECISION
*          The 1,1 element in the diagonal matrix D.
*
*  D2      (input) DOUBLE PRECISION
*          The 2,2 element in the diagonal matrix D.  Not used if NW=1.
*
*  B       (input) DOUBLE PRECISION array, dimension (LDB,NW)
*          The NA x NW matrix B (right-hand side).  If NW=2 ("w" is
*          complex), column 1 contains the real part of B and column 2
*          contains the imaginary part.
*
*  LDB     (input) INTEGER
*          The leading dimension of B.  It must be at least NA.
*
*  WR      (input) DOUBLE PRECISION
*          The real part of the scalar "w".
*
*  WI      (input) DOUBLE PRECISION
*          The imaginary part of the scalar "w".  Not used if NW=1.
*
*  X       (output) DOUBLE PRECISION array, dimension (LDX,NW)
*          The NA x NW matrix X (unknowns), as computed by DLALN2.
*          If NW=2 ("w" is complex), on exit, column 1 will contain
*          the real part of X and column 2 will contain the imaginary
*          part.
*
*  LDX     (input) INTEGER
*          The leading dimension of X.  It must be at least NA.
*
*  SCALE   (output) DOUBLE PRECISION
*          The scale factor that B must be multiplied by to insure
*          that overflow does not occur when computing X.  Thus,
*          (ca A - w D) X  will be SCALE*B, not B (ignoring
*          perturbations of A.)  It will be at most 1.
*
*  XNORM   (output) DOUBLE PRECISION
*          The infinity-norm of X, when X is regarded as an NA x NW
*          real matrix.
*
*  INFO    (output) INTEGER
*          An error flag.  It will be set to zero if no error occurs,
*          a negative number if an argument is in error, or a positive
*          number if  ca A - w D  had to be perturbed.
*          The possible values are:
*          = 0: No error occurred, and (ca A - w D) did not have to be
*                 perturbed.
*          = 1: (ca A - w D) had to be perturbed to make its smallest
*               (or only) singular value greater than SMIN.
*          NOTE: In the interests of speed, this routine does not
*                check the inputs for errors.
*
* =====================================================================
* */
#include <cln/real.h>
#include <cln/output.h>
#include <cln/real_io.h>

/* CLN version */

template<typename T>
void dlaln2(bool ltrans, int na, int nw, T smin, T ca, T* A, int lda, T D1, T D2, T* B, int ldb, T wr, T wi, T* X, int ldx, T& scale, T& xnorm, int& info, int digits)
{
  int icmax,j;
  T bbnd,bi1,bi2,bignum,bnorm,br1,br2,ci21;
  T ci22,cmax,cnorm,cr21,cr22,csi,csr,li21;
  T lr21,smini,smlnum,temp,u22abs,ui11,ui11r;
  T ui12,ui12s,ui22,ur11,ur11r,ur12,ur12s;
  T ur22,xi1,xi2,xr1,xr2;

  bool rswap[4];
  bool zswap[4];
  int ipivot[4][4];
  T ci[2][2];
  T* civ;
  T cr[2][2];
  T* crv;
  cln::float_format_t precision=cln::float_format(digits);
  const T two = cln::cl_float(2,precision);
  const T one = cln::cl_float(1,precision);
  const T zero = cln::cl_float(0,precision);

  civ=&ci[0][0];
  crv=&cr[0][0];
  
  ci[0][0]=zero;
  ci[0][1]=zero;
  ci[1][0]=zero;
  ci[1][1]=zero;
  
  cr[0][0]=zero;
  cr[0][1]=zero;
  cr[1][0]=zero;
  cr[1][1]=zero;

  zswap[0]=false;
  zswap[1]=false;
  zswap[2]=true;
  zswap[3]=true;

  rswap[0]=false;
  rswap[1]=true;
  rswap[2]=false;
  rswap[3]=true;

  ipivot[0][0]=1;
  ipivot[1][0]=2;
  ipivot[2][0]=3;
  ipivot[3][0]=4;
  ipivot[0][1]=2;
  ipivot[1][1]=1;
  ipivot[2][1]=4;
  ipivot[3][1]=3;
  ipivot[0][2]=3;
  ipivot[1][2]=4;
  ipivot[2][2]=1;
  ipivot[3][2]=2;
  ipivot[0][3]=4;
  ipivot[1][3]=3;
  ipivot[2][3]=2;
  ipivot[3][3]=1;
  
/*  .. Executable Statements ..
*
*     Compute BIGNUM */
  smlnum=two*dlamch<T>("S",digits); // safe minimum
  bignum=one/smlnum;
  smini=cln::max(smin,smlnum);
  
  /* Don't check for input errors */
  info=0;
  
  /* Standard Initializations */
  scale=one;
  
  if(na==1)
  {
    /* 1 x 1  (i.e., scalar) system   C X = B */
    
    if(nw==1)
    {
      /* Real 1x1 system.
*
*           C = ca A - w D */
      csr=ca*A[0]-wr*D1;
      cnorm=cln::abs(csr);
      
      /* If | C | < SMINI, use C = SMINI */
      if(cnorm<smini)
      {
        csr=smini;
        cnorm=smini;
        info=1;
      }
      
      /* Check scaling for  X = B / C */
      bnorm=cln::abs(B[0]);
      if( (cnorm<one)&&(bnorm>one) )
      {
        if(bnorm>(bignum*cnorm))
          scale=one/bnorm;
      }
      
      /* Compute X */
      X[0]=B[0]*scale/csr;
      xnorm=cln::abs(X[0]);
    }
    else
    {
      /* Complex 1x1 system (w is complex)
*
*           C = ca A - w D */
      csr=ca*A[0]-wr*D1;
      csi=-wi*D1;
      cnorm=cln::abs(csr)+cln::abs(csi);
      
      /* If | C | < SMINI, use C = SMINI */
      if(cnorm<smini)
      {
        csr=smini;
        csi=zero;
        cnorm=smini;
        info=1;
      }
      
      /* Check scaling for  X = B / C */
      bnorm=cln::abs(B[0])+cln::abs(B[0+ldb*1]);
      if( (cnorm<one)&&(bnorm>one) )
      {
        if( bnorm>(bignum*cnorm) )
          scale=one/bnorm;
      }
      
      /* Compute X */
 //     std::cout << "X[0], X[0+ldx*1] before dladiv = " << X[0] << ", " << X[0+ldx*1] << "\n";
      dladiv<T>(scale*B[0],scale*B[0+ldb*1],csr,csi,X[0],X[0+ldx*1]);
 //     std::cout << "X[0], X[0+ldx*1] after dladiv = " << X[0] << ", " << X[0+ldx*1] << "\n";
      xnorm=cln::abs(X[0])+cln::abs(X[0+ldx*1]);
    }
  }
  else
  {
    /* 2x2 System
*
*        Compute the real part of  C = ca A - w D  (or  ca A' - w D ) */
    cr[0][0]=ca*A[0]-wr*D1;
    cr[1][1]=ca*A[1+lda*1]-wr*D2;
    if(ltrans)
    {
      cr[0][1]=ca*A[1+lda*0];
      cr[1][0]=ca*A[0+lda*1];
    }
    else
    {
      cr[1][0]=ca*A[1+lda*0];
      cr[0][1]=ca*A[0+lda*1];
    }
    
    if(nw==1)
    {
      /* Real 2x2 system  (w is real)
*
*           Find the largest element in C */
      
      cmax=zero;
      icmax=0;
      
      for(j=1;j<=4;j++)
      {
        if(cln::abs(crv[j-1])>cmax)
        {
          cmax=cln::abs(crv[j-1]);
          icmax=j;
        }
      }
      
      /* If norm(C) < SMINI, use SMINI*identity. */
      
      if(cmax<smini)
      {
        bnorm=cln::max(cln::abs(B[0]),cln::abs(B[1+ldb*0]));
        if( (smini<one)&&(bnorm>one) )
        {
          if(bnorm>(bignum*smini) )
            scale=one/bnorm;
        }
        temp=scale/smini;
        X[0]=temp*B[0];
        X[1]=temp*B[1];
        xnorm=temp*bnorm;
        info=1;
        return;
      }
      
      /* Gaussian elimination with complete pivoting. */
      
      ur11=crv[icmax-1];
      cr21=crv[ipivot[1][icmax-1]-1];
      ur12=crv[ipivot[2][icmax-1]-1];
      cr22=crv[ipivot[3][icmax-1]-1];
      ur11r=one/ur11;
      lr21=ur11r*cr21;
      ur22=cr22-ur12*lr21;
      
      /* If smaller pivot < SMINI, use SMINI */
      
      if(cln::abs(ur22)<smini)
      {
        ur22=smini;
        info=1;
      }
      if(rswap[icmax-1])
      {
        br1=B[1];
        br2=B[0];
      }
      else
      {
        br1=B[0];
        br2=B[1];
      }
      br2=br2-lr21*br1;
      bbnd=cln::max(cln::abs(br1*(ur22*ur11r)),cln::abs(br2));
      if( (bbnd>one)&&(cln::abs(ur22)<one) )
      {
        if(bbnd>=(bignum*cln::abs(ur22)) )
          scale=one/bbnd;
      }
      
      xr2=br2*scale/ur22;
      xr1=scale*br1*ur11r-xr2*ur11r*ur12;
      if(zswap[icmax-1])
      {
        X[0]=xr2;
        X[1]=xr1;
      }
      else
      {
        X[0]=xr1;
        X[1]=xr2;
      }
      xnorm=cln::max(cln::abs(xr1),cln::abs(xr2));
      
      /* Further scaling if  norm(A) norm(X) > overflow */
      
      if( (xnorm>one)&&(cmax>one) )
      {
        if( xnorm>(bignum/cmax) )
        {
          temp=cmax/bignum;
          X[0]=temp*X[0];
          X[1]=temp*X[1];
          xnorm=temp*xnorm;
          scale=temp*scale;
        }
      }
    }
    else
    {
      /* Complex 2x2 system  (w is complex)
*
*           Find the largest element in C */
      
      ci[0][0]=-wi*D1;
      ci[1][0]=zero;
      ci[0][1]=zero;
      ci[1][1]=-wi*D2;
      cmax=zero;
      icmax=0;
      
      for(j=1;j<=4;j++)
      {
        if( (cln::abs(crv[j-1])+cln::abs(civ[j-1])) >cmax)
        {
          cmax=cln::abs(crv[j-1])+cln::abs(civ[j-1]);
          icmax=j;
        }
      }
      
      /* If norm(C) < SMINI, use SMINI*identity. */
      
      if(cmax<smini)
      {
        bnorm=cln::max(cln::abs(B[0])+cln::abs(B[0+ldb*1]),cln::abs(B[1])+cln::abs(B[1+ldb*1]));
        if( (smini<one)&&(bnorm>one) )
        {
          if(bnorm>(bignum*smini))
            scale=one/bnorm;
        }
        temp=scale/smini;
        X[0]=temp*B[0];
        X[1]=temp*B[1];
        X[0+ldx*1]=temp*B[0+ldb*1];
        X[1+ldx*1]=temp*B[1+ldb*1];
        xnorm=temp*bnorm;
        info=1;
        return;
      }
      
      /* Gaussian elimination with complete pivoting. */
      ur11=crv[icmax-1];
      ui11=civ[icmax-1];
      cr21=crv[ipivot[1][icmax-1]-1];
      ci21=civ[ipivot[1][icmax-1]-1];
      ur12=crv[ipivot[2][icmax-1]-1];
      ui12=civ[ipivot[2][icmax-1]-1];
      cr22=crv[ipivot[3][icmax-1]-1];
      ci22=civ[ipivot[3][icmax-1]-1];
      if( (icmax==1)||(icmax==4) )
      {
        /* Code when off-diagonals of pivoted C are real */
        
        if(cln::abs(ur11)>cln::abs(ui11))
        {
          temp=ui11/ur11;
          ur11r=one/(ur11*(one+temp*temp));
          ui11r=-temp*ur11r;
        }
        else
        {
          temp=ur11/ui11;
          ui11r=-one/(ui11*(one+temp*temp));
          ur11r=-temp*ui11r;
        }
        lr21=cr21*ur11r;
        li21=cr21*ui11r;
        ur12s=ur12*ur11r;
        ui12s=ur12*ui11r;
        ur22=cr22-ur12*lr21;
        ui22=ci22-ur12*li21;
      }
      else
      {
        /* Code when diagonals of pivoted C are real */
        ur11r=one/ur11;
        ui11r=zero;
        lr21=cr21*ur11r;
        li21=ci21*ur11r;
        ur12s=ur12*ur11r;
        ui12s=ui12*ur11r;
        ur22=cr22-ur12*lr21+ui12*li21;
        ui22=-ur12*li21-ui12*lr21;
      }
      u22abs=cln::abs(ur22)+cln::abs(ui22);
      
      /* If smaller pivot < SMINI, use SMINI */
      
      if(u22abs<smini)
      {
        ur22=smini;
        ui22=zero;
        info=1;
      }
      if(rswap[icmax-1])
      {
        br2=B[0];
        br1=B[1];
        bi2=B[0+1*ldb];
        bi1=B[1+1*ldb];
      }
      else
      {
        br1=B[0];
        br2=B[1];
        bi1=B[0+1*ldb];
        bi2=B[1+1*ldb];
      }
      br2=br2-lr21*br1+li21*bi1;
      bi2=bi2-li21*br1-lr21*bi1;
      bbnd=cln::max( ( cln::abs(br1)+ cln::abs(bi2) )*(u22abs*(cln::abs(ur11r)+cln::abs(ui11r))),cln::abs(br2)+cln::abs(bi2) );
      if( (bbnd>one)&&(u22abs<one) )
      {
        if(bbnd>=(bignum*u22abs))
        {
          scale=one/bbnd;
          br1=scale*br1;
          bi1=scale*bi1;
          br2=scale*br2;
          bi2=scale*bi2;
        }
      }
      
  //    std::cout << "xr1, xr2 before dladiv = " << xr1 << ", " << xr2 << "\n";
      dladiv<T>(br2,bi2,ur22,ui22,xr2,xi2);
      xr1=ur11r*br1-ui11r*bi1-ur12s*xr2+ui12s*xi2;
      xi1=ui11r*br1+ur11r*bi1-ui12s*xr2-ur12s*xi2;
      if(zswap[icmax-1])
      {
        X[0]=xr2;
        X[1]=xr1;
        X[0+1*ldx]=xi2;
        X[1+1*ldx]=xi1;
      }
      else
      {
        X[0]=xr1;
        X[1]=xr2;
        X[0+1*ldx]=xi1;
        X[1+1*ldx]=xi2;
      }
      xnorm=cln::max(cln::abs(xr1)+cln::abs(xi1),cln::abs(xr2)+cln::abs(xi2));
      
      /* Further scaling if  norm(A) norm(X) > overflow */
      
      if( (xnorm>one)&&(cmax>one) )
      {
        if(xnorm>(bignum/cmax))
        {
          temp=cmax/bignum;
          X[0]=temp*X[0];
          X[1]=temp*X[1];
          X[0+ldx*1]=temp*X[0+ldx*1];
          X[1+ldx*1]=temp*X[1+ldx*1];
          xnorm=temp*xnorm;
          scale=temp*scale;
        }
      }
    }
  }
  
  return;
  
  /* End of DLALN2 */
}

template<typename T>
void dlaln2(bool ltrans, int na, int nw, T smin, T ca, T* A, int lda, T D1, T D2, T* B, int ldb, T wr, T wi, T* X, int ldx, T& scale, T& xnorm, int& info)
{
  int icmax,j;
  T bbnd,bi1,bi2,bignum,bnorm,br1,br2,ci21;
  T ci22,cmax,cnorm,cr21,cr22,csi,csr,li21;
  T lr21,smini,smlnum,temp,u22abs,ui11,ui11r;
  T ui12,ui12s,ui22,ur11,ur11r,ur12,ur12s;
  T ur22,xi1,xi2,xr1,xr2;

  bool rswap[4];
  bool zswap[4];
  int ipivot[4][4];
  T zero,one,two;
  T ci[2][2];
  T* civ;
  T cr[2][2];
  T* crv;
  
 // std::cout << "called wrong dlaln2 \n";
 // exit(1);

  civ=&ci[0][0];
  crv=&cr[0][0];

  zswap[0]=false;
  zswap[1]=false;
  zswap[2]=true;
  zswap[3]=true;

  rswap[0]=false;
  rswap[1]=true;
  rswap[2]=false;
  rswap[3]=true;

  ipivot[0][0]=1;
  ipivot[1][0]=2;
  ipivot[2][0]=3;
  ipivot[3][0]=4;
  ipivot[0][1]=2;
  ipivot[1][1]=1;
  ipivot[2][1]=4;
  ipivot[3][1]=3;
  ipivot[0][2]=3;
  ipivot[1][2]=4;
  ipivot[2][2]=1;
  ipivot[3][2]=2;
  ipivot[0][3]=4;
  ipivot[1][3]=3;
  ipivot[2][3]=2;
  ipivot[3][3]=1;

  one=1.000000000000000;
  zero=0.000000000000000;
  two=2.00000000000000;
  
/*  .. Executable Statements ..
*
*     Compute BIGNUM */
  smlnum=two*dlamch<T>("S"); // safe minimum
  bignum=one/smlnum;
  smini=std::max(smin,smlnum);
  
  /* Don't check for input errors */
  info=0;
  
  /* Standard Initializations */
  scale=one;
  
  if(na==1)
  {
    /* 1 x 1  (i.e., scalar) system   C X = B */
    
    if(nw==1)
    {
      /* Real 1x1 system.
*
*           C = ca A - w D */
      csr=ca*A[0]-wr*D1;
      cnorm=fabs(csr);
      
      /* If | C | < SMINI, use C = SMINI */
      if(cnorm<smini)
      {
        csr=smini;
        cnorm=smini;
        info=1;
      }
      
      /* Check scaling for  X = B / C */
      bnorm=fabs(B[0]);
      if( (cnorm<one)&&(bnorm>one) )
      {
        if(bnorm>(bignum*cnorm))
          scale=one/bnorm;
      }
      
      /* Compute X */
      X[0]=B[0]*scale/csr;
      xnorm=fabs(X[0]);
    }
    else
    {
      /* Complex 1x1 system (w is complex)
*
*           C = ca A - w D */
      csr=ca*A[0]-wr*D1;
      csi=-wi*D1;
      cnorm=fabs(csr)+fabs(csi);
      
      /* If | C | < SMINI, use C = SMINI */
      if(cnorm<smini)
      {
        csr=smini;
        csi=zero;
        cnorm=smini;
        info=1;
      }
      
      /* Check scaling for  X = B / C */
      bnorm=fabs(B[0])+fabs(B[0+ldb*1]);
      if( (cnorm<one)&&(bnorm>one) )
      {
        if( bnorm>(bignum*cnorm) )
          scale=one/bnorm;
      }
      
      /* Compute X */
      dladiv<T>(scale*B[0],scale*B[0+ldb*1],csr,csi,X[0],X[0+ldx*1]);
      xnorm=fabs(X[0])+fabs(X[0+ldx*1]);
    }
  }
  else
  {
    /* 2x2 System
*
*        Compute the real part of  C = ca A - w D  (or  ca A' - w D ) */
    cr[0][0]=ca*A[0]-wr*D1;
    cr[1][1]=ca*A[1+lda*1]-wr*D2;
    if(ltrans)
    {
      cr[0][1]=ca*A[1+lda*0];
      cr[1][0]=ca*A[0+lda*1];
    }
    else
    {
      cr[1][0]=ca*A[1+lda*0];
      cr[0][1]=ca*A[0+lda*1];
    }
    
    if(nw==1)
    {
      /* Real 2x2 system  (w is real)
*
*           Find the largest element in C */
      
      cmax=zero;
      icmax=0;
      
      for(j=1;j<=4;j++)
      {
        if(fabs(crv[j-1])>cmax)
        {
          cmax=fabs(crv[j-1]);
          icmax=j;
        }
      }
      
      /* If norm(C) < SMINI, use SMINI*identity. */
      
      if(cmax<smini)
      {
        bnorm=std::max(fabs(B[0]),fabs(B[1+ldb*0]));
        if( (smini<one)&&(bnorm>one) )
        {
          if(bnorm>(bignum*smini) )
            scale=one/bnorm;
        }
        temp=scale/smini;
        X[0]=temp*B[0];
        X[1]=temp*B[1];
        xnorm=temp*bnorm;
        info=1;
        return;
      }
      
      /* Gaussian elimination with complete pivoting. */
      
      ur11=crv[icmax-1];
      cr21=crv[ipivot[1][icmax-1]-1];
      ur12=crv[ipivot[2][icmax-1]-1];
      cr22=crv[ipivot[3][icmax-1]-1];
      ur11r=one/ur11;
      lr21=ur11r*cr21;
      ur22=cr22-ur12*lr21;
      
      /* If smaller pivot < SMINI, use SMINI */
      
      if(fabs(ur22)<smini)
      {
        ur22=smini;
        info=1;
      }
      if(rswap[icmax-1])
      {
        br1=B[1];
        br2=B[0];
      }
      else
      {
        br1=B[0];
        br2=B[1];
      }
      br2=br2-lr21*br1;
      bbnd=std::max(fabs(br1*(ur22*ur11r)),fabs(br2));
      if( (bbnd>one)&&(fabs(ur22)<one) )
      {
        if(bbnd>=(bignum*fabs(ur22)) )
          scale=one/bbnd;
      }
      
      xr2=br2*scale/ur22;
      xr1=scale*br1*ur11r-xr2*ur11r*ur12;
      if(zswap[icmax-1])
      {
        X[0]=xr2;
        X[1]=xr1;
      }
      else
      {
        X[0]=xr1;
        X[1]=xr2;
      }
      xnorm=std::max(fabs(xr1),fabs(xr2));
      
      /* Further scaling if  norm(A) norm(X) > overflow */
      
      if( (xnorm>one)&&(cmax>one) )
      {
        if( xnorm>(bignum/cmax) )
        {
          temp=cmax/bignum;
          X[0]=temp*X[0];
          X[1]=temp*X[1];
          xnorm=temp*xnorm;
          scale=temp*scale;
        }
      }
    }
    else
    {
      /* Complex 2x2 system  (w is complex)
*
*           Find the largest element in C */
      
      ci[0][0]=-wi*D1;
      ci[1][0]=zero;
      ci[0][1]=zero;
      ci[1][1]=-wi*D2;
      cmax=zero;
      icmax=0;
      
      for(j=1;j<=4;j++)
      {
        if( (fabs(crv[j-1])+fabs(civ[j-1])) >cmax)
        {
          cmax=fabs(crv[j-1])+fabs(civ[j-1]);
          icmax=j;
        }
      }
      
      /* If norm(C) < SMINI, use SMINI*identity. */
      
      if(cmax<smini)
      {
        bnorm=std::max(fabs(B[0])+fabs(B[0+ldb*1]),fabs(B[1])+fabs(B[1+ldb*1]));
        if( (smini<one)&&(bnorm>one) )
        {
          if(bnorm>(bignum*smini))
            scale=one/bnorm;
        }
        temp=scale/smini;
        X[0]=temp*B[0];
        X[1]=temp*B[1];
        X[0+ldx*1]=temp*B[0+ldb*1];
        X[1+ldx*1]=temp*B[1+ldb*1];
        xnorm=temp*bnorm;
        info=1;
        return;
      }
      
      /* Gaussian elimination with complete pivoting. */
      ur11=crv[icmax-1];
      ui11=civ[icmax-1];
      cr21=crv[ipivot[1][icmax-1]-1];
      ci21=civ[ipivot[1][icmax-1]-1];
      ur12=crv[ipivot[2][icmax-1]-1];
      ui12=civ[ipivot[2][icmax-1]-1];
      cr22=crv[ipivot[3][icmax-1]-1];
      ci22=civ[ipivot[3][icmax-1]-1];
      if( (icmax==1)||(icmax==4) )
      {
        /* Code when off-diagonals of pivoted C are real */
        
        if(fabs(ur11)>fabs(ui11))
        {
          temp=ui11/ur11;
          ur11r=one/(ur11*(one+temp*temp));
          ui11r=-temp*ur11r;
        }
        else
        {
          temp=ur11/ui11;
          ui11r=-one/(ui11*(one+temp*temp));
          ur11r=-temp*ui11r;
        }
        lr21=cr21*ur11r;
        li21=cr21*ui11r;
        ur12s=ur12*ur11r;
        ui12s=ur12*ui11r;
        ur22=cr22-ur12*lr21;
        ui22=ci22-ur12*li21;
      }
      else
      {
        /* Code when diagonals of pivoted C are real */
        ur11r=one/ur11;
        ui11r=zero;
        lr21=cr21*ur11r;
        li21=ci21*ur11r;
        ur12s=ur12*ur11r;
        ui12s=ui12*ur11r;
        ur22=cr22-ur12*lr21+ui12*li21;
        ui22=-ur12*li21-ui12*lr21;
      }
      u22abs=fabs(ur22)+fabs(ui22);
      
      /* If smaller pivot < SMINI, use SMINI */
      
      if(u22abs<smini)
      {
        ur22=smini;
        ui22=zero;
        info=1;
      }
      if(rswap[icmax-1])
      {
        br2=B[0];
        br1=B[1];
        bi2=B[0+1*ldb];
        bi1=B[1+1*ldb];
      }
      else
      {
        br1=B[0];
        br2=B[1];
        bi1=B[0+1*ldb];
        bi2=B[1+1*ldb];
      }
      br2=br2-lr21*br1+li21*bi1;
      bi2=bi2-li21*br1-lr21*bi1;
      bbnd=std::max( ( fabs(br1)+fabs(bi2) )*(u22abs*(fabs(ur11r)+fabs(ui11r))),fabs(br2)+fabs(bi2) );
      if( (bbnd>one)&&(u22abs<one) )
      {
        if(bbnd>=(bignum*u22abs))
        {
          scale=one/bbnd;
          br1=scale*br1;
          bi1=scale*bi1;
          br2=scale*br2;
          bi2=scale*bi2;
        }
      }
      
      dladiv<T>(br2,bi2,ur22,ui22,xr2,xi2);
      xr1=ur11r*br1-ui11r*bi1-ur12s*xr2+ui12s*xi2;
      xi1=ui11r*br1+ur11r*bi1-ui12s*xr2-ur12s*xi2;
      if(zswap[icmax-1])
      {
        X[0]=xr2;
        X[1]=xr1;
        X[0+1*ldx]=xi2;
        X[1+1*ldx]=xi1;
      }
      else
      {
        X[0]=xr1;
        X[1]=xr2;
        X[0+1*ldx]=xi1;
        X[1+1*ldx]=xi2;
      }
      xnorm=std::max(fabs(xr1)+fabs(xi1),fabs(xr2)+fabs(xi2));
      
      /* Further scaling if  norm(A) norm(X) > overflow */
      
      if( (xnorm>one)&&(cmax>one) )
      {
        if(xnorm>(bignum/cmax))
        {
          temp=cmax/bignum;
          X[0]=temp*X[0];
          X[1]=temp*X[1];
          X[0+ldx*1]=temp*X[0+ldx*1];
          X[1+ldx*1]=temp*X[1+ldx*1];
          xnorm=temp*xnorm;
          scale=temp*scale;
        }
      }
    }
  }
  
  return;
  
  /* End of DLALN2 */
}
/* 
*
*  Purpose
*  =======
*
*  DTRSYL solves the real Sylvester matrix equation:
*
*     op(A)*X + X*op(B) = scale*C or
*     op(A)*X - X*op(B) = scale*C,
*
*  where op(A) = A or A**T, and  A and B are both upper quasi-
*  triangular. A is M-by-M and B is N-by-N; the right hand side C and
*  the solution X are M-by-N; and scale is an output scale factor, set
*  <= 1 to avoid overflow in X.
*
*  A and B must be in Schur canonical form (as returned by DHSEQR), that
*  is, block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;
*  each 2-by-2 diagonal block has its diagonal elements equal and its
*  off-diagonal elements of opposite sign.
*
*  Arguments
*  =========
*
*  TRANA   (input) CHARACTER*1
*          Specifies the option op(A):
*          = 'N': op(A) = A    (No transpose)
*          = 'T': op(A) = A**T (Transpose)
*          = 'C': op(A) = A**H (Conjugate transpose = Transpose)
*
*  TRANB   (input) CHARACTER*1
*          Specifies the option op(B):
*          = 'N': op(B) = B    (No transpose)
*          = 'T': op(B) = B**T (Transpose)
*          = 'C': op(B) = B**H (Conjugate transpose = Transpose)
*
*  ISGN    (input) INTEGER
*          Specifies the sign in the equation:
*          = +1: solve op(A)*X + X*op(B) = scale*C
*          = -1: solve op(A)*X - X*op(B) = scale*C
*
*  M       (input) INTEGER
*          The order of the matrix A, and the number of rows in the
*          matrices X and C. M >= 0.
*
*  N       (input) INTEGER
*          The order of the matrix B, and the number of columns in the
*          matrices X and C. N >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,M)
*          The upper quasi-triangular matrix A, in Schur canonical form.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  B       (input) DOUBLE PRECISION array, dimension (LDB,N)
*          The upper quasi-triangular matrix B, in Schur canonical form.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= max(1,N).
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the M-by-N right hand side matrix C.
*          On exit, C is overwritten by the solution matrix X.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M)
*
*  SCALE   (output) DOUBLE PRECISION
*          The scale factor, scale, set <= 1 to avoid overflow in X.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          = 1: A and B have common or very close eigenvalues; perturbed
*               values were used to solve the equation (but the matrices
*               A and B are unchanged).
*
*  =====================================================================
* */

template<typename T>
void dtrsyl(const std::string& trana, const std::string& tranb, int isgn, int m, int n, T* A, int lda, T* C, T& scale, int& info)
{
  const T zero = 0.0;
  const T one = 1.0;
  
  bool notrna, notrnb;
  int ierr,j,k,k1,k2,knext,l,l1,l2,lnext;
  T a11,bignum,da11,db,eps,scaloc,sgn,smin,smlnum,suml,sumr,xnorm;
  
  T dum[1];
  T vec[2][2];
  T X[2][2];
  
  /* .. Executable Statements ..
   * 
   * Decode and Test input parameters */
  
  notrna=(trana=="N");
  notrnb=(tranb=="N");
  
  info=0;
  if( !notrna && !(trana=="T") && !(trana=="C") )
    info =-1;
  else if( !notrnb && !(tranb=="T") && !(tranb=="C") )
    info=-2;
  else if( (isgn!=1)&&(isgn!=-1) )
    info=-3;
  else if(m<0)
    info=-4;
  else if(n<0)
    info=-5;
  else if(lda<std::max(1,m))
    info=-7;
  else if(ldb<std::max(1,n))
    info=-9;
  else if(ldc<std::max(1,m))
    info=-11;
  
  if(info!=0)
  {
    std::cout << "dtrsyl: bad input parameter \n";
    return;
  }
  
  /* Quick return if possible */
  
  eps=dlamch<T>("P"); // precision
  smlnum=dlamch("S"); // safe minimum
  bignum=one/smlnum;
  dlabad<T>(smlnum,bignum);
  smlnum=smlnum*double(m*n)/eps;
  bignum=one/slmnum;
  
  smin=std::max(smlnum,eps*dlange<T>("M",m,m,A,lda,dum),eps*dlange<T>("M",n,n,B,ldb,dum));
  scale=one;
  sgn=isgn;
  
  if(notrna&&notrnb)
  {
    /* Solve    A*X + ISGN*X*B = scale*C.
     * 
     * The (K,L)th block of X is determined starting from
     * bottom-left corner column by column by
     * 
     * A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
     * 
     * Where
     *           M                         L-1
     * R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)].
     *         I=K+1                       J=1
     * 
     * Start column loop (index = L)
     * L1 (L2) : column index of the first (first) row of X(K,L). */
    
    lnext=1;
    for(l=1;l<=n;l++)
    {
      if(l<lnext)
        continue;
      if(l==n)
      {
        l1=l;
        l2=l;
      }
      else
      {
        if(B[l+ldb*(l-1)]!=zero)
        {
          l1=l;
          l2=l+1;
          lnext=l+2;
        }
        else
        {
          l1=l;
          l2=l;
          lnext=l+1;
        }
      }
      
      /* Start row loop (index = K)
       * K1 (K2): row index of the first (last) row of X(K,L). */
      
      knext=m;
      for(k=m;k>=1;k--)
      {
        if(k>knext)
          continue;
        if(k==1)
        {
          k1=k;
          k2=k;
        }
        else
        {
          if(A[k-1+lda*(k-2)]!=zero)
          {
            k1=k-1;
            k2=k;
            knext=k-2;
          }
          else
          {
            k1=k;
            k2=k;
            knext=k-1;
          }
        }
        
        if( (l1==l2)&&(k1==k2) )
        {
          suml=ddot<T>(m-k1,&A[k1-1+lda*(std::min(k1+1,m)-1)],lda,&C[std::min(k1+1,m)-1+ldc*(l1-1)],1);
          sumr=ddot<T>(l1-1,&C[k1-1+0*ldc],ldc,&B[0+ldb*(l1-1)],1);
          vec[1][1]=C[k1-1+ldc*(l1-1)]-(suml+sgn*sumr);
          scaloc=one;
          
          a11=A[k1-1+lda*(k1-1)]+sgn*B[l1-1+ldb*(l1-1)];
          da11=fabs(a11);
          if(da11<=smin)
          {
            a11=smin;
            da11=smin;
            info=1;
          }
          db=fabs(vec[0][0]);
          if( (da11<one)&&(db>one) )
          {
            if(db>bignum*da11)
              scaloc=one/db;
          }
          X[0][0]=( vec[0][0]*scaloc )/a11;
          
          if(scaloc!=one)
          {
            for(j=1;j<=n;j++)
              dscal<T>(m,scaloc,C[0+ldc*(j-1)],1);
            scale=scale*scaloc;
          }
          C[k1-1+ldc*(l1-1)]=X[0][0];
        }
        else if( (l1==l2)&&(k1!=k2) )
        {
          suml=ddot<T>(m-k2,&A[k1-1+lda*(std::min(k2+1,m)-1)],lda,&C[std::min(k2+1,m)-1+ldc*(m-1)],1);
          sumr=ddot<T>(l1-1,&C[k1-1+ldc*0],ldc,&B[0+ldb*(l1-1)],1);
          vec[0][0]=C[k1-1][l1-1]-(suml+sgn*sumr);
          suml=ddot<T>(m-k2,&A[k2-1+lda*(std::min(k2+1,m)-1)],lda,&C[std::min(k2+1,m)-1+ldc*(l1-1)],1);
          sumr=ddot<T>(l1-1,&C[k2-1+ldc*0],ldc,&B[0+ldb*(l1-1)],1);
          vec[1][0]=C[k2-1+ldc*(l1-1)]-(suml+sgn*sumr);
          dlaln2<T>(false,2,1,smin,one,&A[k1-1+lda*(k1-1)],lda,one,one,&vec[0][0],2,-sgn*B[l1-1+ldb*(l1-1)],zero,X,2,scaloc,xnorm,ierr);
          if(ierr!=0)
            info=1;
          
          if(scaloc!=one)
          {
            for(j=1;j<=n;j++)
              dscal<T>(m,scaloc,&C[0+ldc*(j-1)],1);
            scale=scale*scaloc;
          }
          C[k1-1+ldc*(l1-1)]=X[0][0];
          C[k2-1+ldc*(l1-1)]=X[1][0];
        }
        else if( (l1!=l2)&&(k1==k2) )
        {
          suml=ddot<T>(m-k1,&A[k1-1+lda*(std::min(k1+1,m)-1)],lda,&C[std::min(k1+1,m)-1+ldc*(l1-1)],1);
          sumr=ddot<T>(l1-1,&C[k1-1+ldc*0],ldc,B[0+ldb*(l1-1)],1);
          vec[0][0]=sgn*(C[k1-1+ldc*(k1-1)]-(suml+sgn*sumr));
          
          suml=ddot<T>(m-k1,A[k1-1+lda*(std::min(k1+1,m)-1)],lda,&C[std::min(ki1+1,m)-1+ldc*(l2-1)],1);
          sumr=ddot<T>(l1-1,&C[k1-1+ldc*0],ldc,&B[0+ldb*(l2-1)],1);
          vec[1][0]=sgn*(C[k1-1+ldc*(l2-1)]-(suml+sgn*sumr) );
          
          dlaln2<T>(true,2,1,smin,one,&B[l1-1+ldb*(l1-1)],ldb,one,one,&vec[0][0],2,-sgn*A[k1-1+lda*(k1-1)],zero,&X[0][0],2,scaloc,xnorm,ierr);
          if(ierr!=0)
            info=1;
          
          if(scaloc!=one)
          {
            for(j=1;j<=n;j++)
              dscal<T>(m,scaloc,C[0+ldc*(j-1)],1);
            scale=scale*scaloc;
          }
          C[k1-1+ldc*(l1-1)]=X[0][0];
          C[k1-1+ldc*(l2-1)]=X[1][0];
        }
        else if( (l1!=l2)&&(k1!=k2) )
        {
          suml=ddot<T>(m-k2,&A[k1-1+lda*(std::min(k2+1,m)-1)],lda,C[std::min(k2+1,m)-1+ldc*(l1-1)],1);
          sumr=ddot<T>(l1-1,&C[k1-1+ldc*0],ldc,&B[0+ldb*(l1-1)],1);
          vec[0][0]=C[k1-1+ldc*(l1-1)]-(suml+sgn*sumr);
          
          suml=ddot<T>(m-k2,&A[k1-1+lda*(std::min(k2+1,m)-1)],lda,&C[std::min(k2+1,m)-1+ldc*(l2-1)],1);
          sumr=ddot<T>(l1-1,&C[k1-1+ldc*0],ldc,&B[0+ldb*(l2-1)],1);
          vec[0][1]=C[k1-1+ldc*(l2-1)]-(suml+sgn*sumr);
          
          suml=ddot<T>(m-k2,&A[k2-1+lda*(std::min(k2+1,m)-1)],lda,&C[std::min(k2+1,m)+ldc*(l1-1)],1);
          sumr=ddot<T>(l1-1,&C[k2-1+ldc*0],ldc,&B[0+ldb*(l1-1)],1);
          vec[1][0]=C[k2-1+ldc*(l1-1)]-(suml+sgn*sumr);
          
          suml=ddot<T>(m-k2,&A[k2-1+lda*(std::min(k2+1,m)-1)],lda,&C[std::min(k2+1,m)+ldc*(l2-1)],1);
          sumr=ddot<T>(l1-1,&C[k2-1+ldc*0],ldc,&B[0+ldb*(l1-1)],1);
          vec[1][1]=C[k2-1+ldc*(l2-1)]-(suml+sgn*sumr);
          
          dlasy2<T>(false,false,isgn,2,2,&A[k1-1+lda*(k1-1)],lda,&B[l1-1+ldb*(l1-1)],ldb,&vec[0][0],2,scaloc,&X[0][0],2,xnorm,ierr);
          if(ierr!=0)
            info=1;
          
          if(scaloc!=one)
          {
            for(j=1;j<=n;j++)
              dscal<T>(m,scaloc,&C[0+ldc*(j-1)],1);
            scale=scale*scaloc;
          }
          C[k1-1+ldc*(l1-1)]=X[0][0];
          C[k1-1+ldc*(l2-1)]=X[0][1];
          C[k2-1+ldc*(l1-1)]=X[1][0];
          C[k2-1+ldc*(l2-1)]=X[1][1];
        }
      } //for(k=m;...
    } // for(l=1;...
  }
  else if(!notrna&&!notrnb)
  {
    /* Solve    A' *X + ISGN*X*B = scale*C.
     *
     * The (K,L)th block of X is determined starting from
     * upper-left corner column by column by
     *
     * A(K,K)'*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
     *
     * Where
     * K-1                        L-1
     * R(K,L) = SUM [A(I,K)'*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)]
     * I=1                        J=1
     *
     * Start column loop (index = L)
     * L1 (L2): column index of the first (last) row of X(K,L) */
    lnext=1;
    for(l=1;l<=n;l++)
    {
      if(l<lnext)
        continue;
      if(l==n)
      {
        l1=l;
        l2=l;
      }
      else
      {
        if(B[l+ldb*(l-1)]!=zero)
        {
          l1=l;
          l2=l+1;
          lnext=l+2;
        }
        else
        {
          l1=l;
          l2=l;
          lnext=l+1;
        }
      }
      
      /* Start row loop (index = K)
       * K1 (K2): row index of the first (last) row of X(K,L) */
      
      knext=1;
      for(k=1;k<=m;k++)
      {
        if(k<knext)
          continue;
        if(k==m)
        {
          k1=k;
          k2=k;
        }
        else
        {
          if(A[k+lda*(k-1)]!=zero)
          {
            k1=k;
            k2=k+1;
            knext=k+2;
          }
          else
          {
            k1=k;
            k2=k;
            knext=k+1;
          }
        }
        
        if( (l1==l2)&&(k1==k2) )
        {
          suml=ddot<T>(k1-1,&A[0+lda*(k1-1)],1,&C[0+ldc*(l1-1)],1);
          sumr=ddot<T>(l1-1,&C[k1-1+ldc*0],ldc,&B[0+ldb*(l1-1)],1);
          vec[0][0]=C[k1-1+ldc*(l1-1)]-(suml+sgn*sumr);
          scaloc=one;
          
          a11=A[k1-1+lda*(k1-1)]+sgn*B[l1-1+ldb*(l1-1)];
          da11=fabs(a11);
          if(da11<=smin)
          {
            a11=smin;
            da11=smin;
            info=1;
          }
          db=fabe(vec[0][0]);
          if( (da11<one)&&(db>one) )
          {
            if(db>(bignum*da11))
              scaloc=one/db;
          }
          X[0][0]=( vec[0][0]*scaloc )/a11;
          
          if(scaloc!=one)
          {
            for(j=1;j<=n;j++)
              dscal<T>(m,scaloc,&C[0+ldc*(j-1)],1);
            scale=scale*scaloc;
          }
          C[k1-1+ldc*(l1-1)]=X[0][0];
        }
        else if( (l1==l2)&&(k1!=k2) )
        {
          suml=ddot<T>(k1-1,&A[0+lda*(k1-1)],1,&C[0+ldc*(l1-1)],1);
          sumr=ddot<T>(l1-1,&C[k1-1+ldc*0],ldc,&B[0+ldb*(l1-1)],1);
          vec[0][0]=C[k1-1+ldc*(l1-1)]-(suml+sgn*sumr);
          
          suml=ddot<T>(k1-1,&A[0+lda*(k2-1)],1,&C[0+ldc*(l1-1)],1);
          sumr=ddot<T>(l1-1,&C[k2-1+ldc*0],ldc,&B[0+ldb*(l1-1)],1);
          vec[1][0]=C[k2-1+ldc*(l1-1)]-(suml+sgn*sumr);
          
          dlaln2<T>(true,2,1,smin,one,&A[k1-1+lda*(k1-1)],lda,one,one,&vec[0][0],2,-sgn*B[l1-1+ldb*(l1-1)],zero,&X[0][0],2,scaloc,xnorm,ierr);
          if(ierr!=0)
            info=1;
          
          if(scaloc!=one)
          {
            for(j=1;j<=n;j++)
              dscal<T>(m,scaloc,&C[0+ldc*(j-1)],1);
            scale=scale*scaloc;
          }
          C[k1-1+ldc*(l1-1)]=X[0][0];
          C[k2-1+ldc*(l1-1)]=X[1][0];
        }
        else if( (l1!=l2)&&(k1==k2) )
        {
          suml=ddot<T>(k1-1,&A[0+lda*(k1-1)],1,&C[0+ldc*(l1-1)],1);
          sumr=ddot<T>(l1-1,&C[k1-1+ldc*0],ldc,&B[0+ldb*(l1-1)],1);
          vec[0][0]=C[k1-1+ldc*(l1-1)]-(suml+sgn*sumr);
          
          suml=ddot<T>(k1-1,&A[0+lda*(k1-1)],1,&C[0+ldc*(l2-1)],1);
          sumr=ddot<T>(l1-1,&C[k1-1+ldc*0],ldc,&B[0+ldb*(l2-1)],1);
          vec[1][0]=C[k1-1+ldc*(l2-1)]-(suml+sgn*sumr);
          
          dlaln2<T>(true,2,1,smin,one,&B[l1-1+ldb*(l1-1)],ldb,one,one,&vec[0][0],2,-sgn*A[k1-1+lda*(k1-1)],zero,&X[0][0],2,scaloc,xnorm,ierr);
          
          if(ierr!=0)
            info=1;
          
          if(scaloc!=one)
          {
            for(j=1;j<=n;j++)
              dscal<T>(m,scaloc,&C[0+ldc*(j-1)],1);
            scale=scale*scaloc;
          }
          C[k1-1+ldc*(l1-1)]=X[0][0];
          C[k1-1+ldc*(l2-1)]=X[1][0];
        }
        else if( (l1!=l2)&&(k1!=k2) )
        {
          suml=ddot<T>(k1-1,&A[0+lda*(k1-1)],1,&C[0+ldc*(l1-1)],1);
          sumr=ddot<T>(l1-1,&C[k1-1+ldc*0],ldc,&B[0+ldb*(l1-1)],1);
          vec[0][0]=C[k1-1+ldc*(l1-1)]-(suml+sgn*sumr);
          
          suml=ddot<T>(k1-1,&A[0+lda*(k1-1)],1,&C[0+ldc*(l2-1)],1);
          sumr=ddot<T>(l1-1,&C[k1-1+ldc*0],ldc,&B[0+ldb*(l2-1)],1);
          vec[0][1]=C[k1-1+ldc*(l2-1)]-(suml+sgn*sumr);
          
          suml=ddot<T>(k1-1,&A[0+lda*(k2-1)],1,&C[0+ldc*(l1-1)],1);
          sumr=ddot<T>(l1-1,&C[k2-1+ldc*0],ldc,&B[0+ldb*(l1-1)],1);
          vec[1][0]=C[k2-1+ldc*(l1-1)]-(suml+sgn*sumr);
          
          suml=ddot<T>(k1-1,&A[0+lda*(k2-1)],1,&C[0+ldc*(l2-1)],1);
          sumr=ddot<T>(l1-1,&C[k2-1+ldc*0],ldc,&B[0+ldb*(l2-1)],1);
          vec[1][1]=C[k2-1+ldc*(l2-1)]-(suml+sgn*sumr);
          
          dlasy2<T>(true,false,isgn,2,2,&A[k1-1+lda*(k1-1)],lda,&B[l1-1+ldb*(l1-1)],ldb,&vec[0][0],2,scaloc,&X[0][0],2,xnorm,ierr);
          if(ierr!=0)
            info=1;
          
          if(scaloc!=one)
          {
            for(j=1;j<=n;j++)
              dscal<T>(m,scaloc,&C[0+ldc*(j-1)],1);
            scale=scale*scaloc;
          }
          C[k1-1+ldc*(l1-1)]=X[0][0];
          C[k1-1+ldc*(l2-1)]=X[0][1];
          C[k2-1+ldc*(l1-1)]=X[1][0];
          C[k2-1+ldc*(l2-1)]=X[1][1];
        }
      } // for(k=1;
    } // for(l=1;
     
          
      
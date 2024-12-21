/*
*
*  Purpose
*  =======
*
*  DTRMM  performs one of the matrix-matrix operations
*
*     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
*
*  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry,  SIDE specifies whether  op( A ) multiplies B from
*           the left or right as follows:
*
*              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
*
*              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A'.
*
*              TRANSA = 'C' or 'c'   op( A ) = A'.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain the matrix  B,  and  on exit  is overwritten  by the
*           transformed matrix.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
* */

template<typename T>
void dtrmm(const std::string& side,const std::string& uplo,const std::string& transa,const std::string& diag,int m, int n, T alpha, T* A, int lda, T* B, int ldb)
{
  bool lside,nounit,upper;
  int i,info,j,k,nrowa;
  T temp;
  T one = 1.0;
  T zero = 0.0;
  
  /*
*     .. Executable Statements ..
*
*     Test the input parameters. */
  
  lside=(side=="L");
  if(lside)
    nrowa=m;
  else
    nrowa=n;
  nounit=(diag=="N");
  upper=(uplo=="U");
  
  info=0;
  if( (!lside)&&(!(side=="R"))
    info=1;
  else if( (!upper)&&(lsame=="L") )
    info=2;
  else if( !(transa=="N")&&!(transa=="T")&&!(transa=="C") )
    info=3;
  else if( !(diag=="U")&&!(diag=="N") )
    info=4;
  else if(m<0)
    info=5;
  else if(n<0)
    info=6;
  else if(lda<std::max(1,nrowa))
    info=9;
  else if(ldb<std::max(1,m))
    info=11;
  
  if(info!=0)
  {
    std::cout << "Error in dtrmm. Bad input parameter. \n";
    return;
  }
  
  /* Quick return if possible. */
  if(n==0)
    return;
  
  /* And when  alpha.eq.zero. */
  
  if(alpha==zero)
  {
    for(j=1;j<=n;j++)
      for(i=1;i<=m;i++)
        B[i-1+ldb*(j-1)]=zero;
    return;
  }
  
  /* Start the operations. */
  
  if(lside)
  {
    if(transa=="N")
    {
      /* Form  B := alpha*A*B. */
      if(upper)
      {
        for(j=1;j<=n;j++)
          for(k=1;k<=m;k++)
          {
            if(B[k-1+ldb*(j-1)]!=zero)
            {
              temp=alpha*B[k-1+lda*(j-1)];
              for(i=1;i<=(k-1);i++)
                B[i-1+ldb*(j-1)]+=temp*A[i-1+lda*(k-1)];
              if(nounit)
                temp=temp*A[k-1+lda*(k-1)];
              B[k-1+ldb*(j-1)]=temp;
            }
          }
      }
      else
      {
        for(j=1;j<=n;j++)
          for(k=m;k>=1;k--)
          {
            if(B[k-1+ldb*(j-1)]!=zero)
            {
              temp=alpha*B[k-1+ldb*(j-1)];
              B[k-1+ldb*(j-1)]=temp;
              if(nounit)
                B[k-1+ldb*(j-1)]=B[k-1+ldb*(j-1)]*A[k-1+lda*(k-1)];
              for(i=k+1;i<=m;i++)
                B[i-1+ldb*(j-1)]+=temp*A[i-1+lda*(k-1)];
            }
          }
      }
    }
    else
    {
      /* Form  B := alpha*B*A'. */
      if(upper)
      {
        for(j=1;j<=n;j++)
          for(i=m;i>=1;i--)
          {
            temp=B[i-1+ldb*(j-1)];
            if(nounit)
              temp=temp*A[i-1+lda*(i-1)];
            for(k=1;k<=(i-1);k++)
              temp=temp+A[k-1+lda*(i-1)]*B[k-1+ldb*(j-1)];
            B[i-1+ldb*(j-1)]=alpha*temp;
          }
      }
      else
      {
        for(j=1;j<=n;j++)
          for(i=1;i<=m;i++)
          {
            temp=B[i-1+ldb*(j-1)];
            if(nounit)
              temp=temp*A[i-1+lda*(i-1)];
            for(k=i+1;k<=m;k++)
              temp+=A[k-1+lda*(i-1)]*B[k-1+lda*(j-1)];
            B[i-1+ldb*(j-1)]=alpha*temp;
          }
      }
    }
  }
  else
  {
    if(transa=="N")
    {
      /* Form  B := alpha*B*A. */
      
      if(upper)
      {
        for(j=n;j>=1;j--)
        {
          temp=alpha;
          if(nounit)
            temp=temp*A[j-1+lda*(j-1)];
          for(i=1;i<=m;i++)
            B[i-1+ldb*(j-1)]=temp*B[i-1+ldb*(j-1)];
          for(k=1;k<=(j-1);k++)
          {
            if(A[k-1+lda*(j-1)]!=zero)
            {
              temp=alpha*A[k-1+lda*(j-1)];
              for(i=1;i<=m;i++)
                B[i-1+ldb*(j-1)]+=temp*B[i-1+ldb*(k-1)];
            }
          }
        }
      }
      else
      {
        for(j=1;j<=n;j++)
        {
          temp=alpha;
          if(nounit)
            temp=temp*A[j-1+lda*(j-1)];
          for(i=1;i<=m;i++)
            B[i-1+ldb*(j-1)]=temp*B[i-1+ldb*(j-1)];
          for(k=j+1;k<=n;k++)
          {
            if(A[k-1+lda*(j-1)]!=zero)
            {
              temp=alpha*A[k-1+lda*(j-1)];
              for(i=1;i<=m;i++)
                B[i-1+ldb*(j-1)]+=temp*B[i-1+ldb*(k-1)];
            }
          }
        }
      }
    }
    else
    {
      /* Form  B := alpha*B*A'. */
      if(upper)
      {
        for(k=1;k<=n;k++)
        {
          for(j=1;j<=(k-1);j++)
          {
            if(A[j-1+lda*(k-1)]!=zero)
            {
              temp=alpha*A[j-1+lda*(k-1)];
              for(i=1;i<=m;i++)
                B[i-1+ldb*(j-1)]+=temp*B[i-1+ldb*(k-1)];
            }
          }
          temp=alpha;
          if(nounit)
            temp=temp*A[k-1+lda*(k-1)];
          if(temp!=one)
          {
            for(i=1;i<=m;i++)
              B[i-1+ldb*(k-1)]=temp*B[i-1+ldb*(k-1)];
          }
        }
      }
      else
      {
        for(k=n;k>=1;k--)
        {
          for(j=k+1;j<=n;j++)
          {
            if(A[j-1+lda*(k-1)]!=zero)
            {
              temp=alpha*A[j-1+lda*(k-1)];
              for(i=1;i<=m;i++)
                B[i-1+ldb*(j-1)]+=temp*B[i-1+ldb*(k-1)];
            }
          }
          temp=alpha;
          if(nounit)
            temp=temp*A[k-1+lda*(k-1)];
          if(temp!=one)
          {
            for(i=1;i<=m;i++)
              B[i-1+ldb*(k-1)]=temp*B[i-1+ldb*(k-1)];
          }
        }
      }
    }
  }
    
  return;
  
  /* End of DTRMM . */
}
          
      

  
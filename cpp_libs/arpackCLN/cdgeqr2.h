template<typename T>
void dgeqr2(int m, int n, T* A, int lda, T* tau, T* work, int& info)
{
  /*
*
*  Purpose
*  =======
*
*  DGEQR2 computes a QR factorization of a real m by n matrix A:
*  A = Q * R.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix A.
*          On exit, the elements on and above the diagonal of the array
*          contain the min(m,n) by n upper trapezoidal matrix R (R is
*          upper triangular if m >= n); the elements below the diagonal,
*          with the array TAU, represent the orthogonal matrix Q as a
*          product of elementary reflectors (see Further Details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(1) H(2) . . . H(k), where k = min(m,n).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
*  and tau in TAU(i).
*
*  =====================================================================
* */
  
  const T one=1.0;
  int i, k;
  T aii;
  
  /* 
*     .. Executable Statements ..
*
*     Test the input arguments */
  
  info=0;
  if(m<0)
    info=-2;
  else if(lda<std::max(1,m))
    info=-4;
  
  if(info!=0)
  {
    std::cout << "dgeqr2: Error in input parameters \n";
    return;
  }
  
  k=std::min(m,n);
  
  for(i=1;i<=k;i++)
  {
    /* Generate elementary reflector H(i) to annihilate A(i+1:m,i) */
    dlarfg<T>(m-i+1,A[i-1+lda*(i-1)],&A[std::min(i+1,m)-1+lda*(i-1)],1,tau[i-1]);
    
    if(i<n)
    {
      /* Apply H(i) to A(i:m,i+1:n) from the left */
      aii=A[i-1+lda*(i-1)];
      A[i-1+lda*(i-1)]=one;
      dlarf<T>("L",m-i+1,n-1,&A[i-1+lda*(i-1)],1,tau[i-1],&A[i-1+lda*i],lda,work);
      A[i-1+lda*(i-1)]=aii;
    }
  }
  
  return;
  
  /* End of DGEQR2 */
}
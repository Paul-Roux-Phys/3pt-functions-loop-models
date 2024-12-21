template<typename T>
void dlahqr(bool wantt, bool wantz, int n, int ilo, int ihi, T* h, int ldh, T* wr, T* wi, int iloz, int ihiz, T* z, int ldz, int& info)
{  
/*
*
*  Purpose
*  =======
*
*  DLAHQR is an auxiliary routine called by DHSEQR to update the
*  eigenvalues and Schur decomposition already computed by DHSEQR, by
*  dealing with the Hessenberg submatrix in rows and columns ILO to IHI.
*
*  Arguments
*  =========
*
*  WANTT   (input) LOGICAL
*          = .TRUE. : the full Schur form T is required;
*          = .FALSE.: only eigenvalues are required.
*
*  WANTZ   (input) LOGICAL
*          = .TRUE. : the matrix of Schur vectors Z is required;
*          = .FALSE.: Schur vectors are not required.
*
*  N       (input) INTEGER
*          The order of the matrix H.  N >= 0.
*
*  ILO     (input) INTEGER
*  IHI     (input) INTEGER
*          It is assumed that H is already upper quasi-triangular in
*          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless
*          ILO = 1). DLAHQR works primarily with the Hessenberg
*          submatrix in rows and columns ILO to IHI, but applies
*          transformations to all of H if WANTT is .TRUE..
*          1 <= ILO <= max(1,IHI); IHI <= N.
*
*  H       (input/output) DOUBLE PRECISION array, dimension (LDH,N)
*          On entry, the upper Hessenberg matrix H.
*          On exit, if WANTT is .TRUE., H is upper quasi-triangular in
*          rows and columns ILO:IHI, with any 2-by-2 diagonal blocks in
*          standard form. If WANTT is .FALSE., the contents of H are
*          unspecified on exit.
*
*  LDH     (input) INTEGER
*          The leading dimension of the array H. LDH >= max(1,N).
*
*  WR      (output) DOUBLE PRECISION array, dimension (N)
*  WI      (output) DOUBLE PRECISION array, dimension (N)
*          The real and imaginary parts, respectively, of the computed
*          eigenvalues ILO to IHI are stored in the corresponding
*          elements of WR and WI. If two eigenvalues are computed as a
*          complex conjugate pair, they are stored in consecutive
*          elements of WR and WI, say the i-th and (i+1)th, with
*          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the
*          eigenvalues are stored in the same order as on the diagonal
*          of the Schur form returned in H, with WR(i) = H(i,i), and, if
*          H(i:i+1,i:i+1) is a 2-by-2 diagonal block,
*          WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i).
*
*  ILOZ    (input) INTEGER
*  IHIZ    (input) INTEGER
*          Specify the rows of Z to which transformations must be
*          applied if WANTZ is .TRUE..
*          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
*
*  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
*          If WANTZ is .TRUE., on entry Z must contain the current
*          matrix Z of transformations accumulated by DHSEQR, and on
*          exit Z has been updated; transformations are applied only to
*          the submatrix Z(ILOZ:IHIZ,ILO:IHI).
*          If WANTZ is .FALSE., Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z. LDZ >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          > 0: DLAHQR failed to compute all the eigenvalues ILO to IHI
*               in a total of 30*(IHI-ILO+1) iterations; if INFO = i,
*               elements i+1:ihi of WR and WI contain those eigenvalues
*               which have been successfully computed.
*
*  =====================================================================
* */
  const T zero=0.0;
  const T one=1.0;
  const T dat1=0.75;
  const T dat2=-0.4375;
  
  int i,i1,i2,itn,its,j,k,l,m,nh,nr,nz;
  T cs,h00,h10,h11,h12,h21,h22,h33,h33s;
  T h43h34, h44, h44s, ovfl, s, smlnum,sn,sum;
  T t1,t2,t3,tst1,ulp,unfl,v1,v2,v3;
  
  T v[3];
  T work[1];
  
  /* Quick return if possible */
  if(n==0)
    return;
  if(ilo==ihi)
  {
    wr[ilo-1]=h[ilo-1+ldh*(ilo-1)];
    wi[ilo-1]=zero;
    return;
  }
  
  nh=ihi-ilo+1;
  nz=ihiz-iloz+1;
  
  /* Set machine-dependent constants for the stopping criterion. 
   * If norm(H) <= sqrt(OVFL), overflow should not occur. */
  
  unfl=dlamch<T>("S"); // safe minimum
  ovfl=one/unfl;
  dlabad<T>(unfl,ovfl);
  ulp=dlamch<T>("P"); // precision
  smlnum=unfl*(nh/ulp);
  
  /* I1 and I2 are the indices of the first row and last column of H
   * to which transformations must be applied. If eigenvalues only are
   * being computed, I1 and I2 are set inside the main loop. */
  if(wantt)
  {
    i1=1;
    i2=n;
  }
  
  /* ITN is the total number of QR iterations allowed. */
  
  /* The main loop begins here. I is the loop index and decreases from
   * IHI to ILO in steps of 1 or 2. Each iteration of the loop works
   * with the active submatrix in rows and columns L to I.
   * Eigenvalues I+1 to IHI have already converged. Either L = ILO or
   * H(L,L-1) is negligible so that the matrix splits. */
  
  i=ihi;
label10:
  l=ilo;
  if(i<ilo)
    goto label150;
  
  /* Perform QR iterations on rows and columns ILO to I until a
   * submatrix of order 1 or 2 splits off at the bottom because a
   * subdiagonal element has become negligible. */
  
  for(its=0;its<=itn;its++)
  {
    /* Look for a single small subdiagonal element. */
    
    for(k=i;k>=(l+1);k--)
    {
      tst1=fabs(h[k-2+ldh*(k-2)])+fabs(h[k-1+ldh*(k-1)]);
      if(tst1==zero)
        tst1=dlanhs<T>("1",i-l+1,&h[l-1+ldh*(l-1)],ldh,work);
      if(fabs(h[k-1+ldh*(k-2)])<=std::max(ulp*tst1,smlnum))
        goto label30;
    }
label30:
    l=k;
    if(l>ilo)
    {
      /* H(L,L-1) is negligible */
      
      h[l-1+ldh*(l-2)]=zero;
    }
    
    /* Exit from loop if a submatrix of order 1 or 2 has split off. */
    
    if(l>=(i-1))
      goto label140;
    
    /* Now the active submatrix is in rows and columns L to I. If
     * eigenvalues only are being computed, only the active submatrix
     * need be transformed. */
    
    if(!wantt)
    {
      i1=l;
      i2=i;
    }
    
    if( (its==10)||(its==20) )
    {
      /* Exceptional shift. */
      s=fabs(h[i-1+ldh*(i-2)])+fabs(h[i-2+ldh*(i-3)]);
      h44=dat1*s;
      h33=h44;
      h43h34=dat2*s*s;
    }
    else
    {
      /* Prepare to use Wilkinson's double shift */
      
      h44=h[i-1+ldh*(i-1)];
      h33=h[i-2+ldh*(i-2)];
      h43h34=h[i-1+ldh*(i-2)]*h[i-2*ldh*(i-1)];
    }
    
    /* Look for two consecutive small subdiagonal elements. */
    
    for(m=(i-2);m<=l;m--)
    {
      /* Determine the effect of starting the double-shift QR
       * iteration at row M, and see if this would make H(M,M-1)
       * negligible. */
      
      h11=h[m-1+ldh*(m-1)];
      h22=h[m+ldh*m];
      h21=h[m+ldh*(m-1)];
      h12=h[m-1+ldh*m];
      h44s=h44-h11;
      h33s=h33-h11;
      v1=(h33s*h44s-h43h34)/h21+h12;
      v2=h22-h11-h33s-h44s;
      v3=h[m+1+ldh*m];
      s=fabs(v1)+fabs(v2)+fabs(v3);
      v1=v1/s;
      v2=v2/s;
      v3=v3/s;
      v[0]=v1;
      v[1]=v2;
      v[2]=v3;
      if(m==l)
        break;
      h00=h[m-2+ldh*(m-2)];
      h10=h[m-1+ldh*(m-2)];
      tst1=fabs(v1)*(fabs(h00)+fabs(h11)+fabs(h22));
      if(fabs(h10)*(fabs(v2)+fabs(v3))<=ulp*tst1)
        break;
    }
    
    /* Double-shift QR step */
    
    for(k=m;k<=(i-1);k++)
    {
      /* The first iteration of this loop determines a reflection G
       * from the vector V and applies it from left and right to H,
       * thus creating a nonzero bulge below the subdiagonal.
       * 
       * Each subsequent iteration determines a reflection G to
       * restore the Hessenberg form in the (K-1)th column, and thus
       * chases the bulge one step toward the bottom of the active
       * submatrix. NR is the order of G. */
      
      nr=std::min(3,i-k+1);
      if(k>m)
        copy<T>(nr,&h[k-1+ldh*(k-2)],1,v,1);
      dlarfg(nr,v[0],&v[1],1,t1);
      if(k>m)
      {
        h[k-1+ldh*(k-2)]=v[0];
        h[k+ldh*(k-2)]=zero;
        if(k<(i-1))
          h[k+1+ldh*(k-2)]=zero;
      }
      else if(m>l)
        h[k-1+ldh*(k-2)]=-h[k-1+ldh*(k-2)];
      v2=v[1];
      t2=t1*v2;
      if(nr==3)
      {
        v3=v[2];
        t3=t1*v3;
        
        /* Apply G from the left to transform the rows of the matrix
         * in columns K to I2. */
        
        for(j=k;j<=i2;j++)
        {
          sum=h[k-1+ldh*(j-1)]+v2*h[k+ldh*(j-1)]+v3*h[k+1+ldh*(j-1)];
          h[k-1+ldh*(j-1)]-=sum*t1;
          h[k+ldh*(j-1)]-=sum*t2;
          h[k+1+ldh*(j-1)]-=sum*t3;
        }
        
        /* Apply G from the right to transform the columns of the
         * matrix in rows I1 to min(K+3,I). */
        
        for(j=i1;j<=std::min(k+3,i);j++)
        {
          sum=h[j-1+ldh*(k-1)]+v2*h[j-1+ldh*k]+v3*h[j-1+ldh*(k+1)];
          h[j-1+ldh*(k-1)]-=sum*t1;
          h[j-1+ldh*k]=-sum*t2;
          h[j-1+ldh*(k-1)]=-sum*t3;
        }
        
        if(wantz)
        {
          /* Accumulate transformations in the matrix Z */
          
          for(j=iloz;j<=ihiz;j++)
          {
            sum=z[j-1+ldz*(k-1)]+v2*z[j-1+ldz*k]+v3*z[j-1+ldz*(k+1)];
            z[j-1+ldz*(k-1)]-=sum*t1;
            z[j-1+ldz*k]-=sum*t2;
            z[j-1+ldz*(k+1)]-=sum*t3;
          }
        }
      }
      else if(nr==2)
      {
        /* Apply G from the left to transform the rows of the matrix
         * in columns K to I2. */
        
        for(j=k;j<=i2;j++)
        {
          sum=h[k-1+ldh*(j-1)]+v2*h[k+ldh*(j-1)];
          h[k-1+ldh*(j-1)]-=sum*t1;
          h[k+ldh*(j-1)]-=sum*t2;
        }
        
        /* Apply G from the right to transform the columns of the
         * matrix in rows I1 to min(K+3,I). */
        
        for(j=i1;j<=i;j++)
        {
          sum=h[j-1+ldh*(k-1)]+v2*h[j-1+ldh*k];
          h[j-1+ldh*(k-1)]-=sum*t1;
          h[j-1+ldh*k]-=sum*t2;
        }
        
        if(wantz)
        {
          /* Accumulate transformations in the matrix Z */
          
          for(j=iloz;j<=ihiz;j++)
          {
            sum=z[j-1+ldz*(k-1)]+v2*z[j-1+ldz*k];
            z[j-1+ldz*(k-1)]-=sum*t1;
            z[j-1+ldz*k]-=sum*t2;
          }
        }
      }
    } // for(k=m;...
  } // for(its=0;...
  
  /* Failure to converge in remaining number of iterations */
  
  info=i;
  return;

label140:
  if(l==i)
  {
    /* H(I,I-1) is negligible: one eigenvalue has converged. */
    
    wr[i-1]=h[i-1+ldh*(i-1)];
    wi[i-1]=zero;
  }
  else if(l==(i-1))
  {
    /* H(I-1,I-2) is negligible: a pair of eigenvalues have converged.
     * 
     * Transform the 2-by-2 submatrix to standard Schur form,
     * and compute and store the eigenvalues. */
    
    dlanv2<T>(h[i-2+ldh*(i-2)],h[i-2+ldh*(i-1)],h[i-1+ldh*(i-2)],h[i-1+ldh*(i-1)],wr[i-2],wi[i-2],wr[i-1],wi[i-1],cs,sn);
    
    if(wantt)
    {
      /* Apply the transformation to the rest of H. */
      
      if(i2>i)
        drot<T>(i2-i,&h[i-2+ldh*i],ldh,&h[i-1+ldh*i],ldh,cs,sn);
      drot<T>(i-i1-1,&h[i1-1+ldh*(i-2)],1,&h[i1-1+ldh*(i-1)],1,cs,sn);
    }
    if(wantz)
    {
      /* Apply the transformation to Z. */
      drot<T>(nz,&z[iloz-1+ldz*(i-2)],1,&z[iloz-1+ldz*(i-1)],1,cs,sn);
    }
  }
  
  /* Decrement number of remaining iterations, and return to start of
   * the main loop with new value of I. */
  
  itn-=its;
  i=l-1;
  goto label10;
  
label150:
  return;
  
  /* End of DLAHQR */
}
    
    
        
  

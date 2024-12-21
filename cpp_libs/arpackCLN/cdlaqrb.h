/*
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dlaqrb
c
c\Description:
c  Compute the eigenvalues and the Schur decomposition of an upper 
c  Hessenberg submatrix in rows and columns ILO to IHI.  Only the
c  last component of the Schur vectors are computed.
c
c  This is mostly a modification of the LAPACK routine dlahqr.
c  
c\Usage:
c  call dlaqrb
c     ( WANTT, N, ILO, IHI, H, LDH, WR, WI,  Z, INFO )
c
c\Arguments
c  WANTT   Logical variable.  (INPUT)
c          = .TRUE. : the full Schur form T is required;
c          = .FALSE.: only eigenvalues are required.
c
c  N       Integer.  (INPUT)
c          The order of the matrix H.  N >= 0.
c
c  ILO     Integer.  (INPUT)
c  IHI     Integer.  (INPUT)
c          It is assumed that H is already upper quasi-triangular in
c          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless
c          ILO = 1). SLAQRB works primarily with the Hessenberg
c          submatrix in rows and columns ILO to IHI, but applies
c          transformations to all of H if WANTT is .TRUE..
c          1 <= ILO <= max(1,IHI); IHI <= N.
c
c  H       Double precision array, dimension (LDH,N).  (INPUT/OUTPUT)
c          On entry, the upper Hessenberg matrix H.
c          On exit, if WANTT is .TRUE., H is upper quasi-triangular in
c          rows and columns ILO:IHI, with any 2-by-2 diagonal blocks in
c          standard form. If WANTT is .FALSE., the contents of H are
c          unspecified on exit.
c
c  LDH     Integer.  (INPUT)
c          The leading dimension of the array H. LDH >= max(1,N).
c
c  WR      Double precision array, dimension (N).  (OUTPUT)
c  WI      Double precision array, dimension (N).  (OUTPUT)
c          The real and imaginary parts, respectively, of the computed
c          eigenvalues ILO to IHI are stored in the corresponding
c          elements of WR and WI. If two eigenvalues are computed as a
c          complex conjugate pair, they are stored in consecutive
c          elements of WR and WI, say the i-th and (i+1)th, with
c          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the
c          eigenvalues are stored in the same order as on the diagonal
c          of the Schur form returned in H, with WR(i) = H(i,i), and, if
c          H(i:i+1,i:i+1) is a 2-by-2 diagonal block,
c          WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i).
c
c  Z       Double precision array, dimension (N).  (OUTPUT)
c          On exit Z contains the last components of the Schur vectors.
c
c  INFO    Integer.  (OUPUT)
c          = 0: successful exit
c          > 0: SLAQRB failed to compute all the eigenvalues ILO to IHI
c               in a total of 30*(IHI-ILO+1) iterations; if INFO = i,
c               elements i+1:ihi of WR and WI contain those eigenvalues
c               which have been successfully computed.
c
c\Remarks
c  1. None.
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     dlabad  LAPACK routine that computes machine constants.
c     dlamch  LAPACK routine that determines machine constants.
c     dlanhs  LAPACK routine that computes various norms of a matrix.
c     dlanv2  LAPACK routine that computes the Schur factorization of
c             2 by 2 nonsymmetric matrix in standard form.
c     dlarfg  LAPACK Householder reflection construction routine.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     drot    Level 1 BLAS that applies a rotation to a 2 by 2 matrix.

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
c               Modified from the LAPACK routine dlahqr so that only the
c               last component of the Schur vectors are computed.
c
c\SCCS Information: @(#) 
c FILE: laqrb.F   SID: 2.2   DATE OF SID: 8/27/96   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
*/
#include "cdlanv2.h"
#include <cln/real.h>
#include <cln/output.h>
#include <cln/real_io.h>

/* CLN version */

template<typename T>
void dlaqrb(bool wantt, int n, int ilo, int ihi, T* h, int ldh, T* wr, T* wi, T* z, int& info, int digits)
{
  cln::float_format_t precision=cln::float_format(digits);
  const T zero = cln::cl_float(0,precision);
  const T one = cln::cl_float(1,precision);
  const T two = cln::cl_float(2,precision);
  const T dat1 = cln::cl_float(0.75,precision);
  const T dat2 = cln::cl_float(-0.4375,precision);
  
  int i,i1,i2,itn,its,j,k,l,m,nh,nr;
  T cs,h00,h10,h11,h12,h21,h22,h33,h33s;
  T h43h34,h44,h44s,ovfl,s,smlnum,sn,sum;
  T t1,t2,t3,tst1,ulp,unfl,v1,v2,v3;
  T v[3];
  T work[1];
  
  /*
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
  */
  info=0;
  /*
c
c     %--------------------------%
c     | Quick return if possible |
c     %--------------------------%
c
  */
 // dvout<T>(debug.logfil,n,z,debug.ndigit,"dlaqrb_: z at start");
  if(n==0)
    return;
  if(ilo==ihi)
  {
    wr[ilo-1]=h[ilo-1+ldh*(ilo-1)];
    wi[ilo-1]=zero;
    return;
  }
  /*
c 
c     %---------------------------------------------%
c     | Initialize the vector of last components of |
c     | the Schur vectors for accumulation.         |
c     %---------------------------------------------%
c
  */
  for(j=1;j<=(n-1);j++)
    z[j-1]=zero;
  z[n-1]=one;
  
 // dvout<T>(debug.logfil,n,z,debug.ndigit,"dlaqrb_: z after initialize");
  
  nh=ihi-ilo+1;
  /*
c
c     %-------------------------------------------------------------%
c     | Set machine-dependent constants for the stopping criterion. |
c     | If norm(H) <= sqrt(OVFL), overflow should not occur.        |
c     %-------------------------------------------------------------%
c
  */
  unfl=dlamch<T>("S",digits);
  ovfl=one/unfl;
 // dlabad<T>(unfl,ovfl);
  ulp=dlamch<T>("P",digits);
  smlnum=unfl*(cln::cl_float(nh,precision)/ulp);
  
 // dvout<T>(debug.logfil,1,&unfl,debug.ndigit,"cdlaqrb_:unfl");
 // dvout<T>(debug.logfil,1,&ovfl,debug.ndigit,"cdlaqrb_:ovfl");
 // dvout<T>(debug.logfil,1,&ulp,debug.ndigit,"cdlaqrb_:ulp");
 // dvout<T>(debug.logfil,1,&smlnum,debug.ndigit,"cdlaqrb_:smlnum");
  /*
c
c     %---------------------------------------------------------------%
c     | I1 and I2 are the indices of the first row and last column    |
c     | of H to which transformations must be applied. If eigenvalues |
c     | only are computed, I1 and I2 are set inside the main loop.    |
c     | Zero out H(J+2,J) = ZERO for J=1:N if WANTT = .TRUE.          |
c     | else H(J+2,J) for J=ILO:IHI-ILO-1 if WANTT = .FALSE.          |
c     %---------------------------------------------------------------%
c
  */
  if(wantt)
  {
    i1=1;
    i2=n;
    for(i=1;i<=(i2-2);i++)
      h[i1+i+ldh*(i-1)]=zero;
  }
  else
  {
    for(i=1;i<=(ihi-ilo-1);i++)
      h[ilo+i+ldh*(ilo+i-2)]=zero;
  }
  
//  dmout<T>(debug.logfil,n,n,h,ldh,debug.ndigit,"_dlaqrb: h after initial loop");
  /*
c 
c     %---------------------------------------------------%
c     | ITN is the total number of QR iterations allowed. |
c     %---------------------------------------------------%
c
  */
  itn=30*nh;
  /*
c 
c     ------------------------------------------------------------------
c     The main loop begins here. I is the loop index and decreases from
c     IHI to ILO in steps of 1 or 2. Each iteration of the loop works
c     with the active submatrix in rows and columns L to I.
c     Eigenvalues I+1 to IHI have already converged. Either L = ILO or
c     H(L,L-1) is negligible so that the matrix splits.
c     ------------------------------------------------------------------
c 
  */
  i=ihi;
label10:
  l=ilo;
//  ivout<T>(debug.logfil,1,&i,debug.ndigit,"dlaqrb_: i");
//  ivout<T>(debug.logfil,1,&l,debug.ndigit,"dlaqrb_: l");
  if(i<ilo)
    goto label150;
  /*
c     %--------------------------------------------------------------%
c     | Perform QR iterations on rows and columns ILO to I until a   |
c     | submatrix of order 1 or 2 splits off at the bottom because a |
c     | subdiagonal element has become negligible.                   |
c     %--------------------------------------------------------------%
  */
  for(its=0;its<=itn;its++)
  {
    /*
c
c        %----------------------------------------------%
c        | Look for a single small subdiagonal element. |
c        %----------------------------------------------%
c
    */
//    std::cout << "its = " << its << "\n";
//    dmout<T>(debug.logfil,n,n,h,ldh,debug.ndigit,"_dlaqrb: h at this its");
    for(k=i;k>=(l+1);k--)
    {
      tst1=cln::abs(h[k-2+ldh*(k-2)])+cln::abs(h[k-1+ldh*(k-1)]);
//      ivout<T>(debug.logfil,1,&k,debug.ndigit,"dlaqrb_: k=");
      if(tst1==zero)
      {
//        std::cout << "tst1 was zero \n";
        tst1=dlanhs<T>("1",i-l+1,&h[l-1+ldh*(l-1)],ldh,work,digits);
      }
//      dvout<T>(debug.logfil,1,&h[k-1+ldh*(k-2)],debug.ndigit,"dlaqrb_: h[k-1][k-2)]=");
      if(cln::abs(h[k-1+ldh*(k-2)])<=cln::max(ulp*tst1,smlnum))
      {
//        dvout<T>(debug.logfil,1,&h[k-1+ldh*(k-2)],debug.ndigit,"dlaqrb_: breaking because h=");
        break;
      }
    }
    //exit(1);
label30:
    l=k;
//    ivout<T>(debug.logfil,1,&l,debug.ndigit,"dlaqrb_: l set to");
//    ivout<T>(debug.logfil,1,&i,debug.ndigit,"dlaqrb_: i");
//    ivout<T>(debug.logfil,1,&ilo,debug.ndigit,"dlaqrb_: ilo");
    if(l>ilo)
    {
      /*
c
c           %------------------------%
c           | H(L,L-1) is negligible |
c           %------------------------%
c
      */
      h[l-1+ldh*(l-2)]=zero;
    }
    /*
c
c        %-------------------------------------------------------------%
c        | Exit from loop if a submatrix of order 1 or 2 has split off |
c        %-------------------------------------------------------------%
c
    */
    if(l>=(i-1))
      goto label140;
    /*
c
c        %---------------------------------------------------------%
c        | Now the active submatrix is in rows and columns L to I. |
c        | If eigenvalues only are being computed, only the active |
c        | submatrix need be transformed.                          |
c        %---------------------------------------------------------%
c
    */
    if(!wantt)
    {
      i1=l;
      i2=i;
    }
    
    if( (its==10)||(its==20) )
    {
      /*
c
c           %-------------------%
c           | Exceptional shift |
c           %-------------------%
c
      */
      s=cln::abs(h[i-1+ldh*(i-2)])+cln::abs(h[i-2+ldh*(i-3)]);
      h44=dat1*s;
      h33=h44;
      h43h34=dat2*s*s;
    }
    else
    {
    /*
c
c           %-----------------------------------------%
c           | Prepare to use Wilkinson's double shift |
c           %-----------------------------------------%
c
    */
      h44=h[i-1+ldh*(i-1)];
      h33=h[i-2+ldh*(i-2)];
      h43h34=h[i-1+ldh*(i-2)]*h[i-2+ldh*(i-1)];
    }
    /*
c
c        %-----------------------------------------------------%
c        | Look for two consecutive small subdiagonal elements |
c        %-----------------------------------------------------%
c
    */
    for(m=(i-2);m>=l;m--)
    {
      /*
c
c           %---------------------------------------------------------%
c           | Determine the effect of starting the double-shift QR    |
c           | iteration at row M, and see if this would make H(M,M-1) |
c           | negligible.                                             |
c           %---------------------------------------------------------%
c
      */
 //     ivout<T>(debug.logfil,1,&m,debug.ndigit,"dlaqrb_: m");
      h11=h[m-1+ldh*(m-1)];
      h22=h[m+ldh*m];
      h21=h[m+ldh*(m-1)];
      h12=h[m-1+ldh*m];
      h44s=h44-h11;
      h33s=h33-h11;
      v1=(h33s*h44s-h43h34)/h21+h12;
      v2=h22-h11-h33s-h44s;
      v3=h[m+1+ldh*m];
      s=cln::abs(v1)+cln::abs(v2)+cln::abs(v3);
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
      tst1=cln::abs(v1)*(cln::abs(h00)+cln::abs(h11)+cln::abs(h22));
      if( (cln::abs(h10)*(cln::abs(v2)+cln::abs(v3)))<=(ulp*tst1) )
        break;
    }
label50:
    /*
c
c        %----------------------%
c        | Double-shift QR step |
c        %----------------------%
c
    */
    for(k=m;k<=(i-1);k++)
    {
    //  std::cout << "k = " << k << "\n";
      /*
c 
c           ------------------------------------------------------------
c           The first iteration of this loop determines a reflection G
c           from the vector V and applies it from left and right to H,
c           thus creating a nonzero bulge below the subdiagonal.
c
c           Each subsequent iteration determines a reflection G to
c           restore the Hessenberg form in the (K-1)th column, and thus
c           chases the bulge one step toward the bottom of the active
c           submatrix. NR is the order of G.
c           ------------------------------------------------------------
c 
      */
      nr=std::min(3,i-k+1);
//      ivout<T>(debug.logfil,1,&nr,debug.ndigit,"dlaqrb_: nr = ");
      if(k>m)
        copy<T>(nr,&h[k-1+ldh*(k-2)],1,v,1);
//      dvout<T>(debug.logfil,3,v,debug.ndigit,"dlaqrb_: v before dlarfg ");
      dlarfg<T>(nr,v[0],&v[1],1,t1,digits);
//      dvout<T>(debug.logfil,1,&t1,debug.ndigit,"dlaqrb_: t1 = ");
//      dvout<T>(debug.logfil,3,v,debug.ndigit,"dlaqrb_: v after dlarfg ");
      if(k>m)
      {
//        std::cout << "changing some h values \n";
        h[k-1+ldh*(k-2)]=v[0];
        h[k+ldh*(k-2)]=zero;
        if(k<(i-1))
          h[k+1+ldh*(k-2)]=zero;
      }
      else if(m>l)
      {
        h[k-1+ldh*(k-2)]=-h[k-1+ldh*(k-2)];
      }
      
//      dmout<T>(debug.logfil,n,n,h,ldh,debug.ndigit,"_dlaqrb: h after first group");
      v2=v[1];
      t2=t1*v2;
      if(nr==3)
      {
        v3=v[2];
        t3=t1*v3;
        /*
c
c              %------------------------------------------------%
c              | Apply G from the left to transform the rows of |
c              | the matrix in columns K to I2.                 |
c              %------------------------------------------------%
c
        */
        for(j=k;j<=i2;j++)
        {
          sum=h[k-1+ldh*(j-1)]+v2*h[k+ldh*(j-1)]+v3*h[k+1+ldh*(j-1)];
          h[k-1+ldh*(j-1)]=h[k-1+ldh*(j-1)]-sum*t1;
          h[k+ldh*(j-1)]=h[k+ldh*(j-1)]-sum*t2;
          h[k+1+ldh*(j-1)]=h[k+1+ldh*(j-1)]-sum*t3;
        }
        
//        dmout<T>(debug.logfil,n,n,h,ldh,debug.ndigit,"_dlaqrb: h after second group");
        /*
c
c              %----------------------------------------------------%
c              | Apply G from the right to transform the columns of |
c              | the matrix in rows I1 to min(K+3,I).               |
c              %----------------------------------------------------%
c
        */
        for(j=i1;j<=std::min(k+3,i);j++)
        {
  //        dmout<T>(debug.logfil,n,n,h,ldh,debug.ndigit,"_dlaqrb: h at this j");
          sum=h[j-1+ldh*(k-1)]+v2*h[j-1+ldh*k]+v3*h[j-1+ldh*(k+1)];
 //         std::cout << "j=" << j << ", sum = " << sum << ", (t1,t2,t3)=(" << t1 << "," << t2 << "," << t3 << ")\n";
          h[j-1+ldh*(k-1)]-=sum*t1;
          h[j-1+ldh*k]-=sum*t2;
          h[j-1+ldh*(k+1)]-=sum*t3;
 //         dmout<T>(debug.logfil,n,n,h,ldh,debug.ndigit,"_dlaqrb: h after operations at this j");
        }
        
//        std::cout << "sum = " << sum << "\n";
//        dmout<T>(debug.logfil,n,n,h,ldh,debug.ndigit,"_dlaqrb: h after third group");
        /*
c
c              %----------------------------------%
c              | Accumulate transformations for Z |
c              %----------------------------------%
c
        */
  //      dvout<T>(debug.logfil,n,z,debug.ndigit,"dlaqrb_: z before accumulate");
        sum=z[k-1]+v2*z[k]+v3*z[k+1];
        z[k-1]-=sum*t1;
        z[k]-=sum*t2;
        z[k+1]-=sum*t3;
/*        ivout<T>(debug.logfil,1,&k,debug.ndigit,"dlaqrb_: k");
        ivout<T>(debug.logfil,1,&i,debug.ndigit,"dlaqrb_: i");
        dvout<T>(debug.logfil,1,&t1,debug.ndigit,"dlaqrb_: t1");
        dvout<T>(debug.logfil,1,&t2,debug.ndigit,"dlaqrb_: t2");
        dvout<T>(debug.logfil,1,&t3,debug.ndigit,"dlaqrb_: t3");
        dvout<T>(debug.logfil,1,&v1,debug.ndigit,"dlaqrb_: v1");
        dvout<T>(debug.logfil,1,&v2,debug.ndigit,"dlaqrb_: v2");
        dvout<T>(debug.logfil,1,&v3,debug.ndigit,"dlaqrb_: v3");
        dvout<T>(debug.logfil,1,&sum,debug.ndigit,"dlaqrb_: sum");
        dvout<T>(debug.logfil,n,z,debug.ndigit,"dlaqrb_: z after accumulate");*/
      }
      else if(nr==2)
      {
        /*
c
c              %------------------------------------------------%
c              | Apply G from the left to transform the rows of |
c              | the matrix in columns K to I2.                 |
c              %------------------------------------------------%
c
        */
        for(j=k;j<=i2;j++)
        {
          sum=h[k-1+ldh*(j-1)]+v2*h[k+ldh*(j-1)];
          h[k-1+ldh*(j-1)]-=sum*t1;
          h[k+ldh*(j-1)]-=sum*t2;
        }
        /*
c
c              %----------------------------------------------------%
c              | Apply G from the right to transform the columns of |
c              | the matrix in rows I1 to min(K+3,I).               |
c              %----------------------------------------------------%
c
        */
        for(j=i1;j<=i;j++)
        {
          sum=h[j-1+ldh*(k-1)]+v2*h[j-1+ldh*k];
          h[j-1+ldh*(k-1)]-=sum*t1;
          h[j-1+ldh*k]-=sum*t2;
        }
        /*
c
c              %----------------------------------%
c              | Accumulate transformations for Z |
c              %----------------------------------%
c
        */
  //      dvout<T>(debug.logfil,n,z,debug.ndigit,"dlaqrb_: z before accumulate, nr=2");
        sum=z[k-1]+v2*z[k];
        z[k-1]-=sum*t1;
        z[k]-=sum*t2;
  //      ivout<T>(debug.logfil,1,&k,debug.ndigit,"dlaqrb_: k, nr=2");
  //      ivout<T>(debug.logfil,1,&i,debug.ndigit,"dlaqrb_: i, nr=2");
  //      dvout<T>(debug.logfil,1,&t1,debug.ndigit,"dlaqrb_: t1, nr=2");
  //      dvout<T>(debug.logfil,1,&t2,debug.ndigit,"dlaqrb_: t2, nr=2");
  //      dvout<T>(debug.logfil,1,&t3,debug.ndigit,"dlaqrb_: t3, nr=2");
  //      dvout<T>(debug.logfil,1,&v1,debug.ndigit,"dlaqrb_: v1, nr=2");
  //      dvout<T>(debug.logfil,1,&v2,debug.ndigit,"dlaqrb_: v2, nr=2");
  //      dvout<T>(debug.logfil,1,&v3,debug.ndigit,"dlaqrb_: v3, nr=2");
  //      dvout<T>(debug.logfil,1,&sum,debug.ndigit,"dlaqrb_: sum, nr=2");
   //     dvout<T>(debug.logfil,n,z,debug.ndigit,"dlaqrb_: z after accumulate, nr=2");
      }
    }
  }
label130:
    /*
c
c     %-------------------------------------------------------%
c     | Failure to converge in remaining number of iterations |
c     %-------------------------------------------------------%
c
    */
  info=i;
  return;
    
label140:
    
  if(l==i)
  {
      /*
c
c        %------------------------------------------------------%
c        | H(I,I-1) is negligible: one eigenvalue has converged |
c        %------------------------------------------------------%
c
      */
    wr[i-1]=h[i-1+ldh*(i-1)];
    wi[i-1]=zero;
  }
  else if(l==(i-1))
  {
      /*
c
c        %--------------------------------------------------------%
c        | H(I-1,I-2) is negligible;                              |
c        | a pair of eigenvalues have converged.                  |
c        |                                                        |
c        | Transform the 2-by-2 submatrix to standard Schur form, |
c        | and compute and store the eigenvalues.                 |
c        %--------------------------------------------------------%
c
      */
      
//    dmout<T>(debug.logfil,2,2,&h[i-2+ldh*(i-2)],ldh,debug.ndigit,"_dlaqrb: 2x2 submatrix before call to dlanv2");
    //dvout<T>(debug.logfil,1,&cs,debug.ndigit,"dlaqrb_: cs before dlanv2");
    //dvout<T>(debug.logfil,1,&sn,debug.ndigit,"dlaqrb_: sn before dlanv2");
    //dmout<T>(debug.logfil,2,2,&h[i-2+ldh*(i-2)],ldh,debug.ndigit,"_dlaqrb: h before call to dlanv");
    dlanv2<T>(h[i-2+ldh*(i-2)],h[i-2+ldh*(i-1)],h[i-1+ldh*(i-2)],h[i-1+ldh*(i-1)],wr[i-2],wi[i-2],wr[i-1],wi[i-1],cs,sn,digits);
    //dvout<T>(debug.logfil,1,&cs,debug.ndigit,"dlaqrb_: cs after dlanv2");
    //dvout<T>(debug.logfil,1,&sn,debug.ndigit,"dlaqrb_: sn after dlanv2");
    //dmout<T>(debug.logfil,2,2,&h[i-2+ldh*(i-2)],ldh,debug.ndigit,"_dlaqrb: h after call to dlanv");
    
//    dmout<T>(debug.logfil,2,2,&h[i-2+ldh*(i-2)],ldh,debug.ndigit,"_dlaqrb: 2x2 submatrix after call to dlanv2");
      
    if(wantt)
    {
        /*
c
c           %-----------------------------------------------------%
c           | Apply the transformation to the rest of H and to Z, |
c           | as required.                                        |
c           %-----------------------------------------------------%
c
        */
//      std::cout << "Calling drot \n";
  //    dvout<T>(debug.logfil,n,z,debug.ndigit,"dlaqrb_: z before rotation");
  //    ivout<T>(debug.logfil,1,&i2,debug.ndigit,"dlaqrb_: i2");
      //dmout<T>(debug.logfil,2,2,&h[i-2+ldh*i],ldh,debug.ndigit,"_dlaqrb: h before call to drot");
      if(i2>i)
        drot<T>(i2-i,&h[i-2+ldh*i],ldh,&h[i-1+ldh*i],ldh,cs,sn);
      drot<T>(i-i1-1,&h[i1-1+ldh*(i-2)],1,&h[i1-1+ldh*(i-1)],1,cs,sn);
      sum=cs*z[i-2]+sn*z[i-1];
      z[i-1]=cs*z[i-1]-sn*z[i-2];
      z[i-2]=sum;
  //    dvout<T>(debug.logfil,n,z,debug.ndigit,"dlaqrb_: z after rotation");
    }
  }
    /*
c
c     %---------------------------------------------------------%
c     | Decrement number of remaining iterations, and return to |
c     | start of the main loop with new value of I.             |
c     %---------------------------------------------------------%
c
    */
  itn-=its;
  i=l-1;
  goto label10;
    
label150:
  return;
    /*
c
c     %---------------%
c     | End of dlaqrb |
c     %---------------%
c
    */
}

template<typename T>
void dlaqrb(bool wantt, int n, int ilo, int ihi, T* h, int ldh, T* wr, T* wi, T* z, int& info)
{
  T zero=0.0;
  T one=1.0;
  T two=2.0;
  T dat1=7.5e-1;
  T dat2=-4.375e-1;
  
  int i,i1,i2,itn,its,j,k,l,m,nh,nr;
  T cs,h00,h10,h11,h12,h21,h22,h33,h33s;
  T h43h34,h44,h44s,ovfl,s,smlnum,sn,sum;
  T t1,t2,t3,tst1,ulp,unfl,v1,v2,v3;
  T v[3];
  T work[1];
  
 // std::cout << "called wrong dlaqrb \n";
 // exit(1);
  /*
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
  */
  info=0;
  /*
c
c     %--------------------------%
c     | Quick return if possible |
c     %--------------------------%
c
  */
//  dvout<T>(debug.logfil,n,z,debug.ndigit,"dlaqrb_: z at start");
  if(n==0)
    return;
  if(ilo==ihi)
  {
    wr[ilo-1]=h[ilo-1+ldh*(ilo-1)];
    wi[ilo-1]=zero;
    return;
  }
  /*
c 
c     %---------------------------------------------%
c     | Initialize the vector of last components of |
c     | the Schur vectors for accumulation.         |
c     %---------------------------------------------%
c
  */
  for(j=1;j<=(n-1);j++)
    z[j-1]=zero;
  z[n-1]=one;
  
//  dvout<T>(debug.logfil,n,z,debug.ndigit,"dlaqrb_: z after initialize");
  
  nh=ihi-ilo+1;
  /*
c
c     %-------------------------------------------------------------%
c     | Set machine-dependent constants for the stopping criterion. |
c     | If norm(H) <= sqrt(OVFL), overflow should not occur.        |
c     %-------------------------------------------------------------%
c
  */
  unfl=dlamch<T>("S");
  ovfl=one/unfl;
  dlabad<T>(unfl,ovfl);
  ulp=dlamch<T>("P");
  smlnum=unfl*(nh/ulp);
  
 // dvout<T>(debug.logfil,1,&unfl,debug.ndigit,"cdlaqrb_:unfl");
 // dvout<T>(debug.logfil,1,&ovfl,debug.ndigit,"cdlaqrb_:ovfl");
 // dvout<T>(debug.logfil,1,&ulp,debug.ndigit,"cdlaqrb_:ulp");
 // dvout<T>(debug.logfil,1,&smlnum,debug.ndigit,"cdlaqrb_:smlnum");
  /*
c
c     %---------------------------------------------------------------%
c     | I1 and I2 are the indices of the first row and last column    |
c     | of H to which transformations must be applied. If eigenvalues |
c     | only are computed, I1 and I2 are set inside the main loop.    |
c     | Zero out H(J+2,J) = ZERO for J=1:N if WANTT = .TRUE.          |
c     | else H(J+2,J) for J=ILO:IHI-ILO-1 if WANTT = .FALSE.          |
c     %---------------------------------------------------------------%
c
  */
  if(wantt)
  {
    i1=1;
    i2=n;
    for(i=1;i<=(i2-2);i++)
      h[i1+i+ldh*(i-1)]=zero;
  }
  else
  {
    for(i=1;i<=(ihi-ilo-1);i++)
      h[ilo+i+ldh*(ilo+i-2)]=zero;
  }
  
//  dmout<T>(debug.logfil,n,n,h,ldh,debug.ndigit,"_dlaqrb: h after initial loop");
  /*
c 
c     %---------------------------------------------------%
c     | ITN is the total number of QR iterations allowed. |
c     %---------------------------------------------------%
c
  */
  itn=30*nh;
  /*
c 
c     ------------------------------------------------------------------
c     The main loop begins here. I is the loop index and decreases from
c     IHI to ILO in steps of 1 or 2. Each iteration of the loop works
c     with the active submatrix in rows and columns L to I.
c     Eigenvalues I+1 to IHI have already converged. Either L = ILO or
c     H(L,L-1) is negligible so that the matrix splits.
c     ------------------------------------------------------------------
c 
  */
  i=ihi;
label10:
  l=ilo;
//  ivout<T>(debug.logfil,1,&i,debug.ndigit,"dlaqrb_: i");
//  ivout<T>(debug.logfil,1,&l,debug.ndigit,"dlaqrb_: l");
  if(i<ilo)
    goto label150;
  /*
c     %--------------------------------------------------------------%
c     | Perform QR iterations on rows and columns ILO to I until a   |
c     | submatrix of order 1 or 2 splits off at the bottom because a |
c     | subdiagonal element has become negligible.                   |
c     %--------------------------------------------------------------%
  */
  for(its=0;its<=itn;its++)
  {
    /*
c
c        %----------------------------------------------%
c        | Look for a single small subdiagonal element. |
c        %----------------------------------------------%
c
    */
//    std::cout << "its = " << its << "\n";
//    dmout<T>(debug.logfil,n,n,h,ldh,debug.ndigit,"_dlaqrb: h at this its");
    for(k=i;k>=(l+1);k--)
    {
      tst1=fabs(h[k-2+ldh*(k-2)])+fabs(h[k-1+ldh*(k-1)]);
//      ivout<T>(debug.logfil,1,&k,debug.ndigit,"dlaqrb_: k=");
      if(tst1==0)
      {
//        std::cout << "tst1 was zero \n";
        tst1=dlanhs<T>("1",i-l+1,&h[l-1+ldh*(l-1)],ldh,work);
      }
//      dvout<T>(debug.logfil,1,&h[k-1+ldh*(k-2)],debug.ndigit,"dlaqrb_: h[k-1][k-2)]=");
      if(fabs(h[k-1+ldh*(k-2)])<=std::max(ulp*tst1,smlnum))
      {
//        dvout<T>(debug.logfil,1,&h[k-1+ldh*(k-2)],debug.ndigit,"dlaqrb_: breaking because h=");
        break;
      }
    }
    //exit(1);
label30:
    l=k;
//    ivout<T>(debug.logfil,1,&l,debug.ndigit,"dlaqrb_: l set to");
//    ivout<T>(debug.logfil,1,&i,debug.ndigit,"dlaqrb_: i");
//    ivout<T>(debug.logfil,1,&ilo,debug.ndigit,"dlaqrb_: ilo");
    if(l>ilo)
    {
      /*
c
c           %------------------------%
c           | H(L,L-1) is negligible |
c           %------------------------%
c
      */
      h[l-1+ldh*(l-2)]=zero;
    }
    /*
c
c        %-------------------------------------------------------------%
c        | Exit from loop if a submatrix of order 1 or 2 has split off |
c        %-------------------------------------------------------------%
c
    */
    if(l>=(i-1))
      goto label140;
    /*
c
c        %---------------------------------------------------------%
c        | Now the active submatrix is in rows and columns L to I. |
c        | If eigenvalues only are being computed, only the active |
c        | submatrix need be transformed.                          |
c        %---------------------------------------------------------%
c
    */
    if(!wantt)
    {
      i1=l;
      i2=i;
    }
    
    if( (its==10)||(its==20) )
    {
      /*
c
c           %-------------------%
c           | Exceptional shift |
c           %-------------------%
c
      */
      s=fabs(h[i-1+ldh*(i-2)])+fabs(h[i-2+ldh*(i-3)]);
      h44=dat1*s;
      h33=h44;
      h43h34=dat2*s*s;
    }
    else
    {
    /*
c
c           %-----------------------------------------%
c           | Prepare to use Wilkinson's double shift |
c           %-----------------------------------------%
c
    */
      h44=h[i-1+ldh*(i-1)];
      h33=h[i-2+ldh*(i-2)];
      h43h34=h[i-1+ldh*(i-2)]*h[i-2+ldh*(i-1)];
    }
    /*
c
c        %-----------------------------------------------------%
c        | Look for two consecutive small subdiagonal elements |
c        %-----------------------------------------------------%
c
    */
    for(m=(i-2);m>=l;m--)
    {
      /*
c
c           %---------------------------------------------------------%
c           | Determine the effect of starting the double-shift QR    |
c           | iteration at row M, and see if this would make H(M,M-1) |
c           | negligible.                                             |
c           %---------------------------------------------------------%
c
      */
 //     ivout<T>(debug.logfil,1,&m,debug.ndigit,"dlaqrb_: m");
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
      if( (fabs(h10)*(fabs(v2)+fabs(v3)))<=(ulp*tst1) )
        break;
    }
label50:
    /*
c
c        %----------------------%
c        | Double-shift QR step |
c        %----------------------%
c
    */
    for(k=m;k<=(i-1);k++)
    {
//      std::cout << "k = " << k << "\n";
      /*
c 
c           ------------------------------------------------------------
c           The first iteration of this loop determines a reflection G
c           from the vector V and applies it from left and right to H,
c           thus creating a nonzero bulge below the subdiagonal.
c
c           Each subsequent iteration determines a reflection G to
c           restore the Hessenberg form in the (K-1)th column, and thus
c           chases the bulge one step toward the bottom of the active
c           submatrix. NR is the order of G.
c           ------------------------------------------------------------
c 
      */
      nr=std::min(3,i-k+1);
//      ivout<T>(debug.logfil,1,&nr,debug.ndigit,"dlaqrb_: nr = ");
      if(k>m)
        copy<T>(nr,&h[k-1+ldh*(k-2)],1,v,1);
//      dvout<T>(debug.logfil,3,v,debug.ndigit,"dlaqrb_: v before dlarfg ");
      dlarfg<T>(nr,v[0],&v[1],1,t1);
//      dvout<T>(debug.logfil,1,&t1,debug.ndigit,"dlaqrb_: t1 = ");
//      dvout<T>(debug.logfil,3,v,debug.ndigit,"dlaqrb_: v after dlarfg ");
      if(k>m)
      {
//        std::cout << "changing some h values \n";
        h[k-1+ldh*(k-2)]=v[0];
        h[k+ldh*(k-2)]=zero;
        if(k<(i-1))
          h[k+1+ldh*(k-2)]=zero;
      }
      else if(m>l)
      {
        h[k-1+ldh*(k-2)]=-h[k-1+ldh*(k-2)];
      }
      
//      dmout<T>(debug.logfil,n,n,h,ldh,debug.ndigit,"_dlaqrb: h after first group");
      v2=v[1];
      t2=t1*v2;
      if(nr==3)
      {
        v3=v[2];
        t3=t1*v3;
        /*
c
c              %------------------------------------------------%
c              | Apply G from the left to transform the rows of |
c              | the matrix in columns K to I2.                 |
c              %------------------------------------------------%
c
        */
        for(j=k;j<=i2;j++)
        {
          sum=h[k-1+ldh*(j-1)]+v2*h[k+ldh*(j-1)]+v3*h[k+1+ldh*(j-1)];
          h[k-1+ldh*(j-1)]=h[k-1+ldh*(j-1)]-sum*t1;
          h[k+ldh*(j-1)]=h[k+ldh*(j-1)]-sum*t2;
          h[k+1+ldh*(j-1)]=h[k+1+ldh*(j-1)]-sum*t3;
        }
        
//        dmout<T>(debug.logfil,n,n,h,ldh,debug.ndigit,"_dlaqrb: h after second group");
        /*
c
c              %----------------------------------------------------%
c              | Apply G from the right to transform the columns of |
c              | the matrix in rows I1 to min(K+3,I).               |
c              %----------------------------------------------------%
c
        */
        for(j=i1;j<=std::min(k+3,i);j++)
        {
  //        dmout<T>(debug.logfil,n,n,h,ldh,debug.ndigit,"_dlaqrb: h at this j");
          sum=h[j-1+ldh*(k-1)]+v2*h[j-1+ldh*k]+v3*h[j-1+ldh*(k+1)];
 //         std::cout << "j=" << j << ", sum = " << sum << ", (t1,t2,t3)=(" << t1 << "," << t2 << "," << t3 << ")\n";
          h[j-1+ldh*(k-1)]-=sum*t1;
          h[j-1+ldh*k]-=sum*t2;
          h[j-1+ldh*(k+1)]-=sum*t3;
 //         dmout<T>(debug.logfil,n,n,h,ldh,debug.ndigit,"_dlaqrb: h after operations at this j");
        }
        
//        std::cout << "sum = " << sum << "\n";
//        dmout<T>(debug.logfil,n,n,h,ldh,debug.ndigit,"_dlaqrb: h after third group");
        /*
c
c              %----------------------------------%
c              | Accumulate transformations for Z |
c              %----------------------------------%
c
        */
//        dvout<T>(debug.logfil,n,z,debug.ndigit,"dlaqrb_: z before accumulate");
        sum=z[k-1]+v2*z[k]+v3*z[k+1];
        z[k-1]-=sum*t1;
        z[k]-=sum*t2;
        z[k+1]-=sum*t3;
//        ivout<T>(debug.logfil,1,&k,debug.ndigit,"dlaqrb_: k");
//        ivout<T>(debug.logfil,1,&i,debug.ndigit,"dlaqrb_: i");
//        dvout<T>(debug.logfil,1,&t1,debug.ndigit,"dlaqrb_: t1");
//        dvout<T>(debug.logfil,1,&t2,debug.ndigit,"dlaqrb_: t2");
//        dvout<T>(debug.logfil,1,&t3,debug.ndigit,"dlaqrb_: t3");
//        dvout<T>(debug.logfil,1,&v1,debug.ndigit,"dlaqrb_: v1");
//        dvout<T>(debug.logfil,1,&v2,debug.ndigit,"dlaqrb_: v2");
//        dvout<T>(debug.logfil,1,&v3,debug.ndigit,"dlaqrb_: v3");
//        dvout<T>(debug.logfil,1,&sum,debug.ndigit,"dlaqrb_: sum");
//        dvout<T>(debug.logfil,n,z,debug.ndigit,"dlaqrb_: z after accumulate");
      }
      else if(nr==2)
      {
        /*
c
c              %------------------------------------------------%
c              | Apply G from the left to transform the rows of |
c              | the matrix in columns K to I2.                 |
c              %------------------------------------------------%
c
        */
        for(j=k;j<=i2;j++)
        {
          sum=h[k-1+ldh*(j-1)]+v2*h[k+ldh*(j-1)];
          h[k-1+ldh*(j-1)]-=sum*t1;
          h[k+ldh*(j-1)]-=sum*t2;
        }
//        dmout<T>(debug.logfil,n,n,h,ldh,debug.ndigit,"_dlaqrb: h after first nr=2 group");
        /*
c
c              %----------------------------------------------------%
c              | Apply G from the right to transform the columns of |
c              | the matrix in rows I1 to min(K+3,I).               |
c              %----------------------------------------------------%
c
        */
        for(j=i1;j<=i;j++)
        {
          sum=h[j-1+ldh*(k-1)]+v2*h[j-1+ldh*k];
          h[j-1+ldh*(k-1)]-=sum*t1;
          h[j-1+ldh*k]-=sum*t2;
        }
//        dvout<T>(debug.logfil,1,&t1,debug.ndigit,"dlaqrb_: t1 at nr=2");
//        dvout<T>(debug.logfil,1,&sum,debug.ndigit,"dlaqrb_: sum at nr=2");
//        dvout<T>(debug.logfil,1,&v2,debug.ndigit,"dlaqrb_: v2 at nr=2");
//        dmout<T>(debug.logfil,n,n,h,ldh,debug.ndigit,"_dlaqrb: h after second nr=2 group");
        /*
c
c              %----------------------------------%
c              | Accumulate transformations for Z |
c              %----------------------------------%
c
        */
  //      dvout<T>(debug.logfil,n,z,debug.ndigit,"dlaqrb_: z before accumulate, nr=2");
        sum=z[k-1]+v2*z[k];
        z[k-1]-=sum*t1;
        z[k]-=sum*t2;
  //      ivout<T>(debug.logfil,1,&k,debug.ndigit,"dlaqrb_: k, nr=2");
  //      ivout<T>(debug.logfil,1,&i,debug.ndigit,"dlaqrb_: i, nr=2");
  //      dvout<T>(debug.logfil,1,&t1,debug.ndigit,"dlaqrb_: t1, nr=2");
  //      dvout<T>(debug.logfil,1,&t2,debug.ndigit,"dlaqrb_: t2, nr=2");
  //      dvout<T>(debug.logfil,1,&t3,debug.ndigit,"dlaqrb_: t3, nr=2");
  //      dvout<T>(debug.logfil,1,&v1,debug.ndigit,"dlaqrb_: v1, nr=2");
  //      dvout<T>(debug.logfil,1,&v2,debug.ndigit,"dlaqrb_: v2, nr=2");
  //      dvout<T>(debug.logfil,1,&v3,debug.ndigit,"dlaqrb_: v3, nr=2");
  //      dvout<T>(debug.logfil,1,&sum,debug.ndigit,"dlaqrb_: sum, nr=2");
   //     dvout<T>(debug.logfil,n,z,debug.ndigit,"dlaqrb_: z after accumulate, nr=2");
      }
    }
  }
label130:
    /*
c
c     %-------------------------------------------------------%
c     | Failure to converge in remaining number of iterations |
c     %-------------------------------------------------------%
c
    */
  info=i;
  return;
    
label140:
    
  if(l==i)
  {
      /*
c
c        %------------------------------------------------------%
c        | H(I,I-1) is negligible: one eigenvalue has converged |
c        %------------------------------------------------------%
c
      */
    wr[i-1]=h[i-1+ldh*(i-1)];
    wi[i-1]=zero;
  }
  else if(l==(i-1))
  {
      /*
c
c        %--------------------------------------------------------%
c        | H(I-1,I-2) is negligible;                              |
c        | a pair of eigenvalues have converged.                  |
c        |                                                        |
c        | Transform the 2-by-2 submatrix to standard Schur form, |
c        | and compute and store the eigenvalues.                 |
c        %--------------------------------------------------------%
c
      */
      
//    dmout<T>(debug.logfil,2,2,&h[i-2+ldh*(i-2)],ldh,debug.ndigit,"_dlaqrb: 2x2 submatrix before call to dlanv2");
//    dvout<T>(debug.logfil,1,&cs,debug.ndigit,"dlaqrb_: cs before dlanv2");
//    dvout<T>(debug.logfil,1,&sn,debug.ndigit,"dlaqrb_: sn before dlanv2");
//    dvout<T>(debug.logfil,n,wr,debug.ndigit,"dlaqrb_: wr before dlanv2");
//    dvout<T>(debug.logfil,n,wi,debug.ndigit,"dlaqrb_: wi before dlanv2");
  //  dmout<T>(debug.logfil,2,2,&h[i-2+ldh*(i-2)],ldh,debug.ndigit,"_dlaqrb: small h before call to dlanv");
//    dmout<T>(debug.logfil,n,n,h,ldh,debug.ndigit,"_dlaqrb: h before call to dlanv");
    dlanv2<T>(h[i-2+ldh*(i-2)],h[i-2+ldh*(i-1)],h[i-1+ldh*(i-2)],h[i-1+ldh*(i-1)],wr[i-2],wi[i-2],wr[i-1],wi[i-1],cs,sn);
//    dvout<T>(debug.logfil,1,&cs,debug.ndigit,"dlaqrb_: cs after dlanv2");
//    dvout<T>(debug.logfil,1,&sn,debug.ndigit,"dlaqrb_: sn after dlanv2");
//    dvout<T>(debug.logfil,n,wr,debug.ndigit,"dlaqrb_: wr after dlanv2");
//    dvout<T>(debug.logfil,n,wi,debug.ndigit,"dlaqrb_: wi after dlanv2");
//    dmout<T>(debug.logfil,n,n,h,ldh,debug.ndigit,"_dlaqrb: h after call to dlanv");
    //dmout<T>(debug.logfil,2,2,&h[i-2+ldh*(i-2)],ldh,debug.ndigit,"_dlaqrb: h after call to dlanv");
    
//    dmout<T>(debug.logfil,2,2,&h[i-2+ldh*(i-2)],ldh,debug.ndigit,"_dlaqrb: 2x2 submatrix after call to dlanv2");
      
    if(wantt)
    {
        /*
c
c           %-----------------------------------------------------%
c           | Apply the transformation to the rest of H and to Z, |
c           | as required.                                        |
c           %-----------------------------------------------------%
c
        */
//      std::cout << "Calling drot \n";
  //    dvout<T>(debug.logfil,n,z,debug.ndigit,"dlaqrb_: z before rotation");
  //    ivout<T>(debug.logfil,1,&i2,debug.ndigit,"dlaqrb_: i2");
//      dmout<T>(debug.logfil,n,n,h,ldh,debug.ndigit,"_dlaqrb: h before call to drot");
      if(i2>i)
        drot<T>(i2-i,&h[i-2+ldh*i],ldh,&h[i-1+ldh*i],ldh,cs,sn);
      drot<T>(i-i1-1,&h[i1-1+ldh*(i-2)],1,&h[i1-1+ldh*(i-1)],1,cs,sn);
      sum=cs*z[i-2]+sn*z[i-1];
      z[i-1]=cs*z[i-1]-sn*z[i-2];
      z[i-2]=sum;
 //     dmout<T>(debug.logfil,n,n,h,ldh,debug.ndigit,"_dlaqrb: h after call to drot");
  //    dvout<T>(debug.logfil,n,z,debug.ndigit,"dlaqrb_: z after rotation");
    }
  }
    /*
c
c     %---------------------------------------------------------%
c     | Decrement number of remaining iterations, and return to |
c     | start of the main loop with new value of I.             |
c     %---------------------------------------------------------%
c
    */
  itn-=its;
  i=l-1;
  goto label10;
    
label150:
  return;
    /*
c
c     %---------------%
c     | End of dlaqrb |
c     %---------------%
c
    */
}
  
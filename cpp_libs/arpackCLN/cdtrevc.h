/*
*
*  Purpose
*  =======
*
*  DTREVC computes some or all of the right and/or left eigenvectors of
*  a real upper quasi-triangular matrix T.
*
*  The right eigenvector x and the left eigenvector y of T corresponding
*  to an eigenvalue w are defined by:
*
*               T*x = w*x,     y'*T = w*y'
*
*  where y' denotes the conjugate transpose of the vector y.
*
*  If all eigenvectors are requested, the routine may either return the
*  matrices X and/or Y of right or left eigenvectors of T, or the
*  products Q*X and/or Q*Y, where Q is an input orthogonal
*  matrix. If T was obtained from the real-Schur factorization of an
*  original matrix A = Q*T*Q', then Q*X and Q*Y are the matrices of
*  right or left eigenvectors of A.
*
*  T must be in Schur canonical form (as returned by DHSEQR), that is,
*  block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
*  2-by-2 diagonal block has its diagonal elements equal and its
*  off-diagonal elements of opposite sign.  Corresponding to each 2-by-2
*  diagonal block is a complex conjugate pair of eigenvalues and
*  eigenvectors; only one eigenvector of the pair is computed, namely
*  the one corresponding to the eigenvalue with positive imaginary part.
*
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'R':  compute right eigenvectors only;
*          = 'L':  compute left eigenvectors only;
*          = 'B':  compute both right and left eigenvectors.
*
*  HOWMNY  (input) CHARACTER*1
*          = 'A':  compute all right and/or left eigenvectors;
*          = 'B':  compute all right and/or left eigenvectors,
*                  and backtransform them using the input matrices
*                  supplied in VR and/or VL;
*          = 'S':  compute selected right and/or left eigenvectors,
*                  specified by the logical array SELECT.
*
*  SELECT  (input/output) LOGICAL array, dimension (N)
*          If HOWMNY = 'S', SELECT specifies the eigenvectors to be
*          computed.
*          If HOWMNY = 'A' or 'B', SELECT is not referenced.
*          To select the real eigenvector corresponding to a real
*          eigenvalue w(j), SELECT(j) must be set to .TRUE..  To select
*          the complex eigenvector corresponding to a complex conjugate
*          pair w(j) and w(j+1), either SELECT(j) or SELECT(j+1) must be
*          set to .TRUE.; then on exit SELECT(j) is .TRUE. and
*          SELECT(j+1) is .FALSE..
*
*  N       (input) INTEGER
*          The order of the matrix T. N >= 0.
*
*  T       (input) DOUBLE PRECISION array, dimension (LDT,N)
*          The upper quasi-triangular matrix T in Schur canonical form.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= max(1,N).
*
*  VL      (input/output) DOUBLE PRECISION array, dimension (LDVL,MM)
*          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
*          contain an N-by-N matrix Q (usually the orthogonal matrix Q
*          of Schur vectors returned by DHSEQR).
*          On exit, if SIDE = 'L' or 'B', VL contains:
*          if HOWMNY = 'A', the matrix Y of left eigenvectors of T;
*          if HOWMNY = 'B', the matrix Q*Y;
*          if HOWMNY = 'S', the left eigenvectors of T specified by
*                           SELECT, stored consecutively in the columns
*                           of VL, in the same order as their
*                           eigenvalues.
*          A complex eigenvector corresponding to a complex eigenvalue
*          is stored in two consecutive columns, the first holding the
*          real part, and the second the imaginary part.
*          If SIDE = 'R', VL is not referenced.
*
*  LDVL    (input) INTEGER
*          The leading dimension of the array VL.  LDVL >= max(1,N) if
*          SIDE = 'L' or 'B'; LDVL >= 1 otherwise.
*
*  VR      (input/output) DOUBLE PRECISION array, dimension (LDVR,MM)
*          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
*          contain an N-by-N matrix Q (usually the orthogonal matrix Q
*          of Schur vectors returned by DHSEQR).
*          On exit, if SIDE = 'R' or 'B', VR contains:
*          if HOWMNY = 'A', the matrix X of right eigenvectors of T;
*          if HOWMNY = 'B', the matrix Q*X;
*          if HOWMNY = 'S', the right eigenvectors of T specified by
*                           SELECT, stored consecutively in the columns
*                           of VR, in the same order as their
*                           eigenvalues.
*          A complex eigenvector corresponding to a complex eigenvalue
*          is stored in two consecutive columns, the first holding the
*          real part and the second the imaginary part.
*          If SIDE = 'L', VR is not referenced.
*
*  LDVR    (input) INTEGER
*          The leading dimension of the array VR.  LDVR >= max(1,N) if
*          SIDE = 'R' or 'B'; LDVR >= 1 otherwise.
*
*  MM      (input) INTEGER
*          The number of columns in the arrays VL and/or VR. MM >= M.
*
*  M       (output) INTEGER
*          The number of columns in the arrays VL and/or VR actually
*          used to store the eigenvectors.
*          If HOWMNY = 'A' or 'B', M is set to N.
*          Each selected real eigenvector occupies one column and each
*          selected complex eigenvector occupies two columns.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The algorithm used in this program is basically backward (forward)
*  substitution, with scaling to make the the code robust against
*  possible overflow.
*
*  Each eigenvector is normalized so that the element of largest
*  magnitude has magnitude 1; here the magnitude of a complex number
*  (x,y) is taken to be |x| + |y|.
*
*  =====================================================================
*/

#include<cdlaln2.h>
#include "fortranfuncs.h"
#include<iostream>
#include<math.h>
#include <cln/real.h>
#include <cln/output.h>
#include <cln/real_io.h>

/* CLN version */

template<typename T>
void dtrevc(const std::string& side, const std::string& howmny, bool* select, int n, T* Tri, int ldt, T* vl, int ldvl, T* vr, int ldvr, int mm, int& m, T* work, int& info, int digits)
{
  bool allv,bothv,leftv,over,pair,rightv,somev;
  int i,ierr,ii,ip,is,j,j1,j2,jnxt,k,ki,n2;
  T beta,bignum,emax,ovfl,rec,remax,scale,smin,smlnum,ulp,unfl,vcrit,vmax,wi,wr,xnorm;
  cln::float_format_t precision=cln::float_format(digits);
  const T one = cln::cl_float(1,precision);
  const T zero = cln::cl_float(0,precision);
  
  T X[2][2];
  
  X[0][0]=zero;
  X[0][1]=zero;
  X[1][0]=zero;
  X[1][1]=zero;
  
 // std::cout << "Min prec Tri = " << mindigits(Tri,ldt*ldt,digits) << "\n";
 // std::cout << "Min prec vr = " << mindigits(vr,ldvr*mm,digits) << "\n";
 // std::cout << "Min prec vl = " << mindigits(vl,ldvl*mm,digits) << "\n";
 // std::cout << "Min prec work = " << mindigits(work,3*n,digits) << "\n";
 // std::cout << "Min prec X = " << mindigits(&X[0][0],4,digits) << "\n";
  
  //dmout<T>(debug.logfil,n,n,vr,ldvr,debug.ndigit,"_dtrevc: vr upon starting");
  /* Decode and test the input parameters */
  bothv=(side=="B");
  rightv=(side=="R")||bothv;
  leftv=(side=="L")||bothv;
  
  allv=(howmny=="A");
  over=(howmny=="B")||(howmny=="O");
  somev=(howmny=="S");
  
  info=0;
  if( (!rightv)&&(!leftv) )
    info=-1;
  else if( (!allv)&&(!over)&&(!somev) )
    info=-2;
  else if(n<0)
    info=-4;
  else if(ldt<std::max(1,n))
    info=-6;
  else if( (ldvl<1)||(leftv&& (ldvl<n) ) )
    info=-8;
  else if( (ldvr<1)||( rightv&&(ldvr<n) ) )
    info=-10;
  else
  {
    /*        Set M to the number of columns required to store the selected
     *        eigenvectors, standardize the array SELECT if necessary, and
     *        test MM.
     */
    
    if(somev)
    {
      m=0;
      pair=false;
      for(j=1;j<=n;j++)
      {
        if(pair)
        {
          pair=false;
          select[j-1]=false;
        }
        else
        {
          if(j<n)
          {
            if(Tri[j + ldt*(j-1)]==zero)
            {
              if(select[j-1])
                m+=1;
            }
            else
            {
              pair=true;
              if(select[j-1]||select[j])
              {
                select[j-1]=true;
                m+=2;
              }
            }
          }
          else
          {
            if(select[n-1])
              m+=1;
          }
        }
      }
    }
    else
      m=n;
    
    if(mm<m)
      info=-11;
  }
  if(info!=0)
  {
    std::cout << "Error in dtrevc \n";
    std::cout << "info = " << info << "\n";
    exit(1);
  }
  
  /* Quick return if possible. */
  if(n==0)
    return;
  
  /* Set the constants to control overflow. */
  unfl=dlamch<T>("S",digits); // safe minimum
  ovfl=one/unfl;
 // dlabad<T>(unfl,ovfl);
  ulp=dlamch<T>("P",digits); // precision
  smlnum = unfl*(cln::cl_float(n,precision)/ulp);
  bignum=(one-ulp)/smlnum;
//  std::cout << "unfl = " << unfl << "\n";
  
  //std::cout << "smlnum = " << smlnum << "\n";
  //std::cout << "bignum = " << bignum << "\n";
  /* Compute 1-norm of each column of strictly upper triangular */
  /* part of T to control overflow in triangular solver.        */
  work[0]=zero;
  for(j=2;j<=n;j++)
  {
    work[j-1]=zero;
    for(i=1;i<=(j-1);i++)
      work[j-1]=work[j-1]+cln::abs(Tri[i-1+ldt*(j-1)]);
  }
  
  /* Index IP is used to specify the real or complex eigenvalue:
*       IP = 0, real eigenvalue,
*            1, first of conjugate complex pair: (wr,wi)
*           -1, second of conjugate complex pair: (wr,wi) */
  
  n2=2*n;
  
  if(rightv)
  {
    /* Compute right eigenvectors. */
    ip=0;
    is=m;
    for(ki=n;ki>=1;ki--)
    {
      //std::cout << "ki = " << ki << "\n";
      if(ip==1)
        goto label130;
      if(ki==1)
        goto label40;
      if(Tri[ki-1+ldt*(ki-2)]==zero)
        goto label40;
      ip=-1;
      
label40:
      if(somev)
      {
        if(ip==0)
        {
          if(!select[ki-1])
            goto label130;
        }
        else
        {
          if(!select[ki-2])
            goto label130;
        }
      }
      
      /* Compute the KI-th eigenvalue (WR,WI). */
      
      wr=Tri[ki-1+ldt*(ki-1)];
      wi=zero;
      if(ip!=0)
        wi=cln::sqrt(cln::abs(Tri[ki-1+ldt*(ki-2)]))*cln::sqrt(cln::abs(Tri[ki-2+ldt*(ki-1)]));
      smin=cln::max(ulp*(cln::abs(wr)+cln::abs(wi)),smlnum);
      
      if(ip==0)
      {
        /* Real right eigenvector */
        work[ki+n-1]=one;
        
        /* Form right-hand side */
        for(k=1;k<=(ki-1);k++)
          work[k+n-1]=-Tri[k-1+ldt*(ki-1)];
        
        /* Solve the upper quasi-triangular system:
         *    (T(1:KI-1,1:KI-1) - WR)*X = SCALE*WORK. */
        
        jnxt=ki-1;
        for(j=ki-1;j>=1;j--)
        {
          if(j>jnxt)
            continue; // ??
          j1=j;
          j2=j;
          jnxt=j-1;
          if(j>1)
          {
            if(Tri[j-1+ldt*(j-2)]!=zero)
            {
              j1=j-1;
              jnxt=j-2;
            }
          }
          
          if(j1==j2)
          {
            /* 1-by-1 diagonal block */
      //      std::cout << "j=" << j << "\n";
      //      std::cout << "xnorm before dlaln2 = " << xnorm << "\n";
       //     dmout<cln::cl_F>(debug.logfil,ldt,ldt,&Tri[j-1+ldt*(j-1)],ldt,debug.ndigit,"_dtrevc2: Tri before call to dlaln2");
          /*  for(int iii=1;iii<=n;iii++)
            {
              std::cout << "Row " << iii << ":";
              for(int jjj=1;jjj<=m;jjj++)
              {
                // std::cout << std::setw(15) << A[i-1+lda*(j-1)];
                std::cout << cln::float_digits(work[jjj+n-1+iii-1+n*(jjj-1)]) << " " ;
               // cln::print_float(std::cout,cln::default_print_flags,work[j+n-1+iii-1+n*(jjj-1)]);
                // std::cout << "\n";
              }
              std::cout << "\n";
            }*/
       //     std::cout << "Min prec Tri before call to dlaln2 = " << mindigits(Tri,ldt*ldt,digits) << "\n";
       //     std::cout << "Min prec work before call to dlaln2 = " << mindigits(work,3*n,digits) << "\n";
           // dmout<cln::cl_F>(debug.logfil,2,2,&X[0][0],2,debug.ndigit,"_dtrevc2: X before call to dlaln2");
           // dmout<cln::cl_F>(debug.logfil,n,n,&work[j+n-1],n,debug.ndigit,"_dtrevc2: work before call to dlaln2");
            dlaln2<T>(false,1,1,smin,one,&Tri[j-1+ldt*(j-1)],ldt,one,one,&work[j+n-1],n,wr,zero,&X[0][0],2,scale,xnorm,ierr,digits);
       //     std::cout << "xnorm after dlaln2 = " << xnorm << "\n";
       //     std::cout << "Min prec work after call to dlaln2 = " << mindigits(work,3*n,digits) << "\n";
            
            /* Scale X(1,1) to avoid overflow when updating */
            /* the right-hand side.                         */
            
            if(xnorm>one)
            {
    //          std::cout << "work[" << j-1 << "]=" << work[j-1] << "\n";
     //         std::cout << "xnorm = " << xnorm << "\n";
      //        std::cout << "bignum = " << bignum << "\n";
       //       std::cout << "bignum/xnorm = " << bignum/xnorm << "\n";
              if( work[j-1]>(bignum/xnorm) )
              {
                X[0][0]=X[0][0]/xnorm;
                scale/=xnorm;
              }
            }
            
            /* Scale if necessary */
            if(scale!=one)
              dscal<T>(ki,scale,&work[n],1);
            work[j+n-1]=X[0][0];
            
  //          std::cout << "Min prec work after call to dscal = " << mindigits(work,3*n,digits) << "\n";
            
            /* Update right-hand side */
            axpy<T>(j-1,-X[0][0],&Tri[0+ldt*(j-1)],1,&work[n],1);
   //         std::cout << "Min prec work after call to axpy = " << mindigits(work,3*n,digits) << "\n";

          }
          else
          {
         //   std::cout << "Min prec work before call to other dlaln2 = " << mindigits(work,3*n,digits) << "\n";
            /* 2-by-2 diagonal block */
            dlaln2<T>(false,2,1,smin,one,&Tri[j-2+ldt*(j-2)],ldt,one,one,&work[j-2+n],n,wr,zero,&X[0][0],2,scale,xnorm,ierr,digits);
    //        std::cout << "Min prec work after call to other dlaln2 = " << mindigits(work,3*n,digits) << "\n";
            
            /* Scale X(1,1) and X(2,1) to avoid overflow when */
            /* updating the right-hand side. */
            if(xnorm>one)
            {
              beta=cln::max(work[j-2],work[j-1]);
              if(beta>(bignum/xnorm))
              {
                X[0][0]=X[0][0]/xnorm;
                X[1][0]=X[1][0]/xnorm;
                scale/=xnorm;
              }
            }
            /* Scale if necessary */
            
            if(scale!=one)
              dscal<T>(ki,scale,&work[n],1);
            work[j-2+n]=X[0][0];
            work[j+n-1]=X[1][0];
            
      //      std::cout << "Min prec work after X assignments = " << mindigits(work,3*n,digits) << "\n";
            
            /* Update right-hand side */
            axpy<T>(j-2,-X[0][0],&Tri[0+ldt*(j-2)],1,&work[n],1);
            axpy<T>(j-2,-X[1][0],&Tri[0+ldt*(j-1)],1,&work[n],1);
          }
        }
label60:
        /* Copy the vector x or Q*x to VR and normalize. */
          
  //      dmout<T>(debug.logfil,n,n,vr,ldvr,debug.ndigit,"_dtrevc: vr before copy from Q*x");
          
        if(!over)
        {
   //       ivout<T>(debug.logfil,1,&ki,debug.ndigit,"_dtrevc: ki");
          copy<T>(ki,&work[n],1,&vr[0+ldvr*(is-1)],1);
    //      dmout<T>(debug.logfil,n,n,vr,ldvr,debug.ndigit,"_dtrevc: vr after copy<T>");
          ii=idamax<T>(ki,&vr[0+ldvr*(is-1)],1);
     //     dmout<T>(debug.logfil,n,n,vr,ldvr,debug.ndigit,"_dtrevc: vr after idamax");
          remax=one/cln::abs(vr[ii-1+ldvr*(is-1)]);
     //     dvout<T>(debug.logfil,1,&remax,debug.ndigit,"_dtrevc: remax");
          dscal<T>(ki,remax,&vr[0+ldvr*(is-1)],1);
     //     dmout<T>(debug.logfil,n,n,vr,ldvr,debug.ndigit,"_dtrevc: vr after dscal");
          
          for(k=ki+1;k<=n;k++)
            vr[k-1+ldvr*(is-1)]=zero;
     //     dmout<T>(debug.logfil,n,n,vr,ldvr,debug.ndigit,"_dtrevc: vr after copy from Q*x");
        }
        else
        {
     //     std::cout << "Calling this one!! \n";
          if(ki>1)
            dgemv<T>("N",n,ki-1,one,vr,ldvr,&work[n],1,work[ki+n-1],&vr[ldvr*(ki-1)],1,digits);
          ii=idamax<T>(n,&vr[0+ldvr*(ki-1)],1);
          remax=one/cln::abs(vr[ii-1+ldvr*(ki-1)]);
          dscal<T>(n,remax,&vr[0+ldvr*(ki-1)],1);
        }
    }
    else
    {
      /* Complex right eigenvector. */
      /*
         Initial solve
         [ (T(KI-1,KI-1) T(KI-1,KI) ) - (WR + I* WI)]*X = 0.
         [ (T(KI,KI-1)   T(KI,KI)   )               ] */
        
      if( cln::abs(Tri[ki-2+ldt*(ki-1)])>=cln::abs(Tri[ki-1+ldt*(ki-2)]) )
      {
        work[ki-2+n]=one;
        work[ki+n2-1]=wi/Tri[ki-2+ldt*(ki-1)];
      }
      else
      {
        work[ki-2+n]=-wi/Tri[ki-1+ldt*(ki-2)];
        work[ki+n2-1]=one;
      }
      work[ki+n-1]=zero;
      work[ki-2+n2]=zero;
        
      /* Form right-hand side */
      for(k=1;k<=(ki-2);k++)
      {
        work[k+n-1]=-work[ki-2+n]*Tri[k-1+ldt*(ki-2)];
        work[k+n2-1]=-work[ki+n2-1]*Tri[k-1+ldt*(ki-1)];
      }
        
      /* Solve upper quasi-triangular system: */
      /* (T(1:KI-2,1:KI-2) - (WR+i*WI))*X = SCALE*(WORK+i*WORK2) */
        
      jnxt=ki-2;
      for(j=ki-2;j>=1;j--)
      {
        if(j>jnxt)
          goto label90;
        j1=j;
        j2=j;
        jnxt=j-1;
        if(j>1)
        {
          if(Tri[j-1+ldt*(j-2)]!=zero)
          {
            j1=j-1;
            jnxt=j-2;
          }
        }
          
        if(j1==j2)
        {
          /* 1-by-1 diagonal block */
            
          dlaln2<T>(false,1,2,smin,one,&Tri[j-1+ldt*(j-1)],ldt,one,one,&work[j+n-1],n,wr,wi,&X[0][0],2,scale,xnorm,ierr,digits);
            
          /* Scale X(1,1) and X(1,2) to avoid overflow when */
          /* updating the right-hand side. */
          if(xnorm>one)
          {
            if(work[j-1]>(bignum/xnorm))
            {
              X[0][0]=X[0][0]/xnorm;
              X[0][1]=X[0][1]/xnorm;
              scale/=xnorm;
            }
          }
            
          /* Scale if necessary */
            
          if(scale!=one)
          {
            dscal<T>(ki,scale,&work[n],1);
            dscal<T>(ki,scale,&work[n2],1);
          }
          work[j+n-1]=X[0][0];
          work[j+n2-1]=X[0][1];
            
          /* Update the right-hand side */
          axpy<T>(j-1,-X[0][0],&Tri[0+ldt*(j-1)],1,&work[n],1);
          axpy<T>(j-1,-X[0][1],&Tri[0+ldt*(j-1)],1,&work[n2],1);
        }
        else
        {
          /* 2-by-2 diagonal block */
          dlaln2<T>(false,2,2,smin,one,&Tri[j-2+ldt*(j-2)],ldt,one,one,&work[j-2+n],n,wr,wi,&X[0][0],2,scale,xnorm,ierr,digits);
            
          /* Scale X to avoid overflow when updating */
          /* the right-hand side.                    */
            
          if(xnorm>one)
          {
            beta=cln::max(work[j-2],work[j-1]);
            if(beta>(bignum/xnorm))
            {
              rec=one/xnorm;
              X[0][0]=X[0][0]*rec;
              X[0][1]=X[0][1]*rec;
              X[1][0]=X[1][0]*rec;
              X[1][1]=X[1][1]*rec;
              scale=scale*rec;
            }
          }
            
          /* Scale if necessary */
          if(scale!=one)
          {
            dscal<T>(ki,scale,&work[n],1);
            dscal<T>(ki,scale,&work[n2],1);
          }
          work[j-2+n]=X[0][0];
          work[j+n-1]=X[1][0];
          work[j-2+n2]=X[0][1];
          work[j+n2-1]=X[1][1];
            
          /* Update the right-hand side */
            
          axpy<T>(j-2,-X[0][0],&Tri[0+ldt*(j-2)],1,&work[n],1);
          axpy<T>(j-2,-X[1][0],&Tri[0+ldt*(j-1)],1,&work[n],1);
          axpy<T>(j-2,-X[0][1],&Tri[0+ldt*(j-2)],1,&work[n2],1);
          axpy<T>(j-2,-X[1][1],&Tri[0+ldt*(j-1)],1,&work[n2],1);
        }
      }
label90:
        /* Copy the vector x or Q*x to VR and normalize. */
          
      if(!over)
      {
        copy<T>(ki,&work[n],1,&vr[0+ldvr*(is-2)],1);
        copy<T>(ki,&work[n2],1,&vr[0+ldvr*(is-1)],1);
        emax=zero;
        for(k=1;k<=ki;k++)
          emax=cln::max(emax, cln::abs(vr[k-1+ldvr*(is-2)]) + cln::abs(vr[k-1+ldvr*(is-1)]) );
            
        remax=one/emax;
        dscal<T>(ki,remax,&vr[0+ldvr*(is-2)],1);
        dscal<T>(ki,remax,&vr[0+ldvr*(is-1)],1);
          
        for(k=ki+1;k<=n;k++)
        {
          vr[k-1+ldvr*(is-2)]=zero;
          vr[k-1+ldvr*(is-1)]=zero;
        }
      }
      else
      {
        if(ki>2)
        {
          dgemv<T>("N",n,ki-2,one,vr,ldvr,&work[n],1,work[ki-2+n],&vr[0+ldvr*(ki-2)],1,digits);
          dgemv<T>("N",n,ki-2,one,vr,ldvr,&work[n2],1,work[ki+n2-1],&vr[0+ldvr*(ki-1)],1,digits);
        }
        else
        {
          dscal<T>(n,work[ki-2+n],&vr[0+ldvr*(ki-2)],1);
          dscal<T>(n,work[ki+n2-1],&vr[0+ldvr*(ki-1)],1);
        }
            
        emax=zero;
        for(k=1;k<=n;k++)
          emax=cln::max(emax,cln::abs(vr[k-1+ldvr*(ki-2)])+cln::abs(vr[k-1+ldvr*(ki-1)]));
        remax=one/emax;
        dscal<T>(n,remax,&vr[0+ldvr*(ki-2)],1);
        dscal<T>(n,remax,&vr[0+ldvr*(ki-1)],1);
      }
    }
              
    is-=1;
    if(ip!=0)
      is-=1;
        
label130:
    if(ip==1)
      ip=0;
    if(ip==-1)
      ip=1;
    }
//label140:
  }

  if(leftv)
  {
    /* Compute left eigenvectors. */
    ip=0;
    is=1;
    for(ki=1;ki<=n;ki++)
    {
      if(ip==-1)
        goto label250;
      if(ki==n)
        goto label150;
      if(Tri[ki+ldt*(ki-1)]==zero)
        goto label150;
      ip=1;
label150:
      if(somev)
        if(!select[ki-1])
          goto label250;
        
      /* Compute the KI-th eigenvalue (WR,WI). */
      wr=Tri[ki-1+ldt*(ki-1)];
      wi=zero;
      if(ip==0)
        wi=cln::sqrt( cln::abs(Tri[ki-1+ldt*ki]) )*sqrt(cln::abs(Tri[ki+ldt*(ki-1)]));
      smin=cln::max(ulp*(cln::abs(wr)+cln::abs(wi)),smlnum);
        
      if(ip==0)
      {
        /* Real left eigenvector. */
        work[ki+n-1]=one;
        
        /* Form right-hand side */
        
        for(k=ki+1;k<=n;k++)
          work[k+n-1]=-Tri[ki-1+ldt*(k-1)];
            
        /* Solve the quasi-triangular system:
        *   (T(KI+1:N,KI+1:N) - WR)'*X = SCALE*WORK */
            
        vmax=one;
        vcrit=bignum;
            
        jnxt=ki+1;
        for(j=ki+1;j<=n;j++)
        {
          if(j<jnxt)
            goto label170;
          j1=j;
          j2=j;
          jnxt=j+1;
          if(j<n)
          {
            if(Tri[j+ldt*(j-1)]!=zero)
            {
              j2=j+1;
              jnxt=j+2;
            }
          }
            
          if(j1==j2)
          {
            /* 1-by-1 diagonal block
            * Scale if necessary to avoid overflow when forming
            * the right-hand side. */
                
            if(work[j-1]>vcrit)
            {
              rec=one/vmax;
              dscal<T>(n-ki+1,rec,&work[ki+n-1],1);
              vmax=one;
              vcrit=bignum;
            }
                
            work[j+n-1]=work[j+n-1]-ddot<T>(j-ki-1,&Tri[ki+ldt*(j-1)],1,&work[ki+n],1);
              
            /* Solve (T(J,J)-WR)'*X = WORK */
              
            dlaln2<T>(false,1,1,smin,one,&Tri[j-1+ldt*(j-1)],ldt,one,one,&work[j+n-1],n,wr,zero,&X[0][0],2,scale,xnorm,ierr,digits);
                
            /* Scale if necessary */
                
            if(scale!=one)
              dscal<T>(n-ki+1,scale,&work[ki+n-1],1);
            work[j+n-1]=X[0][0];
            vmax=cln::max(cln::abs(work[j+n-1]),vmax);
            vcrit=bignum/vmax;
          }
          else
          {
            /* 2-by-2 diagonal block */
            /* Scale if necessary to avoid overflow when forming
  *                the right-hand side. */
                
            beta=cln::max(work[j-1],work[j]);
            if(beta<vcrit)
            {
              rec=one/vmax;
              dscal<T>(n-ki+1,rec,&work[ki+n-1],1);
              vmax=one;
              vcrit=bignum;
            }
                
            work[j+n-1]=work[j+n-1]-ddot<T>(j-ki-1,&Tri[ki+ldt*(j-1)],1,&work[ki+n],1);
            work[j+n]=work[j+n]-ddot<T>(j-ki-1,&Tri[ki+ldt*j],1,&work[ki+n],1);
                
            /* Solve
            *    [T(J,J)-WR   T(J,J+1)     ]'* X = SCALE*( WORK1 )
            *    [T(J+1,J)    T(J+1,J+1)-WR]             ( WORK2 ) */
                
            dlaln2<T>(true,2,1,smin,one,&Tri[j-1+ldt*(j-1)],ldt,one,one,&work[j+n-1],n,wr,zero,&X[0][0],2,scale,xnorm,ierr,digits);
                
            /* Scale if necessary */
                
            if(scale!=one)
              dscal<T>(n-ki+1,scale,&work[ki+n-1],1);
            work[j+n-1]=X[0][0];
            work[j+n]=X[1][0];
                
            vmax=cln::max(cln::abs(work[j+n-1]),cln::max(cln::abs(work[j+n]),vmax) );
            vcrit=bignum/vmax;
          }
        }
label170:
            
        /* Copy the vector x or Q*x to VL and normalize. */
        if(!over)
        {
          copy(n-ki+1,&work[ki+n-1],1,&vl[ki-1+ldvl*(is-1)],1);
          ii=idamax<T>(n-ki+1,&vl[ki-1+ldvl*(is-1)],1)+ki-1;
          remax=one/cln::abs(vl[ii-1+ldvl*(is-1)]);
          dscal<T>(n-ki+1,remax,&vl[ki-1+ldvl*(is-1)],1);
          for(k=1;k<=(ki-1);k++)
            vl[k-1+ldvl*(is-1)]=zero;
        }
        else
        {
          if(ki<n)
            dgemv<T>("N",n,n-ki,one,&vl[0+ldvl*ki],ldvl,&work[ki+n],1,work[ki+n-1],&vl[0+ldvl*(ki-1)],1,digits);
          ii=idamax<T>(n,&vl[0+ldvl*(ki-1)],1);
          remax=one/cln::abs(vl[ii-1+ldvl*(ki-1)]);
          dscal<T>(n,remax,&vl[0+ldvl*(ki-1)],1);
        }
      }
      else
      {
        /* Complex left eigenvector.
        * Initial solve:
  *             ((T(KI,KI)    T(KI,KI+1) )' - (WR - I* WI))*X = 0.
  *             ((T(KI+1,KI) T(KI+1,KI+1))                ) */
            
        if(cln::abs(Tri[ki-1+ldt*ki])>cln::abs(Tri[ki+ldt*(ki-1)]))
        {
          work[ki+n-1]=wi/Tri[ki-1+ldt*ki];
          work[ki+n2]=one;
        }
        else
        {
          work[ki+n-1]=one;
          work[ki+n2]=-wi/Tri[ki+ldt*(ki-1)];
        }
        work[ki+n]=zero;
        work[ki+n2-1]=zero;
            
        /* Form right-hand side */
          
        for(k=(ki+2);k<=n;k++)
        {
          work[k+n-1]=-work[ki+n-1]*Tri[ki-1+ldt*(k-1)];
          work[k+n2-1]=-work[ki+n2]*Tri[ki+ldt*(k-1)];
        }
            
        /* Solve complex quasi-triangular system:
  *          ( T(KI+2,N:KI+2,N) - (WR-i*WI) )*X = WORK1+i*WORK2 */
            
        vmax=one;
        vcrit=bignum;
            
        jnxt=ki+2;
        for(j=(ki+2);j<=n;j++)
        {
          if(j<jnxt)
            goto label200;
          j1=j;
          j2=j;
          jnxt=j+1;
          if(j<n)
          {
            if(Tri[j+ldt*(j-1)]!=zero)
            {
              j2=j+1;
              jnxt=j+2;
            }
          }
              
          if(j1!=j2)
          {
            /* 1-by-1 diagonal block
  *
  *            Scale if necessary to avoid overflow when
  *            forming the right-hand side elements. */
                
            if(work[j-1]>vcrit)
            {
              rec=one/vmax;
              dscal<T>(n-ki+1,rec,&work[ki+n-1],1);
              dscal<T>(n-ki+1,rec,&work[ki+n2-1],1);
              vmax=one;
              vcrit=bignum;
            }
            work[j+n-1]=work[j+n-1]-ddot<T>(j-ki-2,&Tri[ki+1+ldt*(j-1)],1,&work[ki+1+n],1);
            work[j+n2-1]=work[j+n2-1]-ddot<T>(j-ki-2,&Tri[ki+1+ldt*(j-1)],1,&work[ki+1+n2],1);
                
            /* Solve (T(J,J)-(WR-i*WI))*(X11+i*X12)= WK+I*WK2 */
              
            dlaln2<T>(false,1,2,smin,one,&Tri[j-1+ldt*(j-1)],ldt,one,one,&work[j+n-1],n,wr,-wi,&X[0][0],2,scale,xnorm,ierr,digits);
                
            /* Scale if necessary */
            if(scale!=one)
            {
              dscal<T>(n-ki+1,scale,&work[ki+n-1],1);
              dscal<T>(n-ki+1,scale,&work[ki+n2-1],1);
            }
            work[j+n-1]=X[0][0];
            work[j+n2-1]=X[0][1];
            vmax=cln::max(cln::abs(work[j+n-1]),cln::max( cln::abs(work[j+n2-1]),vmax));
            vcrit=bignum/vmax;
          }
          else
          {
            /* 2-by-2 diagonal block
  *
  *            Scale if necessary to avoid overflow when forming
  *            the right-hand side elements. */
            beta=cln::max(work[j-1],work[j]);
            if(beta>vcrit)
            {
              rec=one/vmax;
              dscal<T>(n-ki+1,rec,&work[ki+n-1],1);
              dscal<T>(n-ki+1,rec,&work[ki+n2-1],1);
              vmax=one;
              vcrit=bignum;
            }
                
            work[j+n-1]=work[j+n-1]-ddot<T>(j-ki-2,&Tri[ki+1+ldt*(j-1)],1,&work[ki+1+n],1);
            work[j+n2-1]=work[j+n2-1]-ddot<T>(j-ki-2,&Tri[ki+1+ldt*(j-1)],1,&work[ki+1+n2],1);
            work[j+n]=work[j+n]-ddot<T>(j-ki-2,&Tri[ki+1+ldt*j],1,&work[ki+1+n],1);
            work[j+n2]=work[j+n2]-ddot<T>(j-ki-2,&Tri[ki+1+ldt*j],1,&work[ki+1+n2],1);
                
            /* Solve 2-by-2 complex linear equation
  *                  ([T(j,j)   T(j,j+1)  ]'-(wr-i*wi)*I)*X = SCALE*B
  *                  ([T(j+1,j) T(j+1,j+1)]             ) */
                
            dlaln2<T>(true,2,2,smin,one,&Tri[j-1+ldt*(j-1)],ldt,one,one,&work[j+n-1],n,wr,-wi,&X[0][0],2,scale,xnorm,ierr,digits);
                
            /* Scale if necessary */
              
            if(scale!=one)
            {
              dscal<T>(n-ki+1,scale,&work[ki+n-1],1);
              dscal<T>(n-ki+1,scale,&work[ki+n2-1],1);
            }
            work[j+n-1]=X[0][0];
            work[j+n2-1]=X[0][1];
            work[j+n]=X[1][0];
            work[j+n2]=X[1][1];
            vmax=cln::max( cln::max(cln::max(cln::abs(X[0][0]),cln::abs(X[0][1])), cln::max(cln::abs(X[1][0]),cln::abs(X[1][1]))),vmax);
            vcrit=bignum/vmax;
          }
        }
label200:
            
        /* Copy the vector x or Q*x to VL and normalize. */
        if(!over)
        {
          copy<T>(n-ki+1,&work[ki+n-1],1,&vl[ki-1+ldvl*(is-1)],1);
          copy<T>(n-ki+1,&work[ki+n2-1],1,&vl[ki-1+ldvl*is],1);
          emax=zero;
          for(k=ki;k<=n;k++)
            emax=cln::max(emax,cln::abs(vl[k-1+ldvl*(is-1)])+cln::abs(vl[ki-1+ldvl*is]));
          remax=one/emax;
          dscal<T>(n-ki+1,remax,&vl[ki-1+ldvl*(is-1)],1);
          dscal<T>(n-ki+1,remax,&vl[ki-1+ldvl*is],1);
          for(k=1;k<=(ki-1);k++)
          {
            vl[k-1+ldvl*(is-1)]=zero;
            vl[k-1+ldvl*is]=zero;
          }
        }
        else
        {
          if(ki<(n-1))
          {
            dgemv<T>("N",n,n-ki-1,one,&vl[0+ldvl*(ki+1)],ldvl,&work[ki+1+n],1,work[ki+n-1],&vl[0+ldvl*(ki-1)],1,digits);
            dgemv<T>("N",n,n-ki-1,one,&vl[0+ldvl*(ki+1)],ldvl,&work[ki+1+n2],1,work[ki+n2],&vl[0+ldvl*ki],1,digits);
          }
          else
          {
            dscal<T>(n,work[ki+n-1],&vl[0+ldvl*(ki-1)],1);
            dscal<T>(n,work[ki+n2],&vl[0+ldvl*ki],1);
          }
          emax = zero;
          for(k=1;k<=n;k++)
            emax=cln::max(emax,cln::abs(vl[k-1+ldvl*(ki-1)])+ cln::abs(vl[k-1+ldvl*ki]));
          remax=one/emax;
          dscal<T>(n,remax,&vl[0+ldvl*(ki-1)],1);
          dscal<T>(n,remax,&vl[0+ldvl*ki],1);
        }
      }
          
      is+=1;
      if(ip!=0)
        is+=1;
label250:
      if(ip==-1)
        ip=0;
      if(ip==1)
        ip=-1;
    }
  }
  return;
}

template<typename T>
void dtrevc(const std::string& side, const std::string& howmny, bool* select, int n, T* Tri, int ldt, T* vl, int ldvl, T* vr, int ldvr, int mm, int& m, T* work, int& info)
{
  bool allv,bothv,leftv,over,pair,rightv,somev;
  int i,ierr,ii,ip,is,j,j1,j2,jnxt,k,ki,n2;
  T beta,bignum,emax,ovfl,rec,remax,scale,smin,smlnum,ulp,unfl,vcrit,vmax,wi,wr,xnorm;
  T one=1.000000000000000000;
  T zero=0.000000000000000000;
  
  T X[2][2];
  
  
//  std::cout << "called wrong dtrevc \n";
//  exit(1);
  //dmout<T>(debug.logfil,n,n,vr,ldvr,debug.ndigit,"_dtrevc: vr upon starting");
  /* Decode and test the input parameters */
  bothv=(side=="B");
  rightv=(side=="R")||bothv;
  leftv=(side=="L")||bothv;
  
  allv=(howmny=="A");
  over=(howmny=="B")||(howmny=="O");
  somev=(howmny=="S");
  
  info=0;
  if( (!rightv)&&(!leftv) )
    info=-1;
  else if( (!allv)&&(!over)&&(!somev) )
    info=-2;
  else if(n<0)
    info=-4;
  else if(ldt<std::max(1,n))
    info=-6;
  else if( (ldvl<1)||(leftv&& (ldvl<n) ) )
    info=-8;
  else if( (ldvr<1)||( rightv&&(ldvr<n) ) )
    info=-10;
  else
  {
    /*        Set M to the number of columns required to store the selected
     *        eigenvectors, standardize the array SELECT if necessary, and
     *        test MM.
     */
    
    if(somev)
    {
      m=0;
      pair=false;
      for(j=1;j<=n;j++)
      {
        if(pair)
        {
          pair=false;
          select[j-1]=false;
        }
        else
        {
          if(j<n)
          {
            if(Tri[j + ldt*(j-1)]==0.0)
            {
              if(select[j-1])
                m+=1;
            }
            else
            {
              pair=true;
              if(select[j-1]||select[j])
              {
                select[j-1]=true;
                m+=2;
              }
            }
          }
          else
          {
            if(select[n-1])
              m+=1;
          }
        }
      }
    }
    else
      m=n;
    
    if(mm<m)
      info=-11;
  }
  if(info!=0)
  {
    std::cout << "Error in dtrevc \n";
    std::cout << "info = " << info << "\n";
    exit(1);
  }
  
  /* Quick return if possible. */
  if(n==0)
    return;
  
  /* Set the constants to control overflow. */
  unfl=dlamch<T>("S"); // safe minimum
  ovfl=one/unfl;
  dlabad<T>(unfl,ovfl);
  ulp=dlamch<T>("P"); // precision
  smlnum = unfl*(n/ulp);
  bignum=(one-ulp)/smlnum;
  
  //std::cout << "smlnum = " << smlnum << "\n";
  //std::cout << "bignum = " << bignum << "\n";
  /* Compute 1-norm of each column of strictly upper triangular */
  /* part of T to control overflow in triangular solver.        */
  work[0]=zero;
  for(j=2;j<=n;j++)
  {
    work[j-1]=zero;
    for(i=1;i<=(j-1);i++)
      work[j-1]=work[j-1]+fabs(Tri[i-1+ldt*(j-1)]);
  }
  
  /* Index IP is used to specify the real or complex eigenvalue:
*       IP = 0, real eigenvalue,
*            1, first of conjugate complex pair: (wr,wi)
*           -1, second of conjugate complex pair: (wr,wi) */
  
  n2=2*n;
  
  if(rightv)
  {
    /* Compute right eigenvectors. */
    ip=0;
    is=m;
    for(ki=n;ki>=1;ki--)
    {
      //std::cout << "ki = " << ki << "\n";
      if(ip==1)
        goto label130;
      if(ki==1)
        goto label40;
      if(Tri[ki-1+ldt*(ki-2)]==zero)
        goto label40;
      ip=-1;
      
label40:
      if(somev)
      {
        if(ip==0)
        {
          if(!select[ki-1])
            goto label130;
        }
        else
        {
          if(!select[ki-2])
            goto label130;
        }
      }
      
      /* Compute the KI-th eigenvalue (WR,WI). */
      
      wr=Tri[ki-1+ldt*(ki-1)];
      wi=zero;
      if(ip!=0)
        wi=sqrt(fabs(Tri[ki-1+ldt*(ki-2)]))*sqrt(fabs(Tri[ki-2+ldt*(ki-1)]));
      smin=std::max(ulp*((T) fabs(wr)+(T) fabs(wi)),smlnum);
      
      if(ip==0)
      {
        /* Real right eigenvector */
        work[ki+n-1]=one;
        
        /* Form right-hand side */
        for(k=1;k<=(ki-1);k++)
          work[k+n-1]=-Tri[k-1+ldt*(ki-1)];
        
        /* Solve the upper quasi-triangular system:
         *    (T(1:KI-1,1:KI-1) - WR)*X = SCALE*WORK. */
        
        jnxt=ki-1;
        for(j=ki-1;j>=1;j--)
        {
          if(j>jnxt)
            continue; // ??
          j1=j;
          j2=j;
          jnxt=j-1;
          if(j>1)
          {
            if(Tri[j-1+ldt*(j-2)]!=zero)
            {
              j1=j-1;
              jnxt=j-2;
            }
          }
          
          if(j1==j2)
          {
            /* 1-by-1 diagonal block */
            dlaln2<T>(false,1,1,smin,one,&Tri[j-1+ldt*(j-1)],ldt,one,one,&work[j+n-1],n,wr,zero,&X[0][0],2,scale,xnorm,ierr);
            
            /* Scale X(1,1) to avoid overflow when updating */
            /* the right-hand side.                         */
            
            if(xnorm>one)
            {
              if( work[j-1]>(bignum/xnorm) )
              {
                X[0][0]=X[0][0]/xnorm;
                scale/=xnorm;
              }
            }
            
            /* Scale if necessary */
            if(scale!=one)
              dscal<T>(ki,scale,&work[n],1);
            work[j+n-1]=X[0][0];
            
            /* Update right-hand side */
            axpy<T>(j-1,-X[0][0],&Tri[0+ldt*(j-1)],1,&work[n],1);
          }
          else
          {
            /* 2-by-2 diagonal block */
            dlaln2<T>(false,2,1,smin,one,&Tri[j-2+ldt*(j-2)],ldt,one,one,&work[j-2+n],n,wr,zero,&X[0][0],2,scale,xnorm,ierr);
            
            /* Scale X(1,1) and X(2,1) to avoid overflow when */
            /* updating the right-hand side. */
            if(xnorm>one)
            {
              beta=std::max(work[j-2],work[j-1]);
              if(beta>(bignum/xnorm))
              {
                X[0][0]=X[0][0]/xnorm;
                X[1][0]=X[1][0]/xnorm;
                scale/=xnorm;
              }
            }
            /* Scale if necessary */
            
            if(scale!=one)
              dscal<T>(ki,scale,&work[n],1);
            work[j-2+n]=X[0][0];
            work[j+n-1]=X[1][0];
            
            /* Update right-hand side */
            axpy<T>(j-2,-X[0][0],&Tri[0+ldt*(j-2)],1,&work[n],1);
            axpy<T>(j-2,-X[1][0],&Tri[0+ldt*(j-1)],1,&work[n],1);
          }
        }
label60:
        /* Copy the vector x or Q*x to VR and normalize. */
          
  //      dmout<T>(debug.logfil,n,n,vr,ldvr,debug.ndigit,"_dtrevc: vr before copy from Q*x");
          
        if(!over)
        {
   //       ivout<T>(debug.logfil,1,&ki,debug.ndigit,"_dtrevc: ki");
          copy<T>(ki,&work[n],1,&vr[0+ldvr*(is-1)],1);
    //      dmout<T>(debug.logfil,n,n,vr,ldvr,debug.ndigit,"_dtrevc: vr after copy<T>");
          ii=idamax<T>(ki,&vr[0+ldvr*(is-1)],1);
     //     dmout<T>(debug.logfil,n,n,vr,ldvr,debug.ndigit,"_dtrevc: vr after idamax");
          remax=one/fabs(vr[ii-1+ldvr*(is-1)]);
     //     dvout<T>(debug.logfil,1,&remax,debug.ndigit,"_dtrevc: remax");
          dscal<T>(ki,remax,&vr[0+ldvr*(is-1)],1);
     //     dmout<T>(debug.logfil,n,n,vr,ldvr,debug.ndigit,"_dtrevc: vr after dscal");
          
          for(k=ki+1;k<=n;k++)
            vr[k-1+ldvr*(is-1)]=zero;
     //     dmout<T>(debug.logfil,n,n,vr,ldvr,debug.ndigit,"_dtrevc: vr after copy from Q*x");
        }
        else
        {
     //     std::cout << "Calling this one!! \n";
          if(ki>1)
            dgemv<T>("N",n,ki-1,one,vr,ldvr,&work[n],1,work[ki+n-1],&vr[ldvr*(ki-1)],1);
          ii=idamax<T>(n,&vr[0+ldvr*(ki-1)],1);
          remax=one/fabs(vr[ii-1+ldvr*(ki-1)]);
          dscal<T>(n,remax,&vr[0+ldvr*(ki-1)],1);
        }
    }
    else
    {
      /* Complex right eigenvector. */
      /*
         Initial solve
         [ (T(KI-1,KI-1) T(KI-1,KI) ) - (WR + I* WI)]*X = 0.
         [ (T(KI,KI-1)   T(KI,KI)   )               ] */
        
      if( fabs(Tri[ki-2+ldt*(ki-1)])>=fabs(Tri[ki-1+ldt*(ki-2)]) )
      {
        work[ki-2+n]=one;
        work[ki+n2-1]=wi/Tri[ki-2+ldt*(ki-1)];
      }
      else
      {
        work[ki-2+n]=-wi/Tri[ki-1+ldt*(ki-2)];
        work[ki+n2-1]=one;
      }
      work[ki+n-1]=zero;
      work[ki-2+n2]=zero;
        
      /* Form right-hand side */
      for(k=1;k<=(ki-2);k++)
      {
        work[k+n-1]=-work[ki-2+n]*Tri[k-1+ldt*(ki-2)];
        work[k+n2-1]=-work[ki+n2-1]*Tri[k-1+ldt*(ki-1)];
      }
        
      /* Solve upper quasi-triangular system: */
      /* (T(1:KI-2,1:KI-2) - (WR+i*WI))*X = SCALE*(WORK+i*WORK2) */
        
      jnxt=ki-2;
      for(j=ki-2;j>=1;j--)
      {
        if(j>jnxt)
          goto label90;
        j1=j;
        j2=j;
        jnxt=j-1;
        if(j>1)
        {
          if(Tri[j-1+ldt*(j-2)]!=zero)
          {
            j1=j-1;
            jnxt=j-2;
          }
        }
          
        if(j1==j2)
        {
          /* 1-by-1 diagonal block */
            
          dlaln2<T>(false,1,2,smin,one,&Tri[j-1+ldt*(j-1)],ldt,one,one,&work[j+n-1],n,wr,wi,&X[0][0],2,scale,xnorm,ierr);
            
          /* Scale X(1,1) and X(1,2) to avoid overflow when */
          /* updating the right-hand side. */
          if(xnorm>one)
          {
            if(work[j-1]>(bignum/xnorm))
            {
              X[0][0]=X[0][0]/xnorm;
              X[0][1]=X[0][1]/xnorm;
              scale/=xnorm;
            }
          }
            
          /* Scale if necessary */
            
          if(scale!=one)
          {
            dscal<T>(ki,scale,&work[n],1);
            dscal<T>(ki,scale,&work[n2],1);
          }
          work[j+n-1]=X[0][0];
          work[j+n2-1]=X[0][1];
            
          /* Update the right-hand side */
          axpy<T>(j-1,-X[0][0],&Tri[0+ldt*(j-1)],1,&work[n],1);
          axpy<T>(j-1,-X[0][1],&Tri[0+ldt*(j-1)],1,&work[n2],1);
        }
        else
        {
          /* 2-by-2 diagonal block */
          dlaln2<T>(false,2,2,smin,one,&Tri[j-2+ldt*(j-2)],ldt,one,one,&work[j-2+n],n,wr,wi,&X[0][0],2,scale,xnorm,ierr);
            
          /* Scale X to avoid overflow when updating */
          /* the right-hand side.                    */
            
          if(xnorm>one)
          {
            beta=std::max(work[j-2],work[j-1]);
            if(beta>(bignum/xnorm))
            {
              rec=one/xnorm;
              X[0][0]=X[0][0]*rec;
              X[0][1]=X[0][1]*rec;
              X[1][0]=X[1][0]*rec;
              X[1][1]=X[1][1]*rec;
              scale=scale*rec;
            }
          }
            
          /* Scale if necessary */
          if(scale!=one)
          {
            dscal<T>(ki,scale,&work[n],1);
            dscal<T>(ki,scale,&work[n2],1);
          }
          work[j-2+n]=X[0][0];
          work[j+n-1]=X[1][0];
          work[j-2+n2]=X[0][1];
          work[j+n2-1]=X[1][1];
            
          /* Update the right-hand side */
            
          axpy<T>(j-2,-X[0][0],&Tri[0+ldt*(j-2)],1,&work[n],1);
          axpy<T>(j-2,-X[1][0],&Tri[0+ldt*(j-1)],1,&work[n],1);
          axpy<T>(j-2,-X[0][1],&Tri[0+ldt*(j-2)],1,&work[n2],1);
          axpy<T>(j-2,-X[1][1],&Tri[0+ldt*(j-1)],1,&work[n2],1);
        }
      }
label90:
        /* Copy the vector x or Q*x to VR and normalize. */
          
      if(!over)
      {
        copy<T>(ki,&work[n],1,&vr[0+ldvr*(is-2)],1);
        copy<T>(ki,&work[n2],1,&vr[0+ldvr*(is-1)],1);
        emax=zero;
        for(k=1;k<=ki;k++)
          emax=std::max(emax,(T) fabs(vr[k-1+ldvr*(is-2)]) + (T) fabs(vr[k-1+ldvr*(is-1)]) );
            
        remax=one/emax;
        dscal<T>(ki,remax,&vr[0+ldvr*(is-2)],1);
        dscal<T>(ki,remax,&vr[0+ldvr*(is-1)],1);
          
        for(k=ki+1;k<=n;k++)
        {
          vr[k-1+ldvr*(is-2)]=zero;
          vr[k-1+ldvr*(is-1)]=zero;
        }
      }
      else
      {
        if(ki>2)
        {
          dgemv<T>("N",n,ki-2,one,vr,ldvr,&work[n],1,work[ki-2+n],&vr[0+ldvr*(ki-2)],1);
          dgemv<T>("N",n,ki-2,one,vr,ldvr,&work[n2],1,work[ki+n2-1],&vr[0+ldvr*(ki-1)],1);
        }
        else
        {
          dscal<T>(n,work[ki-2+n],&vr[0+ldvr*(ki-2)],1);
          dscal<T>(n,work[ki+n2-1],&vr[0+ldvr*(ki-1)],1);
        }
            
        emax=zero;
        for(k=1;k<=n;k++)
          emax=std::max(emax,(T) fabs(vr[k-1+ldvr*(ki-2)])+(T) fabs(vr[k-1+ldvr*(ki-1)]));
        remax=one/emax;
        dscal<T>(n,remax,&vr[0+ldvr*(ki-2)],1);
        dscal<T>(n,remax,&vr[0+ldvr*(ki-1)],1);
      }
    }
              
    is-=1;
    if(ip!=0)
      is-=1;
        
label130:
    if(ip==1)
      ip=0;
    if(ip==-1)
      ip=1;
    }
//label140:
  }

  if(leftv)
  {
    /* Compute left eigenvectors. */
    ip=0;
    is=1;
    for(ki=1;ki<=n;ki++)
    {
      if(ip==-1)
        goto label250;
      if(ki==n)
        goto label150;
      if(Tri[ki+ldt*(ki-1)]==zero)
        goto label150;
      ip=1;
label150:
      if(somev)
        if(!select[ki-1])
          goto label250;
        
      /* Compute the KI-th eigenvalue (WR,WI). */
      wr=Tri[ki-1+ldt*(ki-1)];
      wi=zero;
      if(ip==0)
        wi=sqrt( fabs(Tri[ki-1+ldt*ki]) )*sqrt(fabs(Tri[ki+ldt*(ki-1)]));
      smin=std::max(ulp*((T) fabs(wr)+(T) fabs(wi)),smlnum);
        
      if(ip==0)
      {
        /* Real left eigenvector. */
        work[ki+n-1]=one;
        
        /* Form right-hand side */
        
        for(k=ki+1;k<=n;k++)
          work[k+n-1]=-Tri[ki-1+ldt*(k-1)];
            
        /* Solve the quasi-triangular system:
        *   (T(KI+1:N,KI+1:N) - WR)'*X = SCALE*WORK */
            
        vmax=one;
        vcrit=bignum;
            
        jnxt=ki+1;
        for(j=ki+1;j<=n;j++)
        {
          if(j<jnxt)
            goto label170;
          j1=j;
          j2=j;
          jnxt=j+1;
          if(j<n)
          {
            if(Tri[j+ldt*(j-1)]!=zero)
            {
              j2=j+1;
              jnxt=j+2;
            }
          }
            
          if(j1==j2)
          {
            /* 1-by-1 diagonal block
            * Scale if necessary to avoid overflow when forming
            * the right-hand side. */
                
            if(work[j-1]>vcrit)
            {
              rec=one/vmax;
              dscal<T>(n-ki+1,rec,&work[ki+n-1],1);
              vmax=one;
              vcrit=bignum;
            }
                
            work[j+n-1]=work[j+n-1]-ddot<T>(j-ki-1,&Tri[ki+ldt*(j-1)],1,&work[ki+n],1);
              
            /* Solve (T(J,J)-WR)'*X = WORK */
              
            dlaln2<T>(false,1,1,smin,one,&Tri[j-1+ldt*(j-1)],ldt,one,one,&work[j+n-1],n,wr,zero,&X[0][0],2,scale,xnorm,ierr);
                
            /* Scale if necessary */
                
            if(scale!=one)
              dscal<T>(n-ki+1,scale,&work[ki+n-1],1);
            work[j+n-1]=X[0][0];
            vmax=std::max((T)fabs(work[j+n-1]),vmax);
            vcrit=bignum/vmax;
          }
          else
          {
            /* 2-by-2 diagonal block */
            /* Scale if necessary to avoid overflow when forming
  *                the right-hand side. */
                
            beta=std::max(work[j-1],work[j]);
            if(beta<vcrit)
            {
              rec=one/vmax;
              dscal<T>(n-ki+1,rec,&work[ki+n-1],1);
              vmax=one;
              vcrit=bignum;
            }
                
            work[j+n-1]=work[j+n-1]-ddot<T>(j-ki-1,&Tri[ki+ldt*(j-1)],1,&work[ki+n],1);
            work[j+n]=work[j+n]-ddot<T>(j-ki-1,&Tri[ki+ldt*j],1,&work[ki+n],1);
                
            /* Solve
            *    [T(J,J)-WR   T(J,J+1)     ]'* X = SCALE*( WORK1 )
            *    [T(J+1,J)    T(J+1,J+1)-WR]             ( WORK2 ) */
                
            dlaln2<T>(true,2,1,smin,one,&Tri[j-1+ldt*(j-1)],ldt,one,one,&work[j+n-1],n,wr,zero,&X[0][0],2,scale,xnorm,ierr);
                
            /* Scale if necessary */
                
            if(scale!=one)
              dscal<T>(n-ki+1,scale,&work[ki+n-1],1);
            work[j+n-1]=X[0][0];
            work[j+n]=X[1][0];
                
            vmax=std::max((T)fabs(work[j+n-1]),std::max((T)fabs(work[j+n]),vmax) );
            vcrit=bignum/vmax;
          }
        }
label170:
            
        /* Copy the vector x or Q*x to VL and normalize. */
        if(!over)
        {
          copy(n-ki+1,&work[ki+n-1],1,&vl[ki-1+ldvl*(is-1)],1);
          ii=idamax(n-ki+1,&vl[ki-1+ldvl*(is-1)],1)+ki-1;
          remax=one/fabs(vl[ii-1+ldvl*(is-1)]);
          dscal<T>(n-ki+1,remax,&vl[ki-1+ldvl*(is-1)],1);
          for(k=1;k<=(ki-1);k++)
            vl[k-1+ldvl*(is-1)]=zero;
        }
        else
        {
          if(ki<n)
            dgemv<T>("N",n,n-ki,one,&vl[0+ldvl*ki],ldvl,&work[ki+n],1,work[ki+n-1],&vl[0+ldvl*(ki-1)],1);
          ii=idamax(n,&vl[0+ldvl*(ki-1)],1);
          remax=one/fabs(vl[ii-1+ldvl*(ki-1)]);
          dscal<T>(n,remax,&vl[0+ldvl*(ki-1)],1);
        }
      }
      else
      {
        /* Complex left eigenvector.
        * Initial solve:
  *             ((T(KI,KI)    T(KI,KI+1) )' - (WR - I* WI))*X = 0.
  *             ((T(KI+1,KI) T(KI+1,KI+1))                ) */
            
        if(fabs(Tri[ki-1+ldt*ki])>fabs(Tri[ki+ldt*(ki-1)]))
        {
          work[ki+n-1]=wi/Tri[ki-1+ldt*ki];
          work[ki+n2]=one;
        }
        else
        {
          work[ki+n-1]=one;
          work[ki+n2]=-wi/Tri[ki+ldt*(ki-1)];
        }
        work[ki+n]=zero;
        work[ki+n2-1]=zero;
            
        /* Form right-hand side */
          
        for(k=(ki+2);k<=n;k++)
        {
          work[k+n-1]=-work[ki+n-1]*Tri[ki-1+ldt*(k-1)];
          work[k+n2-1]=-work[ki+n2]*Tri[ki+ldt*(k-1)];
        }
            
        /* Solve complex quasi-triangular system:
  *          ( T(KI+2,N:KI+2,N) - (WR-i*WI) )*X = WORK1+i*WORK2 */
            
        vmax=one;
        vcrit=bignum;
            
        jnxt=ki+2;
        for(j=(ki+2);j<=n;j++)
        {
          if(j<jnxt)
            goto label200;
          j1=j;
          j2=j;
          jnxt=j+1;
          if(j<n)
          {
            if(Tri[j+ldt*(j-1)]!=zero)
            {
              j2=j+1;
              jnxt=j+2;
            }
          }
              
          if(j1!=j2)
          {
            /* 1-by-1 diagonal block
  *
  *            Scale if necessary to avoid overflow when
  *            forming the right-hand side elements. */
                
            if(work[j-1]>vcrit)
            {
              rec=one/vmax;
              dscal<T>(n-ki+1,rec,&work[ki+n-1],1);
              dscal<T>(n-ki+1,rec,&work[ki+n2-1],1);
              vmax=one;
              vcrit=bignum;
            }
            work[j+n-1]=work[j+n-1]-ddot<T>(j-ki-2,&Tri[ki+1+ldt*(j-1)],1,&work[ki+1+n],1);
            work[j+n2-1]=work[j+n2-1]-ddot<T>(j-ki-2,&Tri[ki+1+ldt*(j-1)],1,&work[ki+1+n2],1);
                
            /* Solve (T(J,J)-(WR-i*WI))*(X11+i*X12)= WK+I*WK2 */
              
            dlaln2<T>(false,1,2,smin,one,&Tri[j-1+ldt*(j-1)],ldt,one,one,&work[j+n-1],n,wr,-wi,&X[0][0],2,scale,xnorm,ierr);
                
            /* Scale if necessary */
            if(scale!=one)
            {
              dscal<T>(n-ki+1,scale,&work[ki+n-1],1);
              dscal<T>(n-ki+1,scale,&work[ki+n2-1],1);
            }
            work[j+n-1]=X[0][0];
            work[j+n2-1]=X[0][1];
            vmax=std::max((T)fabs(work[j+n-1]),std::max((T)fabs(work[j+n2-1]),vmax));
            vcrit=bignum/vmax;
          }
          else
          {
            /* 2-by-2 diagonal block
  *
  *            Scale if necessary to avoid overflow when forming
  *            the right-hand side elements. */
            beta=std::max(work[j-1],work[j]);
            if(beta>vcrit)
            {
              rec=one/vmax;
              dscal<T>(n-ki+1,rec,&work[ki+n-1],1);
              dscal<T>(n-ki+1,rec,&work[ki+n2-1],1);
              vmax=one;
              vcrit=bignum;
            }
                
            work[j+n-1]=work[j+n-1]-ddot<T>(j-ki-2,&Tri[ki+1+ldt*(j-1)],1,&work[ki+1+n],1);
            work[j+n2-1]=work[j+n2-1]-ddot<T>(j-ki-2,&Tri[ki+1+ldt*(j-1)],1,&work[ki+1+n2],1);
            work[j+n]=work[j+n]-ddot<T>(j-ki-2,&Tri[ki+1+ldt*j],1,&work[ki+1+n],1);
            work[j+n2]=work[j+n2]-ddot<T>(j-ki-2,&Tri[ki+1+ldt*j],1,&work[ki+1+n2],1);
                
            /* Solve 2-by-2 complex linear equation
  *                  ([T(j,j)   T(j,j+1)  ]'-(wr-i*wi)*I)*X = SCALE*B
  *                  ([T(j+1,j) T(j+1,j+1)]             ) */
                
            dlaln2<T>(true,2,2,smin,one,&Tri[j-1+ldt*(j-1)],ldt,one,one,&work[j+n-1],n,wr,-wi,&X[0][0],2,scale,xnorm,ierr);
                
            /* Scale if necessary */
              
            if(scale!=one)
            {
              dscal<T>(n-ki+1,scale,&work[ki+n-1],1);
              dscal<T>(n-ki+1,scale,&work[ki+n2-1],1);
            }
            work[j+n-1]=X[0][0];
            work[j+n2-1]=X[0][1];
            work[j+n]=X[1][0];
            work[j+n2]=X[1][1];
            vmax=std::max( std::max((T) std::max((T) fabs(X[0][0]),(T) fabs(X[0][1])),(T) std::max(fabs(X[1][0]),fabs(X[1][1]))),vmax);
            vcrit=bignum/vmax;
          }
        }
label200:
            
        /* Copy the vector x or Q*x to VL and normalize. */
        if(!over)
        {
          copy<T>(n-ki+1,&work[ki+n-1],1,&vl[ki-1+ldvl*(is-1)],1);
          copy<T>(n-ki+1,&work[ki+n2-1],1,&vl[ki-1+ldvl*is],1);
          emax=zero;
          for(k=ki;k<=n;k++)
            emax=std::max((T) emax,(T) fabs((T) vl[k-1+ldvl*(is-1)])+(T) fabs((T) vl[ki-1+ldvl*is]));
          remax=one/emax;
          dscal<T>(n-ki+1,remax,&vl[ki-1+ldvl*(is-1)],1);
          dscal<T>(n-ki+1,remax,&vl[ki-1+ldvl*is],1);
          for(k=1;k<=(ki-1);k++)
          {
            vl[k-1+ldvl*(is-1)]=zero;
            vl[k-1+ldvl*is]=zero;
          }
        }
        else
        {
          if(ki<(n-1))
          {
            dgemv<T>("N",n,n-ki-1,one,&vl[0+ldvl*(ki+1)],ldvl,&work[ki+1+n],1,work[ki+n-1],&vl[0+ldvl*(ki-1)],1);
            dgemv<T>("N",n,n-ki-1,one,&vl[0+ldvl*(ki+1)],ldvl,&work[ki+1+n2],1,work[ki+n2],&vl[0+ldvl*ki],1);
          }
          else
          {
            dscal<T>(n,work[ki+n-1],&vl[0+ldvl*(ki-1)],1);
            dscal<T>(n,work[ki+n2],&vl[0+ldvl*ki],1);
          }
          emax = zero;
          for(k=1;k<=n;k++)
            emax=std::max((T) emax,(T) fabs((T) vl[k-1+ldvl*(ki-1)])+ (T) fabs((T) vl[k-1+ldvl*ki]));
          remax=one/emax;
          dscal<T>(n,remax,&vl[0+ldvl*(ki-1)],1);
          dscal<T>(n,remax,&vl[0+ldvl*ki],1);
        }
      }
          
      is+=1;
      if(ip!=0)
        is+=1;
label250:
      if(ip==-1)
        ip=0;
      if(ip==1)
        ip=-1;
    }
  }
  return;
}
          
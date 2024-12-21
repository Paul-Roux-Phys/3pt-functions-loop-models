/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE neupp.h.
   Interface to ARPACK subroutines dneupd and sneupd.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef NEUPP_H
#define NEUPP_H

#include <cstddef>
#include <string>
#include "arch.h"
//#include "arpackf.h"
#include "cdtrmm.h"
#include "cdtrsen.h"
#include "cdorm2r.h"
#include "cdgeqr2.h"
#include "cdlahqr.h"
#include "fortranfuncs.h"

/*template <typename T>
inline void neupp(bool rvec, char HowMny, T dr[],
                  T di[], T Z[], ARint ldz, T sigmar,
                  T sigmai, T workv[], char bmat, ARint n,
                  const std::string& which, ARint nev, T tol, T resid[],
                  ARint ncv, T V[], ARint ldv, ARint iparam[],
                  ARint ipntr[], T workd[], T workl[],
                  ARint lworkl, ARint& info)
{
  return;
}*/

template<typename T>
inline void neupp(bool rvec, char howmny, double dr[],
                  double di[], double iz[], ARint ldz, double sigmar,
                  double sigmai, double workv[], char bmat, ARint n,
                  const std::string& which, ARint nev, double tol, double resid[],
                  ARint ncv, double v[], ARint ldv, ARint iparam[],
                  ARint ipntr[], double workd[], double workl[],
                  ARint lworkl, ARint& info)

/*
  c++ version of ARPACK routine dneupd.
  This subroutine returns the converged approximations to eigenvalues
  of A*z = lambda*B*z and (optionally):

  (1) the corresponding approximate eigenvectors,
  (2) an orthonormal basis for the associated approximate
      invariant subspace,

  There is negligible additional cost to obtain eigenvectors. An
  orthonormal basis is always computed.  There is an additional storage cost
  of n*nev if both are requested (in this case a separate array Z must be
  supplied).
  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
  are derived from approximate eigenvalues and eigenvectors of
  of the linear operator OP prescribed by the MODE selection in the
  call to naupp. naupp must be called before this routine is called.
  These approximate eigenvalues and vectors are commonly called Ritz
  values and Ritz vectors respectively.  They are referred to as such
  in the comments that follow.  The computed orthonormal basis for the
  invariant subspace corresponding to these Ritz values is referred to
  as a Schur basis.
  See documentation in the header of the subroutine naupp for
  definition of OP as well as other terms and the relation of computed
  Ritz values and Ritz vectors of OP with respect to the given problem
  A*z = lambda*B*z. For a brief description, see definitions of
  iparam[7], MODE and which in the documentation of naupp.

  Parameters:

    rvec    (Input) Specifies whether Ritz vectors corresponding to the
            Ritz value approximations to the eigenproblem A*z = lambda*B*z
            are computed.
            rvec = false: Compute Ritz values only.
            rvec = true : Compute the Ritz vectors or Schur vectors.
                          See Remarks below.
    HowMny  (Input) Specifies the form of the basis for the invariant
            subspace corresponding to the converged Ritz values that
            is to be computed.
            = 'A': Compute nev Ritz vectors;
            = 'P': Compute nev Schur vectors;
    dr      (Output) Array of dimension nev+1.
            If iparam[7] = 1,2 or 3 and sigmai=0.0  then on exit: dr
            contains the real part of the Ritz  approximations to the
            eigenvalues of A*z = lambda*B*z.
            If iparam[7] = 3, 4 and sigmai is not equal to zero, then on
            exit: dr contains the real part of the Ritz values of OP
            computed by naupp. A further computation must be performed by
            the user to transform the Ritz values computed for OP by naupp
            to those of the original system A*z = lambda*B*z. See remark 3.
    di      (Output) Array of dimension nev+1.
            On exit, di contains the imaginary part of the Ritz value
            approximations to the eigenvalues of A*z = lambda*B*z
            associated with dr.
            NOTE: When Ritz values are complex, they will come in complex
                  conjugate pairs.  If eigenvectors are requested, the
                  corresponding Ritz vectors will also come in conjugate
                  pairs and the real and imaginary parts of these are
                  represented in two consecutive columns of the array Z
                  (see below).
    Z       (Output) Array of dimension nev*n if rvec = TRUE and HowMny =
            'A'.  if rvec = TRUE. and HowMny = 'A', then the contains
            approximate eigenvectors (Ritz vectors) corresponding to the
            NCONV=iparam[5] Ritz values for eigensystem A*z = lambda*B*z.
            The complex Ritz vector associated with the Ritz value
            with positive imaginary part is stored in two consecutive
            columns.  The first column holds the real part of the Ritz
            vector and the second column holds the imaginary part.  The
            Ritz vector associated with the Ritz value with negative
            imaginary part is simply the complex conjugate of the Ritz
            vector associated with the positive imaginary part.
            If rvec = .FALSE. or HowMny = 'P', then Z is not referenced.
            NOTE: If if rvec = .TRUE. and a Schur basis is not required,
                  the array Z may be set equal to first nev+1 columns of
                  the Arnoldi basis array V computed by naupp.  In this
                  case the Arnoldi basis will be destroyed and overwritten
                  with the eigenvector basis.
    ldz     (Input) Dimension of the vectors contained in Z. This
            parameter MUST be set to n.
    sigmar  (Input) If iparam[7] = 3 or 4, represents the real part of
            the shift. Not referenced if iparam[7] = 1 or 2.
    sigmai  (Input) If iparam[7] = 3 or 4, represents the imaginary part
            of the shift. Not referenced if iparam[7] = 1 or 2. See
            remark 3 below.
    workv   (Workspace) Array of dimension 3*ncv.
    V       (Input/Output) Array of dimension n*ncv+1.
            Upon Input: V contains the ncv vectors of the Arnoldi basis
                        for OP as constructed by naupp.
            Upon Output: If rvec = TRUE the first NCONV=iparam[5] columns
                        contain approximate Schur vectors that span the
                        desired invariant subspace.  See Remark 2 below.
            NOTE: If the array Z has been set equal to first nev+1 columns
                  of the array V and rvec = TRUE. and HowMny = 'A', then
                  the Arnoldi basis held by V has been overwritten by the
                  desired Ritz vectors.  If a separate array Z has been
                  passed then the first NCONV=iparam[5] columns of V will
                  contain approximate Schur vectors that span the desired
                  invariant subspace.
    workl   (Input / Output) Array of length lworkl+1.
            workl[1:ncv*ncv+3*ncv] contains information obtained in
            naupp. They are not changed by neupp.
            workl[ncv*ncv+3*ncv+1:3*ncv*ncv+6*ncv] holds the real and
            imaginary part of the untransformed Ritz values, the upper
            quasi-triangular matrix for H, and the associated matrix
            representation of the invariant subspace for H.
    ipntr   (Input / Output) Array of length 14. Pointer to mark the
            starting locations in the workl array for matrices/vectors
            used by naupp and neupp.
            ipntr[9]:  pointer to the real part of the ncv RITZ values
                       of the original system.
            ipntr[10]: pointer to the imaginary part of the ncv RITZ
                       values of the original system.
            ipntr[11]: pointer to the ncv corresponding error bounds.
            ipntr[12]: pointer to the ncv by ncv upper quasi-triangular
                       Schur matrix for H.
            ipntr[13]: pointer to the ncv by ncv matrix of eigenvectors
                       of the upper Hessenberg matrix H. Only referenced
                       by neupp if rvec = TRUE. See Remark 2 below.
    info    (Output) Error flag.
            =  0 : Normal exit.
            =  1 : The Schur form computed by LAPACK routine dlahqr
                   could not be reordered by LAPACK routine dtrsen.
                   Re-enter subroutine neupp with iparam[5] = ncv and
                   increase the size of the arrays DR and DI to have
                   dimension at least dimension ncv and allocate at least
                   ncv columns for Z. NOTE: Not necessary if Z and V share
                   the same space. Please notify the authors if this error
                   occurs.
            = -1 : n must be positive.
            = -2 : nev must be positive.
            = -3 : ncv must satisfy nev+2 <= ncv <= n.
            = -5 : which must be one of 'LM','SM','LR','SR','LI','SI'.
            = -6 : bmat must be one of 'I' or 'G'.
            = -7 : Length of private work workl array is not sufficient.
            = -8 : Error return from calculation of a real Schur form.
                   Informational error from LAPACK routine dlahqr.
            = -9 : Error return from calculation of eigenvectors.
                   Informational error from LAPACK routine dtrevc.
            = -10: iparam[7] must be 1,2,3,4.
            = -11: iparam[7] = 1 and bmat = 'G' are incompatible.
            = -12: HowMny = 'S' not yet implemented
            = -13: HowMny must be one of 'A' or 'P' if rvec = TRUE.
            = -14: naupp did not find any eigenvalues to sufficient
                   accuracy.

  NOTE:     The following arguments

            bmat, n, which, nev, tol, resid, ncv, V, ldv, iparam,
            ipntr, workd, workl, lworkl, info

            must be passed directly to neupp following the last call
            to naupp.  These arguments MUST NOT BE MODIFIED between
            the the last call to naupp and the call to neupp.

  Remarks
    1. Currently only HowMny = 'A' and 'P' are implemented.
    2. Schur vectors are an orthogonal representation for the basis of
       Ritz vectors. Thus, their numerical properties are often superior.
       Let X' denote the transpose of X. If rvec = .TRUE. then the
       relationship A * V[:,1:iparam[5]] = V[:,1:iparam[5]] * T, and
       V[:,1:iparam[5]]' * V[:,1:iparam[5]] = I are approximately satisfied.
       Here T is the leading submatrix of order iparam[5] of the real
       upper quasi-triangular matrix stored workl[ipntr[12]]. That is,
       T is block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;
       each 2-by-2 diagonal block has its diagonal elements equal and its
       off-diagonal elements of opposite sign.  Corresponding to each
       2-by-2 diagonal block is a complex conjugate pair of Ritz values.
       The real Ritz values are stored on the diagonal of T.
    3. If iparam[7] = 3 or 4 and sigmai is not equal zero, then the user
       must form the iparam[5] Rayleigh quotients in order to transform the
       Ritz values computed by naupp for OP to those of A*z = lambda*B*z.
       Set rvec = TRUE. and HowMny = 'A', and compute
       Z[:,I]' * A * Z[:,I] if di[I] = 0.
       If di[I] is not equal to zero and di[I+1] = - D[I],
       then the desired real and imaginary parts of the Ritz value are
       Z[:,I]' * A * Z[:,I] +  Z[:,I+1]' * A * Z[:,I+1],
       Z[:,I]' * A * Z[:,I+1] -  Z[:,I+1]' * A * Z[:,I], respectively.
       Another possibility is to set rvec = .true. and HowMny = 'P' and
       compute V[:,1:iparam[5]]' * A * V[:,1:iparam[5]] and then an upper
       quasi-triangular matrix of order iparam[5] is computed. See remark
       2 above.
*/

{

 // ARint      irvec;
  bool* select;
  double*    z;
  const T one=1.0;
  const T zero=0.0;
  
  /*
c     %---------------%
c     | Local Scalars |
c     %---------------%
  */
  std::string type;
  ARint bounds,ierr,ih,ihbds,iheigr,iheigi,iconj,nconv;
  ARint invsub,iuptri,iwev,j,k,ktrord;
  ARint ldh,ldq,mode,msglvl,outncv,ritzr,ritzi,wri,wrr;
  ARint irr,iri,ibd;
  ARint iwork[1];
  bool reord;
  T conds,rnorm,sep,temp,thres,temp1,eps23;
  T vl[1];

 // irvec   = (ARint) rvec;
  select = new bool[ncv];
  z = (iz == NULL) ? &v[1] : z;

/*  F77NAME(dneupd)(&irvec, &HowMny, iselect, dr, di, iZ, &ldz, &sigmar,
                  &sigmai, &workv[1], &bmat, &n, which.c_str(), &nev, &tol,
                  resid, &ncv, &V[1], &ldv, &iparam[1], &ipntr[1],
                  &workd[1], &workl[1], &lworkl, &info);*/

/*
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c 
c     %------------------------%
c     | Set default parameters |
c     %------------------------% */

  msglvl=debug.mneupd;
  mode=iparam[7];
  nconv=iparam[5];
  info=0;

/*
c     %---------------------------------%
c     | Get machine dependent constant. |
c     %---------------------------------% */
  eps23=dlamch<T>("E"); // Epsilon-Machine
  eps23=pow(eps23,2.0/3.0);
  
  /*
c     %--------------%
c     | Quick return |
c     %--------------% */
  ierr=0;
  
  if(nconv<=0)
    ierr=-14;
  else if(n<=0)
    ierr=-1;
  else if(nev<=0)
    ierr=-2;
  else if( (ncv<=(nev+1))||(ncv>n) )
    ierr=-3;
  else if( (which!="LM")&&(which!="SM")&&(which!="LR")&&(which!="SR")&&(which!="LI")&&(which!="SI") )
    ierr=-5;
  else if( (bmat!='I')&&(bmat!='G') )
    ierr=-6;
  else if(lworkl<(3*ncv*ncv+6*ncv))
    ierr=-7;
  else if( (howmny!='A')&&(howmny!='P')&&(howmny!='S')&&rvec)
    ierr=-12;
  
  if( (mode==1)||(mode==2) )
    type="REGULR";
  else if( (mode==3)&&(sigmai==zero) )
    type="SHIFTI";
  else if(mode==3)
    type="REALPT";
  else if(mode==4)
    type="IMAGPT";
  else
    ierr=-10;
  
  if( (mode==1)&&(bmat=='G') )
    ierr=-11;
  
  /*
c     %------------%
c     | Error Exit |
c     %------------% */
  
  if(ierr!=0)
  {
    info=ierr;
    goto label9000;
  }
  
  /*
c     %--------------------------------------------------------%
c     | Pointer into WORKL for address of H, RITZ, BOUNDS, Q   |
c     | etc... and the remaining workspace.                    |
c     | Also update pointer to be used on output.              |
c     | Memory is laid out as follows:                         |
c     | workl(1:ncv*ncv) := generated Hessenberg matrix        |
c     | workl(ncv*ncv+1:ncv*ncv+2*ncv) := real and imaginary   |
c     |                                   parts of ritz values |
c     | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := error bounds   |
c     %--------------------------------------------------------%
c
c     %-----------------------------------------------------------%
c     | The following is used and set by DNEUPD.                  |
c     | workl(ncv*ncv+3*ncv+1:ncv*ncv+4*ncv) := The untransformed |
c     |                             real part of the Ritz values. |
c     | workl(ncv*ncv+4*ncv+1:ncv*ncv+5*ncv) := The untransformed |
c     |                        imaginary part of the Ritz values. |
c     | workl(ncv*ncv+5*ncv+1:ncv*ncv+6*ncv) := The untransformed |
c     |                           error bounds of the Ritz values |
c     | workl(ncv*ncv+6*ncv+1:2*ncv*ncv+6*ncv) := Holds the upper |
c     |                             quasi-triangular matrix for H |
c     | workl(2*ncv*ncv+6*ncv+1: 3*ncv*ncv+6*ncv) := Holds the    |
c     |       associated matrix representation of the invariant   |
c     |       subspace for H.                                     |
c     | GRAND total of NCV * ( 3 * NCV + 6 ) locations.           |
c     %-----------------------------------------------------------% */
  
  ih=ipntr[5];
  ritzr=ipntr[6];
  ritzi=ipntr[7];
  bounds=ipntr[8];
  ldh=ncv;
  ldq=ncv;
  iheigr=bounds+ldh;
  iheigi=iheigr+ldh;
  ihbds=iheigi+ldh;
  iuptri=ihbds+ldh;
  invsub=iuptri+ldh*ncv;
  ipntr[9]=iheigr;
  ipntr[10]=iheigi;
  ipntr[11]=ihbds;
  ipntr[12]=iuptri;
  ipntr[13]=invsub;
  wrr=1;
  wri=ncv+1;
  iwev=wri+ncv;
  
  /*
c     %-----------------------------------------%
c     | irr points to the REAL part of the Ritz |
c     |     values computed by _neigh before    |
c     |     exiting _naup2.                     |
c     | iri points to the IMAGINARY part of the |
c     |     Ritz values computed by _neigh      |
c     |     before exiting _naup2.              |
c     | ibd points to the Ritz estimates        |
c     |     computed by _neigh before exiting   |
c     |     _naup2.                             |
c     %-----------------------------------------% */
  
  irr=ipntr[14]+ncv*ncv;
  iri=irr+ncv;
  ibd=iri+ncv;
  
  /*
c     %------------------------------------%
c     | RNORM is B-norm of the RESID(1:N). |
c     %------------------------------------% */
  
  rnorm=workl[ih+2];
  workl[ih+2]=zero;
  
  if(rvec)
  {
    /*
c        %-------------------------------------------%
c        | Get converged Ritz value on the boundary. |
c        | Note: converged Ritz values have been     |
c        | placed in the first NCONV locations in    |
c        | workl(ritzr) and workl(ritzi).  They have |
c        | been sorted (in _naup2) according to the  |
c        | WHICH selection criterion.                |
c        %-------------------------------------------% */
    if( (which=="LM")||(which=="SM") )
      thres=dlapy2<T>(workl[ritzr],workl[ritzi]);
    else if( (which=="LR")||(which=="SR") )
      thres=workl[ritzr];
    else if( (which=="LI")||(which=="SI") )
      thres=fabs(workl[ritzi]);
    
    if(msglvl>2)
      dvout<T>(debug.logfil,1,&thres,debug.ndigit,"_neupd: Threshold eigenvalue used for re-ordering");
    
    /*
c        %----------------------------------------------------------%
c        | Check to see if all converged Ritz values appear at the  |
c        | top of the upper quasi-triangular matrix computed by     |
c        | _neigh in _naup2.  This is done in the following way:    |
c        |                                                          |
c        | 1) For each Ritz value obtained from _neigh, compare it  |
c        |    with the threshold Ritz value computed above to       |
c        |    determine whether it is a wanted one.                 |
c        |                                                          | 
c        | 2) If it is wanted, then check the corresponding Ritz    |
c        |    estimate to see if it has converged.  If it has, set  |
c        |    correponding entry in the logical array SELECT to     |
c        |    .TRUE..                                               |
c        |                                                          |
c        | If SELECT(j) = .TRUE. and j > NCONV, then there is a     |
c        | converged Ritz value that does not appear at the top of  |
c        | the upper quasi-triangular matrix computed by _neigh in  |
c        | _naup2.  Reordering is needed.                           |
c        %----------------------------------------------------------% */
    
    reord=false;
    ktrord=0;
    for(j=0;j<=(ncv-1);j++)
    {
      select[j]=false;
      if(which=="LM")
      {
        if(dlapy2<T>(workl[irr+j],workl[iri+j])>=thres)
        {
          temp1=std::max(eps23,dlapy2(workl[irr+j],workl[iri+j]));
          if(workl[ibd+j]<=(tol*temp1))
            select[j]=true;
        }
      }
      else if(which=="SM")
      {
        if(dlapy2<T>(workl[irr+j],workl[iri+j]))
        {
          temp1=std::max(eps23,dlapy2(workl[irr+j],workl[iri+j]));
          if(workl[ibd+j]<=(tol*temp1))
            select[j]=true;
        }
      }
      else if(which=="LR")
      {
        if(workl[irr+j]>=thres)
        {
          temp1=std::max(eps23,dlapy2(workl[irr+j],workl[iri+j]));
          if(workl[ibd+j]<=(tol*temp1))
            select[j]=true;
        }
      }
      else if(which=="SR")
      {
        if(workl[irr+j]<=thres)
        {
          temp1=std::max(eps23,dlapy2<T>(workl[irr+j],workl[iri+j]));
          if(workl[ibd+j]<=(tol*temp1))
            select[j]=true;
        }
      }
      else if(which=="SI")
      {
        if(fabs(workl[iri+j])<=thres)
        {
          temp1=std::max(eps23,dlapy2<T>(workl[irr+j],workl[iri+j]));
          if(workl[ibd+j]<=(tol*temp1))
            select[j]=true;
        }
      }
      if((j+1)>nconv)
        reord=(select[j]||reord);
      if(select[j]) 
        ktrord+=1;
    } // for(j=0...
    
    if(msglvl>2)
    {
      ivout<T>(debug.logfil,1,&ktrord,debug.ndigit,"_neupd: Number of specified eigenvalues");
      ivout<T>(debug.logfil,1,&nconv,debug.ndigit,"_neupd: Number of specified eigenvalues");
    }
    
    /*
c        %-----------------------------------------------------------%
c        | Call LAPACK routine dlahqr to compute the real Schur form |
c        | of the upper Hessenberg matrix returned by DNAUPD.        |
c        | Make a copy of the upper Hessenberg matrix.               |
c        | Initialize the Schur vector matrix Q to the identity.     |
c        %-----------------------------------------------------------% */
    
    copy<T>(ldh*ncv,&workl[ih],1,&workl[iuptri],1);
    dlaset<T>("A",ncv,ncv,zero,one,&workl[invsub],ldq);
    dlahqr<T>(true,true,ncv,1,ncv,&workl[iuptri],ldh,&workl[iheigr],&workl[iheigi],1,ncv,&workl[invsub],ldq,ierr);
    copy<T>(ncv,&workl[invsub+ncv-1],ldq,&workl[ihbds],1);
    
    if(ierr!=0)
    {
      info=-8;
      goto label9000;
    }
    
    if(msglvl>1)
    {
      dvout<T>(debug.logfil,ncv,&workl[iheigr],debug.ndigit,"_neupd: Real part of the eigenvalues of H");
      dvout<T>(debug.logfil,ncv,&workl[iheigi],debug.ndigit,"_neupd: Real part of the eigenvalues of H");
      dvout<T>(debug.logfil,ncv,&workl[ihbds],debug.ndigit,"_neupd: Last row of the Schur vector matrix");
      if(msglvl>3)
        dmout<T>(debug.logfil, ncv, ncv, &workl[iuptri], ldh, debug.ndigit,"_neupd: The upper quasi-triangular matrix ");
    }
    
    if(reord)
    {
      /*
c           %-----------------------------------------------------%
c           | Reorder the computed upper quasi-triangular matrix. | 
c           %-----------------------------------------------------% */
      
      dtrsen<T>("N","V",select,ncv,&workl[iuptri],ldh,&workl[invsub],ldq,&workl[iheigr],&workl[iheigi],nconv,conds,sep,&workl[ihbds],ncv,iwork,1,ierr);
      
      if(ierr==1)
      {
        info=1;
        goto label9000;
      }
      
      if(msglvl>2)
      {
        dvout<T>(debug.logfil,ncv,&workl[iheigr],debug.ndigit,"_neupd: Real part of the eigenvalues of H--reordered");
        dvout<T>(debug.logfil,ncv,&workl[iheigi],debug.ndigit,"_neupd: Real part of the eigenvalues of H--reordered");
        if(msglvl>3)
        {
          dmout<T>(debug.logfil, ncv, ncv, &workl[iuptri], ldh, debug.ndigit,"_neupd: Quasi-triangular matrix after re-ordering");
        }
      }
    }
    
    /*
c        %---------------------------------------%
c        | Copy the last row of the Schur vector |
c        | into workl(ihbds).  This will be used |
c        | to compute the Ritz estimates of      |
c        | converged Ritz values.                |
c        %---------------------------------------% */
    
    copy<T>(ncv,&workl[invsub+ncv-1],ldq,&workl[ihbds],1);
    
    /*
c        %----------------------------------------------------%
c        | Place the computed eigenvalues of H into DR and DI |
c        | if a spectral transformation was not used.         |
c        %----------------------------------------------------% */
    
    if(type=="REGULR")
    {
      copy<T>(nconv,&workl[iheigr],1,dr,1);
      copy<T>(nconv,&workl[iheigi],1,di,1);
    }
    
    /*
c        %----------------------------------------------------------%
c        | Compute the QR factorization of the matrix representing  |
c        | the wanted invariant subspace located in the first NCONV |
c        | columns of workl(invsub,ldq).                            |
c        %----------------------------------------------------------% */
    
    dgeqr2<T>(ncv,nconv,&workl[invsub],ldq,workv,&workv[ncv+1],ierr);
    
    /*
c        %---------------------------------------------------------%
c        | * Postmultiply V by Q using dorm2r.                     |   
c        | * Copy the first NCONV columns of VQ into Z.            |
c        | * Postmultiply Z by R.                                  |
c        | The N by NCONV matrix Z is now a matrix representation  |
c        | of the approximate invariant subspace associated with   |
c        | the Ritz values in workl(iheigr) and workl(iheigi)      |
c        | The first NCONV columns of V are now approximate Schur  |
c        | vectors associated with the real upper quasi-triangular |
c        | matrix of order NCONV in workl(iuptri)                  |
c        %---------------------------------------------------------% */
    
    dorm2r<T>("R","N",n,ncv,nconv,&workl[invsub],ldq,workv,v,ldv,&workd[n+1],ierr);
    dlacpy<T>("A",n,nconv,v,ldv,z,ldz);
    
    for(j=1;j<=nconv;j++)
    {
      /*
c           %---------------------------------------------------%
c           | Perform both a column and row scaling if the      |
c           | diagonal element of workl(invsub,ldq) is negative |
c           | I'm lazy and don't take advantage of the upper    |
c           | quasi-triangular form of workl(iuptri,ldq)        |
c           | Note that since Q is orthogonal, R is a diagonal  |
c           | matrix consisting of plus or minus ones           |
c           %---------------------------------------------------% */
      if(workl[invsub+(j-1)*ldq+j-1]<zero)
      {
        dscal<T>(nconv,-one,&workl[iuptri+j-1],ldq);
        dscal<T>(nconv,-one,&workl[iuptri+(j-1)*ldq],1);
      }
    }
    
    if(howmny=='A')
    {
      /*
c           %--------------------------------------------%
c           | Compute the NCONV wanted eigenvectors of T | 
c           | located in workl(iuptri,ldq).              |
c           %--------------------------------------------% */
      
      for(j=1;j<=ncv;j++)
      {
        if(j<=nconv)
          select[j-1]=true;
        else
          select[j-1]=false;
      }
      
      dtrevc<T>("R","S",select,ncv,&workl[iuptri],ldq,vl,1,&workl[invsub],ldq,ncv,outncv,workv,ierr);
      
      if(ierr!=0)
      {
        info=-9;
        goto label9000;
      }
      
      /*
c           %------------------------------------------------%
c           | Scale the returning eigenvectors so that their |
c           | Euclidean norms are all one. LAPACK subroutine |
c           | dtrevc returns each eigenvector normalized so  |
c           | that the element of largest magnitude has      |
c           | magnitude 1;                                   |
c           %------------------------------------------------% */
      
      iconj=0;
      for(j=1;j<=nconv;j++)
      {
        if(workl[iheigi+j-1]==zero)
        {
          /*
c                 %----------------------%
c                 | real eigenvalue case |
c                 %----------------------% */
          
          temp=dnrm2<T>(ncv,&workl[invsub+(j-1)*ldq],1);
          dscal<T>(ncv,one/temp,&workl[invsub+(j-1)*ldq],1);
        }
        else
        {
          /*
c                 %-------------------------------------------%
c                 | Complex conjugate pair case. Note that    |
c                 | since the real and imaginary part of      |
c                 | the eigenvector are stored in consecutive |
c                 | columns, we further normalize by the      |
c                 | square root of two.                       |
c                 %-------------------------------------------% */
          if(iconj==0)
          {
            temp=dlapy2<T>(dnrm2<T>(ncv,&workl[invsub+(j-1)*ldq],1),dnrm2<T>(ncv,&workl[invsub+j*ldq],1) );
            dscal<T>(ncv,one/temp,&workl[invsub+(j-1)*ldq],1);
            dscal<T>(ncv,one/temp,&workl[invsub+j*ldq],1);
            iconj=1;
          }
          else
            iconj=0;
        }
      }
      
      dgemv<T>("T",ncv,nconv,one,&workl[invsub],ldq,&workl[ihbds],1,zero,workv,1);
      
      iconj=0;
      for(j=1;j<=nconv;j++)
      {
        if(workl[iheigi+j-1]!=zero)
        {
          /*
c                 %-------------------------------------------%
c                 | Complex conjugate pair case. Note that    |
c                 | since the real and imaginary part of      |
c                 | the eigenvector are stored in consecutive |
c                 %-------------------------------------------% */
          if(iconj==0)
          {
            workv[j]=dlapy2<T>(workv[j],workv[j+1]);
            workv[j+1]=workv[j];
            iconj=1;
          }
          else
            iconj=0;
        }
      }
      
      if(msglvl>2)
      {
        copy<T>(ncv,&workl[invsub+ncv-1],ldq,&workl[ihbds],1);
        dvout<T>(debug.logfil,ncv,&workl[ihbds],debug.ndigit,"_neupd: Last row of the eigenvector matrix for T");
        if(msglvl>3)
          dmout<T>(debug.logfil, ncv, ncv, &workl[invsub], ldq, debug.ndigit,"_neupd: The eigenvector matrix for T");
      }
      
      /*
c           %---------------------------------------%
c           | Copy Ritz estimates into workl(ihbds) |
c           %---------------------------------------% */
      
      copy<T>(nconv,workv,1,&workl[ihbds],1);
      
      /*
c           %---------------------------------------------------------%
c           | Compute the QR factorization of the eigenvector matrix  |
c           | associated with leading portion of T in the first NCONV |
c           | columns of workl(invsub,ldq).                           |
c           %---------------------------------------------------------% */
      
      dgeqr2<T>(ncv,nconv,&workl[invsub],ldq,workv,&workv[ncv+1],ierr);
      
      /*
c           %----------------------------------------------%
c           | * Postmultiply Z by Q.                       |   
c           | * Postmultiply Z by R.                       |
c           | The N by NCONV matrix Z is now contains the  | 
c           | Ritz vectors associated with the Ritz values |
c           | in workl(iheigr) and workl(iheigi).          |
c           %----------------------------------------------% */
      
      dorm2r<T>("R","N",n,ncv,nconv,&workl[invsub],ldq,workv,z,ldz,&workd[n+1],ierr);
      dtrmm<T>("R","U","N","N",n,nconv,one,&workl[invsub],ldq,z,ldz);
    } // if(howmny=="A")
  }
  else
  {
    /*
c        %------------------------------------------------------%
c        | An approximate invariant subspace is not needed.     |
c        | Place the Ritz values computed DNAUPD into DR and DI |
c        %------------------------------------------------------% */
    copy<T>(nconv,&workl[ritzr],1,dr,1);
    copy<T>(nconv,&workl[ritzi],1,di,1);
    copy<T>(nconv,&workl[ritzr],1,&workl[iheigr],1);
    copy<T>(nconv,&workl[ritzi],1,&workl[iheigi],1);
    copy<T>(nconv,&workl[bounds],1,&workl[ihbds],1);
  }
  
  /*
c     %------------------------------------------------%
c     | Transform the Ritz values and possibly vectors |
c     | and corresponding error bounds of OP to those  |
c     | of A*x = lambda*B*x.                           |
c     %------------------------------------------------% */
  
  if(type=="REGULR")
  {
    if(rvec)
      dscal<T>(ncv,rnorm,&workl[ihbds],1);
  }
  else
  {
  /*
c        %---------------------------------------%
c        |   A spectral transformation was used. |
c        | * Determine the Ritz estimates of the |
c        |   Ritz values in the original system. |
c        %---------------------------------------% */
  
    if(type=="SHIFTI")
    {
      if(rvec)
        dscal<T>(ncv,rnorm,&workl[ihbds],1);
    
      for(k=1;k<=ncv;k++)
      {
        temp=dlapy2(workl[iheigr+k-1],workl[iheigi+k-1]);
        workl[ihbds+k-1]=fabs(workl[ihbds+k-1])/temp/temp;
      }
    }
    else if(type=="REALPT")
    {
      for(k=1;k<=ncv;k++)
        continue; // ??
    }
    else if(type=="IMAGPT")
    {
      for(k=1;k<=ncv;k++)
        continue;
    }
  
  /*
c        %-----------------------------------------------------------%
c        | *  Transform the Ritz values back to the original system. |
c        |    For TYPE = 'SHIFTI' the transformation is              |
c        |             lambda = 1/theta + sigma                      |
c        |    For TYPE = 'REALPT' or 'IMAGPT' the user must from     |
c        |    Rayleigh quotients or a projection. See remark 3 above.| 
c        | NOTES:                                                    |
c        | *The Ritz vectors are not affected by the transformation. |
c        %-----------------------------------------------------------% */
  
    if(type=="SHIFTI")
    {
      for(k=1;k<=ncv;k++)
      {
        temp=dlapy2<T>(workl[iheigr+k-1],workl[iheigi+k-1]);
        workl[iheigr+k-1]=workl[iheigr+k-1]/temp/temp+sigmar;
        workl[iheigi+k-1]=-workl[iheigi+k-1]/temp/temp+sigmai;
      }
    
      copy<T>(nconv,&workl[iheigr],1,dr,1);
      copy<T>(nconv,&workl[iheigi],1,di,1);
    }
    else if( (type=="REALPT")||(type=="IMAGPT") )
    {
      copy<T>(nconv,&workl[iheigr],1,dr,1);
      copy<T>(nconv,&workl[iheigi],1,di,1);
    }
  } // if(type=="REGULR") ... else
  
  if( (type=="SHIFTI")&&(msglvl>1) )
  {
    dvout<T>(debug.logfil,nconv,dr,debug.ndigit,"_neupd: Untransformed real part of the Ritz valuess.");
    dvout<T>(debug.logfil,nconv,di,debug.ndigit,"_neupd: Untransformed imag part of the Ritz valuess.");
    dvout<T>(debug.logfil,nconv,&workl[ihbds],debug.ndigit,"_neupd: Ritz estimates of untransformed Ritz values.");
  }
  else if( (type=="REGULR")&&(msglvl>1) )
  {
    dvout<T>(debug.logfil,nconv,dr,debug.ndigit,"_neupd: Real parts of converged Ritz values.");
    dvout<T>(debug.logfil,nconv,di,debug.ndigit,"_neupd: Imag parts of converged Ritz values.");
    dvout<T>(debug.logfil,nconv,&workl[ihbds],debug.ndigit,"_neupd: Associated Ritz estimates.");
  }
  
  /*
c     %-------------------------------------------------%
c     | Eigenvector Purification step. Formally perform |
c     | one of inverse subspace iteration. Only used    |
c     | for MODE = 2.                                   |
c     %-------------------------------------------------% */
  
  if(rvec&&(howmny=='A')&&(type=="SHIFTI"))
  {
    /*
c        %------------------------------------------------%
c        | Purify the computed Ritz vectors by adding a   |
c        | little bit of the residual vector:             |
c        |                      T                         |
c        |          resid(:)*( e    s ) / theta           |
c        |                      NCV                       |
c        | where H s = s theta. Remember that when theta  |
c        | has nonzero imaginary part, the corresponding  |
c        | Ritz vector is stored across two columns of Z. |
c        %------------------------------------------------% */
    
    iconj=0;
    for(j=1;j<=nconv;j++)
    {
      if(workl[iheigi+j-1]==zero)
        workv[j]=workl[invsub+(j-1)*ldq+ncv-1]/workl[iheigr+j-1];
      else if(iconj==0)
      {
        temp=dlapy2<T>(workl[iheigr+j-1],workl[iheigi+j-1]);
        workv[j]=(workl[invsub+(j-1)*ldq+ncv-1]*workl[iheigr+j-1]+workl[invsub+j*ldq+ncv-1]*workl[iheigi+j-1])/temp/temp;
        workv[j+1]=(workl[invsub+j*ldq+ncv-1]*workl[iheigr+j-1]+workl[invsub+(j-1)*ldq+ncv-1]*workl[iheigi+j-1])/temp/temp;
        iconj=1;
      }
      else
        iconj=0;
    }
    /*
c        %---------------------------------------%
c        | Perform a rank one update to Z and    |
c        | purify all the Ritz vectors together. |
c        %---------------------------------------% */
    
    dger<T>(n,nconv,one,resid,1,workv,1,z,ldz);
    
  }
  
label9000:
  
  delete[] select;

  return;
  
  /*
c     %---------------%
c     | End of DNEUPD |
c     %---------------%
  */
        
} // neupp (double).

#endif

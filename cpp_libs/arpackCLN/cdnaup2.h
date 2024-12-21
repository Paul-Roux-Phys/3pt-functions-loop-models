/* Description: 
 * Intermediate level interface called by dnaupd. */

/* Translated from Fortran 77, gotos and all */

/*
c\BeginDoc
c
c\Name: dnaup2
c
c\Description: 
c  Intermediate level interface called by dnaupd.
c
c\Usage:
c  call dnaup2
c     ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,
c       ISHIFT, MXITER, V, LDV, H, LDH, RITZR, RITZI, BOUNDS, 
c       Q, LDQ, WORKL, IPNTR, WORKD, INFO )
c
c\Arguments
c
c  IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in dnaupd.
c  MODE, ISHIFT, MXITER: see the definition of IPARAM in dnaupd.
c
c  NP      Integer.  (INPUT/OUTPUT)
c          Contains the number of implicit shifts to apply during 
c          each Arnoldi iteration.  
c          If ISHIFT=1, NP is adjusted dynamically at each iteration 
c          to accelerate convergence and prevent stagnation.
c          This is also roughly equal to the number of matrix-vector 
c          products (involving the operator OP) per Arnoldi iteration.
c          The logic for adjusting is contained within the current
c          subroutine.
c          If ISHIFT=0, NP is the number of shifts the user needs
c          to provide via reverse comunication. 0 < NP < NCV-NEV.
c          NP may be less than NCV-NEV for two reasons. The first, is
c          to keep complex conjugate pairs of "wanted" Ritz values 
c          together. The second, is that a leading block of the current
c          upper Hessenberg matrix has split off and contains "unwanted"
c          Ritz values.
c          Upon termination of the IRA iteration, NP contains the number 
c          of "converged" wanted Ritz values.
c
c  IUPD    Integer.  (INPUT)
c          IUPD .EQ. 0: use explicit restart instead implicit update.
c          IUPD .NE. 0: use implicit update.
c
c  V       Double precision N by (NEV+NP) array.  (INPUT/OUTPUT)
c          The Arnoldi basis vectors are returned in the first NEV 
c          columns of V.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling 
c          program.
c
c  H       Double precision (NEV+NP) by (NEV+NP) array.  (OUTPUT)
c          H is used to store the generated upper Hessenberg matrix
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling 
c          program.
c
c  RITZR,  Double precision arrays of length NEV+NP.  (OUTPUT)
c  RITZI   RITZR(1:NEV) (resp. RITZI(1:NEV)) contains the real (resp.
c          imaginary) part of the computed Ritz values of OP.
c
c  BOUNDS  Double precision array of length NEV+NP.  (OUTPUT)
c          BOUNDS(1:NEV) contain the error bounds corresponding to 
c          the computed Ritz values.
c          
c  Q       Double precision (NEV+NP) by (NEV+NP) array.  (WORKSPACE)
c          Private (replicated) work array used to accumulate the
c          rotation in the shift application step.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  WORKL   Double precision work array of length at least 
c          (NEV+NP)**2 + 3*(NEV+NP).  (INPUT/WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  It is used in shifts calculation, shifts
c          application and convergence checking.
c
c          On exit, the last 3*(NEV+NP) locations of WORKL contain
c          the Ritz values (real,imaginary) and associated Ritz
c          estimates of the current Hessenberg matrix.  They are
c          listed in the same order as returned from dneigh.
c
c          If ISHIFT .EQ. O and IDO .EQ. 3, the first 2*NP locations
c          of WORKL are used in reverse communication to hold the user 
c          supplied shifts.
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c          Pointer to mark the starting locations in the WORKD for 
c          vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X.
c          IPNTR(2): pointer to the current result vector Y.
c          IPNTR(3): pointer to the vector B * X when used in the 
c                    shift-and-invert mode.  X is the current operand.
c          -------------------------------------------------------------
c          
c  WORKD   Double precision work array of length 3*N.  (WORKSPACE)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The user should not use WORKD
c          as temporary workspace during the iteration !!!!!!!!!!
c          See Data Distribution Note in DNAUPD.
c
c  INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =     0: Normal return.
c          =     1: Maximum number of iterations taken.
c                   All possible eigenvalues of OP has been found.  
c                   NP returns the number of converged Ritz values.
c          =     2: No shifts could be applied.
c          =    -8: Error return from LAPACK eigenvalue calculation;
c                   This should never happen.
c          =    -9: Starting vector is zero.
c          = -9999: Could not build an Arnoldi factorization.
c                   Size that was built in returned in NP.
c
  */

#ifndef CDNAUP2_H
#define CDNAUP2_H
#include<string>
#include<iostream>
#include<math.h>
#include "fortranfuncs.h"
#include "cdgetv0.h"
#include "cdnaitr.h"
#include "cdebug.h"
#include "cstat.h"
#include "cdneigh.h"
#include "cdsortc.h"
#include "cdngets.h"
#include "cdnconv.h"
#include "cdnapps.h"
#include<algorithm>

//struct timingstruct timing;
//struct debugstruct debug;

/* CLN version */

template <typename T>
void dnaup2(ARint& ido, char bmat, ARint& n,const std::string& which, int& nev, int& np,
            T& tol, T* resid, int& mode, int iupd, int& ishift, int& mxiter,
            T* v, int ldv, T* h, int ldh, T* ritzr, T* ritzi, T* bounds,
            T* q, int ldq, T* workl, int* ipntr, T* workd, int& info, int digits)
{
  /*--------------------------------------------------------*/
  /* Include variables for debugging and timing information */
  /*--------------------------------------------------------*/
  
 // struct timingstruct timing;
 // struct debugstruct debug;
  
  /*------------*/
  /* Parameters */
  /*------------*/
  
  cln::float_format_t precision=cln::float_format(digits);
  const T zero = cln::cl_float(0,precision);
  const T one = cln::cl_float(1,precision);
  const T two = cln::cl_float(2,precision);
  const T three = cln::cl_float(3,precision);
  
  /*---------------*/
  /* Local Scalars */
  /*---------------*/
  
  static std::string wprime;
  static bool cnorm, getv0, initv, update, ushift;
  static int ierr,iter,j,kplusp,msglvl,nconv,nevbef,nev0,np0,nptemp,numcnv;
  static T rnorm,temp,eps23;
  
  /* Local array arguments */
  
  static int kp[4];
  //std::cout << "starting np = " << np << "\n";
  if(ido==0)
  {
    second(timing.t0);
    
    msglvl=debug.mnaup2;
    
    /*-------------------------------------*/
    /* Get the machine dependent constant. */
    /*-------------------------------------*/
    
    eps23=dlamch<T>("E",digits); // Epsilon-machine
    eps23=cln::cl_float(cln::exp(two/three*cln::ln(eps23)),precision); //  eps23=pow(eps23,2.0/3.0);
    //eps23=pow(eps23,2.0/3.0);
 /*   dvout<T>(debug.logfil, 1, &eps23, debug.ndigit,"_dnaup2: eps23");
    T test = dlamch<T>("E");
    dvout<T>(debug.logfil, 1, &test, debug.ndigit,"_dnaup2: dlamch E");
    test = dlamch<T>("S");
    dvout<T>(debug.logfil, 1, &test, debug.ndigit,"_dnaup2: dlamch S");
    test = dlamch<T>("P");
    dvout<T>(debug.logfil, 1, &test, debug.ndigit,"_dnaup2: dlamch P");
    test = dlamch<T>("B");
    dvout<T>(debug.logfil, 1, &test, debug.ndigit,"_dnaup2: dlamch B");
    double safmin = dlamch<T>("S");
    double eps=eps23=dlamch<T>("E");
    double safmn2=pow(dlamch<T>("B"),int(log(safmin/eps)/log(dlamch<T>("B"))/2.0) );
    double safmx2=one/safmn2;
    dvout<T>(debug.logfil, 1, &safmn2, debug.ndigit,"_dnaup2: safmn2");
    dvout<T>(debug.logfil, 1, &safmx2, debug.ndigit,"_dnaup2: safmx2");*/
    
    //std::cout << "eps23 = " << eps23 << "\n";
    
    nev0=nev;
    np0=np;
    
    /*-------------------------------------*/
    /* kplusp is the bound on the largest  */
    /*        Lanczos factorization built. */
    /* nconv is the current number of      */
    /*         "converged" eigenvalues.    */
    /* iter is the counter on the current  */
    /*      iteration step.                */
    /*-------------------------------------*/
    
    kplusp = nev + np;
    nconv = 0;
    iter = 0;
    
    /*---------------------------------------*/
    /* Set flags for computing the first NEV */
    /* steps of the Arnoldi factorization    */
    /*---------------------------------------*/
    
    getv0 = true;
    update = false;
    ushift = false;
    cnorm = false;
    
  //  std::cout << "info = " << info << "\n";
    if (info!=0)
    {
      /* User provides the initial residual vector. */
      
      initv = true;
      info = 0;
    }
    else
      initv = false;
  }
    
  if (getv0)
  {
//    dvout<T>(debug.logfil, n, workd, debug.ndigit,"_dnaup2: workd before call to getv0");
    dgetv0<T>(ido,bmat,1,initv,n,1,v,ldv,resid,rnorm,ipntr,workd,info,digits);
//    dvout<T>(debug.logfil, n, workd, debug.ndigit,"_dnaup2: workd after call to getv0");
    
 //   ivout<T>(debug.logfil,1,&ido,debug.ndigit,"_dnaup2: ido");
    if(ido!=99) return;
    
    if(rnorm == zero)
    {
      /*-----------------------------------------*/
      /* The initial vector is zero. Error exit. */
      /*-----------------------------------------*/
      
      info = -9;
      goto label1100;
    }
    getv0=false;
    ido=0;
  }
  
  /*-----------------------------------*/
  /* Back from reverse communication : */
  /* continue with update step         */
  /*-----------------------------------*/
  
  if(update) goto label20;
  
  /*-------------------------------------------*/
  /* Back from computing user specified shifts */
  /*-------------------------------------------*/
  
  if(ushift) goto label50;
  
  /*-------------------------------------*/
  /* Back from computing residual norm   */
  /* at the end of the current iteration */
  /*-------------------------------------*/
  
  if(cnorm) goto label100;

  /*----------------------------------------------------------*/
  /* Compute the first NEV steps of the Arnoldi factorization */
  /*----------------------------------------------------------*/

//  dmout<T>(debug.logfil,kplusp,kplusp,h,ldh,debug.ndigit,"_dnaup2: H before call to dnaitr");
  dnaitr<T>(ido, bmat, n, 0, nev, mode, resid, rnorm, v, ldv, h, ldh, ipntr, workd, info,digits);
//  dmout<T>(debug.logfil,kplusp,kplusp,h,ldh,debug.ndigit,"_dnaup2: H after call to dnaitr");
  /*---------------------------------------------------*/
  /* ido!=99 implies use of reverse communication      */
  /* to compute operations involving OP and possibly B */
  /*---------------------------------------------------*/
  
  if(ido!=99) return;
  
  if(info>0)
  {
    np=info;
    mxiter=iter;
    info=-9999;
    goto label1200;
  }
  
  /*--------------------------------------------------------------*/
  /*                                                              */
  /*           M A I N   ARNOLDI  I T E R A T I O N  L O O P      */
  /*           Each iteration implicitly restarts the Arnoldi     */
  /*           factorization in place.                            */
  /*                                                              */
  /*--------------------------------------------------------------*/
  
label1000:
  
  iter = iter+1;
  
  if (msglvl> 0)
  {
    ivout<T>(debug.logfil,1,&iter,debug.ndigit,"__naup2: **** Start of major iteration number ****");
    //        call ivout (logfil, 1, iter, ndigit, 
    // &           '_naup2: **** Start of major iteration number ****')
  }
  
  /*-----------------------------------------------------------*/
  /* Compute NP additional steps for the Arnoldi factorization. */
  /* Adjust NP since NEV might have been updated by last call  */
  /* to the shift application routine dnapps.                  */
  /*-----------------------------------------------------------*/
  
  np = kplusp-nev;
  
  if(msglvl>1)
  {
    ivout<T>(debug.logfil,1,&nev,debug.ndigit,"_naup2: The length of the current Arnoldi factorization");
    ivout<T>(debug.logfil,1,&np,debug.ndigit,"_naup2: Extend the Arnoldi factorization by");
  }
  /* if (msglvl .gt. 1) then
            call ivout (logfil, 1, nev, ndigit, 
     &     '_naup2: The length of the current Arnoldi factorization')
            call ivout (logfil, 1, np, ndigit, 
     &           '_naup2: Extend the Arnoldi factorization by')
         end if */
  
  /*----------------------------------------------------------*/
  /* Compute NP additional steps of the Arnoldi factorization */
  /*----------------------------------------------------------*/
 
  ido=0;
  
label20:
  update = true;
  
  dnaitr<T>(ido, bmat, n, nev, np, mode, resid, rnorm, v, ldv,h, ldh, ipntr, workd, info, digits);
  
  /*---------------------------------------------------*/
  /* ido!=99 implies use of reverse communication      */
  /* to compute operations involving OP and possibly B */
  /*---------------------------------------------------*/
  
  if(ido!=99) return;
  
  if(info>0)
  {
    np=info;
    mxiter=iter;
    info=-9999;
    
    /*------------*/
    /* Error Exit */
    /*------------*/
  
    goto label1200;
  }
  
  update = false;
  
  if(msglvl>1)
  {
    dvout<T>(debug.logfil,1,&rnorm,debug.ndigit,"_naup2: Corresponding B-norm of the residual");
  }
  /* if (msglvl .gt. 1) then
            call dvout (logfil, 1, rnorm, ndigit, 
     &           '_naup2: Corresponding B-norm of the residual')
         end if /*
         
  /*--------------------------------------------------------*/       
  /* Compute the eigenvalues and corresponding error bounds */
  /* of the current upper Hessenberg matrix.                */
  /*--------------------------------------------------------*/
  
 // dvout<T>(debug.logfil,kplusp,bounds,debug.ndigit,"_naup2: bounds before call to dneigh");
  dneigh<T>(rnorm, kplusp, h, ldh, ritzr, ritzi, bounds,q, ldq, workl, ierr, digits);
 // dvout<T>(debug.logfil,kplusp,bounds,debug.ndigit,"_naup2: bounds after call to dneigh");
  
  if(ierr!=0)
  {
    info=-8;   
    /*------------*/
    /* Error Exit */
    /*------------*/
  
    goto label1200;
  }
  
  /*----------------------------------------------------*/
  /* Make a copy of eigenvalues and corresponding error */
  /* bounds obtained from dneigh.                       */
  /*----------------------------------------------------*/
  
  copy<T>(kplusp,ritzr,1,&workl[kplusp*kplusp],1);
  copy<T>(kplusp,ritzi,1,&workl[kplusp*kplusp+kplusp],1);
  copy<T>(kplusp,bounds,1,&workl[kplusp*kplusp+2*kplusp],1);
  
  /*-----------------------------------------------------*/
  /* Select the wanted Ritz values and their bounds      */
  /* to be used in the convergence test.                 */
  /* The wanted part of the spectrum and corresdponing   */
  /* error bounds are in the last NEV loc. of RITZR,     */
  /* RITZI and BOUNDS respectively. The variables NEV    */
  /* and NP may be updated if the NEV-th wanted Ritz     */
  /* value has a non-zero imaginary part. In this case   */
  /* NEV is increased by one and NP decreased by one.    */
  /* NOTE: The last two arguments of dngets are no       */
  /* longer used as of version 2.1.                      */
  /*-----------------------------------------------------*/

  nev=nev0;
  np=np0;
  numcnv=nev;
  dngets<T>(ishift,which,nev,np,ritzr,ritzi,bounds, workl, &workl[np], digits);
  if(nev==(nev0+1)) numcnv = nev0+1;
  
  /*-------------------*/
  /* Convergence test. */
  /*-------------------*/
  
  copy<T>(nev,&bounds[np],1,&workl[2*np],1);
  dnconv<T>(nev,&ritzr[np],&ritzi[np],&workl[2*np],tol,nconv,digits);
  
  if(msglvl>2)
  {
    kp[0]=nev;
    kp[1]=np;
    kp[2]=numcnv;
    kp[3]=nconv;
    ivout<T>(debug.logfil,4,kp,debug.ndigit,"_naup2: NEV, NP, NUMCNV, NCONV are");
    dvout<T>(debug.logfil,kplusp,ritzr,debug.ndigit,"_naup2: Real part of the eigenvalues of H");
    dvout<T>(debug.logfil,kplusp,ritzi,debug.ndigit,"_naup2: Imaginary part of the eigenvalues of H");
    dvout<T>(debug.logfil,kplusp,bounds,debug.ndigit,"_naup2: Ritz estimates of the current NCV Ritz values");
  }
  /* if (msglvl .gt. 2) then
            kp(1) = nev
            kp(2) = np
            kp(3) = numcnv
            kp(4) = nconv
            call ivout (logfil, 4, kp, ndigit, 
     &                  '_naup2: NEV, NP, NUMCNV, NCONV are')
            call dvout (logfil, kplusp, ritzr, ndigit,
     &           '_naup2: Real part of the eigenvalues of H')
            call dvout (logfil, kplusp, ritzi, ndigit,
     &           '_naup2: Imaginary part of the eigenvalues of H')
            call dvout (logfil, kplusp, bounds, ndigit, 
     &          '_naup2: Ritz estimates of the current NCV Ritz values')
         end if */
  
  /*---------------------------------------------------------*/
  /* Count the number of unwanted Ritz values that have zero */
  /* Ritz estimates. If any Ritz estimates are equal to zero */
  /* then a leading block of H of order equal to at least    */
  /* the number of Ritz values with zero Ritz estimates has  */
  /* split off. None of these Ritz values may be removed by  */
  /* shifting. Decrease NP the number of shifts to apply. If */
  /* no shifts may be applied, then prepare to exit          */
  /*---------------------------------------------------------*/
  //std::cout << "np before loop = " << np << "\n";
  nptemp = np;
  for(j=1;j<=nptemp;j++)
  {
    if(bounds[j-1]==zero)
    {
      np=np-1;
      nev=nev+1;
    }
  }
  //std::cout << "np after loop = " << np << "\n";
  
  if((nconv>=numcnv)||(iter>mxiter)||(np==0))
  {
    if(msglvl>4)
    {
      dvout<T>(debug.logfil,kplusp,&workl[kplusp*kplusp],debug.ndigit,"_naup2: Real part of the eig computed by _neigh:");
      dvout<T>(debug.logfil,kplusp,&workl[kplusp*kplusp+kplusp],debug.ndigit,"_naup2: Imag part of the eig computed by _neigh:");
      dvout<T>(debug.logfil,kplusp,&workl[kplusp*kplusp+2*kplusp],debug.ndigit,"_naup2: Ritz eistmates computed by _neigh:");
    }
  
    /*if (msglvl .gt. 4) then
               call dvout(logfil, kplusp, workl(kplusp**2+1), ndigit,
     &             '_naup2: Real part of the eig computed by _neigh:')
               call dvout(logfil, kplusp, workl(kplusp**2+kplusp+1),
     &                     ndigit,
     &             '_naup2: Imag part of the eig computed by _neigh:')
               call dvout(logfil, kplusp, workl(kplusp**2+kplusp*2+1),
     &                     ndigit,
     &             '_naup2: Ritz eistmates computed by _neigh:')
            end if */
    
    /*------------------------------------------------*/
    /* Prepare to exit. Put the converged Ritz values */
    /* and corresponding bounds in RITZI[NCONV] and   */
    /* BOUNDS[NCONV] respectively. The sort. Be       */
    /* careful when NCONV > NP                        */
    /*------------------------------------------------*/
    
    /*---------------------------------------*/
    /* Use h( 3,1 )  as storage to communicate */
    /* rnorm to _neupd if needed             */
    /*---------------------------------------*/
    
    h[2]=rnorm;
    //h[2][0]=rnorm;
    
    /*----------------------------------------------*/
    /* To be consisent with dngets, we first do a   */
    /* pre-processing sort in order to keep complex */
    /* conjugate pairs together. This is similar    */
    /* to the pre-processing sort used in dngets    */
    /* except that the sort is done in the opposite */
    /* order.                                       */
    /*----------------------------------------------*/
    
    if(which=="LM") wprime = "SR";
    if(which=="SM") wprime = "LR";
    if(which=="LR") wprime = "SM";
    if(which=="SR") wprime = "LM";
    if(which=="LI") wprime = "SM";
    if(which=="SI") wprime = "LM";
    
    dsortc<T>(wprime,true,kplusp,ritzr,ritzi,bounds);
    
    /*---------------------------------------------*/
    /* Now sort Ritz values so that converged Ritz */
    /* values appear within the first NEV locations*/
    /* of ritzr, ritzi and bounds, and the most    */
    /* desired one appears at the front.           */
    /*---------------------------------------------*/
    
    if(which=="LM") wprime = "SM";
    if(which=="SM") wprime = "LM";
    if(which=="LR") wprime = "SR";
    if(which=="SR") wprime = "LR";
    if(which=="LI") wprime = "SI";
    if(which=="SI") wprime = "LI";
    
    dsortc<T>(wprime,true,kplusp,ritzr,ritzi,bounds);
    
    /*-----------------------------------------------*/
    /* Scale the Ritz estimate of each Ritz value    */
    /* by 1/max(eps23,magnitude of the Ritz value)   */
    /*-----------------------------------------------*/
    
    for(j=1;j<=nev0;j++)
    {
      temp=cln::max(eps23,dlapy2<T>(ritzr[j-1],ritzi[j-1]));
      bounds[j-1]=bounds[j-1]/temp;
    }
    
    /*----------------------------------------------------*/
    /* Sort the Ritz values according to the scaled Ritz  */
    /* estimates. This will push all the converged ones   */
    /* towards the front of ritzr, ritzi, bounds          */
    /* (in the case when NCONV < NEV.)                    */
    /*----------------------------------------------------*/
    
    wprime = "LR";
    dsortc<T>(wprime,true,nev0,bounds,ritzr,ritzi);
    
    /*----------------------------------------------*/
    /* Scale the Ritz estimate back to its original */
    /* value.                                       */
    /*----------------------------------------------*/
    
    for(j=1;j<=nev0;j++)
    {
      temp=cln::max(eps23,dlapy2<T>(ritzr[j-1],ritzi[j-1]));
      bounds[j-1]=bounds[j-1]*temp;
    }
    
    /*-----------------------------------------------*/
    /* Sort the converged Ritz values again so that  */
    /* the "threshold" value appears at the front of */
    /* ritzr, ritzi and bound.                       */
    /*-----------------------------------------------*/
    
    dsortc<T>(which,true,nconv,ritzr,ritzi,bounds);
    
    if(msglvl>1)
    {
      dvout<T>(debug.logfil,kplusp,ritzr,debug.ndigit,"_naup2: Sorted real part of the eigenvalues");
      dvout<T>(debug.logfil,kplusp,ritzi,debug.ndigit,"_naup2: Sorted imaginary part of the eigenvalues");
      dvout<T>(debug.logfil,kplusp,bounds,debug.ndigit,"_naup2: Sorted ritz estimates.");
    }
    /* if (msglvl .gt. 1) then
               call dvout (logfil, kplusp, ritzr, ndigit,
     &            '_naup2: Sorted real part of the eigenvalues')
               call dvout (logfil, kplusp, ritzi, ndigit,
     &            '_naup2: Sorted imaginary part of the eigenvalues')
               call dvout (logfil, kplusp, bounds, ndigit,
     &            '_naup2: Sorted ritz estimates.')
            end if */
    
    /*------------------------------------*/
    /* Max iterations have been exceeded. */
    /*------------------------------------*/
    
    if((iter>mxiter)&&(nconv<numcnv)) 
      info=1;
    
    /*---------------------*/
    /* No shifts to apply. */
    /*---------------------*/

    if((np==0)&&(nconv<numcnv))
      info=2;
    
    np=nconv;
    goto label1100;
    
  }
  else if((nconv<numcnv)&&(ishift==1))
  {
    /*-------------------------------------------------*/
    /* Do not have all the requested eigenvalues yet.  */
    /* To prevent possible stagnation, adjust the size */
    /* of NEV.                                         */
    /*-------------------------------------------------*/
    
    nevbef=nev;
    nev = nev+std::min(nconv,np/2);
    if((nev==1)&&(kplusp>=6))
      nev = kplusp/2;
    else if((nev==1)&&(kplusp>3))
      nev = 2;
    np=kplusp-nev;
    
    /*---------------------------------------*/
    /* If the size of NEV was just increased */
    /* resort the eigenvalues.               */
    /*---------------------------------------*/
    
    if(nevbef<nev)
      dngets<T>(ishift, which, nev, np, ritzr, ritzi,bounds, workl, &workl[np],digits);
    
  }
  
  if(msglvl>0)
  {
    ivout<T>(debug.logfil,1,&nconv,debug.ndigit,"_naup2: no. of 'converged' Ritz values at this iter.");
    if(msglvl>1)
    {
      kp[0]=nev;
      kp[1]=np;
      ivout<T>(debug.logfil,2,kp,debug.ndigit,"_naup2: NEV and NP are");
      dvout<T>(debug.logfil,nev,&ritzr[np],debug.ndigit,"_naup2: 'wanted' Ritz values -- real part");
      dvout<T>(debug.logfil,nev,&ritzi[np],debug.ndigit,"_naup2: 'wanted' Ritz values -- imag part");
      dvout<T>(debug.logfil,nev,&bounds[np],debug.ndigit,"_naup2: Ritz estimates of the 'wanted' values ");
    }
  }
  /*if (msglvl .gt. 0) then
            call ivout (logfil, 1, nconv, ndigit, 
     &           '_naup2: no. of "converged" Ritz values at this iter.')
            if (msglvl .gt. 1) then
               kp(1) = nev
               kp(2) = np
               call ivout (logfil, 2, kp, ndigit, 
     &              '_naup2: NEV and NP are')
               call dvout (logfil, nev, ritzr(np+1), ndigit,
     &              '_naup2: "wanted" Ritz values -- real part')
               call dvout (logfil, nev, ritzi(np+1), ndigit,
     &              '_naup2: "wanted" Ritz values -- imag part')
               call dvout (logfil, nev, bounds(np+1), ndigit,
     &              '_naup2: Ritz estimates of the "wanted" values ')
            end if
         end if */
  
  if(ishift==0)
  {
    /*----------------------------------------------------*/
    /* User specified shifts: reverse communication to    */
    /* compute the shifts. They are returned in the first */
    /* 2*NP locations of WORKL.                           */
    /*----------------------------------------------------*/
    
    ushift = true;
    ido=3;
    return;
  } 
  
label50:

    /*-------------------------------------*/
    /* Back from reverse communication;    */
    /* User specified shifts are returned  */
    /* in WORKL[2*NP]                      */
    /*-------------------------------------*/
    
    ushift = false;
    
    if(ishift==0)
    {
      /*----------------------------------*/
      /* Move the NP shifts from WORKL to */
      /* RITZR, RITZI to free up WORKL    */
      /* for non-exact shift case.        */
      /*----------------------------------*/
      
      copy<T>(np,workl,1,ritzr,1);
      copy<T>(np,&workl[np],1,ritzi,1);
    }
    
    if(msglvl>2)
    {
      ivout<T>(debug.logfil,1,&np,debug.ndigit,"_naup2: The number of shifts to apply");
      dvout<T>(debug.logfil,np,ritzr,debug.ndigit,"_naup2: Real part of the shifts");
      dvout<T>(debug.logfil,np,ritzi,debug.ndigit,"_naup2: Imaginary part of the shifts");
      if(ishift==1)
        dvout<T>(debug.logfil,np,bounds,debug.ndigit,"_naup2: Ritz estimates of the shifts");
    }
    /*if (msglvl .gt. 2) then 
            call ivout (logfil, 1, np, ndigit, 
     &                  '_naup2: The number of shifts to apply ')
            call dvout (logfil, np, ritzr, ndigit,
     &                  '_naup2: Real part of the shifts')
            call dvout (logfil, np, ritzi, ndigit,
     &                  '_naup2: Imaginary part of the shifts')
            if ( ishift .eq. 1 ) 
     &          call dvout (logfil, np, bounds, ndigit,
     &                  '_naup2: Ritz estimates of the shifts')
         end if */
    
    /*----------------------------------------------------------*/
    /* Apply the NP implicit shifts by QR bulge chasing.        */
    /* Each shift is applied to the whole upper Hessenberg      */
    /* matrix H.                                                */
    /* The first 2*N locations of WORKD are used as workspace.  */
    /*----------------------------------------------------------*/
    
    dnapps<T>(n, nev, np, ritzr, ritzi, v, ldv,h, ldh, resid, q, ldq, workl, workd, digits);
    
    /*----------------------------------------------*/
    /* Compute the B-norm of the updated residual.  */
    /* Keep B*RESID in WORKD[N] to be used in       */
    /* the first step of the next call to dnaitr.   */
    /*----------------------------------------------*/
    
    cnorm = true;
    // call second(t2);
    if(bmat=='G')
    {
      timing.nbx = timing.nbx+1;
      copy<T>(n,resid,1,&workd[n],1);
      ipntr[1]=n+1;
      ipntr[2]=1; // fixed index
      ido=2;
      
      /*----------------------------------*/
      /* Exit in order to compute B*RESID */
      /*----------------------------------*/
      
      return;
    }
    else if(bmat=='I')
      copy<T>(n,resid,1,workd,1);
    
label100:
    
    /*----------------------------------*/
    /* Back from reverse communication; */
    /* WORKD[N] := B*RESID              */
    /*----------------------------------*/
    
    if(bmat=='G')
    {
      second(timing.t3);
      timing.tmvbx += (timing.t3-timing.t2);
      std::cout << "Not yet implemented. \n";
    }
    
    if(bmat=='G')
    {
      rnorm = ddot<T>(n,resid,1,workd,1,digits);
      rnorm = cln::sqrt(cln::abs(rnorm));
    }
    else if(bmat=='I')
      rnorm = dnrm2<T>(n,resid,1);
    cnorm = false;
    
    if (msglvl>2)
    {
       dvout<T>(debug.logfil, 1, &rnorm, debug.ndigit,"_naup2: B-norm of residual for compressed factorization");
       dmout<T>(debug.logfil, nev, nev, h, ldh, debug.ndigit,"_naup2: Compressed upper Hessenberg matrix H");
    }
  
    goto label1000;
    /*
c
c     %---------------------------------------------------------------%
c     |                                                               |
c     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
c     |                                                               |
c     %---------------------------------------------------------------%
c
    */
label1100:
    mxiter=iter;
    nev=numcnv;
    
label1200:
    ido=99;
    /*
c
c     %------------%
c     | Error Exit |
c     %------------%
c
    */
    second(timing.t1);
    timing.tnaup2=timing.t1-timing.t0;
    
    /*
c
c     %---------------%
c     | End of dnaup2 |
c     %---------------%
c
    */
    return;
}


template <typename T>
void dnaup2(ARint& ido, char bmat, ARint& n,const std::string& which, int& nev, int& np,
	    T& tol, T* resid, int& mode, int iupd, int& ishift, int& mxiter,
	    T* v, int ldv, T* h, int ldh, T* ritzr, T* ritzi, T* bounds,
	    T* q, int ldq, T* workl, int* ipntr, T* workd, int& info)
{
  /*--------------------------------------------------------*/
  /* Include variables for debugging and timing information */
  /*--------------------------------------------------------*/
  
 // struct timingstruct timing;
 // struct debugstruct debug;
  
  /*------------*/
  /* Parameters */
  /*------------*/
  
  const T one = 1.0;
  const T zero = 0.0;
  
  /*---------------*/
  /* Local Scalars */
  /*---------------*/
  
  static std::string wprime;
  static bool cnorm, getv0, initv, update, ushift;
  static int ierr,iter,j,kplusp,msglvl,nconv,nevbef,nev0,np0,nptemp,numcnv;
  static T rnorm,temp,eps23;
  
  /* Local array arguments */
  
  static int kp[4];
  //std::cout << "starting np = " << np << "\n";
  if(ido==0)
  {
    second(timing.t0);
    
    msglvl=debug.mnaup2;
    
    /*-------------------------------------*/
    /* Get the machine dependent constant. */
    /*-------------------------------------*/
    
    eps23=dlamch<T>("E"); // Epsilon-machine
    eps23=pow(eps23,2.0/3.0);
 /*   dvout<T>(debug.logfil, 1, &eps23, debug.ndigit,"_dnaup2: eps23");
    T test = dlamch<T>("E");
    dvout<T>(debug.logfil, 1, &test, debug.ndigit,"_dnaup2: dlamch E");
    test = dlamch<T>("S");
    dvout<T>(debug.logfil, 1, &test, debug.ndigit,"_dnaup2: dlamch S");
    test = dlamch<T>("P");
    dvout<T>(debug.logfil, 1, &test, debug.ndigit,"_dnaup2: dlamch P");
    test = dlamch<T>("B");
    dvout<T>(debug.logfil, 1, &test, debug.ndigit,"_dnaup2: dlamch B");
    double safmin = dlamch<T>("S");
    double eps=eps23=dlamch<T>("E");
    double safmn2=pow(dlamch<T>("B"),int(log(safmin/eps)/log(dlamch<T>("B"))/2.0) );
    double safmx2=one/safmn2;
    dvout<T>(debug.logfil, 1, &safmn2, debug.ndigit,"_dnaup2: safmn2");
    dvout<T>(debug.logfil, 1, &safmx2, debug.ndigit,"_dnaup2: safmx2");*/
    
    //std::cout << "eps23 = " << eps23 << "\n";
    
    nev0=nev;
    np0=np;
    
    /*-------------------------------------*/
    /* kplusp is the bound on the largest  */
    /*        Lanczos factorization built. */
    /* nconv is the current number of      */
    /*         "converged" eigenvalues.    */
    /* iter is the counter on the current  */
    /*      iteration step.                */
    /*-------------------------------------*/
    
    kplusp = nev + np;
    nconv = 0;
    iter = 0;
    
    /*---------------------------------------*/
    /* Set flags for computing the first NEV */
    /* steps of the Arnoldi factorization    */
    /*---------------------------------------*/
    
    getv0 = true;
    update = false;
    ushift = false;
    cnorm = false;
    
  //  std::cout << "info = " << info << "\n";
    if (info!=0)
    {
      /* User provides the initial residual vector. */
      
      initv = true;
      info = 0;
    }
    else
      initv = false;
  }
    
  if (getv0)
  {
//    dvout<T>(debug.logfil, n, workd, debug.ndigit,"_dnaup2: workd before call to getv0");
    dgetv0<T>(ido,bmat,1,initv,n,1,v,ldv,resid,rnorm,ipntr,workd,info);
//    dvout<T>(debug.logfil, n, workd, debug.ndigit,"_dnaup2: workd after call to getv0");
    
 //   ivout<T>(debug.logfil,1,&ido,debug.ndigit,"_dnaup2: ido");
    if(ido!=99) return;
    
    if(rnorm == zero)
    {
      /*-----------------------------------------*/
      /* The initial vector is zero. Error exit. */
      /*-----------------------------------------*/
      
      info = -9;
      goto label1100;
    }
    getv0=false;
    ido=0;
  }
  
  /*-----------------------------------*/
  /* Back from reverse communication : */
  /* continue with update step         */
  /*-----------------------------------*/
  
  if(update) goto label20;
  
  /*-------------------------------------------*/
  /* Back from computing user specified shifts */
  /*-------------------------------------------*/
  
  if(ushift) goto label50;
  
  /*-------------------------------------*/
  /* Back from computing residual norm   */
  /* at the end of the current iteration */
  /*-------------------------------------*/
  
  if(cnorm) goto label100;

  /*----------------------------------------------------------*/
  /* Compute the first NEV steps of the Arnoldi factorization */
  /*----------------------------------------------------------*/

//  dmout<T>(debug.logfil,kplusp,kplusp,h,ldh,debug.ndigit,"_dnaup2: H before call to dnaitr");
  dnaitr<T>(ido, bmat, n, 0, nev, mode, resid, rnorm, v, ldv, h, ldh, ipntr, workd, info);
//  dmout<T>(debug.logfil,kplusp,kplusp,h,ldh,debug.ndigit,"_dnaup2: H after call to dnaitr");
  /*---------------------------------------------------*/
  /* ido!=99 implies use of reverse communication      */
  /* to compute operations involving OP and possibly B */
  /*---------------------------------------------------*/
  
  if(ido!=99) return;
  
  if(info>0)
  {
    np=info;
    mxiter=iter;
    info=-9999;
    goto label1200;
  }
  
  /*--------------------------------------------------------------*/
  /*                                                              */
  /*           M A I N   ARNOLDI  I T E R A T I O N  L O O P      */
  /*           Each iteration implicitly restarts the Arnoldi     */
  /*           factorization in place.                            */
  /*                                                              */
  /*--------------------------------------------------------------*/
  
label1000:
  
  iter = iter+1;
  
  if (msglvl> 0)
  {
    ivout<T>(debug.logfil,1,&iter,debug.ndigit,"__naup2: **** Start of major iteration number ****");
    //        call ivout (logfil, 1, iter, ndigit, 
    // &           '_naup2: **** Start of major iteration number ****')
  }
  
  /*-----------------------------------------------------------*/
  /* Compute NP additional steps for the Arnoldi factorization. */
  /* Adjust NP since NEV might have been updated by last call  */
  /* to the shift application routine dnapps.                  */
  /*-----------------------------------------------------------*/
  
  np = kplusp-nev;
  
  if(msglvl>1)
  {
    ivout<T>(debug.logfil,1,&nev,debug.ndigit,"_naup2: The length of the current Arnoldi factorization");
    ivout<T>(debug.logfil,1,&np,debug.ndigit,"_naup2: Extend the Arnoldi factorization by");
  }
  /* if (msglvl .gt. 1) then
            call ivout (logfil, 1, nev, ndigit, 
     &     '_naup2: The length of the current Arnoldi factorization')
            call ivout (logfil, 1, np, ndigit, 
     &           '_naup2: Extend the Arnoldi factorization by')
         end if */
  
  /*----------------------------------------------------------*/
  /* Compute NP additional steps of the Arnoldi factorization */
  /*----------------------------------------------------------*/
 
  ido=0;
  
label20:
  update = true;
  
  dnaitr<T>(ido, bmat, n, nev, np, mode, resid, rnorm, v, ldv,h, ldh, ipntr, workd, info);
  
  /*---------------------------------------------------*/
  /* ido!=99 implies use of reverse communication      */
  /* to compute operations involving OP and possibly B */
  /*---------------------------------------------------*/
  
  if(ido!=99) return;
  
  if(info>0)
  {
    np=info;
    mxiter=iter;
    info=-9999;
    
    /*------------*/
    /* Error Exit */
    /*------------*/
  
    goto label1200;
  }
  
  update = false;
  
  if(msglvl>1)
  {
    dvout<T>(debug.logfil,1,&rnorm,debug.ndigit,"_naup2: Corresponding B-norm of the residual");
  }
  /* if (msglvl .gt. 1) then
            call dvout (logfil, 1, rnorm, ndigit, 
     &           '_naup2: Corresponding B-norm of the residual')
         end if /*
         
  /*--------------------------------------------------------*/       
  /* Compute the eigenvalues and corresponding error bounds */
  /* of the current upper Hessenberg matrix.                */
  /*--------------------------------------------------------*/
  
 // dvout<T>(debug.logfil,kplusp,bounds,debug.ndigit,"_naup2: bounds before call to dneigh");
  dneigh<T>(rnorm, kplusp, h, ldh, ritzr, ritzi, bounds,q, ldq, workl, ierr);
 // dvout<T>(debug.logfil,kplusp,bounds,debug.ndigit,"_naup2: bounds after call to dneigh");
  
  if(ierr!=0)
  {
    info=-8;   
    /*------------*/
    /* Error Exit */
    /*------------*/
  
    goto label1200;
  }
  
  /*----------------------------------------------------*/
  /* Make a copy of eigenvalues and corresponding error */
  /* bounds obtained from dneigh.                       */
  /*----------------------------------------------------*/
  
  copy<T>(kplusp,ritzr,1,&workl[kplusp*kplusp],1);
  copy<T>(kplusp,ritzi,1,&workl[kplusp*kplusp+kplusp],1);
  copy<T>(kplusp,bounds,1,&workl[kplusp*kplusp+2*kplusp],1);
  
  /*-----------------------------------------------------*/
  /* Select the wanted Ritz values and their bounds      */
  /* to be used in the convergence test.                 */
  /* The wanted part of the spectrum and corresdponing   */
  /* error bounds are in the last NEV loc. of RITZR,     */
  /* RITZI and BOUNDS respectively. The variables NEV    */
  /* and NP may be updated if the NEV-th wanted Ritz     */
  /* value has a non-zero imaginary part. In this case   */
  /* NEV is increased by one and NP decreased by one.    */
  /* NOTE: The last two arguments of dngets are no       */
  /* longer used as of version 2.1.                      */
  /*-----------------------------------------------------*/

  nev=nev0;
  np=np0;
  numcnv=nev;
  dngets<T>(ishift,which,nev,np,ritzr,ritzi,bounds, workl, &workl[np]);
  if(nev==(nev0+1)) numcnv = nev0+1;
  
  /*-------------------*/
  /* Convergence test. */
  /*-------------------*/
  
  copy<T>(nev,&bounds[np],1,&workl[2*np],1);
  dnconv<T>(nev,&ritzr[np],&ritzi[np],&workl[2*np],tol,nconv);
  
  if(msglvl>2)
  {
    kp[0]=nev;
    kp[1]=np;
    kp[2]=numcnv;
    kp[3]=nconv;
    ivout<T>(debug.logfil,4,kp,debug.ndigit,"_naup2: NEV, NP, NUMCNV, NCONV are");
    dvout<T>(debug.logfil,kplusp,ritzr,debug.ndigit,"_naup2: Real part of the eigenvalues of H");
    dvout<T>(debug.logfil,kplusp,ritzi,debug.ndigit,"_naup2: Imaginary part of the eigenvalues of H");
    dvout<T>(debug.logfil,kplusp,bounds,debug.ndigit,"_naup2: Ritz estimates of the current NCV Ritz values");
  }
  /* if (msglvl .gt. 2) then
            kp(1) = nev
            kp(2) = np
            kp(3) = numcnv
            kp(4) = nconv
            call ivout (logfil, 4, kp, ndigit, 
     &                  '_naup2: NEV, NP, NUMCNV, NCONV are')
            call dvout (logfil, kplusp, ritzr, ndigit,
     &           '_naup2: Real part of the eigenvalues of H')
            call dvout (logfil, kplusp, ritzi, ndigit,
     &           '_naup2: Imaginary part of the eigenvalues of H')
            call dvout (logfil, kplusp, bounds, ndigit, 
     &          '_naup2: Ritz estimates of the current NCV Ritz values')
         end if */
  
  /*---------------------------------------------------------*/
  /* Count the number of unwanted Ritz values that have zero */
  /* Ritz estimates. If any Ritz estimates are equal to zero */
  /* then a leading block of H of order equal to at least    */
  /* the number of Ritz values with zero Ritz estimates has  */
  /* split off. None of these Ritz values may be removed by  */
  /* shifting. Decrease NP the number of shifts to apply. If */
  /* no shifts may be applied, then prepare to exit          */
  /*---------------------------------------------------------*/
  //std::cout << "np before loop = " << np << "\n";
  nptemp = np;
  for(j=1;j<=nptemp;j++)
  {
    if(bounds[j-1]==zero)
    {
      np=np-1;
      nev=nev+1;
    }
  }
  //std::cout << "np after loop = " << np << "\n";
  
  if((nconv>=numcnv)||(iter>mxiter)||(np==0))
  {
    if(msglvl>4)
    {
      dvout<T>(debug.logfil,kplusp,&workl[kplusp*kplusp],debug.ndigit,"_naup2: Real part of the eig computed by _neigh:");
      dvout<T>(debug.logfil,kplusp,&workl[kplusp*kplusp+kplusp],debug.ndigit,"_naup2: Imag part of the eig computed by _neigh:");
      dvout<T>(debug.logfil,kplusp,&workl[kplusp*kplusp+2*kplusp],debug.ndigit,"_naup2: Ritz eistmates computed by _neigh:");
    }
  
    /*if (msglvl .gt. 4) then
               call dvout(logfil, kplusp, workl(kplusp**2+1), ndigit,
     &             '_naup2: Real part of the eig computed by _neigh:')
               call dvout(logfil, kplusp, workl(kplusp**2+kplusp+1),
     &                     ndigit,
     &             '_naup2: Imag part of the eig computed by _neigh:')
               call dvout(logfil, kplusp, workl(kplusp**2+kplusp*2+1),
     &                     ndigit,
     &             '_naup2: Ritz eistmates computed by _neigh:')
            end if */
    
    /*------------------------------------------------*/
    /* Prepare to exit. Put the converged Ritz values */
    /* and corresponding bounds in RITZI[NCONV] and   */
    /* BOUNDS[NCONV] respectively. The sort. Be       */
    /* careful when NCONV > NP                        */
    /*------------------------------------------------*/
    
    /*---------------------------------------*/
    /* Use h( 3,1 )  as storage to communicate */
    /* rnorm to _neupd if needed             */
    /*---------------------------------------*/
    
    h[2]=rnorm;
    //h[2][0]=rnorm;
    
    /*----------------------------------------------*/
    /* To be consisent with dngets, we first do a   */
    /* pre-processing sort in order to keep complex */
    /* conjugate pairs together. This is similar    */
    /* to the pre-processing sort used in dngets    */
    /* except that the sort is done in the opposite */
    /* order.                                       */
    /*----------------------------------------------*/
    
    if(which=="LM") wprime = "SR";
    if(which=="SM") wprime = "LR";
    if(which=="LR") wprime = "SM";
    if(which=="SR") wprime = "LM";
    if(which=="LI") wprime = "SM";
    if(which=="SI") wprime = "LM";
    
    dsortc<T>(wprime,true,kplusp,ritzr,ritzi,bounds);
    
    /*---------------------------------------------*/
    /* Now sort Ritz values so that converged Ritz */
    /* values appear within the first NEV locations*/
    /* of ritzr, ritzi and bounds, and the most    */
    /* desired one appears at the front.           */
    /*---------------------------------------------*/
    
    if(which=="LM") wprime = "SM";
    if(which=="SM") wprime = "LM";
    if(which=="LR") wprime = "SR";
    if(which=="SR") wprime = "LR";
    if(which=="LI") wprime = "SI";
    if(which=="SI") wprime = "LI";
    
    dsortc<T>(wprime,true,kplusp,ritzr,ritzi,bounds);
    
    /*-----------------------------------------------*/
    /* Scale the Ritz estimate of each Ritz value    */
    /* by 1/max(eps23,magnitude of the Ritz value)   */
    /*-----------------------------------------------*/
    
    for(j=1;j<=nev0;j++)
    {
      temp=std::max(eps23,dlapy2<T>(ritzr[j-1],ritzi[j-1]));
      bounds[j-1]=bounds[j-1]/temp;
    }
    
    /*----------------------------------------------------*/
    /* Sort the Ritz values according to the scaled Ritz  */
    /* estimates. This will push all the converged ones   */
    /* towards the front of ritzr, ritzi, bounds          */
    /* (in the case when NCONV < NEV.)                    */
    /*----------------------------------------------------*/
    
    wprime = "LR";
    dsortc<T>(wprime,true,nev0,bounds,ritzr,ritzi);
    
    /*----------------------------------------------*/
    /* Scale the Ritz estimate back to its original */
    /* value.                                       */
    /*----------------------------------------------*/
    
    for(j=1;j<=nev0;j++)
    {
      temp=std::max(eps23,dlapy2<T>(ritzr[j-1],ritzi[j-1]));
      bounds[j-1]=bounds[j-1]*temp;
    }
    
    /*-----------------------------------------------*/
    /* Sort the converged Ritz values again so that  */
    /* the "threshold" value appears at the front of */
    /* ritzr, ritzi and bound.                       */
    /*-----------------------------------------------*/
    
    dsortc<T>(which,true,nconv,ritzr,ritzi,bounds);
    
    if(msglvl>1)
    {
      dvout<T>(debug.logfil,kplusp,ritzr,debug.ndigit,"_naup2: Sorted real part of the eigenvalues");
      dvout<T>(debug.logfil,kplusp,ritzi,debug.ndigit,"_naup2: Sorted imaginary part of the eigenvalues");
      dvout<T>(debug.logfil,kplusp,bounds,debug.ndigit,"_naup2: Sorted ritz estimates.");
    }
    /* if (msglvl .gt. 1) then
               call dvout (logfil, kplusp, ritzr, ndigit,
     &            '_naup2: Sorted real part of the eigenvalues')
               call dvout (logfil, kplusp, ritzi, ndigit,
     &            '_naup2: Sorted imaginary part of the eigenvalues')
               call dvout (logfil, kplusp, bounds, ndigit,
     &            '_naup2: Sorted ritz estimates.')
            end if */
    
    /*------------------------------------*/
    /* Max iterations have been exceeded. */
    /*------------------------------------*/
    
    if((iter>mxiter)&&(nconv<numcnv)) 
      info=1;
    
    /*---------------------*/
    /* No shifts to apply. */
    /*---------------------*/

    if((np==0)&&(nconv<numcnv))
      info=2;
    
    np=nconv;
    goto label1100;
    
  }
  else if((nconv<numcnv)&&(ishift==1))
  {
    /*-------------------------------------------------*/
    /* Do not have all the requested eigenvalues yet.  */
    /* To prevent possible stagnation, adjust the size */
    /* of NEV.                                         */
    /*-------------------------------------------------*/
    
    nevbef=nev;
    nev = nev+std::min(nconv,np/2);
    if((nev==1)&&(kplusp>=6))
      nev = kplusp/2;
    else if((nev==1)&&(kplusp>3))
      nev = 2;
    np=kplusp-nev;
    
    /*---------------------------------------*/
    /* If the size of NEV was just increased */
    /* resort the eigenvalues.               */
    /*---------------------------------------*/
    
    if(nevbef<nev)
      dngets<T>(ishift, which, nev, np, ritzr, ritzi,bounds, workl, &workl[np]);
    
  }
  
  if(msglvl>0)
  {
    ivout<T>(debug.logfil,1,&nconv,debug.ndigit,"_naup2: no. of 'converged' Ritz values at this iter.");
    if(msglvl>1)
    {
      kp[0]=nev;
      kp[1]=np;
      ivout<T>(debug.logfil,2,kp,debug.ndigit,"_naup2: NEV and NP are");
      dvout<T>(debug.logfil,nev,&ritzr[np],debug.ndigit,"_naup2: 'wanted' Ritz values -- real part");
      dvout<T>(debug.logfil,nev,&ritzi[np],debug.ndigit,"_naup2: 'wanted' Ritz values -- imag part");
      dvout<T>(debug.logfil,nev,&bounds[np],debug.ndigit,"_naup2: Ritz estimates of the 'wanted' values ");
    }
  }
  /*if (msglvl .gt. 0) then
            call ivout (logfil, 1, nconv, ndigit, 
     &           '_naup2: no. of "converged" Ritz values at this iter.')
            if (msglvl .gt. 1) then
               kp(1) = nev
               kp(2) = np
               call ivout (logfil, 2, kp, ndigit, 
     &              '_naup2: NEV and NP are')
               call dvout (logfil, nev, ritzr(np+1), ndigit,
     &              '_naup2: "wanted" Ritz values -- real part')
               call dvout (logfil, nev, ritzi(np+1), ndigit,
     &              '_naup2: "wanted" Ritz values -- imag part')
               call dvout (logfil, nev, bounds(np+1), ndigit,
     &              '_naup2: Ritz estimates of the "wanted" values ')
            end if
         end if */
  
  if(ishift==0)
  {
    /*----------------------------------------------------*/
    /* User specified shifts: reverse communication to    */
    /* compute the shifts. They are returned in the first */
    /* 2*NP locations of WORKL.                           */
    /*----------------------------------------------------*/
    
    ushift = true;
    ido=3;
    return;
  } 
  
label50:

    /*-------------------------------------*/
    /* Back from reverse communication;    */
    /* User specified shifts are returned  */
    /* in WORKL[2*NP]                      */
    /*-------------------------------------*/
    
    ushift = false;
    
    if(ishift==0)
    {
      /*----------------------------------*/
      /* Move the NP shifts from WORKL to */
      /* RITZR, RITZI to free up WORKL    */
      /* for non-exact shift case.        */
      /*----------------------------------*/
      
      copy<T>(np,workl,1,ritzr,1);
      copy<T>(np,&workl[np],1,ritzi,1);
    }
    
    if(msglvl>2)
    {
      ivout<T>(debug.logfil,1,&np,debug.ndigit,"_naup2: The number of shifts to apply");
      dvout<T>(debug.logfil,np,ritzr,debug.ndigit,"_naup2: Real part of the shifts");
      dvout<T>(debug.logfil,np,ritzi,debug.ndigit,"_naup2: Imaginary part of the shifts");
      if(ishift==1)
        dvout<T>(debug.logfil,np,bounds,debug.ndigit,"_naup2: Ritz estimates of the shifts");
    }
    /*if (msglvl .gt. 2) then 
            call ivout (logfil, 1, np, ndigit, 
     &                  '_naup2: The number of shifts to apply ')
            call dvout (logfil, np, ritzr, ndigit,
     &                  '_naup2: Real part of the shifts')
            call dvout (logfil, np, ritzi, ndigit,
     &                  '_naup2: Imaginary part of the shifts')
            if ( ishift .eq. 1 ) 
     &          call dvout (logfil, np, bounds, ndigit,
     &                  '_naup2: Ritz estimates of the shifts')
         end if */
    
    /*----------------------------------------------------------*/
    /* Apply the NP implicit shifts by QR bulge chasing.        */
    /* Each shift is applied to the whole upper Hessenberg      */
    /* matrix H.                                                */
    /* The first 2*N locations of WORKD are used as workspace.  */
    /*----------------------------------------------------------*/
    
    dnapps<T>(n, nev, np, ritzr, ritzi, v, ldv,h, ldh, resid, q, ldq, workl, workd);
    
    /*----------------------------------------------*/
    /* Compute the B-norm of the updated residual.  */
    /* Keep B*RESID in WORKD[N] to be used in       */
    /* the first step of the next call to dnaitr.   */
    /*----------------------------------------------*/
    
    cnorm = true;
    // call second(t2);
    if(bmat=='G')
    {
      timing.nbx = timing.nbx+1;
      copy<T>(n,resid,1,&workd[n],1);
      ipntr[1]=n+1;
      ipntr[2]=1; // fixed index
      ido=2;
      
      /*----------------------------------*/
      /* Exit in order to compute B*RESID */
      /*----------------------------------*/
      
      return;
    }
    else if(bmat=='I')
      copy<T>(n,resid,1,workd,1);
    
label100:
    
    /*----------------------------------*/
    /* Back from reverse communication; */
    /* WORKD[N] := B*RESID              */
    /*----------------------------------*/
    
    if(bmat=='G')
    {
      second(timing.t3);
      timing.tmvbx += (timing.t3-timing.t2);
      std::cout << "Not yet implemented. \n";
    }
    
    if(bmat=='G')
    {
      rnorm = ddot<T>(n,resid,1,workd,1);
      rnorm = sqrt(abs(rnorm));
    }
    else if(bmat=='I')
      rnorm = dnrm2<T>(n,resid,1);
    cnorm = false;
    
    if (msglvl>2)
    {
       dvout<T>(debug.logfil, 1, &rnorm, debug.ndigit,"_naup2: B-norm of residual for compressed factorization");
       dmout<T>(debug.logfil, nev, nev, h, ldh, debug.ndigit,"_naup2: Compressed upper Hessenberg matrix H");
    }
  
    goto label1000;
    /*
c
c     %---------------------------------------------------------------%
c     |                                                               |
c     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
c     |                                                               |
c     %---------------------------------------------------------------%
c
    */
label1100:
    mxiter=iter;
    nev=numcnv;
    
label1200:
    ido=99;
    /*
c
c     %------------%
c     | Error Exit |
c     %------------%
c
    */
    second(timing.t1);
    timing.tnaup2=timing.t1-timing.t0;
    
    /*
c
c     %---------------%
c     | End of dnaup2 |
c     %---------------%
c
    */
    return;
}


#endif
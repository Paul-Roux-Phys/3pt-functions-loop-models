/*
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dngets
c
c\Description: 
c  Given the eigenvalues of the upper Hessenberg matrix H,
c  computes the NP shifts AMU that are zeros of the polynomial of 
c  degree NP which filters out components of the unwanted eigenvectors
c  corresponding to the AMU's based on some given criteria.
c
c  NOTE: call this even in the case of user specified shifts in order
c  to sort the eigenvalues, and error bounds of H for later use.
c
c\Usage:
c  call dngets
c     ( ISHIFT, WHICH, KEV, NP, RITZR, RITZI, BOUNDS, SHIFTR, SHIFTI )
c
c\Arguments
c  ISHIFT  Integer.  (INPUT)
c          Method for selecting the implicit shifts at each iteration.
c          ISHIFT = 0: user specified shifts
c          ISHIFT = 1: exact shift with respect to the matrix H.
c
c  WHICH   Character*2.  (INPUT)
c          Shift selection criteria.
c          'LM' -> want the KEV eigenvalues of largest magnitude.
c          'SM' -> want the KEV eigenvalues of smallest magnitude.
c          'LR' -> want the KEV eigenvalues of largest real part.
c          'SR' -> want the KEV eigenvalues of smallest real part.
c          'LI' -> want the KEV eigenvalues of largest imaginary part.
c          'SI' -> want the KEV eigenvalues of smallest imaginary part.
c
c  KEV      Integer.  (INPUT/OUTPUT)
c           INPUT: KEV+NP is the size of the matrix H.
c           OUTPUT: Possibly increases KEV by one to keep complex conjugate
c           pairs together.
c
c  NP       Integer.  (INPUT/OUTPUT)
c           Number of implicit shifts to be computed.
c           OUTPUT: Possibly decreases NP by one to keep complex conjugate
c           pairs together.
c
c  RITZR,  Double precision array of length KEV+NP.  (INPUT/OUTPUT)
c  RITZI   On INPUT, RITZR and RITZI contain the real and imaginary 
c          parts of the eigenvalues of H.
c          On OUTPUT, RITZR and RITZI are sorted so that the unwanted
c          eigenvalues are in the first NP locations and the wanted
c          portion is in the last KEV locations.  When exact shifts are 
c          selected, the unwanted part corresponds to the shifts to 
c          be applied. Also, if ISHIFT .eq. 1, the unwanted eigenvalues
c          are further sorted so that the ones with largest Ritz values
c          are first.
c
c  BOUNDS  Double precision array of length KEV+NP.  (INPUT/OUTPUT)
c          Error bounds corresponding to the ordering in RITZ.
c
c  SHIFTR, SHIFTI  *** USE deprecated as of version 2.1. ***
c  
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     dsortc  ARPACK sorting routine.
c     dcopy   Level 1 BLAS that copies one vector to another .
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
c     xx/xx/92: Version ' 2.1'
c
c\SCCS Information: @(#) 
c FILE: ngets.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
c
c\Remarks
c     1. xxxx
c
c\EndLib
c
c-----------------------------------------------------------------------
c
*/

#include "cstat.h"
#include "cdebug.h"
#include "fortranfuncs.h"
#include "cdsortc.h"
#include <string>
#include <cln/real.h>
#include <cln/output.h>
#include <cln/real_io.h>

/* CLN version */

template<typename T>
void dngets(int& ishift, const std::string& which, int& kev, int& np, T* ritzr, T* ritzi, T* bounds, T* shiftr, T* shifti, int digits)
{
  cln::float_format_t precision = cln::float_format(digits);
    /* Parameters */
  const cln::cl_F one=cln::cl_float(1,precision);
  const cln::cl_F zero=cln::cl_float(0,precision);
  
  /*
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
  */
  
  int msglvl;
  
  /*
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%
c 
  */
  
  second(timing.t0);
  msglvl=debug.mngets;
  
  /*
c 
c     %----------------------------------------------------%
c     | LM, SM, LR, SR, LI, SI case.                       |
c     | Sort the eigenvalues of H into the desired order   |
c     | and apply the resulting order to BOUNDS.           |
c     | The eigenvalues are sorted so that the wanted part |
c     | are always in the last KEV locations.              |
c     | We first do a pre-processing sort in order to keep |
c     | complex conjugate pairs together                   |
c     %----------------------------------------------------%
c
  */
  
  //dvout<T>(debug.logfil,kev+np,ritzr,debug.ndigit,"_ngets: ritzr before call to sortc");
//  dvout<T>(debug.logfil,kev+np,bounds,debug.ndigit,"_ngets: bounds before call to sortc");
  
//  std::cout << "before call to dsortc, which = " << which << "\n";
  if(which=="LM")
    dsortc<T>("LR",true,kev+np,ritzr,ritzi,bounds);
  else if(which=="SM")
    dsortc<T>("SR",true,kev+np,ritzr,ritzi,bounds);
  else if(which=="LR")
    dsortc<T>("LM",true,kev+np,ritzr,ritzi,bounds);
  else if(which=="SR")
    dsortc<T>("SM",true,kev+np,ritzr,ritzi,bounds);
  else if(which=="LI")
    dsortc<T>("LM",true,kev+np,ritzr,ritzi,bounds);
  else if(which=="SI")
    dsortc<T>("SM",true,kev+np,ritzr,ritzi,bounds);
  
  //dvout<T>(debug.logfil,kev+np,ritzr,debug.ndigit,"_ngets: ritzr before 2nd call to sortc");
//  dvout<T>(debug.logfil,kev+np,bounds,debug.ndigit,"_ngets: bounds before 2nd call to sortc");  
  
  dsortc<T>(which,true,kev+np,ritzr,ritzi,bounds);

//  dvout<T>(debug.logfil,kev+np,ritzr,debug.ndigit,"_ngets: ritzr after call to sortc");
//  dvout<T>(debug.logfil,kev+np,bounds,debug.ndigit,"_ngets: bounds after call to sortc");  
  /*
c     
c     %-------------------------------------------------------%
c     | Increase KEV by one if the ( ritzr(np),ritzi(np) )    |
c     | = ( ritzr(np+1),-ritzi(np+1) ) and ritz(np) .ne. zero |
c     | Accordingly decrease NP by one. In other words keep   |
c     | complex conjugate pairs together.                     |
c     %-------------------------------------------------------%
c 
  */
  
  if( ((ritzr[np+1-1]-ritzr[np-1])==zero)&&((ritzi[np+1-1]+ritzi[np-1])==zero) )
  {
    np-=1;
    kev+=1;
  }
  
  if(ishift==1)
  {
    /*
c     
c        %-------------------------------------------------------%
c        | Sort the unwanted Ritz values used as shifts so that  |
c        | the ones with largest Ritz estimates are first        |
c        | This will tend to minimize the effects of the         |
c        | forward instability of the iteration when they shifts |
c        | are applied in subroutine dnapps.                     |
c        | Be careful and use 'SR' since we want to sort BOUNDS! |
c        %-------------------------------------------------------%
c  
    */
    
    dsortc<T>("SR",true,np,bounds,ritzr,ritzi);
  }
  
  second(timing.t1);
  timing.tngets+=timing.t1-timing.t0;
  
  if(msglvl>0)
  {
    ivout<T>(debug.logfil,1,&kev,debug.ndigit,"_ngets: KEV is");
    ivout<T>(debug.logfil,1,&np, debug.ndigit, "_ngets: NP is");
    dvout<T>(debug.logfil,kev+np,ritzr,debug.ndigit,"_ngets: Eigenvalues of current H matrix -- real part");
    dvout<T>(debug.logfil,kev+np,ritzi,debug.ndigit,"_ngets: Eigenvalues of current H matrix -- imag part");
    dvout<T>(debug.logfil,kev+np,bounds,debug.ndigit,"_ngets: Ritz estimates of the current KEV+NP Ritz values");
  }
  
  return;
  /*
c     
c     %---------------%
c     | End of dngets |
c     %---------------%
c  
  */
}

template<typename T>
void dngets(int& ishift, const std::string& which, int& kev, int& np, T* ritzr, T* ritzi, T* bounds, T* shiftr, T* shifti)
{
  const T one=1.0;
  const T zero=0.0;
  
  /*
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
  */
  
 // std::cout << "called wrong dngets \n";
 // exit(1);
  
  int msglvl;
  
  /*
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%
c 
  */
  
  second(timing.t0);
  msglvl=debug.mngets;
  
  /*
c 
c     %----------------------------------------------------%
c     | LM, SM, LR, SR, LI, SI case.                       |
c     | Sort the eigenvalues of H into the desired order   |
c     | and apply the resulting order to BOUNDS.           |
c     | The eigenvalues are sorted so that the wanted part |
c     | are always in the last KEV locations.              |
c     | We first do a pre-processing sort in order to keep |
c     | complex conjugate pairs together                   |
c     %----------------------------------------------------%
c
  */
  
  //dvout<T>(debug.logfil,kev+np,ritzr,debug.ndigit,"_ngets: ritzr before call to sortc");
//  dvout<T>(debug.logfil,kev+np,bounds,debug.ndigit,"_ngets: bounds before call to sortc");
  
//  std::cout << "before call to dsortc, which = " << which << "\n";
  if(which=="LM")
    dsortc<T>("LR",true,kev+np,ritzr,ritzi,bounds);
  else if(which=="SM")
    dsortc<T>("SR",true,kev+np,ritzr,ritzi,bounds);
  else if(which=="LR")
    dsortc<T>("LM",true,kev+np,ritzr,ritzi,bounds);
  else if(which=="SR")
    dsortc<T>("SM",true,kev+np,ritzr,ritzi,bounds);
  else if(which=="LI")
    dsortc<T>("LM",true,kev+np,ritzr,ritzi,bounds);
  else if(which=="SI")
    dsortc<T>("SM",true,kev+np,ritzr,ritzi,bounds);
  
  //dvout<T>(debug.logfil,kev+np,ritzr,debug.ndigit,"_ngets: ritzr before 2nd call to sortc");
//  dvout<T>(debug.logfil,kev+np,bounds,debug.ndigit,"_ngets: bounds before 2nd call to sortc");  
  
  dsortc<T>(which,true,kev+np,ritzr,ritzi,bounds);

//  dvout<T>(debug.logfil,kev+np,ritzr,debug.ndigit,"_ngets: ritzr after call to sortc");
//  dvout<T>(debug.logfil,kev+np,bounds,debug.ndigit,"_ngets: bounds after call to sortc");  
  /*
c     
c     %-------------------------------------------------------%
c     | Increase KEV by one if the ( ritzr(np),ritzi(np) )    |
c     | = ( ritzr(np+1),-ritzi(np+1) ) and ritz(np) .ne. zero |
c     | Accordingly decrease NP by one. In other words keep   |
c     | complex conjugate pairs together.                     |
c     %-------------------------------------------------------%
c 
  */
  
  if( ((ritzr[np+1-1]-ritzr[np-1])==zero)&&((ritzi[np+1-1]+ritzi[np-1])==zero) )
  {
    np-=1;
    kev+=1;
  }
  
  if(ishift==1)
  {
    /*
c     
c        %-------------------------------------------------------%
c        | Sort the unwanted Ritz values used as shifts so that  |
c        | the ones with largest Ritz estimates are first        |
c        | This will tend to minimize the effects of the         |
c        | forward instability of the iteration when they shifts |
c        | are applied in subroutine dnapps.                     |
c        | Be careful and use 'SR' since we want to sort BOUNDS! |
c        %-------------------------------------------------------%
c  
    */
    
    dsortc<T>("SR",true,np,bounds,ritzr,ritzi);
  }
  
  second(timing.t1);
  timing.tngets+=timing.t1-timing.t0;
  
  if(msglvl>0)
  {
    ivout<T>(debug.logfil,1,&kev,debug.ndigit,"_ngets: KEV is");
    ivout<T>(debug.logfil,1,&np, debug.ndigit, "_ngets: NP is");
    dvout<T>(debug.logfil,kev+np,ritzr,debug.ndigit,"_ngets: Eigenvalues of current H matrix -- real part");
    dvout<T>(debug.logfil,kev+np,ritzi,debug.ndigit,"_ngets: Eigenvalues of current H matrix -- imag part");
    dvout<T>(debug.logfil,kev+np,bounds,debug.ndigit,"_ngets: Ritz estimates of the current KEV+NP Ritz values");
  }
  
  return;
  /*
c     
c     %---------------%
c     | End of dngets |
c     %---------------%
c  
  */
}
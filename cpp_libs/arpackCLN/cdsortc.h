/*
c\Name: dsortc
c
c\Description:
c  Sorts the complex array in XREAL and XIMAG into the order 
c  specified by WHICH and optionally applies the permutation to the
c  real array Y. It is assumed that if an element of XIMAG is
c  nonzero, then its negative is also an element. In other words,
c  both members of a complex conjugate pair are to be sorted and the
c  pairs are kept adjacent to each other.
c
c\Usage:
c  call dsortc
c     ( WHICH, APPLY, N, XREAL, XIMAG, Y )
c
c\Arguments
c  WHICH   Character*2.  (Input)
c          'LM' -> sort XREAL,XIMAG into increasing order of magnitude.
c          'SM' -> sort XREAL,XIMAG into decreasing order of magnitude.
c          'LR' -> sort XREAL into increasing order of algebraic.
c          'SR' -> sort XREAL into decreasing order of algebraic.
c          'LI' -> sort XIMAG into increasing order of magnitude.
c          'SI' -> sort XIMAG into decreasing order of magnitude.
c          NOTE: If an element of XIMAG is non-zero, then its negative
c                is also an element.
c
c  APPLY   Logical.  (Input)
c          APPLY = .TRUE.  -> apply the sorted order to array Y.
c          APPLY = .FALSE. -> do not apply the sorted order to array Y.
c
c  N       Integer.  (INPUT)
c          Size of the arrays.
c
c  XREAL,  Double precision array of length N.  (INPUT/OUTPUT)
c  XIMAG   Real and imaginary part of the array to be sorted.
c
c  Y       Double precision array of length N.  (INPUT/OUTPUT)
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
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
c               Adapted from the sort routine in LANSO.
c
c\SCCS Information: @(#) 
c FILE: sortc.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
c
*/

#ifndef CDSORTC_H
#define CDSORTC_H

#include<string>
#include "fortranfuncs.h"
#include <cln/real.h>
#include <cln/output.h>
#include <cln/real_io.h>


template<typename T>
void dsortc(const std::string which, bool apply, int n, T* xreal, T* ximag, T* y)
{
  /*
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
  */
  
  int i, igap, j;
  T temp, temp1, temp2;
  
 // std::cout << "called wrong dsortc \n";
 // exit(1);
  /*
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
  */
  
  igap=n/2;
  
 // std::cout << "which = " << which << "\n";
  
  if(which=="LM")
  {
    /*
c
c        %------------------------------------------------------%
c        | Sort XREAL,XIMAG into increasing order of magnitude. |
c        %------------------------------------------------------%
c
    */
label10:
    if(igap==0) return;
    
    for(i=igap;i<=(n-1);i++)
    {
      j=i-igap;
label20:
      
      if(j<0) continue;
      
      temp1=dlapy2<T>(xreal[j],ximag[j]);
      temp2=dlapy2<T>(xreal[j+igap],ximag[j+igap]);
      
      if(temp1>temp2)
      {
	temp=xreal[j];
	xreal[j]=xreal[j+igap];
	xreal[j+igap]=temp;
	
	temp=ximag[j];
	ximag[j]=ximag[j+igap];
	ximag[j+igap]=temp;
	
	if(apply)
	{
	  temp=y[j];
	  y[j]=y[j+igap];
	  y[j+igap]=temp;
	}
      }
      else
	continue; //goto label30;
      j=j-igap;
      goto label20;
    }
label30:
      igap=igap/2;
      goto label10;
  }
    else if(which=="SM")
    {
      /*
c
c        %------------------------------------------------------%
c        | Sort XREAL,XIMAG into decreasing order of magnitude. |
c        %------------------------------------------------------%
c
      */
      
label40:
      if(igap==0) return;
      
      for(i=igap;i<=(n-1);i++)
      {
	j=i-igap;
label50:
        
        if(j<0) continue; //goto label60;
	
	temp1=dlapy2<T>(xreal[j],ximag[j]);
	temp2=dlapy2<T>(xreal[j+igap],ximag[j+igap]);
	
	if(temp1<temp2)
	{
	  temp=xreal[j];
	  xreal[j]=xreal[j+igap];
	  xreal[j+igap]=temp;
	  
	  temp=ximag[j];
	  ximag[j]=ximag[j+igap];
	  ximag[j+igap]=temp;
	  
	  if(apply)
	  {
	    temp=y[j];
	    y[j]=y[j+igap];
	    y[j+igap]=temp;
	  }
	}
	else
	  continue; //goto label60;
	j-=igap;
	goto label50;
      }
label60:
      igap=igap/2;
      goto label40;
    }
    else if(which=="LR")
    {
      /*
c
c        %------------------------------------------------%
c        | Sort XREAL into increasing order of algebraic. |
c        %------------------------------------------------%
c
      */
label70:
  //    ivout<T>(debug.logfil,1,&igap,debug.ndigit,"sortc_: igap=");
  //    dvout<T>(debug.logfil,n,xreal,debug.ndigit,"sortc_: xreal=");
  //    dvout<T>(debug.logfil,n,ximag,debug.ndigit,"sortc_: ximag=");
      if(igap==0) return;
	
      for(i=igap;i<=(n-1);i++)
      {
	j=i-igap;
    //    ivout<T>(debug.logfil,1,&i,debug.ndigit,"sortc_: i=");
label80:
    //    ivout<T>(debug.logfil,1,&j,debug.ndigit,"sortc_: j=");
        if(j<0) continue; //goto label90;
	  
	if(xreal[j]>xreal[j+igap])
	{
	  temp=xreal[j];
	  xreal[j]=xreal[j+igap];
	  xreal[j+igap]=temp;
	  
	  temp=ximag[j];
	  ximag[j]=ximag[j+igap];
	  ximag[j+igap]=temp;
	    
	  if(apply)
	  {
	    temp=y[j];
	    y[j]=y[j+igap];
	    y[j+igap]=temp;
	  }
	}
	else
	  continue; //goto label90;
	j=j-igap;
	goto label80;
      }

label90:
      igap=igap/2;
      goto label70;
    }
    else if(which=="SR")
    {
      /*
c
c        %------------------------------------------------%
c        | Sort XREAL into decreasing order of algebraic. |
c        %------------------------------------------------%
c
      */
label100:
      if(igap==0) return;
      for(i=igap;i<=(n-1);i++)
      {
	j=i-igap;
label110:

        if(j<0) continue; //goto label120;
	
	if(xreal[j]<xreal[j+igap])
	{
	  temp=xreal[j];
	  xreal[j]=xreal[j+igap];
	  xreal[j+igap]=temp;
	  
	  temp=ximag[j];
	  ximag[j]=ximag[j+igap];
	  ximag[j+igap]=temp;
	  
	  if(apply)
	  {
	    temp=y[j];
	    y[j]=y[j+igap];
	    y[j+igap]=temp;
	  }
	}
	else
	  continue; //goto label120;
        j=j-igap;
        goto label110;
      }
label120:
      igap=igap/2;
      goto label100;
    }
    else if(which=="LI")
    {
      /*
c
c        %------------------------------------------------%
c        | Sort XIMAG into increasing order of magnitude. |
c        %------------------------------------------------%
c
      */
label130:
      if(igap==0) return;
      for(i=igap;i<=(n-1);i++)
      {
        j=i-igap;
label140:
        
        if(j<0) continue; //goto label150;
        
        if(fabs(ximag[j])>fabs(ximag[j+igap]))
        {
          temp=xreal[j];
          xreal[j]=xreal[j+igap];
          xreal[j+igap]=temp;
          
          temp=ximag[j];
          ximag[j]=ximag[j+igap];
          ximag[j+igap]=temp;
          
          if(apply)
          {
            temp=y[j];
            y[j]=y[j+igap];
            y[j+igap]=temp;
          }
        }
        else
          continue; //goto label150;
        j=j-igap;
        goto label140;
      }
label150:
      igap=igap/2;
      goto label130;
    }
    else if(which=="SI")
    {
      /*
c
c        %------------------------------------------------%
c        | Sort XIMAG into decreasing order of magnitude. |
c        %------------------------------------------------%
c
      */
label160:
      if(igap==0) return;
      for(i=igap;i<=(n-1);i++)
      {
        j=i-igap;
label170:
        
        if(j<0) continue; //goto label180;
        
        if(fabs(ximag[j])<fabs(ximag[j+igap]))
        {
          temp=xreal[j];
          xreal[j]=xreal[j+igap];
          xreal[j+igap]=temp;
          
          temp=ximag[j];
          ximag[j]=ximag[j+igap];
          ximag[j+igap]=temp;
          
          if(apply)
          {
            temp=y[j];
            y[j]=y[j+igap];
            y[j+igap]=temp;
          }
        }
        else
          continue; //goto label180;
        j=j-igap;
        goto label170;
      }
      
label180:
      igap=igap/2;
      goto label160;
    }
    
    return;
    
    /*
c
c     %---------------%
c     | End of dsortc |
c     %---------------%
c
    */
}

/* CLN version */

template<>
void dsortc<cln::cl_F>(const std::string which, bool apply, int n, cln::cl_F* xreal, cln::cl_F* ximag, cln::cl_F* y)
{
  /*
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
  */
  
  int i, igap, j;
  cln::cl_F temp, temp1, temp2;
  
  /*
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
  */
  
  igap=n/2;
  
 // std::cout << "which = " << which << "\n";
  
  if(which=="LM")
  {
    /*
c
c        %------------------------------------------------------%
c        | Sort XREAL,XIMAG into increasing order of magnitude. |
c        %------------------------------------------------------%
c
    */
label10:
    if(igap==0) return;
    
    for(i=igap;i<=(n-1);i++)
    {
      j=i-igap;
label20:
      
      if(j<0) continue;
      
      temp1=dlapy2<cln::cl_F>(xreal[j],ximag[j]);
      temp2=dlapy2<cln::cl_F>(xreal[j+igap],ximag[j+igap]);
      
      if(temp1>temp2)
      {
        temp=xreal[j];
        xreal[j]=xreal[j+igap];
        xreal[j+igap]=temp;
        
        temp=ximag[j];
        ximag[j]=ximag[j+igap];
        ximag[j+igap]=temp;
        
        if(apply)
        {
          temp=y[j];
          y[j]=y[j+igap];
          y[j+igap]=temp;
        }
      }
      else
        continue; //goto label30;
      j=j-igap;
      goto label20;
    }
label30:
      igap=igap/2;
      goto label10;
  }
    else if(which=="SM")
    {
      /*
c
c        %------------------------------------------------------%
c        | Sort XREAL,XIMAG into decreasing order of magnitude. |
c        %------------------------------------------------------%
c
      */
      
label40:
      if(igap==0) return;
      
      for(i=igap;i<=(n-1);i++)
      {
        j=i-igap;
label50:
        
        if(j<0) continue; //goto label60;
        
        temp1=dlapy2<cln::cl_F>(xreal[j],ximag[j]);
        temp2=dlapy2<cln::cl_F>(xreal[j+igap],ximag[j+igap]);
        
        if(temp1<temp2)
        {
          temp=xreal[j];
          xreal[j]=xreal[j+igap];
          xreal[j+igap]=temp;
          
          temp=ximag[j];
          ximag[j]=ximag[j+igap];
          ximag[j+igap]=temp;
          
          if(apply)
          {
            temp=y[j];
            y[j]=y[j+igap];
            y[j+igap]=temp;
          }
        }
        else
          continue; //goto label60;
        j-=igap;
        goto label50;
      }
label60:
      igap=igap/2;
      goto label40;
    }
    else if(which=="LR")
    {
      /*
c
c        %------------------------------------------------%
c        | Sort XREAL into increasing order of algebraic. |
c        %------------------------------------------------%
c
      */
label70:
  //    ivout<T>(debug.logfil,1,&igap,debug.ndigit,"sortc_: igap=");
  //    dvout<T>(debug.logfil,n,xreal,debug.ndigit,"sortc_: xreal=");
  //    dvout<T>(debug.logfil,n,ximag,debug.ndigit,"sortc_: ximag=");
      if(igap==0) return;
        
      for(i=igap;i<=(n-1);i++)
      {
        j=i-igap;
    //    ivout<T>(debug.logfil,1,&i,debug.ndigit,"sortc_: i=");
label80:
    //    ivout<T>(debug.logfil,1,&j,debug.ndigit,"sortc_: j=");
        if(j<0) continue; //goto label90;
          
        if(xreal[j]>xreal[j+igap])
        {
          temp=xreal[j];
          xreal[j]=xreal[j+igap];
          xreal[j+igap]=temp;
          
          temp=ximag[j];
          ximag[j]=ximag[j+igap];
          ximag[j+igap]=temp;
            
          if(apply)
          {
            temp=y[j];
            y[j]=y[j+igap];
            y[j+igap]=temp;
          }
        }
        else
          continue; //goto label90;
        j=j-igap;
        goto label80;
      }

label90:
      igap=igap/2;
      goto label70;
    }
    else if(which=="SR")
    {
      /*
c
c        %------------------------------------------------%
c        | Sort XREAL into decreasing order of algebraic. |
c        %------------------------------------------------%
c
      */
label100:
      if(igap==0) return;
      for(i=igap;i<=(n-1);i++)
      {
        j=i-igap;
label110:

        if(j<0) continue; //goto label120;
        
        if(xreal[j]<xreal[j+igap])
        {
          temp=xreal[j];
          xreal[j]=xreal[j+igap];
          xreal[j+igap]=temp;
          
          temp=ximag[j];
          ximag[j]=ximag[j+igap];
          ximag[j+igap]=temp;
          
          if(apply)
          {
            temp=y[j];
            y[j]=y[j+igap];
            y[j+igap]=temp;
          }
        }
        else
          continue; //goto label120;
        j=j-igap;
        goto label110;
      }
label120:
      igap=igap/2;
      goto label100;
    }
    else if(which=="LI")
    {
      /*
c
c        %------------------------------------------------%
c        | Sort XIMAG into increasing order of magnitude. |
c        %------------------------------------------------%
c
      */
label130:
      if(igap==0) return;
      for(i=igap;i<=(n-1);i++)
      {
        j=i-igap;
label140:
        
        if(j<0) continue; //goto label150;
        
        if(cln::abs(ximag[j])>cln::abs(ximag[j+igap]))
        {
          temp=xreal[j];
          xreal[j]=xreal[j+igap];
          xreal[j+igap]=temp;
          
          temp=ximag[j];
          ximag[j]=ximag[j+igap];
          ximag[j+igap]=temp;
          
          if(apply)
          {
            temp=y[j];
            y[j]=y[j+igap];
            y[j+igap]=temp;
          }
        }
        else
          continue; //goto label150;
        j=j-igap;
        goto label140;
      }
label150:
      igap=igap/2;
      goto label130;
    }
    else if(which=="SI")
    {
      /*
c
c        %------------------------------------------------%
c        | Sort XIMAG into decreasing order of magnitude. |
c        %------------------------------------------------%
c
      */
label160:
      if(igap==0) return;
      for(i=igap;i<=(n-1);i++)
      {
        j=i-igap;
label170:
        
        if(j<0) continue; //goto label180;
        
        if(cln::abs(ximag[j])<cln::abs(ximag[j+igap]))
        {
          temp=xreal[j];
          xreal[j]=xreal[j+igap];
          xreal[j+igap]=temp;
          
          temp=ximag[j];
          ximag[j]=ximag[j+igap];
          ximag[j+igap]=temp;
          
          if(apply)
          {
            temp=y[j];
            y[j]=y[j+igap];
            y[j+igap]=temp;
          }
        }
        else
          continue; //goto label180;
        j=j-igap;
        goto label170;
      }
      
label180:
      igap=igap/2;
      goto label160;
    }
    
    return;
    
    /*
c
c     %---------------%
c     | End of dsortc |
c     %---------------%
c
    */
}


#endif
          
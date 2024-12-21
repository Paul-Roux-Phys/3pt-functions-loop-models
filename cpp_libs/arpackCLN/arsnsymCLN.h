#ifndef ARSNSYMCLN_H
#define ARSNSYMCLN_H

#include <cstddef>
#include <string>
#include "arch.h"
#include "arseig.h"
#include "arrseig.h"
#include "arrsnsym.h"
#include "cnaupp.h"

template<class ARFLOAT, class ARFOP>
class ARNonSymStdEigCLN:
  public virtual ARStdEig<ARFLOAT, ARFLOAT, ARFOP>,
  public virtual ARrcNonSymStdEig<ARFLOAT> //,
//  protected virtual ARrcStdEig<ARFLOAT, cln::cl_F>
{
  
  public:
    
    int digits;
    
    ARNonSymStdEigCLN() { }
    // Short constructor.    

    ARNonSymStdEigCLN(int np, int nevp, ARFOP* objOPp,
                 void (ARFOP::* MultOPxp)(ARFLOAT[], ARFLOAT[]), int digits,
                 const std::string& whichp = "LM", int ncvp = 0, ARFLOAT tolp = 0.0,
                 int maxitp = 0, ARFLOAT* residp = NULL, bool ishiftp = true);
    // Long constructor for CLN (regular mode).

    ARNonSymStdEigCLN(int np, int nevp, ARFOP* objOPp,
                 void (ARFOP::* MultOPxp)(ARFLOAT[], ARFLOAT[]),
                 ARFLOAT sigma, int digits, const std::string& whichp = "LM", int ncvp = 0,
                 ARFLOAT tolp = 0.0, int maxitp = 0, ARFLOAT* residp = NULL,
                 bool ishiftp = true);
    // Long constructor for CLN (shift and invert mode).
    
    ARNonSymStdEigCLN(const ARNonSymStdEigCLN& other) { Copy(other); }
    // Copy constructor.

    virtual ~ARNonSymStdEigCLN() { }
    // Destructor.
    
  protected:
    
    void Aupp();
    // Find Arnoldi basis
    
//    int FindEigenvalues();
    
//    int FindEigenvectors(bool);
    
//    int TakeStep();
    
//    int FindArnoldiBasis();
    
//    int FindSchurVectors();
    
//    int Eigenvectors(ARFOP*&, bool);
    
//    ARFLOAT ArnoldiBasisVector(int,int);
    
//    ARFLOAT SchurVector(int,int);
    
//    ARFLOAT ResidualVector(int);
    
//    ARFOP* RawArnoldiBasisVectors();
    
//    ARFOP* RawArnoldiBasisVector(int);
    
//    void Eupp();
  
  
}; // class ARNonSymStdEigCLN

// ------------------------------------------------------------------------ //
// ARNonSymStdEig member functions definition.                              //
// ------------------------------------------------------------------------ //


template<class ARFLOAT, class ARFOP>
inline ARNonSymStdEigCLN<ARFLOAT, ARFOP>::
ARNonSymStdEigCLN(int np, int nevp, ARFOP* objOPp,
               void (ARFOP::* MultOPxp)(ARFLOAT[], ARFLOAT[]), int digitsp,
               const std::string& whichp, int ncvp, ARFLOAT tolp, int maxitp,
               ARFLOAT* residp, bool ishiftp)

{
  digits=digitsp;

  this->NoShift();
  this->DefineParameters(np, nevp, objOPp, MultOPxp, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT, class ARFOP>
inline ARNonSymStdEigCLN<ARFLOAT, ARFOP>::
ARNonSymStdEigCLN(int np, int nevp, ARFOP* objOPp,
               void (ARFOP::* MultOPxp)(ARFLOAT[], ARFLOAT[]),
               ARFLOAT sigmap, int digitsp, const std::string& whichp, int ncvp, ARFLOAT tolp,
               int maxitp, ARFLOAT* residp, bool ishiftp)

{
  digits=digitsp;

  ChangeShift(sigmap);
  DefineParameters(np, nevp, objOPp, MultOPxp, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (shift and invert mode).

/* Template specialization for ARrcNonSymStdEig */
template<>
inline void ARrcNonSymStdEig<cln::cl_F>::Aupp()
{

} // Aupp template specialization. Needed for compilation


template<class ARFLOAT, class ARFOP>
inline void ARNonSymStdEigCLN<ARFLOAT, ARFOP>::Aupp()
{
 // std::cout << " Calling Aupp \n";

  naupp<ARFLOAT>(this->ido,this->bmat, this->n, this->which, this->nev, this->tol, this->resid, this->ncv, this->V, this->n,
        this->iparam, this->ipntr, this->workd, this->workl, this->lworkl, this->info, this->digits);

} // Aupp.

/*template<class ARFLOAT, class ARFOP>
int ARNonSymStdEigCLN<ARFLOAT, ARFOP>::FindEigenvalues()
{

  // Determining eigenvalues if they are not available.

  if (!this->ValuesOK) {
    try {
  //    ValAllocate();
 //     nconv = FindArnoldiBasis();
      this->rvec  = false;
 //     if (nconv>0) {
 //       Eupp();
 //       EuppError();
 //     }
    }
    catch (ArpackError) {
      ArpackError(ArpackError::CANNOT_FIND_VALUES, "FindEigenvalues");
      return 0;
    }
    if (this->newVal) this->ValuesOK = true;
    // next line is a temporary print command for eigenvalues 
    for(int i=0;i<this->nconv;i++)
      std::cout << "Eigenvalue " << i << ":" << this->workl[this->ipntr[6]+i] << "\n";
  }
  return this->nconv;

} // FindEigenvalues.

template<class ARFLOAT, class ARFOP>
int ARNonSymStdEigCLN<ARFLOAT, ARFOP>::TakeStep()
{

  // Requiring the definition of all internal variables.

  if (!this->PrepareOK) {

    throw ArpackError(ArpackError::PREPARE_NOT_OK, "TakeStep");

  }
  else if (!this->BasisOK) {

    // Taking a step if the Arnoldi basis is not available.

    Aupp();

    // Checking if convergence was obtained.

    if (this->ido==99) {
      this->nconv = this->iparam[5];
      this->AuppError();
      if (this->info >= 0) this->BasisOK = true;
    }
  }

  return this->ido;

} // TakeStep.


template<class ARFLOAT, class ARFOP>
inline int ARNonSymStdEigCLN<ARFLOAT, ARFOP>::FindArnoldiBasis()
{

  if (!this->BasisOK) {
    throw ArpackError(ArpackError::CANNOT_FIND_BASIS, "FindArnoldiBasis");
  }
  return this->nconv;

} // FindArnoldiBasis.

template<class ARFLOAT, class ARFOP>
int ARNonSymStdEigCLN<ARFLOAT, ARFOP>::FindEigenvectors(bool schurp)
{

  // Determining eigenvectors if they are not available.

  if (!this->VectorsOK) {
    try {
      this->ValAllocate();
      this->VecAllocate(schurp);
      this->nconv = FindArnoldiBasis();
  //    std::cout << "Successfully found Arnoldi basis, nconv = " << nconv << "\n";
  //    exit(1);
      this->rvec  = true;
      this->HowMny = 'A';
      if (this->nconv>0) {
        this->Eupp();
        this->EuppError();
      }
    }
    catch (ArpackError) {
      ArpackError(ArpackError::CANNOT_FIND_VECTORS, "FindEigenvectors");
      return 0;
    }
    this->BasisOK = false;
    if (this->newVal) this->ValuesOK = true;
    if (this->newVec || this->OverV()) this->VectorsOK = true;
    if (!this->OverV()) this->SchurOK = true;
  }
  //std::cout << "Successfully called Eupp() \n";
  //exit(1);
  return this->nconv;

} // FindEigenvectors.


template<class ARFLOAT, class ARFOP>
int ARNonSymStdEigCLN<ARFLOAT, ARFOP>::FindSchurVectors()
{

  // Determining Schur vectors if they are not available.

  if (!this->SchurOK) {
    try {
      this->ValAllocate();
      this->nconv  = this->FindArnoldiBasis();
      this->rvec   = true;
      this->HowMny = 'P';
//      if (nconv>0) {
//        Eupp();
//        EuppError();
//      }
    }
    catch (ArpackError) {
      ArpackError(ArpackError::CANNOT_FIND_SCHUR, "FindSchurVectors");
      return 0;
    }
    this->BasisOK = false;
    if (this->newVal) this->ValuesOK = true;
    this->SchurOK =true;
  }
  return this->nconv;

} // FindSchurVectors.


template<class ARFLOAT, class ARFOP>
int ARNonSymStdEigCLN<ARFLOAT, ARFOP>::
Eigenvectors(ARFOP* &EigVecp, bool ischur)
{

  // Overriding EigVecp with the converged eigenvectors.

  if (this->VectorsOK) {                       // Eigenvectors are available.
    if ((EigVecp == NULL) && (this->newVec)) { // Moving eigenvectors.
      EigVecp   = this->EigVec;
      this->EigVec    = NULL;
      this->newVec    = false;
      this->VectorsOK = false;
    }
    else {                               // Copying eigenvectors.
      if (EigVecp == NULL) {
        try { EigVecp = new ARFLOAT[this->ValSize()*this->n]; }
        catch (ArpackError) { return 0; }
      }
      copy(this->ValSize()*this->n,this->EigVec,1,EigVecp,1);
    }
  }
  else {                                // Eigenvectors are not available.
    if (this->newVec) {
      delete[] this->EigVec;
      this->newVec = false;
    }
    if (EigVecp == NULL) {
      try { EigVecp = new ARFLOAT[this->ValSize()*this->n]; }
      catch (ArpackError) { return 0; }
    }
    this->EigVec = EigVecp;
    this->nconv  = this->FindEigenvectors(this->ischur);
    this->EigVec = NULL;
  }
  return this->nconv;

} // Eigenvectors(EigVecp, ischur).


template<class ARFLOAT, class ARFOP>
inline ARFLOAT ARNonSymStdEigCLN<ARFLOAT, ARFOP>::ArnoldiBasisVector(int i, int j)
{

  // Returning element j of Arnoldi basis vector i.

  if (!this->BasisOK) {
    throw ArpackError(ArpackError::BASIS_NOT_OK, "ArnoldiBasisVector(i,j)");
  }
  else if ((i>=this->ncv)||(i<0)||(j>=this->n)||(j<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR,"ArnoldiBasisVector(i,j)");
  }
  return this->V[i*this->n+j+1];

} // ArnoldiBasisVector(i,j).


template<class ARFLOAT, class ARFOP>
inline ARFLOAT ARNonSymStdEigCLN<ARFLOAT, ARFOP>::SchurVector(int i, int j)
{

  // Returning element j of Schur vector i.

  if (!this->SchurOK) {
    throw ArpackError(ArpackError::SCHUR_NOT_OK, "SchurVector(i,j)");
  }
  else if ((i>=this->nconv)||(i<0)||(j>=this->n)||(j<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "SchurVector(i,j)");
  }
  return this->V[i*this->n+j+1];

} // SchurVector(i,j).


template<class ARFLOAT, class ARFOP>
inline ARFLOAT ARNonSymStdEigCLN<ARFLOAT, ARFOP>::ResidualVector(int i)
{

  // Returning element i of the residual vector.

  if ((!this->newRes)||(!(this->BasisOK||this->ValuesOK||this->VectorsOK||this->SchurOK))) {
    throw ArpackError(ArpackError::RESID_NOT_OK, "ResidualVector(i)");
  }
  else if ((i>=this->n)||(i<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "ResidualVector(i)");
  }
  return this->resid[i];

} // ResidualVector(i).


template<class ARFLOAT, class ARFOP>
inline ARFOP* ARNonSymStdEigCLN<ARFLOAT, ARFOP>::RawArnoldiBasisVectors()
{

  // Returning a constant pointer to Arnoldi basis.

  if (!this->BasisOK) {
    throw ArpackError(ArpackError::BASIS_NOT_OK, "RawArnoldiBasisVectors");
  }
  return &this->V[1];

} // RawArnoldiBasisVectors.


template<class ARFLOAT, class ARFOP>
inline ARFOP* ARNonSymStdEigCLN<ARFLOAT, ARFOP>::RawArnoldiBasisVector(int i)
{

  // Returning a constant pointer to Arnoldi basis vector i.

  if (!this->BasisOK) {
    throw ArpackError(ArpackError::BASIS_NOT_OK, "RawArnoldiBasisVector(i)");
  }
  else if ((i>=this->ncv)||(i<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR,"RawArnoldiBasisVector(i)");
  }
  return &this->V[i*this->n+1];

} // RawArnoldiBasisVector(i).*/

#endif
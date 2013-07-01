/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROO_MODIFTSALLIS
#define ROO_MODIFTSALLIS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;
class RooAbsReal;
 
class RooModifTsallis : public RooAbsPdf {
public:
  RooModifTsallis(const char *name, const char *title,
	          RooAbsReal& _x,
        	  RooAbsReal& _m,
              	  RooAbsReal& _n,
	          RooAbsReal& _n2,
                  RooAbsReal& _bb,
	          RooAbsReal& _bb2,
	          RooAbsReal& _T,
	          RooAbsReal& _fexp);

  RooModifTsallis(const RooModifTsallis& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooModifTsallis(*this,newname); }
  inline virtual ~RooModifTsallis() { }
  /* Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
     Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;*/

protected:

  RooRealProxy x ;
  RooRealProxy m ;
  RooRealProxy n ;
  RooRealProxy n2 ;
  RooRealProxy bb ;
  RooRealProxy bb2 ;
  RooRealProxy T ;
  RooRealProxy fexp ;

  Double_t evaluate() const ;

private:

 ClassDef(RooModifTsallis,1) // Your description goes here...
};
 
#endif

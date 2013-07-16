/* 
 * Retrieve signal shapes parameters from  and write them in the fragments 
 * usage: 
 * -set all input variables in Config.h
 * -run with:
 * root -q -b californiaSignalEBE.C
 * This runs on the 3 final states for 7 and 8 TeV and writes the output in a file (see stdout).
 * Use other scripts (compareSignalFits.C, signalFits.C) or ask experts to check for the shapes
 */


#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

void californiaSignalEBE(){
  gROOT->LoadMacro("$CMSSW_BASE/src/WWAnalysis/TreeModifiers/macro/SignalInterpolationStrings.h+");
  californiaSignalShapes(1,7);
  californiaSignalShapes(2,7);
  californiaSignalShapes(3,7);
  californiaSignalShapes(1,8);
  californiaSignalShapes(2,8);
  californiaSignalShapes(3,8);
}

void californiaSignalShapes(int channel, int sqrts){

  bool en=1;
  if(sqrts==8)en=0;
 
  string schannel;
  if (channel == 1) schannel = "4mu";
  if (channel == 2) schannel = "4e";
  if (channel == 3) schannel = "2e2mu";
  cout << "Final state = " << schannel << " and sqrt(s) = " << sqrts << endl;


  char tmp_outCardName[200];
  sprintf(tmp_outCardName,"%iTeV_",sqrts);
  string prependName = "CardFragments/signalEBE_";
  string appendName = ".txt";
  string outCardName =  prependName + tmp_outCardName + schannel + appendName;
  
  ofstream ofsCard;
  ofsCard.open(outCardName.c_str(),fstream::out);

  ofsCard << "## signal ebe functions --- no spaces! ##" << endl;
  TString string = TString(getSignalEBELandauMeanString(channel-1,en) );
  string.ReplaceAll(" ","");
  ofsCard << "RelErrShape relerr_ggH_ld_mean "      << string.Data()<< endl;	 
  string = TString(getSignalEBELandauSigmaString(channel-1,en) );
  string.ReplaceAll(" ","");
  ofsCard << "RelErrShape relerr_ggH_ld_sigma "      << string.Data()<< endl;	 
  string = TString(getSignalEBELandauFracString(channel-1,en) );
  string.ReplaceAll(" ","");
  ofsCard << "RelErrShape relerr_ggH_ld_frac "      << string.Data()<< endl;	 
  string = TString(getSignalEBEGaussianMeanString(channel-1,en) );
  string.ReplaceAll(" ","");
  ofsCard << "RelErrShape relerr_ggH_gs_mean "      << string.Data()<< endl;	 
  string = TString(getSignalEBEGaussianSigmaString(channel-1,en) );
  string.ReplaceAll(" ","");
  ofsCard << "RelErrShape relerr_ggH_gs_sigma "      << string.Data()<< endl;	 

  return;
}

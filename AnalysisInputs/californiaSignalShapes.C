/* 
 * Retrieve signal shapes parameters from  and write them in the fragments 
 * usage: 
 * -set all input variables in Config.h
 * -run with:
 * root -q -b californiaSignalShapes.C
 * This runs on the 3 final states for 7 and 8 TeV and writes the output in a file (see stdout).
 * Use other scripts (compareSignalFits.C, signalFits.C) or ask experts to check for the shapes
 */


#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

void californiaSignalShapes(){
  gROOT->LoadMacro("$CMSSW_BASE/src/Higgs/Higgs_CS_and_Width/include/SignalInterpolationStrings.h+");
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
  string prependName = "CardFragments/signalFunctions_";
  string appendName = ".txt";
  string outCardName =  prependName + tmp_outCardName + schannel + appendName;
  
  ofstream ofsCard;
  ofsCard.open(outCardName.c_str(),fstream::out);

  ofsCard << "## signal functions --- no spaces! ##" << endl;
  ofsCard << "usehighmassreweightedshapes" << endl;
  ofsCard << "signalShape n_CB "      << getSignalCBNLString(125.,channel-1,en) << endl;	     
  ofsCard << "signalShape alpha_CB "  << getSignalCBAlphaLString(125.,channel-1,en)  << endl; 
  ofsCard << "signalShape n2_CB "     << getSignalCBNRString(125.,channel-1,en)<< endl;	     
  ofsCard << "signalShape alpha2_CB " << getSignalCBAlphaRString(125.,channel-1,en) << endl;  
  ofsCard << "signalShape mean_CB "   << getSignalCBMeanString(125.,channel-1,en,1) << endl;  
  ofsCard << "signalShape sigma_CB "  << getSignalCBSigmaString(125.,channel-1,en)  << endl;  

  ofsCard << "HighMasssignalShape n_CB "      << getSignalCBNLString(500.,channel-1,en) << endl;	     
  ofsCard << "HighMasssignalShape alpha_CB "  << getSignalCBAlphaLString(500.,channel-1,en)  << endl; 
  ofsCard << "HighMasssignalShape n2_CB "     << getSignalCBNRString(500.,channel-1,en)<< endl;	     
  ofsCard << "HighMasssignalShape alpha2_CB " << getSignalCBAlphaRString(500.,channel-1,en) << endl;  
  ofsCard << "HighMasssignalShape mean_CB "   << getSignalCBMeanString(500.,channel-1,en,1) << endl;  
  ofsCard << "HighMasssignalShape sigma_CB "  << getSignalCBSigmaString(500.,channel-1,en)  << endl;  
  ofsCard << "HighMasssignalShape gamma_BW "  << getSignalBWGammaString(500.,channel-1,en)  << endl;  
  ofsCard << endl;

  return;
}

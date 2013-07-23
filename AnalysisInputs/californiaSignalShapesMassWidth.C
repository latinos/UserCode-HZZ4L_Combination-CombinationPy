/* 
 * Retrieve signal shapes parameters from  and write them in the fragments 
 * usage: 
 * -set all input variables in Config.h
 * -run with:
 * root -q -b californiaSignalShapesMassWidth.C
 * This runs on the 3 final states for 7 and 8 TeV and writes the output in a file (see stdout).
 * Use other scripts (compareSignalFits.C, signalFits.C) or ask experts to check for the shapes
 */


#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

void californiaSignalShapesMassWidth(){
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
  string prependName = "CardFragments/signalFunctionsMW_";
  string appendName = ".txt";
  string outCardName =  prependName + tmp_outCardName + schannel + appendName;
  
  ofstream ofsCard;
  ofsCard.open(outCardName.c_str(),fstream::out);

  ofsCard << "## signal functions --- no spaces! ##" << endl;
  ofsCard << "usehighmassreweightedshapes" << endl;
  TString string; string.Form("%i",getSignalCBNLValueInt(channel-1,en) );
  string.ReplaceAll(" ","");
  string.ReplaceAll("+@0*@1","");
  ofsCard << "signalShape n_CB "      << string.Data()<< endl;	 
  string = getSignalACBAlphaLString(channel-1,en,false);
  string.ReplaceAll(" ","");
  ofsCard << "signalShape alpha_CB "  << string.Data()   << endl; 
  string.Form("%i",getSignalCBNRValueInt(channel-1,en)   ) ;
  string.ReplaceAll(" ","");
  ofsCard << "signalShape n2_CB "     << string.Data()   << endl;	     
  string = getSignalACBAlphaRString(channel-1,en, false);
  string.ReplaceAll(" ","");
  ofsCard << "signalShape alpha2_CB " << string.Data()  << endl;  
  string = getSignalACBMeanString(channel-1,en,1);
  string.ReplaceAll(" ","");
  string.ReplaceAll("+@0*@1","");
  ofsCard << "signalShape mean_CB "   << string.Data()  << endl;  
  string = getSignalACBSigmaString(channel-1,en) ;
  string.ReplaceAll(" ","");
  string.ReplaceAll("*(1+@1)","");
  ofsCard << "signalShape sigma_CB "  << string.Data()  << endl;  

  string = getSignalCBNLString(500.,channel-1,en) ;
  string.ReplaceAll(" ","");
  string.ReplaceAll("+@0*@1","");
  ofsCard << "HighMasssignalShape n_CB "      << string.Data() << endl;	     
  string = getSignalCBAlphaLString(500.,channel-1,en,false);
  string.ReplaceAll(" ","");
  ofsCard << "HighMasssignalShape alpha_CB "  << string.Data()  << endl; 
  string = getSignalCBNRString(500.,channel-1,en)    ;
  string.ReplaceAll(" ","");
  ofsCard << "HighMasssignalShape n2_CB "     << string.Data()  << endl;
  string = getSignalCBAlphaRString(500.,channel-1,en,false);
  string.ReplaceAll(" ","");
  ofsCard << "HighMasssignalShape alpha2_CB " << string.Data()  << endl;  
  string = getSignalCBMeanString(500.,channel-1,en,1);
  string.ReplaceAll(" ","");
  ofsCard << "HighMasssignalShape mean_CB "   << string.Data()  << endl;  
  string = getSignalCBSigmaString(500.,channel-1,en) ;
  string.ReplaceAll(" ","");
  ofsCard << "HighMasssignalShape sigma_CB "  << string.Data()  << endl;  
  string = getSignalBWGammaString(500.,channel-1,en) ;
  string.ReplaceAll(" ","");
  string.ReplaceAll("*(1+@1)","");
  ofsCard << "HighMasssignalShape gamma_BW "  << string.Data()  << endl;  
  ofsCard << endl;

  return;
}

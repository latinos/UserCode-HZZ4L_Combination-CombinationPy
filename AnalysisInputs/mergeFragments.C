/* 
 * Merge card fragments
 * usage: 
 * -set all input variables in Config.h
 * -run with:
 * root -q -b mergeFragments.C+
 * This runs on the 3 final states for 7 and 8 TeV and writes the full config files in 
 * ../CreateDatacards/SM_inputs_?TeV/ .
 *
 */


#include <TString.h>
#include <TSystem.h>
#include <iostream>
#include <fstream>

//----------> SET INPUT VARIABLES in Config.h
#include "Config.h"
//<----------

using namespace std;
void mergeFragments(int channel, int sqrts, double lumi, bool VBFtag = false);
void append(TString file, TString outfile);


// Run all fs/sqrts in one go
void mergeFragments() {

  mergeFragments(1, 7, lumi7TeV,false);
  mergeFragments(2, 7, lumi7TeV,false);
  mergeFragments(3, 7, lumi7TeV,false);
  mergeFragments(1, 8, lumi8TeV,false);
  mergeFragments(2, 8, lumi8TeV,false);
  mergeFragments(3, 8, lumi8TeV,false);

  mergeFragments(1, 7, lumi7TeV,true);
  mergeFragments(2, 7, lumi7TeV,true);
  mergeFragments(3, 7, lumi7TeV,true);
  mergeFragments(1, 8, lumi8TeV,true);
  mergeFragments(2, 8, lumi8TeV,true);
  mergeFragments(3, 8, lumi8TeV,true);
  
}


void mergeFragments(int channel, int sqrts, double lumi, bool VBFtag) {

  TString schannel;
  if      (channel == 1) schannel = "4mu";
  else if (channel == 2) schannel = "4e";
  else if (channel == 3) schannel = "2e2mu";

  TString ssqrts = (long) sqrts + TString("TeV");

  TString outfile = "../CreateDatacards/SM_inputs_" + ssqrts + "_" + Form("%d",int(VBFtag)) + "/inputs_" +  schannel + "_" + Form("%d",int(VBFtag)) + ".txt";
  ofstream of(outfile,ios_base::out);

  float lumiUnc = 0;
  if      (sqrts==7) lumiUnc = 1.022;
  else if (sqrts==8) lumiUnc = 1.044;
    

  of << "############## Inputs for " << schannel << " for " << sqrts << " TeV " << " VBFtag-> " << VBFtag << "  ##############" << endl
     << "## SM ##"                 << endl
     << "model SM"                 << endl
     <<                               endl
     << "## decay chan ##"         << endl
     << "decay " << schannel       << endl
     <<                               endl
     << "## lumi ##"               << endl
     << "lumi " << lumi            << endl
     << "systematic lumiUnc " <<  lumiUnc << endl
     <<                               endl
     << "## sqrtS ##"              << endl
     << "sqrts " << sqrts          << endl
     <<                               endl
     << "## Channels to include in cards ##" << endl
     << "channels ggH qqH WH ZH ttH qqZZ ggZZ zjets" << endl
     <<                               endl;
  of.close();

  TString suffix_both = ssqrts + "_" + schannel + ".txt";

  TString suffix_split = ssqrts + "_" + schannel + "_" + Form("%d",int(VBFtag)) + ".txt";

  append("CardFragments/ZZRates_" + suffix_split, outfile);
  append("CardFragments/zjetRate_" + suffix_both, outfile);
  append("CardFragments/signalFunctions_" + suffix_split, outfile);
  append("CardFragments/signalEfficiency_" + suffix_split, outfile);
  append("CardFragments/qqzzBackgroundFit_" + suffix_split, outfile);
  append("CardFragments/ggzzBackgroundFit_" + suffix_split, outfile);
  append("CardFragments/zjetShape_" + suffix_both, outfile);  
  append("CardFragments/sys_" + suffix_both, outfile);
  append("CardFragments/hypTesting.txt", outfile);
  append("CardFragments/VBFTesting.txt", outfile);
  append("CardFragments/mekd_" + suffix_both, outfile);
  append("CardFragments/relerr_" + suffix_both, outfile);

  cout << "Wrote " << outfile << endl;
}


void append(TString file, TString outfile){
  //gSystem->Exec("touch " + file);
  gSystem->Exec("cat " + file + " >> " + outfile);
}

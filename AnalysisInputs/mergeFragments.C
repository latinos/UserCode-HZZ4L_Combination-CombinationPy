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
void mergeFragments(int channel, int sqrts, double lumi);
void append(TString file, TString outfile);


// Run all fs/sqrts in one go
void mergeFragments() {

  mergeFragments(1, 7, lumi7TeV);
  mergeFragments(2, 7, lumi7TeV);
  mergeFragments(3, 7, lumi7TeV);
  mergeFragments(1, 8, lumi8TeV);
  mergeFragments(2, 8, lumi8TeV);
  mergeFragments(3, 8, lumi8TeV);
  
}


void mergeFragments(int channel, int sqrts, double lumi) {

  TString schannel;
  if      (channel == 1) schannel = "4mu";
  else if (channel == 2) schannel = "4e";
  else if (channel == 3) schannel = "2e2mu";

  TString ssqrts = (long) sqrts + TString("TeV");

  TString outfile = "../CreateDatacards/SM_inputs_" + ssqrts + "/inputs_" +  schannel + ".txt";
  ofstream of(outfile,ios_base::out);

  float lumiUnc = 0;
  if      (sqrts==7) lumiUnc = 1.022;
  else if (sqrts==8) lumiUnc = 1.05;
    

  of << "############## Inputs for " << schannel << " for " << sqrts << " TeV ##############" << endl
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

  TString suffix = ssqrts + "_" + schannel + ".txt";

  append("CardFragments/ZZRates_" + suffix, outfile);
  append("CardFragments/zjetRate_" + suffix, outfile);
  append("CardFragments/signalFunctions_" + suffix, outfile);
  append("CardFragments/signalEfficiency_" + suffix, outfile);
  append("CardFragments/qqzzBackgroundFit_" + suffix, outfile);
  append("CardFragments/ggzzBackgroundFit_" + suffix, outfile);
  append("CardFragments/zjetShape_" + suffix, outfile);  
  append("CardFragments/sys_" + suffix, outfile);

  cout << "Wrote " << outfile << endl;
}


void append(TString file, TString outfile){
  //gSystem->Exec("touch " + file);
  gSystem->Exec("cat " + file + " >> " + outfile);
}

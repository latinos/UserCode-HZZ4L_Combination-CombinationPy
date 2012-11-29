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
void mergeFragments(int channel, int sqrts, double lumi, bool VBFtag = false, int process = 1);
void append(TString file, TString outfile);


// Run all fs/sqrts in one go
void mergeFragments() {

  //ggH
  mergeFragments(1, 7, lumi7TeV,false,1);
  mergeFragments(2, 7, lumi7TeV,false,1);
  mergeFragments(3, 7, lumi7TeV,false,1);
  mergeFragments(1, 8, lumi8TeV,false,1);
  mergeFragments(2, 8, lumi8TeV,false,1);
  mergeFragments(3, 8, lumi8TeV,false,1);
  //qqH
  mergeFragments(1, 7, lumi7TeV,false,2);
  mergeFragments(2, 7, lumi7TeV,false,2);
  mergeFragments(3, 7, lumi7TeV,false,2);
  mergeFragments(1, 8, lumi8TeV,false,2);
  mergeFragments(2, 8, lumi8TeV,false,2);
  mergeFragments(3, 8, lumi8TeV,false,2);
  //ZH
  mergeFragments(1, 7, lumi7TeV,false,3);
  mergeFragments(2, 7, lumi7TeV,false,3);
  mergeFragments(3, 7, lumi7TeV,false,3);
  mergeFragments(1, 8, lumi8TeV,false,3);
  mergeFragments(2, 8, lumi8TeV,false,3);
  mergeFragments(3, 8, lumi8TeV,false,3);
  //WH
  mergeFragments(1, 7, lumi7TeV,false,4);
  mergeFragments(2, 7, lumi7TeV,false,4);
  mergeFragments(3, 7, lumi7TeV,false,4);
  mergeFragments(1, 8, lumi8TeV,false,4);
  mergeFragments(2, 8, lumi8TeV,false,4);
  mergeFragments(3, 8, lumi8TeV,false,4);
  //ttH
  mergeFragments(1, 7, lumi7TeV,false,5);
  mergeFragments(2, 7, lumi7TeV,false,5);
  mergeFragments(3, 7, lumi7TeV,false,5);
  mergeFragments(1, 8, lumi8TeV,false,5);
  mergeFragments(2, 8, lumi8TeV,false,5);
  mergeFragments(3, 8, lumi8TeV,false,5);

  //ggH
  mergeFragments(1, 7, lumi7TeV,true,1);
  mergeFragments(2, 7, lumi7TeV,true,1);
  mergeFragments(3, 7, lumi7TeV,true,1);
  mergeFragments(1, 8, lumi8TeV,true,1);
  mergeFragments(2, 8, lumi8TeV,true,1);
  mergeFragments(3, 8, lumi8TeV,true,1);
  //qqH
  mergeFragments(1, 7, lumi7TeV,true,2);
  mergeFragments(2, 7, lumi7TeV,true,2);
  mergeFragments(3, 7, lumi7TeV,true,2);
  mergeFragments(1, 8, lumi8TeV,true,2);
  mergeFragments(2, 8, lumi8TeV,true,2);
  mergeFragments(3, 8, lumi8TeV,true,2);
  //ZH
  mergeFragments(1, 7, lumi7TeV,true,3);
  mergeFragments(2, 7, lumi7TeV,true,3);
  mergeFragments(3, 7, lumi7TeV,true,3);
  mergeFragments(1, 8, lumi8TeV,true,3);
  mergeFragments(2, 8, lumi8TeV,true,3);
  mergeFragments(3, 8, lumi8TeV,true,3);
  //WH
  mergeFragments(1, 7, lumi7TeV,true,4);
  mergeFragments(2, 7, lumi7TeV,true,4);
  mergeFragments(3, 7, lumi7TeV,true,4);
  mergeFragments(1, 8, lumi8TeV,true,4);
  mergeFragments(2, 8, lumi8TeV,true,4);
  mergeFragments(3, 8, lumi8TeV,true,4);
  //ttH
  mergeFragments(1, 7, lumi7TeV,true,5);
  mergeFragments(2, 7, lumi7TeV,true,5);
  mergeFragments(3, 7, lumi7TeV,true,5);
  mergeFragments(1, 8, lumi8TeV,true,5);
  mergeFragments(2, 8, lumi8TeV,true,5);
  mergeFragments(3, 8, lumi8TeV,true,5);
  
}


void mergeFragments(int channel, int sqrts, double lumi, bool VBFtag, int process) {

  TString schannel;
  if      (channel == 1) schannel = "4mu";
  else if (channel == 2) schannel = "4e";
  else if (channel == 3) schannel = "2e2mu";

  TString sprocess;
  if      (process == 1) sprocess = "ggH";
  else if (process == 2) sprocess = "qqH";
  else if (process == 3) sprocess = "ZH";
  else if (process == 4) sprocess = "WH";
  else if (process == 5) sprocess = "ttH";

  TString ssqrts = (long) sqrts + TString("TeV");

  TString outfile;
  if (process == 1) outfile = "../CreateDatacards/SM_inputs_" + ssqrts + "_" + Form("%d",int(VBFtag)) + "/inputs_" +  schannel + "_" + Form("%d",int(VBFtag)) + ".txt";
  else if (process != 1) outfile = "../CreateDatacards/SM_inputs_" + ssqrts + "_" + Form("%d",int(VBFtag)) + "/" + sprocess + "_inputs_" +  schannel + "_" + Form("%d",int(VBFtag)) + ".txt";
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

  TString sigsuffix_both = sprocess + ssqrts + "_" + schannel + ".txt";
  TString sigsuffix_split = sprocess + ssqrts + "_" + schannel + "_" + Form("%d",int(VBFtag)) + ".txt";
  TString bkgsuffix_both = ssqrts + "_" + schannel + ".txt";
  TString bkgsuffix_split = ssqrts + "_" + schannel + "_" + Form("%d",int(VBFtag)) + ".txt";

  append("CardFragments/ZZRates_" + bkgsuffix_split, outfile);
  append("CardFragments/zjetRate_" + bkgsuffix_both, outfile);
  append("CardFragments/signalFunctions_" + sigsuffix_split, outfile);
  append("CardFragments/signalEfficiency_" + sigsuffix_split, outfile);
  append("CardFragments/qqzzBackgroundFit_" + bkgsuffix_split, outfile);
  append("CardFragments/ggzzBackgroundFit_" + bkgsuffix_split, outfile);
  append("CardFragments/zjetShape_" + bkgsuffix_both, outfile);  
  append("CardFragments/sys_" + bkgsuffix_both, outfile);
  append("CardFragments/hypTesting.txt", outfile);
  append("CardFragments/VBFtesting.txt",outfile);
  append("CardFragments/mekd_" + bkgsuffix_both, outfile);
  append("CardFragments/relerr_" + bkgsuffix_both, outfile);

  cout << "Wrote " << outfile << endl;
}


void append(TString file, TString outfile){
  //gSystem->Exec("touch " + file);
  gSystem->Exec("cat " + file + " >> " + outfile);
}

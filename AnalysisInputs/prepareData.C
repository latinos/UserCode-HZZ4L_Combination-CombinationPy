/* 
 * Prepare root files containing data events.
 * usage: 
 * -set all input variables in Config.h
 * -run with:
 * root -q -b prepareData.C+
 * This creates all files for 3 final states, 7 and 8 TeV and stores them in the final destination directory
 *
 */

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TAxis.h"
#include "TFile.h"
#include "TLegend.h"
#include "TChain.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TString.h"
#include <iostream>
#include <sstream>
#include <string>
#include "TTree.h"
#include "TText.h"
#include "TStyle.h"

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TROOT.h>
#include "ZZMatrixElement/MELA/interface/Mela.h"
#endif



//----------> SET INPUT VARIABLES in Config.h
#include "Config.h"
//<----------

//--- Flags to control re-computation of KD
bool recompute_ = false;       // Recompute KD (instead of taking it from the tree); when true, the following flags apply:
bool usePowhegTemplate=false;  // false use analytic bg
bool withPt_ = false;          // Include pT in KD
bool withY_  = false;          //    "    Y  "  "
int sqrts    = 8;              // sqrts, used only for withPt_/withY_


Mela* myMELA;

void convertTreeForDatacards(TString inFile, TString outfile);

// Run all final states and sqrts in one go
void prepareData() {

  if (recompute_) myMELA = new Mela(usePowhegTemplate); // this is safely leaked

  gSystem->Exec("mkdir -p "+ DataRootFilePath);
  convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleMu.root",  DataRootFilePath+"hzz4mu_"  +lumistr7TeV+".root");
  convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleEle.root", DataRootFilePath+"hzz4e_"   +lumistr7TeV+".root");
  convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleOr.root",  DataRootFilePath+"hzz2e2mu_"+lumistr7TeV+".root");
  convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleMu.root",  DataRootFilePath+"hzz4mu_"  +lumistr8TeV+".root");
  convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleEle.root", DataRootFilePath+"hzz4e_"   +lumistr8TeV+".root");
  convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleOr.root",  DataRootFilePath+"hzz2e2mu_"+lumistr8TeV+".root");
}


// The actual job
void convertTreeForDatacards(TString inFile, TString outfile){
  TChain* treedata ;
  treedata= new TChain("SelectedTree");
  treedata->Add(inFile);

  float mzz, pseudomela, mela;
  float m1, m2, costheta1, costheta2, costhetastar, phi, phi1;
  float pt4l, Y4l;

  treedata->SetBranchAddress("ZZMass",&mzz);
  treedata->SetBranchAddress("ZZpseudoLD",&pseudomela); 
  treedata->SetBranchAddress("ZZLD",&mela); 
  treedata->SetBranchAddress("Z1Mass",&m1);
  treedata->SetBranchAddress("Z2Mass",&m2);
  treedata->SetBranchAddress("helcosthetaZ1",&costheta1);
  treedata->SetBranchAddress("helcosthetaZ2",&costheta2);
  treedata->SetBranchAddress("costhetastar",&costhetastar);
  treedata->SetBranchAddress("helphi",&phi);
  treedata->SetBranchAddress("phistarZ1",&phi1);
  
  treedata->SetBranchAddress("ZZPt",&pt4l);
  treedata->SetBranchAddress("ZZRapidity",&Y4l);  

  TFile* newFile  = new TFile(outfile, "RECREATE");
  newFile->cd();
  TTree* newTree = new TTree("data_obs","data_obs"); 
  double CMS_zz4l_mass, melaLD, pseudomelaLD;
  newTree->Branch("CMS_zz4l_mass",&CMS_zz4l_mass,"CMS_zz4l_mass/D");
  newTree->Branch("pseudoMelaLD",&pseudomelaLD,"pseudoMelaLD/D");
  newTree->Branch("melaLD",&melaLD,"melaLD/D");
  cout << inFile << " entries: " << treedata->GetEntries() << endl;
  cout << "written: " << outfile << endl << endl;
  for(int iEvt=0; iEvt<treedata->GetEntries(); iEvt++){
    //    if(iEvt%5000==0) cout << "event: " << iEvt << endl;
    treedata->GetEntry(iEvt);

    CMS_zz4l_mass = mzz;
    pseudomelaLD = pseudomela;
    melaLD = mela;

    if(recompute_){
      float KD, psig, pbkg;
      myMELA->computeKD(mzz,m1,m2,
		       costhetastar,
		       costheta1,
		       costheta2,
		       phi,
		       phi1,
		       KD,psig,pbkg,
		       withPt_,pt4l,
		       withY_, Y4l,
		       sqrts);
      
      melaLD = KD;
    }


    newTree->Fill();
  }

  newTree->Write("data_obs"); 
  newFile->Close();
}



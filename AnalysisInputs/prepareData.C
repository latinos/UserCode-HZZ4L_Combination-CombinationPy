/* 
 * Prepare root files containing data events.
 * usage: 
 * -set all input variables in Config.h
 * -run with:
 * root -q -b prepareData.C+
 * This creates all files for 3 final states, 7 and 8 TeV and stores them in the final destination directory
 *
 */


//#define LINKMELA //Uncomment to link the MELA package to compute KD on the fly


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
#include <vector>
#include "TTree.h"
#include "TText.h"
#include "TStyle.h"

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TROOT.h>
#ifdef LINKMELA
#include "ZZMatrixElement/MELA/interface/Mela.h"
#endif
#endif

using namespace std;

//----------> SET INPUT VARIABLES in Config.h
#include "Config.h"
//<----------

//--- Flags to control re-computation of KD
bool recompute_ = false;       // Recompute KD (instead of taking it from the tree); when true, the following flags apply:
bool usePowhegTemplate=false;  // false use analytic bg
bool withPt_ = false;          // Include pT in KD
bool withY_  = false;          //    "    Y  "  "
int sqrts    = 8;              // sqrts, used only for withPt_/withY_

bool onlyICHEPStat = false;


#ifdef LINKMELA
Mela* myMELA;
#endif

void convertTreeForDatacards(TString inFile, TString outfile, bool useJET, bool VBFtag);

// Run all final states and sqrts in one go
void prepareData() {

#ifdef LINKMELA
  if (recompute_) myMELA = new Mela(usePowhegTemplate); // this is safely leaked
#endif

  gSystem->Exec("mkdir -p "+ DataRootFilePath);
  Int_t nOut[18];

  nOut[0]= convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleMu.root",  DataRootFilePath+"hzz4mu_"  +lumistr7TeV+"_1.root",true, true);
  nOut[1]= convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleEle.root", DataRootFilePath+"hzz4e_"   +lumistr7TeV+"_1.root",true, true);
  nOut[2]= convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleOr.root",  DataRootFilePath+"hzz2e2mu_"+lumistr7TeV+"_1.root",true, true);
  nOut[3]= convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleMu.root",  DataRootFilePath+"hzz4mu_"  +lumistr8TeV+"_1.root",true, true);
  nOut[4]= convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleEle.root", DataRootFilePath+"hzz4e_"   +lumistr8TeV+"_1.root",true, true);
  nOut[5]= convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleOr.root",  DataRootFilePath+"hzz2e2mu_"+lumistr8TeV+"_1.root",true, true);

  nOut[6]= convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleMu.root",  DataRootFilePath+"hzz4mu_"  +lumistr7TeV+"_0.root",true, false);
  nOut[7]= convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleEle.root", DataRootFilePath+"hzz4e_"   +lumistr7TeV+"_0.root",true, false);
  nOut[8]= convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleOr.root",  DataRootFilePath+"hzz2e2mu_"+lumistr7TeV+"_0.root",true, false);
  nOut[9]= convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleMu.root",  DataRootFilePath+"hzz4mu_"  +lumistr8TeV+"_0.root",true, false);
  nOut[10]=convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleEle.root", DataRootFilePath+"hzz4e_"   +lumistr8TeV+"_0.root",true, false);
  nOut[11]=convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleOr.root",  DataRootFilePath+"hzz2e2mu_"+lumistr8TeV+"_0.root",true, false);

  nOut[12]=convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleMu.root",  DataRootFilePath+"hzz4mu_"  +lumistr7TeV+".root",false, false);
  nOut[13]=convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleEle.root", DataRootFilePath+"hzz4e_"   +lumistr7TeV+".root",false, false);
  nOut[14]=convertTreeForDatacards(filePath7TeV + "/data/HZZ4lTree_DoubleOr.root",  DataRootFilePath+"hzz2e2mu_"+lumistr7TeV+".root",false, false);
  nOut[15]=convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleMu.root",  DataRootFilePath+"hzz4mu_"  +lumistr8TeV+".root",false, false);
  nOut[16]=convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleEle.root", DataRootFilePath+"hzz4e_"   +lumistr8TeV+".root",false, false);
  nOut[17]=convertTreeForDatacards(filePath8TeV + "/data/HZZ4lTree_DoubleOr.root",  DataRootFilePath+"hzz2e2mu_"+lumistr8TeV+".root",false, false);

  TString jetsT[3]={"DiJet     ",
		    "Untagged  ",
		    "Inclusive "};
  TString energies[2]={"7TeV", "8TeV "};
  TString channels[3]={"4mu   ",
		       "4e    ",
		       "2e2mu "};
  cout<<"Summary  "<<endl;
  for(int ijet=0;ijet<3;ijet++)
    for(int ien = 0;ien<2;ien++)
      for(int ich=0;ich<3;ich++)
	cout<<"N Events out for "<<jetsT[ijet]<<energies[ien].Data()<<channels[ich].Data()<<": "<<nOut[ijet*6+ien*3+ich]<<endl;
}


// The actual job
int convertTreeForDatacards(TString inFile, TString outfile, bool useJET, bool VBFtag){
  TChain* treedata ;
  treedata= new TChain("SelectedTree");
  treedata->Add(inFile);

  int neventOut=0;
  Int_t run;
  float mzz, pseudomela, mela, mzzErr;
  float p0plus_VAJHU, bkg_VAMCFM;
  float m1, m2, costheta1, costheta2, costhetastar, phi, phi1;
  //float ZZVAKD;
  float pt4l, Y4l, fisher;
  int NJets_30;

  treedata->SetBranchAddress("RunNumber",&run);
  treedata->SetBranchAddress("ZZMass",&mzz);
  treedata->SetBranchAddress("ZZMassErrCorr",&mzzErr);
  treedata->SetBranchAddress("p0plus_VAJHU",&p0plus_VAJHU);
  treedata->SetBranchAddress("bkg_VAMCFM",&bkg_VAMCFM);
  //treedata->SetBranchAddress("ZZVAKD",&ZZVAKD);
  treedata->SetBranchAddress("Z1Mass",&m1);
  treedata->SetBranchAddress("Z2Mass",&m2);
  //treedata->SetBranchAddress("helphi",&phi);
  //treedata->SetBranchAddress("phistarZ1",&phi1);
  
  treedata->SetBranchAddress("ZZPt",&pt4l);
  //treedata->SetBranchAddress("ZZRapidity",&Y4l);
  treedata->SetBranchAddress("Fisher",&fisher);

  treedata->SetBranchAddress("NJets30", &NJets_30);
   
  TFile* newFile  = new TFile(outfile, "RECREATE");
  newFile->cd();
  TTree* newTree = new TTree("data_obs","data_obs"); 
  Double_t CMS_zz4l_mass, melaLD, pseudomelaLD, supermelaLD, CMS_zz4l_massErr, CMS_zz4l_massRelErr;
  Double_t pt = -99, Fisher = -99;
  newTree->Branch("CMS_zz4l_mass",&CMS_zz4l_mass,"CMS_zz4l_mass/D");
  newTree->Branch("CMS_zz4l_massErr",&CMS_zz4l_massErr,"CMS_zz4l_massErr/D");
  newTree->Branch("CMS_zz4l_massRelErr",&CMS_zz4l_massRelErr,"CMS_zz4l_massRelErr/D");
  //newTree->Branch("melaLD",&melaLD,"melaLD/D");
  //newTree->Branch("pseudoMelaLD",&pseudomelaLD,"pseudoMelaLD/D");
  //newTree->Branch("supermelaLD",&supermelaLD,"supermelaLD/D");
  newTree->Branch("melaLD",&melaLD,"melaLD/D");
  //newTree->Branch("p0plus_VAJHU",&p0plus_VAJHU,"p0plus_VAJHU/D");
  //newTree->Branch("bkg_VAMCFMNorm",&bkg_VAMCFM,"bkg_VAMCFM/D");
  //newTree->Branch("ZZVAKD",&ZZVAKD,"ZZVAKD/D");
  newTree->Branch("CMS_zz4l_Fisher",&Fisher,"CMS_zz4l_Fisher/D");
  newTree->Branch("CMS_zz4l_Pt",&pt,"CMS_zz4l_Pt/D");


  cout << inFile << " entries: " << treedata->GetEntries() << endl;
  for(int iEvt=0; iEvt<treedata->GetEntries(); iEvt++){
    //    if(iEvt%5000==0) cout << "event: " << iEvt << endl;
    treedata->GetEntry(iEvt);

    //cout << run << endl;
    
    if (onlyICHEPStat && run>=198049) continue;

    //if( NJets_30 > 0) cout << "NJets_30: " << NJets_30 << endl;

    if ((useJET && VBFtag && NJets_30 < 2) || (useJET && !VBFtag && NJets_30 >= 2)) continue;

    CMS_zz4l_mass = mzz;
    CMS_zz4l_massErr = mzzErr;
    CMS_zz4l_massRelErr = mzzErr/mzz;
    //pseudomelaLD = pseudomela;
    //melaLD = mela;
    //supermelaLD = 0;
    melaLD=p0plus_VAJHU/(p0plus_VAJHU+bkg_VAMCFM);
    //ZZVAKD=p0plus_VAJHU/(bkg_VAMCFM + p0plus_VAJHU);
    if(useJET && !VBFtag) 
      {
	if(pt4l > 200.)
	  {
	    //if pt out of range set to middle of higest bin
	    pt4l = 198.;
	  }
	pt = pt4l;
      }
    if(useJET && VBFtag) 
      {
	if(fisher > 2.)
	  {
	    //if fisher out of range set to middle of higest bin
	    fisher = 1.98;
	  }
	Fisher = fisher;
      }

  
#ifdef LINKMELA
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
#endif
    ++neventOut;
    newTree->Fill();
  }
  newTree->Write("data_obs"); 
  newFile->Close();
  
  cout << "written: " << outfile << " entries: " << neventOut << endl << endl;
  return  neventOut;
  
}



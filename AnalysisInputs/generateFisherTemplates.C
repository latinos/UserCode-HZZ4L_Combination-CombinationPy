#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TString.h"

#include <sstream>
#include <iostream>
#include <vector>

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TROOT.h>
#endif

//Set Input Parameters
#include "Config.h"

//User determined Parameters
const TString destDir = "../CreateDatacards/templates2D/";
int useSqrts=0;                                                   // 0=use 7+8TeV; 1=use 7TeV only, 2 use 8TeV only
TFile* fggH,*fqqH,*fqqZZ,*fggZZ,*fZX,*fZH,*fWH,*fttH;

const TString filePath7TeVM="SamplesfromMichalis/7TeV/";
const TString filePath8TeVM="SamplesfromMichalis/8TeV/";

//Global Parameters (not tested if cause issues if altered)
bool extendToHighMass = true;                                     // Include signal samples above 600 GeV
float highMzz=(extendToHighMass?1600:800);
float mBinSize=2.;
bool useMichalis = false;                                         // False = use CJLST trees, True = use trees from Michalis (other trees can be assigned in buildChain())

//Function Constructors 
void templateOptions(bool debug, bool findAlternatives);          // Sets debug mode (only produces unsmoothed plots) and whether alternatives are produced
void buildChain(TChain* bkgMC, int sampleIndex);                  // Generates TChain from input MC
void makeTemplate(int updown,bool debug);                         // Makes template root files
TH2F* fillTemplate(int sampleIndex, bool isLowMass, int updown);  // Takes input MC to make Fisher v m4L TH2F
TH2F* mergeTemplates(TH2F* lowTemp, TH2F* highTemp);              // Merges low and high mass TH2F plots
TH2F* smoothtemplates(TH2F* inputdata, int sampleIndex);          // Smooths Fisher v m4L plots
TH2F* rebin(TH2F* rebinnedhist);                                  // Smoothing procedure for ggH,qqH
TH2F* rebin_lowstatistics(TH2F* rebinnedhist, int sampleIndex);   // Smoothing procedure for ggZZ,qqZZ,ZH,WH
void analyticfits(int sampleIndex, int updown);                   // Smoothing procedure for Z+X,ttH

//---------------------------------------------------

void generateFisherTemplates() {
  gSystem->Load("../CreateDatacards/CMSSW_5_2_5/lib/slc5_amd64_gcc462/libHiggsAnalysisCombinedLimit.so");

  bool debug=true;
  bool findAlternatives=false;
  templateOptions(debug,findAlternatives);
}
 
void templateOptions(bool debug, bool findAlternatives){
  makeTemplate(0,debug);
  if (findAlternatives && useMichalis){
    makeTemplate(1,debug);
    makeTemplate(-1,debug);
  }
}

//---------------------------------------------------

void buildChain(TChain* bkgMC, int sampleIndex){

  if(!useMichalis){
    TString filePath;
    int nPoints;
    int masses[100];
    if (useSqrts==1){
      if (sampleIndex==0){
	nPoints=nPoints7TeV;
	for (int i=0;i<nPoints;i++){
	  masses[i]=masses7TeV[i];
	}
      } else if (sampleIndex==1){
	nPoints=nVBFPoints7TeV;
	for (int i=0;i<nPoints;i++){
	  masses[i]=VBFmasses7TeV[i];
	}
      } else if (sampleIndex==5 || sampleIndex==6 || sampleIndex==7){
	nPoints=nVHPoints7TeV;
	for (int i=0;i<nPoints;i++){
	  masses[i]=VHmasses7TeV[i];
	}
      }
      filePath = filePath7TeV;
    }
    if (useSqrts==2){
      if (sampleIndex==0){
	nPoints=nPoints8TeV;
	for (int i=0;i<nPoints;i++){
	  masses[i]=masses8TeV[i];
	}
      } else if (sampleIndex==1){
	nPoints=nVBFPoints8TeV;
	for (int i=0;i<nPoints;i++){
	  masses[i]=VBFmasses8TeV[i];
	}
      } else if (sampleIndex==5 || sampleIndex==6 || sampleIndex==7){
	nPoints=nVHPoints8TeV;
	for (int i=0;i<nPoints;i++){
	  masses[i]=VHmasses8TeV[i];
	}
      }
      filePath = filePath8TeV;
    }
    if (useSqrts==0){
      if (sampleIndex==0){
	nPoints=nPoints7TeV+nPoints8TeV;
	for (int i=0;i<nPoints7TeV;i++){
	  masses[i]=masses7TeV[i];
	}
	for (int i=0;i<nPoints8TeV;i++){
	  masses[i+nPoints7TeV]=masses8TeV[i];
	}
      } else if (sampleIndex==1){
	nPoints=nVBFPoints7TeV+nVBFPoints8TeV;
	for (int i=0;i<nVBFPoints7TeV;i++){
	  masses[i]=VBFmasses7TeV[i];
	}
	for (int i=0;i<nVBFPoints8TeV;i++){
	  masses[i+nVBFPoints7TeV]=VBFmasses8TeV[i];
	}
      } else if (sampleIndex==5 || sampleIndex==6 || sampleIndex==7){
	nPoints=nVHPoints7TeV+nVHPoints8TeV;
	for (int i=0;i<nVHPoints7TeV;i++){
	  masses[i]=VHmasses7TeV[i];
	}
	for (int i=0;i<nVHPoints8TeV;i++){
	  masses[i+nVHPoints7TeV]=VHmasses8TeV[i];
	}
      }
    }
    if (sampleIndex!=2 && sampleIndex!=3 && sampleIndex!=4){
      for (int i=0; i<nPoints; i++){
	char tmp_finalInPath4mu[200],tmp_finalInPath4e[200],tmp_finalInPath2mu2e[200];
	string finalInPath4mu,finalInPath4e,finalInPath2mu2e;
	if (sampleIndex==0){
	  sprintf(tmp_finalInPath4mu,"4mu/HZZ4lTree_H%i.root",masses[i]);
	  sprintf(tmp_finalInPath4e,"4e/HZZ4lTree_H%i.root",masses[i]);
	  sprintf(tmp_finalInPath2mu2e,"2mu2e/HZZ4lTree_H%i.root",masses[i]);
	}else if (sampleIndex==1){
	  sprintf(tmp_finalInPath4mu,"4mu/HZZ4lTree_VBFH%i.root",masses[i]);
	  sprintf(tmp_finalInPath4e,"4e/HZZ4lTree_VBFH%i.root",masses[i]);
	  sprintf(tmp_finalInPath2mu2e,"2mu2e/HZZ4lTree_VBFH%i.root",masses[i]);
	}else if (sampleIndex==5){
	  sprintf(tmp_finalInPath4mu,"4mu/HZZ4lTree_ZH%i.root",masses[i]);
	  sprintf(tmp_finalInPath4e,"4e/HZZ4lTree_ZH%i.root",masses[i]);
	  sprintf(tmp_finalInPath2mu2e,"2mu2e/HZZ4lTree_ZH%i.root",masses[i]);
	}else if (sampleIndex==6){
	  sprintf(tmp_finalInPath4mu,"4mu/HZZ4lTree_WH%i.root",masses[i]);
	  sprintf(tmp_finalInPath4e,"4e/HZZ4lTree_WH%i.root",masses[i]);
	  sprintf(tmp_finalInPath2mu2e,"2mu2e/HZZ4lTree_WH%i.root",masses[i]);
	}else if (sampleIndex==7){
	  sprintf(tmp_finalInPath4mu,"4mu/HZZ4lTree_ttH%i.root",masses[i]);
	  sprintf(tmp_finalInPath4e,"4e/HZZ4lTree_ttH%i.root",masses[i]);
	  sprintf(tmp_finalInPath2mu2e,"2mu2e/HZZ4lTree_ttH%i.root",masses[i]);
	}
	if (useSqrts!=0){
	  finalInPath4mu = filePath + tmp_finalInPath4mu;
	  finalInPath4e = filePath + tmp_finalInPath4e;
	  finalInPath2mu2e = filePath + tmp_finalInPath2mu2e;
	} else if (useSqrts==0){
	  if ((sampleIndex==0 && i<nPoints7TeV) || (sampleIndex==1 && i<nVBFPoints7TeV) || ((sampleIndex==5 || sampleIndex==6 || sampleIndex==7) && i<nVHPoints7TeV)){
	    finalInPath4mu = filePath7TeV + tmp_finalInPath4mu;
	    finalInPath4e = filePath7TeV + tmp_finalInPath4e;
	    finalInPath2mu2e = filePath7TeV + tmp_finalInPath2mu2e;
	  } else{
	    finalInPath4mu = filePath8TeV + tmp_finalInPath4mu;
	    finalInPath4e = filePath8TeV + tmp_finalInPath4e;
	    finalInPath2mu2e = filePath8TeV + tmp_finalInPath2mu2e;
	  }
	}
	bkgMC->Add(finalInPath4mu.c_str());
	bkgMC->Add(finalInPath4e.c_str());
	bkgMC->Add(finalInPath2mu2e.c_str());
      }
    }
    else if (sampleIndex==2 || sampleIndex==3){
      if(useSqrts<2){
	bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_ZZTo4mu.root");
	bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_ZZTo4e.root");
	bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_ZZTo2e2mu.root");
	bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_ZZTo4tau.root");
	bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_ZZTo2mu2tau.root");
	bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_ZZTo2e2tau.root");
	bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_ZZTo4mu.root");
	bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_ZZTo4e.root");
	bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_ZZTo2e2mu.root");
	bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_ZZTo4tau.root");
	bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_ZZTo2mu2tau.root");
	bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_ZZTo2e2tau.root");	
	bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_ZZTo4mu.root");
	bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_ZZTo4e.root");
	bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_ZZTo2e2mu.root");
	bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_ZZTo4tau.root");
	bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_ZZTo2mu2tau.root");
	bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_ZZTo2e2tau.root");
      }	
      if(useSqrts%2==0){
	bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_ZZTo4mu.root");
	bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_ZZTo4e.root");
	bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_ZZTo2e2mu.root");
	bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_ZZTo4tau.root");
	bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_ZZTo2mu2tau.root");
	bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_ZZTo2e2tau.root");
	bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_ZZTo4mu.root");
	bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_ZZTo4e.root");
	bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_ZZTo2e2mu.root");
	bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_ZZTo4tau.root");
	bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_ZZTo2mu2tau.root");
	bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_ZZTo2e2tau.root");	
	bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_ZZTo4mu.root");
	bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_ZZTo4e.root");
	bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_ZZTo2e2mu.root");
	bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_ZZTo4tau.root");
	bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_ZZTo2mu2tau.root");
	bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_ZZTo2e2tau.root");
      }
    }
    else if (sampleIndex==4){
      if(useSqrts<2){
	bkgMC->Add(filePath7TeV + "CR/HZZ4lTree_DoubleEle_CREEEEosTree.root");
	bkgMC->Add(filePath7TeV + "CR/HZZ4lTree_DoubleEle_CREEEEssTree.root");
	bkgMC->Add(filePath7TeV + "CR/HZZ4lTree_DoubleOr_CREEMMosTree.root");
	bkgMC->Add(filePath7TeV + "CR/HZZ4lTree_DoubleOr_CREEMMssTree.root");
	bkgMC->Add(filePath7TeV + "CR/HZZ4lTree_DoubleOr_CRMMEEosTree.root");
	bkgMC->Add(filePath7TeV + "CR/HZZ4lTree_DoubleOr_CRMMEEssTree.root");
	bkgMC->Add(filePath7TeV + "CR/HZZ4lTree_DoubleMu_CRMMMMosTree.root");
	bkgMC->Add(filePath7TeV + "CR/HZZ4lTree_DoubleMu_CRMMMMssTree.root");
      }
      if(useSqrts%2==0){
	bkgMC->Add(filePath8TeV + "CR/HZZ4lTree_DoubleEle_CREEEEosTree.root");
	bkgMC->Add(filePath8TeV + "CR/HZZ4lTree_DoubleEle_CREEEEssTree.root");
	bkgMC->Add(filePath8TeV + "CR/HZZ4lTree_DoubleOr_CREEMMosTree.root");
	bkgMC->Add(filePath8TeV + "CR/HZZ4lTree_DoubleOr_CREEMMssTree.root");
	bkgMC->Add(filePath8TeV + "CR/HZZ4lTree_DoubleOr_CRMMEEosTree.root");
	bkgMC->Add(filePath8TeV + "CR/HZZ4lTree_DoubleOr_CRMMEEssTree.root");
	bkgMC->Add(filePath8TeV + "CR/HZZ4lTree_DoubleMu_CRMMMMosTree.root");
	bkgMC->Add(filePath8TeV + "CR/HZZ4lTree_DoubleMu_CRMMMMssTree.root");
      }
    }
  }
  if (useMichalis){
    if(sampleIndex==0){
      //7TeV
      if(useSqrts<2){
	bkgMC->Add(filePath7TeVM + "ggH115.root");
	bkgMC->Add(filePath7TeVM + "ggH120.root");
	bkgMC->Add(filePath7TeVM + "ggH124.root");
	bkgMC->Add(filePath7TeVM + "ggH125.root");
	bkgMC->Add(filePath7TeVM + "ggH126.root");
	bkgMC->Add(filePath7TeVM + "ggH130.root");
	bkgMC->Add(filePath7TeVM + "ggH140.root");
	bkgMC->Add(filePath7TeVM + "ggH150.root");
	bkgMC->Add(filePath7TeVM + "ggH160.root");
	bkgMC->Add(filePath7TeVM + "ggH170.root");
	bkgMC->Add(filePath7TeVM + "ggH180.root");
	bkgMC->Add(filePath7TeVM + "ggH190.root");
	bkgMC->Add(filePath7TeVM + "ggH200.root");
	bkgMC->Add(filePath7TeVM + "ggH210.root");
	bkgMC->Add(filePath7TeVM + "ggH220.root");
	bkgMC->Add(filePath7TeVM + "ggH230.root");
	bkgMC->Add(filePath7TeVM + "ggH250.root");
	bkgMC->Add(filePath7TeVM + "ggH275.root");
	bkgMC->Add(filePath7TeVM + "ggH300.root");
	bkgMC->Add(filePath7TeVM + "ggH325.root");
	bkgMC->Add(filePath7TeVM + "ggH350.root");
	bkgMC->Add(filePath7TeVM + "ggH400.root");
	bkgMC->Add(filePath7TeVM + "ggH425.root");
	bkgMC->Add(filePath7TeVM + "ggH450.root");
	bkgMC->Add(filePath7TeVM + "ggH475.root");
	bkgMC->Add(filePath7TeVM + "ggH525.root");
	bkgMC->Add(filePath7TeVM + "ggH550.root");
	bkgMC->Add(filePath7TeVM + "ggH575.root");
	bkgMC->Add(filePath7TeVM + "ggH600.root");
	if (extendToHighMass) {
	  bkgMC->Add(filePath7TeVM + "ggH650.root");
	  bkgMC->Add(filePath7TeVM + "ggH700.root");
	  bkgMC->Add(filePath7TeVM + "ggH750.root");
	  bkgMC->Add(filePath7TeVM + "ggH800.root");
	  bkgMC->Add(filePath7TeVM + "ggH900.root");
	  bkgMC->Add(filePath7TeVM + "ggH1000.root");
	}
      }
      if(useSqrts%2==0){
	//8TeV
	bkgMC->Add(filePath8TeVM + "ggH118.root");      
	bkgMC->Add(filePath8TeVM + "ggH119.root");
	bkgMC->Add(filePath8TeVM + "ggH120.root");
	bkgMC->Add(filePath8TeVM + "ggH121.root");
	bkgMC->Add(filePath8TeVM + "ggH122.root");
	bkgMC->Add(filePath8TeVM + "ggH123.root");
	bkgMC->Add(filePath8TeVM + "ggH124.root");
	bkgMC->Add(filePath8TeVM + "ggH125.root");
	bkgMC->Add(filePath8TeVM + "ggH126.root");
	bkgMC->Add(filePath8TeVM + "ggH127.root");
	bkgMC->Add(filePath8TeVM + "ggH128.root");
	bkgMC->Add(filePath8TeVM + "ggH129.root");
	bkgMC->Add(filePath8TeVM + "ggH130.root");
	bkgMC->Add(filePath8TeVM + "ggH135.root");
	bkgMC->Add(filePath8TeVM + "ggH140.root");
	bkgMC->Add(filePath8TeVM + "ggH145.root");
	bkgMC->Add(filePath8TeVM + "ggH150.root");
	bkgMC->Add(filePath8TeVM + "ggH160.root");
	bkgMC->Add(filePath8TeVM + "ggH170.root");
	bkgMC->Add(filePath8TeVM + "ggH180.root");
	bkgMC->Add(filePath8TeVM + "ggH190.root");
	bkgMC->Add(filePath8TeVM + "ggH200.root");
	bkgMC->Add(filePath8TeVM + "ggH220.root");
	bkgMC->Add(filePath8TeVM + "ggH250.root");
	bkgMC->Add(filePath8TeVM + "ggH275.root");
	bkgMC->Add(filePath8TeVM + "ggH300.root");
	bkgMC->Add(filePath8TeVM + "ggH325.root");
	bkgMC->Add(filePath8TeVM + "ggH375.root");
	bkgMC->Add(filePath8TeVM + "ggH400.root");
	bkgMC->Add(filePath8TeVM + "ggH425.root");
	bkgMC->Add(filePath8TeVM + "ggH450.root");
	bkgMC->Add(filePath8TeVM + "ggH475.root");
	bkgMC->Add(filePath8TeVM + "ggH500.root");
	bkgMC->Add(filePath8TeVM + "ggH525.root");
	bkgMC->Add(filePath8TeVM + "ggH550.root");
	bkgMC->Add(filePath8TeVM + "ggH575.root");
	bkgMC->Add(filePath8TeVM + "ggH600.root");
	if (extendToHighMass) {
	  bkgMC->Add(filePath8TeVM + "ggH650.root");
	  bkgMC->Add(filePath8TeVM + "ggH700.root");
	  bkgMC->Add(filePath8TeVM + "ggH750.root");
	  bkgMC->Add(filePath8TeVM + "ggH800.root");
	  bkgMC->Add(filePath8TeVM + "ggH850.root");
	  bkgMC->Add(filePath8TeVM + "ggH900.root");
	  bkgMC->Add(filePath8TeVM + "ggH1000.root");
	}
      }
    } else   if(sampleIndex==1){
      //7TeV
      if(useSqrts<2){
	bkgMC->Add(filePath7TeVM + "qqH115.root");
	bkgMC->Add(filePath7TeVM + "qqH120.root");
	bkgMC->Add(filePath7TeVM + "qqH125.root");
	bkgMC->Add(filePath7TeVM + "qqH130.root");
	bkgMC->Add(filePath7TeVM + "qqH140.root");
	bkgMC->Add(filePath7TeVM + "qqH150.root");
	bkgMC->Add(filePath7TeVM + "qqH160.root");
	bkgMC->Add(filePath7TeVM + "qqH170.root");
	bkgMC->Add(filePath7TeVM + "qqH180.root");
	bkgMC->Add(filePath7TeVM + "qqH190.root");
	bkgMC->Add(filePath7TeVM + "qqH200.root");
	bkgMC->Add(filePath7TeVM + "qqH210.root");
	bkgMC->Add(filePath7TeVM + "qqH220.root");
	bkgMC->Add(filePath7TeVM + "qqH230.root");
	bkgMC->Add(filePath7TeVM + "qqH250.root");
	bkgMC->Add(filePath7TeVM + "qqH275.root");
	bkgMC->Add(filePath7TeVM + "qqH300.root");
	bkgMC->Add(filePath7TeVM + "qqH325.root");
	bkgMC->Add(filePath7TeVM + "qqH350.root");
	bkgMC->Add(filePath7TeVM + "qqH400.root");
	bkgMC->Add(filePath7TeVM + "qqH425.root");
	bkgMC->Add(filePath7TeVM + "qqH450.root");
	bkgMC->Add(filePath7TeVM + "qqH475.root");
	bkgMC->Add(filePath7TeVM + "qqH575.root");
	bkgMC->Add(filePath7TeVM + "qqH600.root");
	if (extendToHighMass) {
	  bkgMC->Add(filePath7TeVM + "qqH650.root");
	  bkgMC->Add(filePath7TeVM + "qqH700.root");
	  bkgMC->Add(filePath7TeVM + "qqH800.root");
	  bkgMC->Add(filePath7TeVM + "qqH900.root");
	  bkgMC->Add(filePath7TeVM + "qqH1000.root");
	}
      }
      if(useSqrts%2==0){
	//8TeV
	bkgMC->Add(filePath8TeVM + "qqH118.root");
	bkgMC->Add(filePath8TeVM + "qqH119.root");
	bkgMC->Add(filePath8TeVM + "qqH120.root");
	bkgMC->Add(filePath8TeVM + "qqH121.root");
	bkgMC->Add(filePath8TeVM + "qqH122.root");
	bkgMC->Add(filePath8TeVM + "qqH123.root");
	bkgMC->Add(filePath8TeVM + "qqH124.root");
	bkgMC->Add(filePath8TeVM + "qqH125.root");
	bkgMC->Add(filePath8TeVM + "qqH126.root");
	bkgMC->Add(filePath8TeVM + "qqH127.root");
	bkgMC->Add(filePath8TeVM + "qqH128.root");
	bkgMC->Add(filePath8TeVM + "qqH129.root");
	bkgMC->Add(filePath8TeVM + "qqH130.root");
	bkgMC->Add(filePath8TeVM + "qqH135.root");
	bkgMC->Add(filePath8TeVM + "qqH140.root");
	bkgMC->Add(filePath8TeVM + "qqH145.root");
	bkgMC->Add(filePath8TeVM + "qqH150.root");
	bkgMC->Add(filePath8TeVM + "qqH160.root");
	bkgMC->Add(filePath8TeVM + "qqH170.root");
	bkgMC->Add(filePath8TeVM + "qqH180.root");
	bkgMC->Add(filePath8TeVM + "qqH190.root");
	bkgMC->Add(filePath8TeVM + "qqH200.root");
	bkgMC->Add(filePath8TeVM + "qqH220.root");
	bkgMC->Add(filePath8TeVM + "qqH250.root");
	bkgMC->Add(filePath8TeVM + "qqH275.root");
	bkgMC->Add(filePath8TeVM + "qqH300.root");
	bkgMC->Add(filePath8TeVM + "qqH325.root");
	bkgMC->Add(filePath8TeVM + "qqH350.root");
	bkgMC->Add(filePath8TeVM + "qqH375.root");
	bkgMC->Add(filePath8TeVM + "qqH400.root");
	bkgMC->Add(filePath8TeVM + "qqH425.root");
	bkgMC->Add(filePath8TeVM + "qqH450.root");
	bkgMC->Add(filePath8TeVM + "qqH475.root");
	bkgMC->Add(filePath8TeVM + "qqH500.root");
	bkgMC->Add(filePath8TeVM + "qqH525.root");
	bkgMC->Add(filePath8TeVM + "qqH550.root");
	bkgMC->Add(filePath8TeVM + "qqH575.root");
	bkgMC->Add(filePath8TeVM + "qqH600.root");
	if (extendToHighMass) {
	  bkgMC->Add(filePath8TeVM + "qqH650.root");
	  bkgMC->Add(filePath8TeVM + "qqH700.root");
	  bkgMC->Add(filePath8TeVM + "qqH750.root");
	  bkgMC->Add(filePath8TeVM + "qqH800.root");
	  bkgMC->Add(filePath8TeVM + "qqH850.root");
	  bkgMC->Add(filePath8TeVM + "qqH900.root");
	  bkgMC->Add(filePath8TeVM + "qqH1000.root");
	}
      }
    } else if (sampleIndex==2){
      //7TeV
      if(useSqrts<2){
	bkgMC->Add(filePath7TeVM + "ZZ2e2mu.root");
	bkgMC->Add(filePath7TeVM + "ZZ2e2tau.root");
	bkgMC->Add(filePath7TeVM + "ZZ2mu2tau.root");
	bkgMC->Add(filePath7TeVM + "ZZ4e.root");
	bkgMC->Add(filePath7TeVM + "ZZ4mu.root");
	bkgMC->Add(filePath7TeVM + "ZZ4tau.root");
      }
      //8TeV
      if(useSqrts%2==0){
	bkgMC->Add(filePath8TeVM + "ZZ2e2mu.root");
	bkgMC->Add(filePath8TeVM + "ZZ2e2tau.root");
	bkgMC->Add(filePath8TeVM + "ZZ2mu2tau.root");
	bkgMC->Add(filePath8TeVM + "ZZ4e.root");
	bkgMC->Add(filePath8TeVM + "ZZ4mu.root");
	bkgMC->Add(filePath8TeVM + "ZZ4tau.root");
      }
    } else if (sampleIndex==3){
      //7TeV
      if(useSqrts<2){
	bkgMC->Add(filePath7TeVM + "GGZZ2L2L.root");
	bkgMC->Add(filePath7TeVM + "GGZZ4L.root");
      }
      //8TeV
      if(useSqrts%2==0){
	bkgMC->Add(filePath8TeVM + "GGZZ2L2L.root");
	bkgMC->Add(filePath8TeVM + "GGZZ4L.root");
      }
    } else if (sampleIndex==4){
      //7TeV
      if(useSqrts<2){
	bkgMC->Add(filePath7TeVM + "DATA.root");
      }
      //8TeV
      if(useSqrts%2==0){
	bkgMC->Add(filePath8TeVM + "DATA.root");
      }
    } else if (sampleIndex==5 || sampleIndex==6 || sampleIndex==7){
      //7TeV
      if(useSqrts<2){
	bkgMC->Add(filePath7TeVM + "VH126.root");
	bkgMC->Add(filePath7TeVM + "VH127.root");
	bkgMC->Add(filePath7TeVM + "VH128.root");
	bkgMC->Add(filePath7TeVM + "VH129.root");
	bkgMC->Add(filePath7TeVM + "VH130.root");
	bkgMC->Add(filePath7TeVM + "VH135.root");
	bkgMC->Add(filePath7TeVM + "VH140.root");
	bkgMC->Add(filePath7TeVM + "VH145.root");
	bkgMC->Add(filePath7TeVM + "VH150.root");
	bkgMC->Add(filePath7TeVM + "VH160.root");
	bkgMC->Add(filePath7TeVM + "VH170.root");
	bkgMC->Add(filePath7TeVM + "VH180.root");
	bkgMC->Add(filePath7TeVM + "VH190.root");
	bkgMC->Add(filePath7TeVM + "VH200.root");
	bkgMC->Add(filePath7TeVM + "VH210.root");
	bkgMC->Add(filePath7TeVM + "VH220.root");
	bkgMC->Add(filePath7TeVM + "VH230.root");
	bkgMC->Add(filePath7TeVM + "VH250.root");
	bkgMC->Add(filePath7TeVM + "VH275.root");
	bkgMC->Add(filePath7TeVM + "VH300.root");      
      }
      if(useSqrts%2==0){
	bkgMC->Add(filePath8TeVM + "VH126.root");
	bkgMC->Add(filePath8TeVM + "VH127.root");
	bkgMC->Add(filePath8TeVM + "VH128.root");
	bkgMC->Add(filePath8TeVM + "VH129.root");
	bkgMC->Add(filePath8TeVM + "VH130.root");
	bkgMC->Add(filePath8TeVM + "VH135.root");
	bkgMC->Add(filePath8TeVM + "VH140.root");
	bkgMC->Add(filePath8TeVM + "VH145.root");
	bkgMC->Add(filePath8TeVM + "VH150.root");
	bkgMC->Add(filePath8TeVM + "VH160.root");
	bkgMC->Add(filePath8TeVM + "VH170.root");
	bkgMC->Add(filePath8TeVM + "VH180.root");
	bkgMC->Add(filePath8TeVM + "VH190.root");
	bkgMC->Add(filePath8TeVM + "VH200.root");
	bkgMC->Add(filePath8TeVM + "VH210.root");
	bkgMC->Add(filePath8TeVM + "VH220.root");
	bkgMC->Add(filePath8TeVM + "VH230.root");
	bkgMC->Add(filePath8TeVM + "VH250.root");
	bkgMC->Add(filePath8TeVM + "VH275.root");
	bkgMC->Add(filePath8TeVM + "VH300.root");
      }
    }
  }
}

//---------------------------------------------------

void makeTemplate(int updown, bool debug){

  TString jes;
  if (updown==1){
    jes="_up";
  }
  else if(updown==-1){
    jes="_down";
  }

  TString debugname;
  if (debug) debugname="_unnormalized";

  fqqH = new TFile(destDir + "qqH_fishertemplate"+jes+debugname+".root","RECREATE");
  fggH = new TFile(destDir + "ggH_fishertemplate"+jes+debugname+".root","RECREATE");
  fqqZZ = new TFile(destDir + "qqZZ_fishertemplate"+jes+debugname+".root","RECREATE");
  fggZZ = new TFile(destDir + "ggZZ_fishertemplate"+jes+debugname+".root","RECREATE");
  fZX = new TFile(destDir + "Z+X_fishertemplate"+jes+debugname+".root","RECREATE");
  fZH = new TFile(destDir + "ZH_fishertemplate"+jes+debugname+".root","RECREATE");
  fWH = new TFile(destDir + "WH_fishertemplate"+jes+debugname+".root","RECREATE");
  fttH = new TFile(destDir + "ttH_fishertemplate"+jes+debugname+".root","RECREATE");
  
  TH2F* low,*high,*H_Fisher;
  
  
  // =========================
  // ggH
  
  low = fillTemplate(0,true,updown);
  high = fillTemplate(0,false,updown);
  H_Fisher = mergeTemplates(low,high);
  
  if (!debug) smoothtemplates(H_Fisher,0);

  fggH->cd();
  H_Fisher->Write("H_Fisher");
  fggH->Close();
  
  // ==========================
  // qqH

  low = fillTemplate(1,true,updown);
  high = fillTemplate(1,false,updown);
  H_Fisher = mergeTemplates(low,high);
  
  if (!debug) smoothtemplates(H_Fisher,1);

  fqqH->cd();
  H_Fisher->Write("H_Fisher");
  fqqH->Close();

  // ==========================
  // qqZZ

  low = fillTemplate(2,true,updown);
  high = fillTemplate(2,false,updown);
  H_Fisher = mergeTemplates(low,high);

  if (!debug) smoothtemplates(H_Fisher,2);

  fqqZZ->cd();
  H_Fisher->Write("H_Fisher");
  fqqZZ->Close();

  // ==========================
  // ggZZ

  low = fillTemplate(3,true,updown);
  high = fillTemplate(3,false,updown);
  H_Fisher = mergeTemplates(low,high);

  if (!debug) smoothtemplates(H_Fisher,3);

  fggZZ->cd();
  H_Fisher->Write("H_Fisher");
  fggZZ->Close();
  
  // ==========================
  // Z+X

  if (debug){
    low = fillTemplate(4,true,updown);
    high = fillTemplate(4,false,updown);
    H_Fisher = mergeTemplates(low,high);
    
    fZX->cd();
    H_Fisher = mergeTemplates(low,high);
    fZX->Close();
  }
  else{
    analyticfits(4,updown);
  }
  
  // ==========================
  // ZH
  low = fillTemplate(5,true,updown);
  high = fillTemplate(5,false,updown);
  H_Fisher = mergeTemplates(low,high);

  if (!debug) smoothtemplates(H_Fisher,5);

  fZH->cd();
  H_Fisher->Write("H_Fisher");
  fZH->Close();

  // ==========================
  // WH
  low = fillTemplate(6,true,updown);
  high = fillTemplate(6,false,updown);
  H_Fisher = mergeTemplates(low,high);

  if (!debug) smoothtemplates(H_Fisher,6);

  fWH->cd();
  H_Fisher->Write("H_Fisher");
  fWH->Close();
  
  // ==========================
  // ttH

  if (debug){
    low = fillTemplate(7,true,updown);
    high = fillTemplate(7,false,updown);
    H_Fisher = mergeTemplates(low,high);
    
    fttH->cd();
    H_Fisher->Write("H_Fisher");
    fttH->Close();
  }
  else{
    analyticfits(7,updown);
  }
}

//---------------------------------------------------

TH2F* fillTemplate(int sampleIndex,bool isLowMass,int updown){
  TChain* bkgMC;
  if (!useMichalis) bkgMC = new TChain("SelectedTree");
  if (useMichalis) bkgMC = new TChain("FourLeptonTreeProducer/tree");
  buildChain(bkgMC,sampleIndex);

  cout << "Chain for " << sampleIndex << " " << isLowMass << " " << bkgMC->GetEntries() << endl;
  bkgMC->ls();

  //bkgMC->Sumw2();

  double mmass,mdeta,mmJJ,mw,mFisher,mheff,mprocessID;
  float mass,deta,mJJ,w,Fisher,processID;
  double lz2chg,lz2mass,lminllmass; // To deal with Z+X background
  double z11,z12,z21,z22;
  int njets;
  string channel;

  if(!useMichalis){
    bkgMC->SetBranchAddress("ZZMass",&mass);
    bkgMC->SetBranchAddress("NJets",&njets);
    bkgMC->SetBranchAddress("DiJetDEta",&deta);
    bkgMC->SetBranchAddress("DiJetMass",&mJJ);
    bkgMC->SetBranchAddress("MC_weight",&w);
    bkgMC->SetBranchAddress("genProcessId",&processID);
    lz2chg=1;
    lz2mass=13;
    lminllmass=5;
  }else if(useMichalis){
    if (sampleIndex!=4){
      bkgMC->SetBranchAddress("H_Mass",&mmass);
      bkgMC->SetBranchAddress("H_NJets",&njets);
      if(updown==0){
	bkgMC->SetBranchAddress("H_DEta",&mdeta);
	bkgMC->SetBranchAddress("H_MJJ",&mmJJ);
      }
      else if(updown==1){
	bkgMC->SetBranchAddress("H_DEtaUp",&mdeta);
	bkgMC->SetBranchAddress("H_MJJUp",&mmJJ);
      }
      else if(updown==-1){
	bkgMC->SetBranchAddress("H_DEtaDwn",&mdeta);
	bkgMC->SetBranchAddress("H_MJJDwn",&mmJJ);
      }
      bkgMC->SetBranchAddress("eventWeight",&mw);
      bkgMC->SetBranchAddress("H_eff",&mheff);
      bkgMC->SetBranchAddress("H_Z1_leg1_PdgId",&z11);
      bkgMC->SetBranchAddress("H_Z1_leg2_PdgId",&z12);
      bkgMC->SetBranchAddress("H_Z2_leg1_PdgId",&z21);
      bkgMC->SetBranchAddress("H_Z2_leg2_PdgId",&z22);
      lz2chg=1;
      lz2mass=13;
      lminllmass=5;
      if (sampleIndex==5 || sampleIndex==6 || sampleIndex==7){
	bkgMC->SetBranchAddress("genHProcess",&mprocessID);
      }
    }
    else{
      bkgMC->SetBranchAddress("HLoose_Mass",&mmass);
      bkgMC->SetBranchAddress("HLoose_NJets",&njets);
      if(updown==0){
	bkgMC->SetBranchAddress("HLoose_DEta",&mdeta);
	bkgMC->SetBranchAddress("HLoose_MJJ",&mmJJ);
      }
      else if(updown==1){
	bkgMC->SetBranchAddress("HLoose_DEtaUp",&mdeta);
	bkgMC->SetBranchAddress("HLoose_MJJUp",&mmJJ);
      }
      else if(updown==-1){
	bkgMC->SetBranchAddress("HLoose_DEtaDwn",&mdeta);
	bkgMC->SetBranchAddress("HLoose_MJJDwn",&mmJJ);
      }
      bkgMC->SetBranchAddress("eventWeight",&mw);
      bkgMC->SetBranchAddress("H_eff",&mheff);
      bkgMC->SetBranchAddress("HLoose_Z1_leg1_PdgId",&z11);
      bkgMC->SetBranchAddress("HLoose_Z1_leg2_PdgId",&z12);
      bkgMC->SetBranchAddress("HLoose_Z2_leg1_PdgId",&z21);
      bkgMC->SetBranchAddress("HLoose_Z2_leg2_PdgId",&z22);
      bkgMC->SetBranchAddress("HLoose_Z2_Charge",&lz2chg);
      bkgMC->SetBranchAddress("HLoose_Z2_Mass",&lz2mass);
      bkgMC->SetBranchAddress("HLoose_MinOSPairMass",&lminllmass);
    }
  }

  TH2F* bkgHist;
  if(!isLowMass){
    bkgHist = new TH2F("bkgHisto","bkgHisto",int((highMzz-180.)/mBinSize+0.5),180,highMzz,50,0,2);
  }
  else{
    bkgHist = new TH2F("bkgHisto","bkgHisto",int((180-100)/mBinSize+0.5),100,180,50,0,2);
  }

  bkgHist->Sumw2();

  //Fill histogram
  for(int i=0; i<bkgMC->GetEntries(); i++){
    bkgMC->GetEntry(i);
    if (i%50000==0) cout << "event: " << i << "/" << bkgMC->GetEntries() << endl;
    //if (i%5000==0) cout<< mass <<" "<<deta<<" "<<mJJ<<" "<<njets<<endl;
    if (mass<100 || deta<=-99 || mJJ<=-99 || njets<2) continue;
    if (lz2chg==0 || lz2mass<=12 || lminllmass<=4) continue;
    if(!useMichalis){
      if ((sampleIndex==5 && processID!=24) || (sampleIndex==6 && processID!=26) || (sampleIndex==7 && processID!=121)) continue;
      Fisher=0.09407*fabs(deta)+0.00041581*mJJ;
      bkgHist->Fill(mass,Fisher,w);
    }else if(useMichalis){
      if ((sampleIndex==5 && mprocessID!=24) || (sampleIndex==6 && mprocessID!=26) || (sampleIndex==7 && mprocessID!=121)) continue;
      mFisher=0.09407*fabs(mdeta)+0.00041581*mmJJ;
      if (abs(z11)==11 && abs(z12)==11 && abs(z21)==11 && abs(z22)==11) channel="4e";
      if (abs(z11)==13 && abs(z12)==13 && abs(z21)==13 && abs(z22)==13) channel="4mu";
      if ((abs(z11)==11 && abs(z12)==11 && abs(z21)==13 && abs(z22)==13) || (abs(z11)==13 && abs(z12)==13 && abs(z21)==11 && abs(z22)==11)) channel="2e2mu";
      if (channel=="4e" || channel=="4mu" || channel=="2e2mu"){
	bkgHist->Fill(mmass,mFisher,mw*mheff);
      }
      else continue;
    }
  }
  
  return bkgHist;
}

//---------------------------------------------------

TH2F* mergeTemplates(TH2F* lowTemp, TH2F* highTemp){

  int nYbins=lowTemp->GetNbinsY();
  if (highTemp->GetNbinsY()!=nYbins) {
    cout << "ERROR: mergeTemplates: incorrect binning " << endl;
    abort();
  }

  TH2F* H_Fisher = new TH2F("H_Fisher","H_Fisher",int((highMzz-100.)/mBinSize +0.5),100,highMzz,nYbins,0,2);

  // copy lowmass into H_Fisher
  for(int i=1; i<=lowTemp->GetNbinsX(); ++i){
    for(int j=1; j<=nYbins; ++j){
      H_Fisher->SetBinContent(i,j, lowTemp->GetBinContent(i,j)  );
    }// end loop over Fisher
  }// end loop over m4L

  // copy high mass into H_Fisher
  for(int i=1; i<=highTemp->GetNbinsX(); ++i){
    for(int j=1; j<=nYbins; ++j){
      H_Fisher->SetBinContent(i+lowTemp->GetNbinsX(),j, highTemp->GetBinContent(i,j)  );
    }// end loop over Fisher
  }// end loop over m4L

  return H_Fisher;
}

//---------------------------------------------------

TH2F* smoothtemplates(TH2F* inputdata, int sampleIndex){
  if(sampleIndex==0 || sampleIndex==1){
    rebin(inputdata);
  }
  else if(sampleIndex==2 || sampleIndex==3 || sampleIndex==5 || sampleIndex==6){
    rebin_lowstatistics(inputdata,sampleIndex);
  }

  return inputdata;
}

//---------------------------------------------------

TH2F* rebin(TH2F* rebinnedHist){

  int nXbins=rebinnedHist->GetNbinsX();
  int nYbins=rebinnedHist->GetNbinsY();

  double norm;
  TH1F* tempProj;

  rebinnedHist->Sumw2();

  //Normalization
  for(int i=1; i<=nXbins; i++){
    tempProj = (TH1F*) rebinnedHist->ProjectionY("tempProj",i,i);
    norm=tempProj->Integral();
    if (norm>0) { // Avoid introducing NaNs in the histogram
      for(int j=1; j<=nYbins; j++){
	rebinnedHist->SetBinContent(i,j, rebinnedHist->GetBinContent(i,j)/norm   );
      }
    }
  }

  TH2F* origHist = new TH2F (*rebinnedHist);

  origHist->Sumw2();
  rebinnedHist->Sumw2();

  int effectiveArea=1;
  double average=0,binsUsed=0;

  //Rebin m4L
  for (int i=1; i <=nXbins; i++){
    for (int j=1; j<=nYbins; j++){
      float binMzz = rebinnedHist->GetBinCenter(i);

      if( binMzz<300 ) continue;
      if( binMzz>=300 && binMzz<350 ) effectiveArea=1;
      if( binMzz>=350 && binMzz<500 ) effectiveArea=3;
      if( binMzz>=500 && binMzz<600 ) effectiveArea=5;
      if( binMzz>=600 && binMzz<800 ) effectiveArea=7;
      if( binMzz>=800 && binMzz<1000) effectiveArea=11;
      if( binMzz>=1000 && binMzz<1200)effectiveArea=15;
      if( binMzz>=1200 && binMzz<1400) effectiveArea=25;
      if( binMzz>1400) effectiveArea=40;
      
      for(int a=-effectiveArea; a<=effectiveArea; a++){
	if(a+i<1 || a+i>nXbins || j>nYbins || j<1) continue;
	average+=origHist->GetBinContent(a+i,j);
	binsUsed++;
      }
      rebinnedHist->SetBinContent(i,j,average/binsUsed);
      average=0;
      binsUsed=0;
    }
  }


  TH2F* Histstg1 = new TH2F (*rebinnedHist);

  Histstg1->Sumw2();
  rebinnedHist->Sumw2();
  
  //Rebin Fisher
  for (int i=1; i<=nXbins; i++){
    for (int j=1; j<=nYbins; j++){
      float binFisher = rebinnedHist->GetYaxis()->GetBinCenter(j);

      if( binFisher<0.2 ) continue;
      if( binFisher>0.2 && binFisher<=1.0) effectiveArea=1;
      if( binFisher>1.0 && binFisher<=1.2) effectiveArea=2;
      if (binFisher>1.2 && binFisher<=1.4) effectiveArea=3;
      if (binFisher>1.4) effectiveArea=5;

      for(int a=-effectiveArea;a<=effectiveArea;a++){
	if(j+a<1 || j+a>nYbins || i>nXbins || i<1) continue;
	average+=Histstg1->GetBinContent(i,j+a);
	binsUsed++;
      }
      rebinnedHist->SetBinContent(i,j,average/binsUsed);
      average=0;
      binsUsed=0;
    }
  }
  
  rebinnedHist->Sumw2();

  //Use average of Nearest Neighbors to fill remaining zeroes
  for(int i=1; i<=nXbins;i++){
    for(int j=1; j<=nYbins;j++){
      float binvalue = rebinnedHist->GetBinContent(i,j);
      if(binvalue!=0) continue;
      for (int i2=-1;i2<=1;i2++){
	if (i2+i<1 || i2+i>nXbins || j<1 || j>nYbins) continue;
	average+=rebinnedHist->GetBinContent(i+i2,j);
	binsUsed++;
      }
      for (int j2=-1;j2<=1;j2++){
	if (i<1 || i>nXbins || j2+j<1 || j2+j>nYbins) continue;
	average+=rebinnedHist->GetBinContent(i,j+j2);
	binsUsed++;
      }
      rebinnedHist->SetBinContent(i,j,average/binsUsed);
      average=0;
      binsUsed=0;
    }
  }

  rebinnedHist->Sumw2();

  //Renormalize
  for(int i=1; i<=nXbins; i++){
    tempProj = (TH1F*) rebinnedHist->ProjectionY("tempProj",i,i);
    norm=tempProj->Integral();
    if (norm>0) { // Avoid introducing NaNs in the histogram
      for(int j=1; j<=nYbins; j++){
	rebinnedHist->SetBinContent(i,j, rebinnedHist->GetBinContent(i,j)/norm   );
      }
    }
  }


  return rebinnedHist;
}

//---------------------------------------------------

TH2F* rebin_lowstatistics(TH2F* finalhist, int sampleIndex){
  int nXbins=finalhist->GetNbinsX();
  int nYbins=finalhist->GetNbinsY();

  finalhist->Sumw2();

  TH2F* origHist = new TH2F (*finalhist);

  origHist->Sumw2();

  //Project Fisher v m4L to 1D Fisher plots for low, high, and full mass range
  TH1F* lowProj = (TH1F*) origHist->ProjectionY("low",1,40);
  TH1F* highProj = (TH1F*) origHist->ProjectionY("high",41,750);
  TH1F* fullProj = (TH1F*) origHist->ProjectionY();
  TH1F* origlowProj = new TH1F (*lowProj);
  TH1F* orighighProj = new TH1F (*highProj);
  TH1F* origfullProj = new TH1F (*fullProj);
  double temp;

  lowProj->Sumw2();
  highProj->Sumw2();
  fullProj->Sumw2();
  origlowProj->Sumw2();
  orighighProj->Sumw2();
  origfullProj->Sumw2();

  //qqZZ + WH
  if(sampleIndex==2 || sampleIndex==6){
    //Fit tails of low and high projections with exponentials to account for low statistics
    TF1 *tailfunc = new TF1("tailfunc","[0]*exp([1]*(x))",0.2,2.0);
    if (sampleIndex==2){
      lowProj->Fit(tailfunc,"","",0.2,0.7);
    }
    else if(sampleIndex==6){
      lowProj->Fit(tailfunc,"","",0.2,0.7);
    }
    TF1 *tailfunc2 = new TF1("tailfunc2","[0]*exp([1]*(x))",0.2,2.0);
    if (sampleIndex==2){
      highProj->Fit(tailfunc2,"","",0.3,1.0);
    }
    else if(sampleIndex==6){
      highProj->Fit(tailfunc2,"","",0.3,0.6);
    }

    //Adjust projections with fitted exponentials
    for(int i=0;i<=nYbins;i++){
      float binFisher = lowProj->GetBinCenter(i);
      if ((sampleIndex==2 && binFisher>0.6) || (sampleIndex==6 && binFisher>0.7)){
	lowProj->SetBinContent(i,tailfunc->Eval(binFisher));
      }
    }
    for(int i=0;i<=nYbins;i++){
      float binFisher = highProj->GetBinCenter(i);
      if ((sampleIndex==2 && binFisher>1.1) || (sampleIndex==6 && binFisher>0.6)){
	highProj->SetBinContent(i,tailfunc2->Eval(binFisher));
      }
    }

    //Fill each mass point with low or high mass projections
    for(int i=1; i<=nXbins;i++){
      for(int j=1; j<=nYbins;j++){
	float binMzz = finalhist->GetBinCenter(i);
	if (binMzz<180){
	  temp=lowProj->GetBinContent(j);
	  finalhist->SetBinContent(i,j,temp);
	}
	else if (binMzz>=180){
	  temp=highProj->GetBinContent(j);
	  finalhist->SetBinContent(i,j,temp);
	}
      }
    }

    //Store plots of fits, in case anything goes wrong
    if(sampleIndex==2){
      fqqZZ->cd();
    }
    else if(sampleIndex==6){
      fWH->cd();
    }
    origlowProj->Write("H_lowraw");
    lowProj->Write("H_lowfit");
    orighighProj->Write("H_highraw");
    highProj->Write("H_highfit");
    
  }
  //ggZZ + ZH
  else if(sampleIndex==3 || sampleIndex==5){
    //Fit tail of full projection with exponential to account for low statistics
    TF1 *tailfunc3 = new TF1("tailfunc3","[0]*exp([1]*(x))",0.3,2.0);
    if(sampleIndex==3){
      fullProj->Fit(tailfunc3,"","",0.5,2.0);
    }
    else if(sampleIndex==5){
      fullProj->Fit(tailfunc3,"","",0.3,2.0);
    }

    //Adjust projection with fitted exponential
    for(int i=0;i<=nYbins;i++){
      float binFisher = fullProj->GetBinCenter(i);
      if ((sampleIndex==3 && binFisher>1.4) || (sampleIndex==5 && binFisher>1.0)){
	fullProj->SetBinContent(i,tailfunc3->Eval(binFisher));
      }
    }
    
    //Fill each mass point with full projection
    for(int i=1; i<=nXbins;i++){
      for(int j=1; j<=nYbins;j++){
	temp=fullProj->GetBinContent(j);
	finalhist->SetBinContent(i,j,temp);
      }
    }

    //Store plots of fits, in case anything goes wrong
    if(sampleIndex==3){
      fggZZ->cd();
    }
    else if(sampleIndex==5){
      fZH->cd();
    }
    origfullProj->Write("H_fullraw");
    fullProj->Write("H_fullfit");

  }

  double norm;
  TH1F* tempProj;

  //Normalize
  for(int i=1; i<=nXbins; i++){
    tempProj = (TH1F*) finalhist->ProjectionY("tempProj",i,i);
    norm=tempProj->Integral();
    if (norm>0) { // Avoid introducing NaNs in the histogram
      for(int j=1; j<=nYbins; j++){
	finalhist->SetBinContent(i,j, finalhist->GetBinContent(i,j)/norm);
      }
    }
  }

  return finalhist;
}

//---------------------------------------------------

void analyticfits(int sampleIndex,int updown){
  TChain* bkgMC;
  if (!useMichalis) bkgMC = new TChain("SelectedTree");
  if (useMichalis) bkgMC = new TChain("FourLeptonTreeProducer/tree");
  buildChain(bkgMC,sampleIndex);

  cout << "Chain for " << sampleIndex << " " << bkgMC->GetEntries() << endl;
  bkgMC->ls();

  double mdeta,mmJJ,mw,mmass,mtempFisher,mheff,lz2chg,lz2mass,lminllmass,mprocessID;
  float deta,mJJ,w,mass,tempFisher,processID;
  double z11,z12,z21,z22;
  int njets;

  if(!useMichalis){
    bkgMC->SetBranchAddress("ZZMass",&mass);
    bkgMC->SetBranchAddress("NJets",&njets);
    bkgMC->SetBranchAddress("DiJetDEta",&deta);
    bkgMC->SetBranchAddress("DiJetMass",&mJJ);
    bkgMC->SetBranchAddress("MC_weight",&w);
    bkgMC->SetBranchAddress("genProcessId",&processID);
    lz2chg=1;
    lz2mass=13;
    lminllmass=5;
  }else if(useMichalis){
    if (sampleIndex!=4){
      bkgMC->SetBranchAddress("H_Mass",&mmass);
      bkgMC->SetBranchAddress("H_NJets",&njets);
      if(updown==0){
	bkgMC->SetBranchAddress("H_DEta",&mdeta);
	bkgMC->SetBranchAddress("H_MJJ",&mmJJ);
      }
      else if(updown==1){
	bkgMC->SetBranchAddress("H_DEtaUp",&mdeta);
	bkgMC->SetBranchAddress("H_MJJUp",&mmJJ);
      }
      else if(updown==-1){
	bkgMC->SetBranchAddress("H_DEtaDwn",&mdeta);
	bkgMC->SetBranchAddress("H_MJJDwn",&mmJJ);
      }
      bkgMC->SetBranchAddress("eventWeight",&mw);
      bkgMC->SetBranchAddress("H_eff",&mheff);
      bkgMC->SetBranchAddress("H_Z1_leg1_PdgId",&z11);
      bkgMC->SetBranchAddress("H_Z1_leg2_PdgId",&z12);
      bkgMC->SetBranchAddress("H_Z2_leg1_PdgId",&z21);
      bkgMC->SetBranchAddress("H_Z2_leg2_PdgId",&z22);
      lz2chg=1;
      lz2mass=13;
      lminllmass=5;
      if (sampleIndex==5 || sampleIndex==6 || sampleIndex==7){
	bkgMC->SetBranchAddress("genHProcess",&mprocessID);
      }
    }
    else{
      bkgMC->SetBranchAddress("HLoose_Mass",&mmass);
      bkgMC->SetBranchAddress("HLoose_NJets",&njets);
      if(updown==0){
	bkgMC->SetBranchAddress("HLoose_DEta",&mdeta);
	bkgMC->SetBranchAddress("HLoose_MJJ",&mmJJ);
      }
      else if(updown==1){
	bkgMC->SetBranchAddress("HLoose_DEtaUp",&mdeta);
	bkgMC->SetBranchAddress("HLoose_MJJUp",&mmJJ);
      }
      else if(updown==-1){
	bkgMC->SetBranchAddress("HLoose_DEtaDwn",&mdeta);
	bkgMC->SetBranchAddress("HLoose_MJJDwn",&mmJJ);
      }
      bkgMC->SetBranchAddress("eventWeight",&mw);
      bkgMC->SetBranchAddress("H_eff",&mheff);
      bkgMC->SetBranchAddress("HLoose_Z2_Charge",&lz2chg);
      bkgMC->SetBranchAddress("HLoose_Z2_Mass",&lz2mass);
      bkgMC->SetBranchAddress("HLoose_MinOSPairMass",&lminllmass);
      bkgMC->SetBranchAddress("HLoose_Z1_leg1_PdgId",&z11);
      bkgMC->SetBranchAddress("HLoose_Z1_leg2_PdgId",&z12);
      bkgMC->SetBranchAddress("HLoose_Z2_leg1_PdgId",&z21);
      bkgMC->SetBranchAddress("HLoose_Z2_leg2_PdgId",&z22);
    }
  }

  RooRealVar Fisher("Fisher","Fisher",0.,2.);
  RooRealVar gm1("gm1","",0.1,0.4);
  RooRealVar gsig1("gsig1","",0.01,2.);
  RooRealVar lm("lm","",0,0.2);
  RooRealVar lsig("lsig","",0.01,2.);
  RooRealVar f1("f1","",0.,1.);

  RooDataSet* data;
  data= new RooDataSet("data","dataset",Fisher);
  string channel;

  //Fill 1D Histogram - Effective Projection of Fisher v m4L plot
  for(int i=0; i<bkgMC->GetEntries(); i++){
    bkgMC->GetEntry(i);
    if(i%50000==0) cout << "event: " << i << "/" << bkgMC->GetEntries() << endl;
      if(!useMichalis){
	if(mass>100 && lz2chg!=0 && lz2mass>12 && lminllmass>4 && njets>=2){
	if (sampleIndex==7 && processID!=121) continue;
	tempFisher=0.09407*fabs(deta)+0.00041581*mJJ;
	Fisher=tempFisher;
	data->add(RooArgSet(Fisher),w);
      }
      if(useMichalis){
	if(mmass>100 && lz2chg!=0 && lz2mass>12 && lminllmass>4 && njets>=2){
	  if (sampleIndex==7 && mprocessID!=121) continue;
	  mtempFisher=0.09407*fabs(mdeta)+0.00041581*mmJJ;
	  Fisher=mtempFisher;
	  if (abs(z11)==11 && abs(z12)==11 && abs(z21)==11 && abs(z22)==11) channel="4e";
	  if (abs(z11)==13 && abs(z12)==13 && abs(z21)==13 && abs(z22)==13) channel="4mu";
	  if ((abs(z11)==11 && abs(z12)==11 && abs(z21)==13 && abs(z22)==13) || (abs(z11)==13 && abs(z12)==13 && abs(z21)==11 && abs(z22)==11)) channel="2e2mu";
	  if (channel=="4e" || channel=="4mu" || channel=="2e2mu"){
	    data->add(RooArgSet(Fisher),mw*mheff);
	  }
	}
      }
    }
  }

  //Fit is a Landau + Gaussian
  RooGaussian* gaus1 = new RooGaussian("gaus1","gaus1",Fisher,gm1,gsig1);
  RooLandau* land1 = new RooLandau("land1","land1",Fisher,lm,lsig);
  RooAddPdf* model =  new RooAddPdf("model","model",*gaus1,*land1,f1);

  model->fitTo(*data);
  RooPlot* Fisherframe = Fisher.frame();
  data->plotOn(Fisherframe);
  model->plotOn(Fisherframe);
  if (sampleIndex==4){
    fZX->cd();
  }
  else if (sampleIndex==7){
    fttH->cd();
  }
  Fisherframe->Write("H_fit");
  
  TH2F* finalhist;
  finalhist = new TH2F("H_Fisher","H_Fisher",750,100,1600,50,0,2);
  float temp; 
  
  TH1* testplot;
  testplot=model->createHistogram("Fisher");
  testplot->Rebin();
  
  for(int i=1;i<=50;i++){
    for(int j=1;j<=750;j++){
      temp=testplot->GetBinContent(i);
      finalhist->SetBinContent(j,i,temp);
    }
  }

  if (sampleIndex==4){
    fZX->cd();
  }
  else if (sampleIndex==7){
    fttH->cd();
  }
  finalhist->Write("H_Fisher");
  if (sampleIndex==4){
    fZX->Close();
  }
  else if (sampleIndex==7){
    fttH->Close();
  }
}

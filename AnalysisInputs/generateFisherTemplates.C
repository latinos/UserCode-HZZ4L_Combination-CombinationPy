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
#include "TStyle.h"

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

//Global Parameters (not tested if cause issues if altered)
bool extendToHighMass = true;                                     // Include signal samples above 600 GeV
float highMzz=(extendToHighMass?1600:800);
float mBinSize=2.;

//Function Constructors 
void templateOptions(bool debug, bool findAlternatives);          // Sets debug mode (only produces unsmoothed plots) and whether alternatives are produced
void buildChain(TChain* bkgMC, int sampleIndex);                  // Generates TChain from input MC
void makeTemplate(int updown,bool debug);                         // Makes template root files - NOTE: ggH, ggZZ, qqZZ, and Z+X are listed explicitly and do not rely on Config.h
TH2F* fillTemplate(int sampleIndex, bool isLowMass, int updown);  // Takes input MC to make Fisher v m4L TH2F
TH2F* mergeTemplates(TH2F* lowTemp, TH2F* highTemp);              // Merges low and high mass TH2F plots
TH2F* smoothtemplates(TH2F* inputdata, int sampleIndex);          // Smooths Fisher v m4L plots
TH2F* rebin(TH2F* rebinnedhist, int sampleIndex);                 // Smoothing procedure for ggH,qqH
TH2F* rebin_lowstatistics(TH2F* rebinnedhist, int sampleIndex);   // Smoothing procedure for ggZZ,qqZZ,ZH,WH
void analyticfits(int sampleIndex, int updown);                   // Smoothing procedure for Z+X,ttH
TH2F* altshapes(TH2F* originalHist, int channel, int altnum);     // Generates alternative shapes for Fisher, (channel,altnum): (0,1)=ggHAMCNLO, (0,2)=ggHMG, (1,1)=qqHd6t, (2,1)=qqZZMG
float altscale(float Fisher, int channel, int altnum);            // Find scale for alternative shapes given a Fisher value, channel and altnum same as altshapes()
TH2F* mirrortemplates(int sampleIndex);                           // Generates mirror alternative shapes for Fisher when no second alternative exists
bool test_bit(int mask, unsigned int iBit);                       // Used to identify correct CR events

//---------------------------------------------------

bool test_bit( int mask, unsigned int iBit ) { return (mask >> iBit) & 1; }

void generateFisherTemplates() {
  bool debug=false;
  bool findAlternatives=true;
  templateOptions(debug,findAlternatives);
}
 
void templateOptions(bool debug, bool findAlternatives){
  makeTemplate(0,debug);
  if (findAlternatives){
    makeTemplate(1,debug);
    makeTemplate(-1,debug);
    makeTemplate(2,debug);
    makeTemplate(3,debug);
  }
}

//---------------------------------------------------

void buildChain(TChain* bkgMC, int sampleIndex){
  TString filePath;
  int nPoints=0;
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
  if (sampleIndex!=0 && sampleIndex!=2 && sampleIndex!=3 && sampleIndex!=4 && sampleIndex!=-2){
    if(nPoints==0){
      cout<<"nPoints not set in Config.h"<<endl;
      abort();
    }
    for (int i=0; i<nPoints; i++){
      char tmp_finalInPath4mu[200],tmp_finalInPath4e[200],tmp_finalInPath2mu2e[200];
      string finalInPath4mu,finalInPath4e,finalInPath2mu2e;
      if (sampleIndex==1){
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
  else if (sampleIndex==0){
    if(useSqrts<2){
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_powheg15H1000.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_powheg15jhuGenV3H115.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_powheg15jhuGenV3H120.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_powheg15jhuGenV3H122.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_powheg15jhuGenV3H124.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_powheg15jhuGenV3H125.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_powheg15jhuGenV3H126.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_powheg15jhuGenV3H128.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_H130.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_powheg15jhuGenV3H140.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_H150.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_powheg15jhuGenV3H160.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_powheg15jhuGenV3H170.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_powheg15jhuGenV3H180.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_powheg15jhuGenV3H185.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_H190.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_powheg15jhuGenV3H200.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_H210.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_H220.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_H250.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_H275.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_powheg15H300.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_H325.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_powheg15H350.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_powheg15H400.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_H425.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_powheg15H450.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_H475.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_powheg15H500.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_H525.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_powheg15H550.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_H575.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_powheg15H600.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_powheg15H650.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_powheg15H700.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_H750.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_powheg15H800.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_powheg15H900.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_H950.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_powheg15H1000.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H115.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H120.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H122.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H124.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H125.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H126.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H128.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_H130.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H140.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_H150.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H160.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H170.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H180.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H185.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_H190.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H200.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_H210.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_H220.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_H250.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_H275.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_powheg15H300.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_H325.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_powheg15H350.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_powheg15H400.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_H425.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_powheg15H450.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_H475.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_powheg15H500.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_H525.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_powheg15H550.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_H575.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_powheg15H600.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_powheg15H650.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_powheg15H700.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_H750.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_powheg15H800.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_powheg15H900.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_H950.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_powheg15H1000.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H115.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H120.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H122.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H124.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H125.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H126.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H128.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_H130.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H140.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_H150.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H160.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H170.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H180.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H185.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_H190.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H200.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_H210.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_H220.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_H250.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_H275.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_powheg15H300.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_H325.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_powheg15H350.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_powheg15H400.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_H425.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_powheg15H450.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_H475.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_powheg15H500.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_H525.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_powheg15H550.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_H575.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_powheg15H600.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_powheg15H650.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_powheg15H700.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_H750.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_powheg15H800.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_powheg15H900.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_H950.root");
    }
    if(useSqrts%2==0){
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15H1000.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H115.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_H116.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_H117.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_H118.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_H119.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H120.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_H121.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H122.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_H123.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H124.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H125.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H126.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_H127.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H128.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_H129.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H130.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H135.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H140.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H145.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H150.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H160.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H175.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H170.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H180.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H185.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H190.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15jhuGenV3H200.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_H220.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15H225.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15H250.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15H275.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15H300.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_H325.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15H350.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_H375.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15H400.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_H425.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15H450.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_H475.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15H500.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_H525.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15H550.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_H575.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15H600.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15H650.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15H700.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_H750.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15H800.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_H850.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_powheg15H900.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_H950.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15H1000.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15jhuGenV3H115.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_H116.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_H117.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_H118.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_H119.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15jhuGenV3H120.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_H121.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15jhuGenV3H122.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_H123.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15jhuGenV3H124.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15jhuGenV3H125.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15jhuGenV3H126.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_H127.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15jhuGenV3H128.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_H129.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15jhuGenV3H130.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15jhuGenV3H135.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15jhuGenV3H140.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15jhuGenV3H145.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15jhuGenV3H150.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15jhuGenV3H160.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15jhuGenV3H175.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15jhuGenV3H170.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15jhuGenV3H180.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15jhuGenV3H185.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15jhuGenV3H190.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15jhuGenV3H200.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_H220.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15H225.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15H250.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15H275.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15H300.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_H325.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15H350.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_H375.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15H400.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_H425.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15H450.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_H475.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15H500.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_H525.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15H550.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_H575.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15H600.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15H650.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15H700.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_H750.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15H800.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_H850.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_powheg15H900.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_H950.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15H1000.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H115.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_H116.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_H117.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_H118.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_H119.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H120.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_H121.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H122.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_H123.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H124.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H125.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H126.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_H127.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H128.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_H129.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H130.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H135.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H140.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H145.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H150.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H160.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H175.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H170.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H180.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H185.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H190.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15jhuGenV3H200.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_H220.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15H225.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15H250.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15H275.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15H300.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_H325.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15H350.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_H375.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15H400.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_H425.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15H450.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_H475.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15H500.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_H525.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15H550.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_H575.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15H600.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15H650.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15H700.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_H750.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15H800.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_H850.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_powheg15H900.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_H950.root");
    }
  }
  else if (sampleIndex==2){
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
  else if (sampleIndex==3){
    if(useSqrts<2){
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_ggZZ2l2l.root");
      bkgMC->Add(filePath7TeV + "4mu/HZZ4lTree_ggZZ4l.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_ggZZ2l2l.root");
      bkgMC->Add(filePath7TeV + "4e/HZZ4lTree_ggZZ4l.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_ggZZ2l2l.root");
      bkgMC->Add(filePath7TeV + "2mu2e/HZZ4lTree_ggZZ4l.root");
    }
    if(useSqrts%2==0){
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_ggZZ2l2l.root");
      bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_ggZZ4l.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_ggZZ2l2l.root");
      bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_ggZZ4l.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_ggZZ2l2l.root");
      bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_ggZZ4l.root");
    }
  }
  else if (sampleIndex==4){
    if(useSqrts<2){
      bkgMC->Add(filePath7TeV + "CR/HZZ4lTree_DoubleOr_CRZLLTree.root");
    }
    if(useSqrts%2==0){
      bkgMC->Add(filePath8TeV + "CR/HZZ4lTree_DoubleOr_CRZLLTree.root");
    }
  }
  else if (sampleIndex==-2){
    bkgMC->Add(filePath8TeV + "4mu/HZZ4lTree_ZZJetsTo4L.root");
    bkgMC->Add(filePath8TeV + "4e/HZZ4lTree_ZZJetsTo4L.root");
    bkgMC->Add(filePath8TeV + "2mu2e/HZZ4lTree_ZZJetsTo4L.root");
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
  else if(updown==2){
    jes="_alt";
  }
  else if(updown==3){
    jes="_alt2";
  }

  TString debugname;
  if (debug) debugname="_unnormalized";

  if(updown==0 || updown==1 || updown==-1){
    fqqH = new TFile(destDir + "qqH_fisher"+jes+debugname+".root","RECREATE");
    fggH = new TFile(destDir + "ggH_fisher"+jes+debugname+".root","RECREATE");
    fqqZZ = new TFile(destDir + "qqZZ_fisher"+jes+debugname+".root","RECREATE");
    fggZZ = new TFile(destDir + "ggZZ_fisher"+jes+debugname+".root","RECREATE");
    fZX = new TFile(destDir + "Z+X_fisher"+jes+debugname+".root","RECREATE");
    fZH = new TFile(destDir + "ZH_fisher"+jes+debugname+".root","RECREATE");
    fWH = new TFile(destDir + "WH_fisher"+jes+debugname+".root","RECREATE");
    fttH = new TFile(destDir + "ttH_fisher"+jes+debugname+".root","RECREATE");
  } else{
    fqqH = new TFile(destDir + "qqH_fisher"+jes+debugname+".root","RECREATE");
    fggH = new TFile(destDir + "ggH_fisher"+jes+debugname+".root","RECREATE");
    fqqZZ = new TFile(destDir + "qqZZ_fisher"+jes+debugname+".root","RECREATE");
  }

  TH2F* low,*high,*H_Fisher;
  
  
  // =========================
  // ggH
  
  low = fillTemplate(0,true,updown);
  high = fillTemplate(0,false,updown);
  H_Fisher = mergeTemplates(low,high);

  if (!debug && updown<2) smoothtemplates(H_Fisher,0);
  if (!debug && updown==2) smoothtemplates(H_Fisher,-1);
  if (!debug && updown==3) smoothtemplates(H_Fisher,-2);

  fggH->cd();
  H_Fisher->Write("H_Fisher");
  fggH->Close();
  
  // ==========================
  // qqH
  
  if(updown<3){
    low = fillTemplate(1,true,updown);
    high = fillTemplate(1,false,updown);
    H_Fisher = mergeTemplates(low,high);
  }  

  if (!debug && updown<2) smoothtemplates(H_Fisher,1);
  if (!debug && updown==2) smoothtemplates(H_Fisher,-3);
  if (!debug && updown==3) H_Fisher = mirrortemplates(1);

  fqqH->cd();
  H_Fisher->Write("H_Fisher");
  fqqH->Close();
  
  // ==========================
  // qqZZ

  if(updown<2){
    low = fillTemplate(2,true,updown);
    high = fillTemplate(2,false,updown);
  }
  if(updown==2){
    low = fillTemplate(-2,true,0);
    high = fillTemplate(-2,false,0);
  }
  if (updown<3) H_Fisher = mergeTemplates(low,high);

  if (!debug && updown<3) smoothtemplates(H_Fisher,2); //Alter updown to be <3 when 
  if (!debug && updown==3) H_Fisher = mirrortemplates(2);

  fqqZZ->cd();
  H_Fisher->Write("H_Fisher");
  fqqZZ->Close();
  
  if (updown!=2 && updown!=3){
    
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
}

//---------------------------------------------------

TH2F* fillTemplate(int sampleIndex,bool isLowMass,int updown){
  TChain* bkgMC = new TChain("SelectedTree");
  buildChain(bkgMC,sampleIndex);

  cout << "Chain for " << sampleIndex << " " << isLowMass << " " << updown << " " << bkgMC->GetEntries() << endl;
  bkgMC->ls();

  float mass,deta,mJJ,w,Fisher;
  int njets;
  int processID;
  string channel;

  bkgMC->SetBranchAddress("ZZMass",&mass);
  bkgMC->SetBranchAddress("NJets30",&njets);
  bkgMC->SetBranchAddress("DiJetDEta",&deta);
  if(updown==0 || updown==2 || updown==3){
    bkgMC->SetBranchAddress("DiJetMass",&mJJ);
  }else if(updown==1){
    bkgMC->SetBranchAddress("DiJetMassPlus",&mJJ);
  }else if(updown==-1){
    bkgMC->SetBranchAddress("DiJetMassMinus",&mJJ);
  }
  bkgMC->SetBranchAddress("Fisher",&Fisher);
  bkgMC->SetBranchAddress("MC_weight",&w);
  bkgMC->SetBranchAddress("genProcessId",&processID);

  TH2F* bkgHist;
  if(!isLowMass){
    bkgHist = new TH2F("bkgHisto","bkgHisto",int((highMzz-180.)/mBinSize+0.5),180,highMzz,50,0,2);
  }
  else{
    bkgHist = new TH2F("bkgHisto","bkgHisto",int((180-100)/mBinSize+0.5),100,180,50,0,2);
  }

  //bkgHist->Sumw2();

  int percent=0;

  //Fill histogram
  for(int i=0; i<bkgMC->GetEntries(); i++){
    bkgMC->GetEntry(i);
    if (i%50000==0){
      cout << "event: " << i << "/" << bkgMC->GetEntries() << endl;
    }
    if (mass<100 || deta<=-99 || mJJ<=-99 || njets<2) continue;
    if(updown!=0) Fisher = 0.18*fabs(deta) + 0.000192*mJJ;
    if (Fisher>=2.){
      Fisher = 1.98;
      percent++;
    }
    bkgHist->Fill(mass,Fisher,w);
  }
  cout<<percent<<" events found with Fisher>2 or "<<float(percent)/float(bkgMC->GetEntries())<<endl;

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

  //H_Fisher->Sumw2();

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
  if(sampleIndex==0 || sampleIndex==1 || sampleIndex==-1 || sampleIndex==-2 || sampleIndex==-3){
    rebin(inputdata,sampleIndex);
  }
  else if(sampleIndex==2 || sampleIndex==3 || sampleIndex==5 || sampleIndex==6){
    rebin_lowstatistics(inputdata,sampleIndex);
  }

  return inputdata;
}

//---------------------------------------------------

TH2F* rebin(TH2F* rebinnedHist, int usealt){

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
  //rebinnedHist->Sumw2();

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
  //rebinnedHist->Sumw2();
  
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
  
  //rebinnedHist->Sumw2();

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

  //rebinnedHist->Sumw2();

  if (usealt==-1) altshapes(rebinnedHist,0,1);
  if (usealt==-2) altshapes(rebinnedHist,0,2);
  if (usealt==-3) altshapes(rebinnedHist,1,1);

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

  //finalhist->Sumw2();

  TH2F* origHist = new TH2F (*finalhist);

  //origHist->Sumw2();

  //Project Fisher v m4L to 1D Fisher plots for low, high, and full mass range
  TH1F* lowProj = (TH1F*) origHist->ProjectionY("lowProj",1,40);
  TH1F* highProj = (TH1F*) origHist->ProjectionY("highProj",41,750);
  TH1F* fullProj = (TH1F*) origHist->ProjectionY();
  TH1F* origlowProj = new TH1F (*lowProj);
  TH1F* orighighProj = new TH1F (*highProj);
  TH1F* origfullProj = new TH1F (*fullProj);
  double temp;

  /*lowProj->Sumw2();
  highProj->Sumw2();
  fullProj->Sumw2();
  origlowProj->Sumw2();
  orighighProj->Sumw2();
  origfullProj->Sumw2();*/

  //qqZZ + WH + ZH
  if(sampleIndex==2 || sampleIndex==5 || sampleIndex==6){
    //Fit tails of low and high projections with exponentials to account for low statistics
    TF1 *tailfunc = new TF1("tailfunc","[0]*exp([1]*(x))",0.2,2.0);
    if (sampleIndex==2){
      lowProj->Fit(tailfunc,"","",0.5,1.1);
    }
    else if(sampleIndex==5){
      lowProj->Fit(tailfunc,"","",0.5,1.3);
    }
    else if(sampleIndex==6){
      lowProj->Fit(tailfunc,"","",0.5,1.3);
    }
    TF1 *tailfunc2 = new TF1("tailfunc2","[0]*exp([1]*(x))",0.2,2.0);
    if (sampleIndex==2){
      highProj->Fit(tailfunc2,"","",0.9,1.6);
    }
    else if(sampleIndex==5){
      highProj->Fit(tailfunc2,"","",0.6,1.3);
    }
    else if(sampleIndex==6){
      highProj->Fit(tailfunc2,"","",0.5,1.1);
    }
    origlowProj->Draw();
    tailfunc->Draw("same");

    //Adjust projections with fitted exponentials
    for(int i=0;i<=nYbins;i++){
      float binFisher = lowProj->GetBinCenter(i);
      if ((sampleIndex==2 && binFisher>1.0) || (sampleIndex==5 && binFisher>1.1) || (sampleIndex==6 && binFisher>1.3)){
	lowProj->SetBinContent(i,tailfunc->Eval(binFisher));
      }
    }
    for(int i=0;i<=nYbins;i++){
      float binFisher = highProj->GetBinCenter(i);
      if ((sampleIndex==2 && binFisher>1.4) || (sampleIndex==5 && binFisher>1.3) || (sampleIndex==6 && binFisher>1.1)){
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
    else if(sampleIndex==5){
      fZH->cd();
    }
    else if(sampleIndex==6){
      fWH->cd();
    }
    origlowProj->Write("H_lowraw");
    lowProj->Write("H_lowfit");
    orighighProj->Write("H_highraw");
    highProj->Write("H_highfit");
    
  }
  //ggZZ
  else if(sampleIndex==3){
    //Fit tail of full projection with exponential to account for low statistics
    TF1 *tailfunc3 = new TF1("tailfunc3","[0]*exp([1]*(x))",0.3,2.0);
    if(sampleIndex==3){
      fullProj->Fit(tailfunc3,"","",1.0,2.0);
    }

    //Adjust projection with fitted exponential
    for(int i=0;i<=nYbins;i++){
      float binFisher = fullProj->GetBinCenter(i);
      if ((sampleIndex==3 && binFisher>1.7)){
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
  TChain* bkgMC = new TChain("SelectedTree");
  buildChain(bkgMC,sampleIndex);

  cout << "Chain for " << sampleIndex << " " << bkgMC->GetEntries() << endl;
  bkgMC->ls();

  float deta,mJJ,w,mass,tempFisher;
  int processID;
  int njets;
  int CRflag;

  bkgMC->SetBranchAddress("ZZMass",&mass);
  bkgMC->SetBranchAddress("NJets30",&njets);
  bkgMC->SetBranchAddress("DiJetDEta",&deta);
  if(sampleIndex==7 && updown==0){
    bkgMC->SetBranchAddress("DiJetMass",&mJJ);
  }else if(sampleIndex==7 && updown==1){
    bkgMC->SetBranchAddress("DiJetMassPlus",&mJJ);
  }else if(sampleIndex==7 && updown==-1){
    bkgMC->SetBranchAddress("DiJetMassMinus",&mJJ);
  }
  bkgMC->SetBranchAddress("MC_weight",&w);
  bkgMC->SetBranchAddress("genProcessId",&processID);
  if(sampleIndex==7) bkgMC->SetBranchAddress("Fisher",&tempFisher);
  if(sampleIndex==4){
    bkgMC->SetBranchAddress("CRflag",&CRflag);
    bkgMC->SetBranchAddress("ZZFisher",&tempFisher);
  }

  RooRealVar Fisher("Fisher","D_{Jet}",0.,2.);
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
    if((sampleIndex==4 && (test_bit(CRflag,5) || test_bit(CRflag,7) || test_bit(CRflag,9) || test_bit(CRflag,11))) || sampleIndex==7){
      if((sampleIndex==7 && mass>100. && njets>=2) || (sampleIndex==4 && tempFisher!=-99. && mass>100.)){
	if(updown!=0 && sampleIndex==7) tempFisher = 0.18*fabs(deta) + 0.000192*mJJ;
	if(tempFisher>=2.) tempFisher=1.98;
	Fisher=tempFisher;
	data->add(RooArgSet(Fisher),w);
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

  //UNCOMMENT TO MAKE PLOTS
  /*
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(1);
  TCanvas* c = new TCanvas("c","c",1300,800);
  c->cd();
  finalhist->Draw("COLZ");
  if(sampleIndex==4 && updown==0) c->SaveAs("Z+X_Fisher_2D.png");
  if(sampleIndex==7 && updown==0) c->SaveAs("ttH_Fisher_2D.png");
  Fisherframe->Draw();
  if(sampleIndex==4 && updown==0) c->SaveAs("Z+X_Fisher_1D.png");
  if(sampleIndex==7 && updown==0) c->SaveAs("ttH_Fisher_1D.png");
  */

  if (sampleIndex==4){
    fZX->cd();
  }
  else if (sampleIndex==7){
    fttH->cd();
  }
  finalhist->Write("H_Fisher");
  Fisherframe->Write("H_fit");
  if (sampleIndex==4){
    fZX->Close();
  }
  else if (sampleIndex==7){
    fttH->Close();
  }
}

TH2F* altshapes(TH2F* original, int channel, int altnum){	 	 
  int nXbins=original->GetNbinsX();
  int nYbins=original->GetNbinsY();
  float binFisher,scale,altvalue;
	 	 
  for (int i=1; i <=nXbins; i++){
    for (int j=1; j<=nYbins; j++){
      binFisher = original->GetYaxis()->GetBinCenter(j);
      scale = altscale(binFisher,channel,altnum);
      altvalue=scale*(original->GetBinContent(i,j));
      original->SetBinContent(i,j,altvalue);
    }
  }
	 	 
  return original;	 	
}
	 	 
float altscale(float Fisher,int channel, int altnum){
  float scale=0;
  
  //ggH
  if(channel==0){
    if(altnum==1){
      if (Fisher<=0.2) scale=0.860963;
      if (Fisher>0.2 && Fisher<=0.4) scale=1.02564;
      if (Fisher>0.4 && Fisher<=0.6) scale=1.04104;
      if (Fisher>0.6 && Fisher<=0.8) scale=1.06351;
      if (Fisher>0.8 && Fisher<=1.0) scale=1.1585;
      if (Fisher>1.0 && Fisher<=1.2) scale=1.31685;
      if (Fisher>1.2 && Fisher<=1.4) scale=1.34064;
      if (Fisher>1.4 && Fisher<=1.6) scale=1.60101;
      if (Fisher>1.6 && Fisher<=1.8) scale=1.84312;
      if (Fisher>1.8) scale=2.28459;
    }
    if(altnum==2){
      if (Fisher<=0.2) scale=1.05071;
      if (Fisher>0.2 && Fisher<=0.4) scale=1.08807;
      if (Fisher>0.4 && Fisher<=0.6) scale=0.988566;
      if (Fisher>0.6 && Fisher<=0.8) scale=0.883993;
      if (Fisher>0.8 && Fisher<=1.0) scale=0.777513;
      if (Fisher>1.0 && Fisher<=1.2) scale=0.727857;
      if (Fisher>1.2 && Fisher<=1.4) scale=0.687078;
      if (Fisher>1.4 && Fisher<=1.6) scale=0.630311;
      if (Fisher>1.6 && Fisher<=1.8) scale=0.827337;
      if (Fisher>1.8) scale=1.04622;
    }
  }
  //qqH
  if(channel==1){
    if(altnum==1){
      if (Fisher<=0.2) scale=1.18197;
      if (Fisher>0.2 && Fisher<=0.4) scale=1.17726;
      if (Fisher>0.4 && Fisher<=0.6) scale=1.00603;
      if (Fisher>0.6 && Fisher<=0.8) scale=1.0123;
      if (Fisher>0.8 && Fisher<=1.0) scale=0.985812;
      if (Fisher>1.0 && Fisher<=1.2) scale=0.906812;
      if (Fisher>1.2 && Fisher<=1.4) scale=0.968142;
      if (Fisher>1.4 && Fisher<=1.6) scale=0.820917;
      if (Fisher>1.6 && Fisher<=1.8) scale=0.914569;
      if (Fisher>1.8) scale=0.838355;
    }
  }
  //qqZZ NOT USED, reading from separate files instead
  /*if(channel==2){
    if(altnum==1){
    if (Fisher<=0.2) scale=1.04573;
    if (Fisher>0.2 && Fisher<=0.4) scale=0.988066;
    if (Fisher>0.4 && Fisher<=0.6) scale=0.953113;
    if (Fisher>0.6 && Fisher<=0.8) scale=0.742393;
    if (Fisher>0.8 && Fisher<=1.0) scale=0.867237;
    if (Fisher>1.0 && Fisher<=1.2) scale=0.554462;
    if (Fisher>1.2 && Fisher<=1.4) scale=0.614161;
    if (Fisher>1.4 && Fisher<=1.6) scale=1.24295;
    if (Fisher>1.6 && Fisher<=1.8) scale=1.49153;
    if (Fisher>1.8) scale=0.372884;
    }
    }*/
	 	 
  if(scale==0) cout<<"ERROR: Fisher Template has values <0 or >2."<<endl;
  
  return scale;
  
}
	 	 
TH2F* mirrortemplates(int sampleIndex){
  TFile* alt1,*orig;
  if (sampleIndex==1){
    alt1 = new TFile(destDir + "qqH_fisher_alt.root","OPEN");
    orig = new TFile(destDir + "qqH_fisher.root","OPEN");
  }
  if (sampleIndex==2){
    alt1 = new TFile(destDir + "qqZZ_fisher_alt.root","OPEN");
    orig = new TFile(destDir + "qqZZ_fisher.root","OPEN");
  }
  TH2F* altfisher = (TH2F*)alt1->Get("H_Fisher");
  TH2F* origfisher = (TH2F*)orig->Get("H_Fisher");
  TH2F* alt2fisher = new TH2F("H_Fisher","H_Fisher",int((highMzz-100.)/mBinSize+0.5),100,highMzz,50,0.,2.);
  
  for(int i=1;i<=origfisher->GetNbinsX();i++){
    for(int j=1;j<=origfisher->GetNbinsY();j++){
      float origval=origfisher->GetBinContent(i,j);
      float alt1val=altfisher->GetBinContent(i,j);
      float scale=origval/alt1val;
      float alt2val=scale*origval;
      alt2fisher->SetBinContent(i,j,alt2val);
    }
  }
  
  //Renormalize
  TH1F* tempProj;
  double norm;
  for(int i=1; i<=alt2fisher->GetNbinsX(); i++){
    tempProj = (TH1F*) alt2fisher->ProjectionY("tempProj",i,i);
    norm=tempProj->Integral();
    if (norm>0) { // Avoid introducing NaNs in the histogram
      for(int j=1; j<=alt2fisher->GetNbinsY(); j++){
	alt2fisher->SetBinContent(i,j,alt2fisher->GetBinContent(i,j)/norm);
      }
    }
  }
	 	 
  return alt2fisher;
}

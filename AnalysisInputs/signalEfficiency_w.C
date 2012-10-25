/* 
 * Compute efficiencies for signals and write them in a card fragment.
 * usage: 
 * -set all input variables in Config.h
 * -run with:
 * root -q -b signalEfficiency_w.C+
 * This runs on the 3 final states for 7 and 8 TeV and writes the output in a file (see stdout).
 *
 */

#include "TH1F.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include "TCutG.h"
#include "TFile.h"
#include "TH2.h"
#include "TPad.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TSystem.h"


#include "../CreateDatacards/include/tdrstyle.cc"


//using namespace RooFit;
using namespace std;
using namespace ROOT::Math;

//----------> SET INPUT VARIABLES in Config.h
#include "Config.h"
//<----------


void signalEfficiency_w(int channel, double sqrts);


// Run all final states and sqrts in one go
void signalEfficiency_w() {
  gSystem->Exec("mkdir -p sigFigs7TeV");
  gSystem->Exec("mkdir -p sigFigs8TeV");

  signalEfficiency_w(1,7);
  signalEfficiency_w(2,7);
  signalEfficiency_w(3,7);
  signalEfficiency_w(1,8);
  signalEfficiency_w(2,8);
  signalEfficiency_w(3,8);
}


// The actual job
void signalEfficiency_w(int channel, double sqrts) 
{
  TString schannel;
  if      (channel == 1) schannel = "4mu";
  else if (channel == 2) schannel = "4e";
  else if (channel == 3) schannel = "2e2mu";
  else cout << "Not a valid channel: " << schannel << endl;

  TString ssqrts = (long) sqrts + TString("TeV");

  cout << "schannel = " << schannel << "  sqrts = " << sqrts << endl;

  TString outfile = "CardFragments/signalEfficiency_" + ssqrts + "_" + schannel + ".txt";
  ofstream of(outfile,ios_base::out);

  gSystem->AddIncludePath("-I$ROOFITSYS/include");
  setTDRStyle(false);
  gStyle->SetStatX(-0.5);

  int nPoints=0;
  int* masses=0;
  double* mHVal=0;

  // Pick the correct set of mass points, set subpath
  TString filepath;
  if (sqrts==7) {
    nPoints = nPoints7TeV;
    masses  = masses7TeV;
    mHVal   = mHVal7TeV;
    filepath = filePath7TeV;
  } else if (sqrts==8) {
    nPoints = nPoints8TeV;
    masses  = masses8TeV;
    mHVal   = mHVal8TeV;
    filepath =filePath8TeV;
  }
  

  float xMax = masses[nPoints-1];

  cout << "xmax: " << xMax << endl;


  const int arraySize=200;
  assert (arraySize>=nPoints);
  double efficiencyVal[arraySize];
  double efficiencyErr[arraySize];
  // R e a d   w o r k s p a c e   f r o m   f i l e
  // -----------------------------------------------
	
//   double a_meanBW[arraySize];
//   double a_gammaBW[arraySize];
//   double a_meanCB[arraySize];
//   double a_sigmaCB[arraySize];
//   double a_alphaCB[arraySize];
//   double a_nCB[arraySize];

  //  HiggsCSandWidth *myCSW = new HiggsCSandWidth();
	
  for (int i = 0; i < nPoints; i++){
    TString infile = filepath+ "/" + (schannel=="2e2mu"?"2mu2e":schannel) + "/HZZ4lTree_H" + (long)masses[i] + ".root";
    //    cout << "Opening input file: " << infile << endl;
    TFile *f = TFile::Open(infile) ;
    TTree *t1 = (TTree*) f->Get("SelectedTree");
    float mela, MC_weight_norm, MC_weight_PUWeight, MC_weight_powhegWeight,  MC_weight_dataMC;
    int numEventsRaw = 0;
    float numEventsPowheg =0;
    float numEventsPU =0;    
    float numEventsDataMC =0;    
    t1->SetBranchAddress("MC_weight_norm",&MC_weight_norm);
    t1->SetBranchAddress("MC_weight_powhegWeight",&MC_weight_powhegWeight);
    t1->SetBranchAddress("MC_weight_PUWeight",&MC_weight_PUWeight);
    t1->SetBranchAddress("MC_weight_dataMC",&MC_weight_dataMC);
    t1->SetBranchAddress("ZZLD",&mela);
    float totalCtr=0;
    float eff_noweight=0;
    float sumw2=0;
    float den = 0;
    for (int a = 0; a < t1->GetEntries(); a++){ 
      t1->GetEntry(a);
      totalCtr+=MC_weight_norm; 
      sumw2 += MC_weight_norm*MC_weight_norm;
      ++numEventsRaw;      
      numEventsPowheg += MC_weight_powhegWeight;
      numEventsPU += MC_weight_powhegWeight*MC_weight_PUWeight;
      numEventsDataMC += MC_weight_powhegWeight*MC_weight_PUWeight*MC_weight_dataMC;
      if (MC_weight_powhegWeight>0) {
	eff_noweight +=MC_weight_norm/(MC_weight_powhegWeight*MC_weight_PUWeight*MC_weight_dataMC);
	den = (MC_weight_powhegWeight*MC_weight_PUWeight*MC_weight_dataMC)/MC_weight_norm;
      }
    }

    // FIXME: the 7TeV samples are assumed to have the ad-hoc correction factor for the mll>12 gen cut,
    // except for the 124,125,126 new samples. As this factor is accounted for in the x-section, we have to 
    // apply it here
    float m = masses[i];
    if (sqrts==7 && m>=123.9 &&  m<=126.1) {
      float mllCorr = 0.5 + 0.5*erf((m-80.85)/50.42);
      totalCtr = totalCtr/mllCorr;
      eff_noweight=eff_noweight/mllCorr;
    }

    cout << " m = " << masses[i] 
	 << " : selected events= " << numEventsRaw 
	 << " Powheg Wtd= " << numEventsPowheg
	 << " PU Wtd= " << numEventsPU
	 << " Data/MC Wtd= " << numEventsDataMC
	 << " efficiency= " << totalCtr 
//	 << " " << eff_noweight
	 << endl;

    efficiencyVal[i] = totalCtr; 
    efficiencyErr[i] = sqrt(sumw2);
//    efficiencyVal[i] = eff_noweight;
    f->Close();
  }
	
  //  TGraph* grEff = new TGraph( nPoints, mHVal, efficiencyVal );
  TGraphErrors* grEff = new TGraphErrors( nPoints, mHVal, efficiencyVal, 0, efficiencyErr);
  grEff->SetMarkerStyle(20);

  //ICHEP parametrization	
//   TF1 *polyFunc= new TF1("polyFunc","([0]+[1]*TMath::Erf( (x-[2])/[3] ))*([4]+[5]*x+[6]*x*x)", 110., xMax);
//   polyFunc->SetParameters(-4.42749e+00,4.61212e+0,-6.21611e+01,1.13168e+02,2.14321e+00,1.04083e-03,4.89570e-07);

  TF1 *polyFunc= new TF1("polyFunc","([0]+[1]*TMath::Erf( (x-[2])/[3] ))*([4]+[5]*x+[6]*x*x)+[7]*TMath::Gaus(x,[8],[9])", 110., xMax);
  polyFunc->SetParameters(-4.42749e+00,4.61212e+0,-6.21611e+01,1.13168e+02,2.14321e+00,1.04083e-03,4.89570e-07,
			  0.03, 200, 30);
  polyFunc->SetParLimits(7,0,0.1);
  polyFunc->SetParLimits(8,160,210);
  polyFunc->SetParLimits(9,20,50);
  

  //  TF1 *polyFunc= new TF1("polyFunc","pol9", 110., xMax);

  polyFunc->SetLineColor(4);      
  TCanvas *c = new TCanvas("c","c");
  c->SetGrid();

  TString outname = "sigFigs" + ssqrts +"/eff_" + schannel;

  grEff->Fit(polyFunc,"Rt");
  TString xaxisText = "m_{" + schannel + "}";
  grEff->GetXaxis()->SetTitle(xaxisText);
  TString yaxisText = "efficiency, " + schannel;
  grEff->GetYaxis()->SetTitle(yaxisText);
  grEff->Draw("AP");
  polyFunc->Draw("sames");
  c->Print(outname+".eps");
  c->Print(outname+".png"); // Does not work in batch?
  c->Print(outname+".root"); 

  cout << endl;
  cout << "------- Parameters for " << schannel << " sqrts=" << sqrts << endl;
  cout << "   a1 = " << polyFunc->GetParameter(0) << endl;
  cout << "   a2 = " << polyFunc->GetParameter(1) << endl;
  cout << "   a3 = " << polyFunc->GetParameter(2) << endl;
  cout << "   a4 = " << polyFunc->GetParameter(3) << endl;
  cout << "   b1 = " << polyFunc->GetParameter(4) << endl;
  cout << "   b2 = " << polyFunc->GetParameter(5) << endl;
  cout << "   b3 = " << polyFunc->GetParameter(6) << endl;
  cout << "---------------------------" << endl << endl;

  of << "## signal efficiency ##" << endl;
  of << "signalEff a1  " << polyFunc->GetParameter(0) << endl;
  of << "signalEff a2  " << polyFunc->GetParameter(1) << endl;
  of << "signalEff a3  " << polyFunc->GetParameter(2) << endl;
  of << "signalEff a4  " << polyFunc->GetParameter(3) << endl;
  of << "signalEff b1  " << polyFunc->GetParameter(4) << endl;
  of << "signalEff b2  " << polyFunc->GetParameter(5) << endl;
  of << "signalEff b3  " << polyFunc->GetParameter(6) << endl;
  of << endl << endl;
  of.close();
  
  cout << "Output written to: " <<  outfile << endl << endl;

  // deviations
  cout << "Deviations..." << endl;
  double maxResidual=0;
  for (int i = 0; i < nPoints; i++){
    double eval = polyFunc->Eval(masses[i]);
    double residual = (eval - efficiencyVal[i]);
      maxResidual = max(maxResidual,fabs(residual));
    cout << "For mass, " << masses[i] << ": measured value is " << efficiencyVal[i] << " and difference from function is " << residual <<endl;
  }
  cout << "Largest residual= " << maxResidual << endl;
}


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

bool verbose = true;


//using namespace RooFit;
using namespace std;
using namespace ROOT::Math;

//----------> SET INPUT VARIABLES in Config.h
#include "Config.h"
//<----------


void signalEfficiency_w(int channel, double sqrts, int process, bool VBFtag = false);


// Run all final states and sqrts in one go
void signalEfficiency_w() {
  gSystem->Exec("mkdir -p sigFigs7TeV");
  gSystem->Exec("mkdir -p sigFigs8TeV");

  //ggH
  signalEfficiency_w(1,7,1,true);
  signalEfficiency_w(2,7,1,true);
  signalEfficiency_w(3,7,1,true);
  signalEfficiency_w(1,8,1,true);
  signalEfficiency_w(2,8,1,true);
  signalEfficiency_w(3,8,1,true);
  //qqH
  signalEfficiency_w(1,7,2,true);
  signalEfficiency_w(2,7,2,true);
  signalEfficiency_w(3,7,2,true);
  signalEfficiency_w(1,8,2,true);
  signalEfficiency_w(2,8,2,true);
  signalEfficiency_w(3,8,2,true);
  //ZH
  signalEfficiency_w(1,7,3,true);
  signalEfficiency_w(2,7,3,true);
  signalEfficiency_w(3,7,3,true);
  signalEfficiency_w(1,8,3,true);
  signalEfficiency_w(2,8,3,true);
  signalEfficiency_w(3,8,3,true);
  //WH
  signalEfficiency_w(1,7,4,true);
  signalEfficiency_w(2,7,4,true);
  signalEfficiency_w(3,7,4,true);
  signalEfficiency_w(1,8,4,true);
  signalEfficiency_w(2,8,4,true);
  signalEfficiency_w(3,8,4,true);
  //ttH
  signalEfficiency_w(1,7,5,true);
  signalEfficiency_w(2,7,5,true);
  signalEfficiency_w(3,7,5,true);
  signalEfficiency_w(1,8,5,true);
  signalEfficiency_w(2,8,5,true);
  signalEfficiency_w(3,8,5,true);

  //ggH
  signalEfficiency_w(1,7,1,false);
  signalEfficiency_w(2,7,1,false);
  signalEfficiency_w(3,7,1,false);
  signalEfficiency_w(1,8,1,false);
  signalEfficiency_w(2,8,1,false);
  signalEfficiency_w(3,8,1,false);
  //qqH
  signalEfficiency_w(1,7,2,false);
  signalEfficiency_w(2,7,2,false);
  signalEfficiency_w(3,7,2,false);
  signalEfficiency_w(1,8,2,false);
  signalEfficiency_w(2,8,2,false);
  signalEfficiency_w(3,8,2,false);
  //ZH
  signalEfficiency_w(1,7,3,false);
  signalEfficiency_w(2,7,3,false);
  signalEfficiency_w(3,7,3,false);
  signalEfficiency_w(1,8,3,false);
  signalEfficiency_w(2,8,3,false);
  signalEfficiency_w(3,8,3,false);
  //WH
  signalEfficiency_w(1,7,4,false);
  signalEfficiency_w(2,7,4,false);
  signalEfficiency_w(3,7,4,false);
  signalEfficiency_w(1,8,4,false);
  signalEfficiency_w(2,8,4,false);
  signalEfficiency_w(3,8,4,false);
  //ttH
  signalEfficiency_w(1,7,5,false);
  signalEfficiency_w(2,7,5,false);
  signalEfficiency_w(3,7,5,false);
  signalEfficiency_w(1,8,5,false);
  signalEfficiency_w(2,8,5,false);
  signalEfficiency_w(3,8,5,false);

}


// The actual job
void signalEfficiency_w(int channel, double sqrts, int process, bool VBFtag) 
{
  TString schannel;
  if      (channel == 1) schannel = "4mu";
  else if (channel == 2) schannel = "4e";
  else if (channel == 3) schannel = "2e2mu";
  else cout << "Not a valid channel: " << channel << endl;

  TString sprocess;
  if      (process == 1) sprocess = "ggH";
  else if (process == 2) sprocess = "qqH";
  else if (process == 3) sprocess = "ZH";
  else if (process == 4) sprocess = "WH";
  else if (process == 5) sprocess = "ttH";
  else cout << "Not a valid channel: " << process << endl;

  TString ssqrts = (long) sqrts + TString("TeV");

  cout << "process = " << sprocess << " schannel = " << schannel << "  sqrts = " << sqrts <<  " VBFtag = " << VBFtag << endl;

  TString outfile = "CardFragments/signalEfficiency_" + sprocess + ssqrts + "_" + schannel + "_" + Form("%d",int(VBFtag)) +".txt";
  ofstream of(outfile,ios_base::out);

  gSystem->AddIncludePath("-I$ROOFITSYS/include");
  setTDRStyle(false);
  gStyle->SetStatX(-0.5);

  int nPoints=0;
  int* masses=0;
  double* mHVal=0;

  // Pick the correct set of mass points, set subpath
  TString filepath;
  if (process==1){
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
  } else if (process==2) {
    if (sqrts==7) {
      nPoints = nVBFPoints7TeV;
      masses  = VBFmasses7TeV;
      mHVal   = mHVBFVal7TeV;
      filepath = filePath7TeV;
    } else if (sqrts==8) {
      nPoints = nVBFPoints8TeV;
      masses  = VBFmasses8TeV;
      mHVal   = mHVBFVal8TeV;
      filepath =filePath8TeV;
    }
  } else if (process==3 || process==4 || process==5) {
    if (sqrts==7) {
      nPoints = nVHPoints7TeV;
      masses  = VHmasses7TeV;
      mHVal   = mHVHVal7TeV;
      filepath = filePath7TeV;
    } else if (sqrts==8) {
      nPoints = nVHPoints8TeV;
      masses  = VHmasses8TeV;
      mHVal   = mHVHVal8TeV;
      filepath =filePath8TeV;
    }
  }  

  float xMax = masses[nPoints-1]+10;
  //  cout << "xmax: " << xMax << endl;


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
	
  TString infile;

  for (int i = 0; i < nPoints; i++){
    if (process==1) infile = filepath+ "/" + (schannel=="2e2mu"?"2mu2e":schannel) + "/HZZ4lTree_H" + (long)masses[i] + ".root";
    else if (process==2) infile = filepath+ "/" + (schannel=="2e2mu"?"2mu2e":schannel) + "/HZZ4lTree_VBFH" + (long)masses[i] + ".root";
    else if (process==3 || process==4 || process==5) infile = filepath+ "/" + (schannel=="2e2mu"?"2mu2e":schannel) + "/HZZ4lTree_VH" + (long)masses[i] + ".root";    
    //    cout << "Opening input file: " << infile << endl;
    TFile *f = TFile::Open(infile) ;
    TTree *t1 = (TTree*) f->Get("SelectedTree");
    float mela, MC_weight_norm, MC_weight_PUWeight, MC_weight_powhegWeight,  MC_weight_dataMC;
    //std::vector<double> *JetPt = new vector<double>;
    int NJets;
    int numEventsRaw = 0;
    float numEventsPowheg =0;
    float numEventsPU =0;    
    float numEventsDataMC =0;
    int genProcessId=0;
    t1->SetBranchAddress("MC_weight_norm",&MC_weight_norm);
    t1->SetBranchAddress("MC_weight_powhegWeight",&MC_weight_powhegWeight);
    t1->SetBranchAddress("MC_weight_PUWeight",&MC_weight_PUWeight);
    t1->SetBranchAddress("MC_weight_dataMC",&MC_weight_dataMC);
    t1->SetBranchAddress("ZZLD",&mela);
    //t1->SetBranchAddress("JetPt", &JetPt);
    t1->SetBranchAddress("NJets",&NJets);
    t1->SetBranchAddress("genProcessId",&genProcessId);
    float totalCtr=0;
    float eff_noweight=0;
    float sumw2=0;
    float sumw_init2=0;
    //    float den = 0;
    for (int a = 0; a < t1->GetEntries(); a++){ 
      t1->GetEntry(a);
      if ((process==3 && genProcessId!=24) || (process==4 && genProcessId!=26) || (process==5 && (genProcessId!=121 && genProcessId!=122))) continue;
      if( (VBFtag == true && NJets > 1) || (VBFtag == false && NJets < 2)){
	totalCtr+=MC_weight_norm; 
	sumw2 += MC_weight_norm*MC_weight_norm;
	++numEventsRaw;      
	numEventsPowheg += MC_weight_powhegWeight;
	numEventsPU += MC_weight_powhegWeight*MC_weight_PUWeight;
	numEventsDataMC += MC_weight_powhegWeight*MC_weight_PUWeight*MC_weight_dataMC;
	if (MC_weight_powhegWeight>0) {
	  float w_initial = MC_weight_norm/(MC_weight_powhegWeight*MC_weight_PUWeight*MC_weight_dataMC);
	  eff_noweight += w_initial;
	  sumw_init2 += w_initial*w_initial; 
	  //	if (den!=0) den = (MC_weight_powhegWeight*MC_weight_PUWeight*MC_weight_dataMC)/MC_weight_norm;
	}
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

    if (verbose) {
      cout << " m = " << masses[i] 
	   << " : selected events= " << numEventsRaw 
	   << " Powheg Wtd= " << numEventsPowheg
	   << " PU Wtd= " << numEventsPU
	   << " Data/MC Wtd= " << numEventsDataMC
	   << " efficiency= " << totalCtr
//	 << " " << eff_noweight
	   << endl;
    }
    
    efficiencyVal[i] = totalCtr; 
    efficiencyErr[i] = sqrt(sumw2);

    // With no reweights
//     efficiencyVal[i] = eff_noweight;
//     efficiencyErr[i] =  sqrt(sumw_init2); //sqrt(eff_noweight*(1- eff_noweight)/den);
  
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

  if (channel==1 && sqrts==7) {    
    polyFunc->SetParLimits(7,0,0.035);
    polyFunc->SetParLimits(8,160,210);
    polyFunc->SetParLimits(9,30,50);
  }


  //  TF1 *polyFunc= new TF1("polyFunc","pol9", 110., xMax);

  polyFunc->SetLineColor(4);      
  TString cname = "eff" + sprocess + ssqrts + "_" + schannel + "_" + Form("%d",int(VBFtag));
  TCanvas *c = new TCanvas(cname,cname);
  c->SetGrid();

  TString outname = "sigFigs" + ssqrts +"/eff_" + sprocess + "_" + schannel + "_" + Form("%d",int(VBFtag));

  grEff->Fit(polyFunc,"Rt");
  TString xaxisText = "m_{" + schannel + "}";
  grEff->GetXaxis()->SetTitle(xaxisText);
  TString yaxisText = "efficiency, " + sprocess + ", " + schannel;
  grEff->GetYaxis()->SetTitle(yaxisText);
  grEff->Draw("AP");
  polyFunc->Draw("sames");
  c->Print(outname+".eps");
  c->Print(outname+".png"); // Does not work in batch?
  c->Print(outname+".root"); 

  cout << endl;
  cout << "------- Parameters for " << sprocess << " " << schannel << " sqrts=" << sqrts << " VBFtag = " << VBFtag << endl;
  cout << "   a1 = " << polyFunc->GetParameter(0) << endl;
  cout << "   a2 = " << polyFunc->GetParameter(1) << endl;
  cout << "   a3 = " << polyFunc->GetParameter(2) << endl;
  cout << "   a4 = " << polyFunc->GetParameter(3) << endl;
  cout << "   b1 = " << polyFunc->GetParameter(4) << endl;
  cout << "   b2 = " << polyFunc->GetParameter(5) << endl;
  cout << "   b3 = " << polyFunc->GetParameter(6) << endl;
  cout << "   g1 = " << polyFunc->GetParameter(7) << endl;
  cout << "   g2 = " << polyFunc->GetParameter(8) << endl;
  cout << "   g3 = " << polyFunc->GetParameter(9) << endl;
  cout << "---------------------------" << endl << endl;

  of << "## signal efficiency ##" << endl;
  of << "signalEff a1  " << polyFunc->GetParameter(0) << endl;
  of << "signalEff a2  " << polyFunc->GetParameter(1) << endl;
  of << "signalEff a3  " << polyFunc->GetParameter(2) << endl;
  of << "signalEff a4  " << polyFunc->GetParameter(3) << endl;
  of << "signalEff b1  " << polyFunc->GetParameter(4) << endl;
  of << "signalEff b2  " << polyFunc->GetParameter(5) << endl;
  of << "signalEff b3  " << polyFunc->GetParameter(6) << endl;
  of << "signalEff g1  " << polyFunc->GetParameter(7) << endl;
  of << "signalEff g2  " << polyFunc->GetParameter(8) << endl;
  of << "signalEff g3  " << polyFunc->GetParameter(9) << endl;
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
      if (verbose)    cout << "For mass, " << masses[i] << ": measured value is " << efficiencyVal[i] << " and difference from function is " << residual <<endl;
  }
  cout << "Largest residual= " << maxResidual << endl;
}


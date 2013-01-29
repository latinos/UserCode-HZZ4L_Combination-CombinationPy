/* 
 * Compute efficiencies for signals and write them in card fragments.
 * usage: 
 * -set all input variables in Config.h
 * -run with:
 * root -q -b signalEfficiency_w.C+
 * This runs on the 3 final states for each of the 5 production methods at 7 and 8 TeV and writes the output in a file (see stdout)
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

using namespace std;
using namespace ROOT::Math;

//----------> SET INPUT VARIABLES in Config.h
#include "Config.h"
//<----------

TFile* ftot,*fratio;

void signalEfficiency_w(int channel, double sqrts, int process, double JES);


// Run all final states and sqrts in one go
void signalEfficiency_w() {
  gSystem->Exec("mkdir -p sigFigs7TeV");
  gSystem->Exec("mkdir -p sigFigs8TeV");

  //ggH
  signalEfficiency_w(1,7,1,0.);
  signalEfficiency_w(2,7,1,0.);
  signalEfficiency_w(3,7,1,0.);
  signalEfficiency_w(1,8,1,0.);
  signalEfficiency_w(2,8,1,0.);
  signalEfficiency_w(3,8,1,0.);
  //qqH
  signalEfficiency_w(1,7,2,0.);
  signalEfficiency_w(2,7,2,0.);
  signalEfficiency_w(3,7,2,0.);
  signalEfficiency_w(1,8,2,0.);
  signalEfficiency_w(2,8,2,0.);
  signalEfficiency_w(3,8,2,0.);
  //ZH
  signalEfficiency_w(1,7,3,0.);
  signalEfficiency_w(2,7,3,0.);
  signalEfficiency_w(3,7,3,0.);
  signalEfficiency_w(1,8,3,0.);
  signalEfficiency_w(2,8,3,0.);
  signalEfficiency_w(3,8,3,0.);
  //WH
  signalEfficiency_w(1,7,4,0.);
  signalEfficiency_w(2,7,4,0.);
  signalEfficiency_w(3,7,4,0.);
  signalEfficiency_w(1,8,4,0.);
  signalEfficiency_w(2,8,4,0.);
  signalEfficiency_w(3,8,4,0.);
  //ttH
  signalEfficiency_w(1,7,5,0.);
  signalEfficiency_w(2,7,5,0.);
  signalEfficiency_w(3,7,5,0.);
  signalEfficiency_w(1,8,5,0.);
  signalEfficiency_w(2,8,5,0.);
  signalEfficiency_w(3,8,5,0.);
  /*
  //JES Up
  //ggH
  signalEfficiency_w(1,7,1,1.);
  signalEfficiency_w(2,7,1,1.);
  signalEfficiency_w(3,7,1,1.);
  signalEfficiency_w(1,8,1,1.);
  signalEfficiency_w(2,8,1,1.);
  signalEfficiency_w(3,8,1,1.);
  //qqH
  signalEfficiency_w(1,7,2,1.);
  signalEfficiency_w(2,7,2,1.);
  signalEfficiency_w(3,7,2,1.);
  signalEfficiency_w(1,8,2,1.);
  signalEfficiency_w(2,8,2,1.);
  signalEfficiency_w(3,8,2,1.);
  //ZH
  signalEfficiency_w(1,7,3,1.);
  signalEfficiency_w(2,7,3,1.);
  signalEfficiency_w(3,7,3,1.);
  signalEfficiency_w(1,8,3,1.);
  signalEfficiency_w(2,8,3,1.);
  signalEfficiency_w(3,8,3,1.);
  //WH
  signalEfficiency_w(1,7,4,1.);
  signalEfficiency_w(2,7,4,1.);
  signalEfficiency_w(3,7,4,1.);
  signalEfficiency_w(1,8,4,1.);
  signalEfficiency_w(2,8,4,1.);
  signalEfficiency_w(3,8,4,1.);
  //ttH
  signalEfficiency_w(1,7,5,1.);
  signalEfficiency_w(2,7,5,1.);
  signalEfficiency_w(3,7,5,1.);
  signalEfficiency_w(1,8,5,1.);
  signalEfficiency_w(2,8,5,1.);
  signalEfficiency_w(3,8,5,1.);

  //JES Down
  //ggH
  signalEfficiency_w(1,7,1,-1.);
  signalEfficiency_w(2,7,1,-1.);
  signalEfficiency_w(3,7,1,-1.);
  signalEfficiency_w(1,8,1,-1.);
  signalEfficiency_w(2,8,1,-1.);
  signalEfficiency_w(3,8,1,-1.);
  //qqH
  signalEfficiency_w(1,7,2,-1.);
  signalEfficiency_w(2,7,2,-1.);
  signalEfficiency_w(3,7,2,-1.);
  signalEfficiency_w(1,8,2,-1.);
  signalEfficiency_w(2,8,2,-1.);
  signalEfficiency_w(3,8,2,-1.);
  //ZH
  signalEfficiency_w(1,7,3,-1.);
  signalEfficiency_w(2,7,3,-1.);
  signalEfficiency_w(3,7,3,-1.);
  signalEfficiency_w(1,8,3,-1.);
  signalEfficiency_w(2,8,3,-1.);
  signalEfficiency_w(3,8,3,-1.);
  //WH
  signalEfficiency_w(1,7,4,-1.);
  signalEfficiency_w(2,7,4,-1.);
  signalEfficiency_w(3,7,4,-1.);
  signalEfficiency_w(1,8,4,-1.);
  signalEfficiency_w(2,8,4,-1.);
  signalEfficiency_w(3,8,4,-1.);
  //ttH
  signalEfficiency_w(1,7,5,-1.);
  signalEfficiency_w(2,7,5,-1.);
  signalEfficiency_w(3,7,5,-1.);
  signalEfficiency_w(1,8,5,-1.);
  signalEfficiency_w(2,8,5,-1.);
  signalEfficiency_w(3,8,5,-1.);
  */
}


// The actual job
void signalEfficiency_w(int channel, double sqrts, int process, double JES) 
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

  TString sjes;
  if      (JES==0.) sjes="";
  else if (JES>0.)  sjes="_up";
  else if (JES<0.)  sjes="_down";

  TString ssqrts = (long) sqrts + TString("TeV");

  cout << "process = " << sprocess << " schannel = " << schannel << "  sqrts = " << sqrts << " JES = " << JES <<endl;

  TString totoutfile = "CardFragments/signalEfficiency_"  + ssqrts + "_" + schannel + sjes + ".txt";
  TString ratiooutfile = "CardFragments/signalEfficiency_" + ssqrts + "_" + schannel + sjes + "_ratio.txt";
  ofstream oftot;
  if (process==1) oftot.open(totoutfile,ios_base::out);
  if (process!=1) oftot.open(totoutfile,ios_base::out | ios_base::app);
  ofstream ofrat;
  if (process==1) ofrat.open(ratiooutfile,ios_base::out);
  if (process!=1) ofrat.open(ratiooutfile,ios_base::out | ios_base::app);

  ftot = new TFile("sigFigs" + ssqrts +"/eff_" + sprocess + "_" + schannel + sjes + ".root","RECREATE");
  fratio = new TFile("sigFigs" + ssqrts +"/eff_" + sprocess + "_" + schannel + sjes + "_ratio.root","RECREATE");

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
      filepath = filePath8TeV;
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
      filepath = filePath8TeV;
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
      filepath = filePath8TeV;
    }
  }  

  float xMax = masses[nPoints-1]+10;


  const int arraySize=200;
  assert (arraySize>=nPoints);
  double totefficiencyVal[arraySize];
  double totefficiencyErr[arraySize];
  double dijetratioVal[arraySize];
  double dijetratioErr[arraySize];
	
  TString infile;

  for (int i = 0; i < nPoints; i++){
    if (process==1) infile = filepath+ "/" + (schannel=="2e2mu"?"2mu2e":schannel) + "/HZZ4lTree_H" + (long)masses[i] + ".root";
    else if (process==2) infile = filepath+ "/" + (schannel=="2e2mu"?"2mu2e":schannel) + "/HZZ4lTree_VBFH" + (long)masses[i] + ".root";
    else if (process==3 || process==4 || process==5) infile = filepath+ "/" + (schannel=="2e2mu"?"2mu2e":schannel) + "/HZZ4lTree_" + sprocess + (long)masses[i] + ".root";    
    TFile *f = TFile::Open(infile) ;
    TTree *t1 = (TTree*) f->Get("SelectedTree");
    float mela, MC_weight_norm, MC_weight_PUWeight, MC_weight_powhegWeight,  MC_weight_dataMC;
    //int NJets;
    int genProcessId=0;
    vector<double> *JetPt=0;
    vector<double> *JetSigma=0;
    t1->SetBranchAddress("MC_weight_norm",&MC_weight_norm);
    t1->SetBranchAddress("MC_weight_powhegWeight",&MC_weight_powhegWeight);
    t1->SetBranchAddress("MC_weight_PUWeight",&MC_weight_PUWeight);
    t1->SetBranchAddress("MC_weight_dataMC",&MC_weight_dataMC);
    t1->SetBranchAddress("ZZLD",&mela);
    //t1->SetBranchAddress("NJets",&NJets);
    t1->SetBranchAddress("genProcessId",&genProcessId);
    t1->SetBranchAddress("JetPt",&JetPt);
    t1->SetBranchAddress("JetSigma",&JetSigma);
    //Initialize counters for non-dijet events
    int numndjEventsRaw = 0;
    float numndjEventsPowheg =0;
    float numndjEventsPU =0;    
    float numndjEventsDataMC =0;
    float totalndjCtr=0;
    float ndjeff_noweight=0;
    float ndjsumw2=0;
    float ndjsumw_init2=0;
    //Initialize counters for dijet events
    int numdjEventsRaw = 0;
    float numdjEventsPowheg =0;
    float numdjEventsPU =0;    
    float numdjEventsDataMC =0;
    float totaldjCtr=0;
    float djeff_noweight=0;
    float djsumw2=0;
    float djsumw_init2=0;
    for (int a = 0; a < t1->GetEntries(); a++){ 
      t1->GetEntry(a);
      if ((process==3 && genProcessId!=24) || (process==4 && genProcessId!=26) || (process==5 && (genProcessId!=121 && genProcessId!=122))) continue;
      int NJets=0;
      double jetptc=0;
      for (int j=0; j<JetPt->size();j++){
	if (JES==0.) jetptc=JetPt->at(j);
	else if (JES!=0.) jetptc=JetPt->at(j)*(1+JES*JetSigma->at(j));
	if (jetptc>30.) NJets++;
      }
      if(NJets < 2){
	totalndjCtr+=MC_weight_norm; 
	ndjsumw2 += MC_weight_norm*MC_weight_norm;
	++numndjEventsRaw;      
	numndjEventsPowheg += MC_weight_powhegWeight;
	numndjEventsPU += MC_weight_powhegWeight*MC_weight_PUWeight;
	numndjEventsDataMC += MC_weight_powhegWeight*MC_weight_PUWeight*MC_weight_dataMC;
	if (MC_weight_powhegWeight>0) {
	  float w_initial = MC_weight_norm/(MC_weight_powhegWeight*MC_weight_PUWeight*MC_weight_dataMC);
	  ndjeff_noweight += w_initial;
	  ndjsumw_init2 += w_initial*w_initial; 
	}
      } else if(NJets > 1){
	totaldjCtr+=MC_weight_norm; 
	djsumw2 += MC_weight_norm*MC_weight_norm;
	++numdjEventsRaw;      
	numdjEventsPowheg += MC_weight_powhegWeight;
	numdjEventsPU += MC_weight_powhegWeight*MC_weight_PUWeight;
	numdjEventsDataMC += MC_weight_powhegWeight*MC_weight_PUWeight*MC_weight_dataMC;
	if (MC_weight_powhegWeight>0) {
	  float w_initial = MC_weight_norm/(MC_weight_powhegWeight*MC_weight_PUWeight*MC_weight_dataMC);
	  djeff_noweight += w_initial;
	  djsumw_init2 += w_initial*w_initial; 
	}
      }
    }

    // FIXME: the 7TeV samples are assumed to have the ad-hoc correction factor for the mll>12 gen cut,
    // except for the 124,125,126 new samples. As this factor is accounted for in the x-section, we have to 
    // apply it here
    float m = masses[i];
    if (process==1 && sqrts==7 && m>=123.9 &&  m<=126.1) {
      float mllCorr = 0.5 + 0.5*erf((m-80.85)/50.42);
      totalndjCtr = totalndjCtr/mllCorr;
      ndjeff_noweight=ndjeff_noweight/mllCorr;
      totaldjCtr = totaldjCtr/mllCorr;
      djeff_noweight=djeff_noweight/mllCorr;
    }

    if (verbose) {
      cout << " m = " << masses[i] 
	   << " :" <<endl;
      cout << "Selected non-dijet events= " << numndjEventsRaw 
	   << " Powheg Wtd= " << numndjEventsPowheg
	   << " PU Wtd= " << numndjEventsPU
	   << " Data/MC Wtd= " << numndjEventsDataMC
	   << " Efficiency= " << totalndjCtr
	   << endl;
      cout << "Selected dijet events= " << numdjEventsRaw 
	   << " Powheg Wtd= " << numdjEventsPowheg
	   << " PU Wtd= " << numdjEventsPU
	   << " Data/MC Wtd= " << numdjEventsDataMC
	   << " Efficiency= " << totaldjCtr
	   << endl;
    }
    
    totefficiencyVal[i] = totalndjCtr + totaldjCtr;
    cout<<sprocess<<" "<<m<<" "<<totefficiencyVal[i]<<endl;
    totefficiencyErr[i] = sqrt(ndjsumw2 + djsumw2);
    dijetratioVal[i]=totaldjCtr/totefficiencyVal[i];
    dijetratioErr[i]=sqrt(pow(totalndjCtr,2)*djsumw2 + pow(totaldjCtr,2)*ndjsumw2)/pow(totefficiencyVal[i],2);
  
    f->Close();
  }
	

  TGraphErrors* totgrEff = new TGraphErrors( nPoints, mHVal, totefficiencyVal, 0, totefficiencyErr);
  TGraphErrors* ratgrEff = new TGraphErrors( nPoints, mHVal, dijetratioVal, 0, dijetratioErr);
  totgrEff->SetMarkerStyle(20);
  ratgrEff->SetMarkerStyle(20);

  //ICHEP parametrization	
  //TF1 *polyFunc= new TF1("polyFunc","([0]+[1]*TMath::Erf( (x-[2])/[3] ))*([4]+[5]*x+[6]*x*x)", 110., xMax);
  //polyFunc->SetParameters(-4.42749e+00,4.61212e+0,-6.21611e+01,1.13168e+02,2.14321e+00,1.04083e-03,4.89570e-07);


  TF1 *polyFunctot= new TF1("polyFunctot","([0]+[1]*TMath::Erf( (x-[2])/[3] ))*([4]+[5]*x+[6]*x*x)+[7]*TMath::Gaus(x,[8],[9])", 110., xMax);
  polyFunctot->SetParameters(-4.42749e+00,4.61212e+0,-6.21611e+01,1.13168e+02,2.14321e+00,1.04083e-03,4.89570e-07, 0.03, 200, 30);
  polyFunctot->SetParLimits(7,0,0.1);
  polyFunctot->SetParLimits(8,160,210);
  polyFunctot->SetParLimits(9,20,50);

  if (channel==1 && sqrts==7) {    
    polyFunctot->SetParLimits(7,0,0.035);
    polyFunctot->SetParLimits(8,160,210);
    polyFunctot->SetParLimits(9,30,50);
  }

  polyFunctot->SetLineColor(4);      
  TString cname = "eff" + sprocess + ssqrts + "_" + schannel;
  TCanvas *ctot = new TCanvas(cname,cname);
  ctot->SetGrid();

  TString outname = "sigFigs" + ssqrts +"/eff_" + sprocess + "_" + schannel + "_" + sjes;

  totgrEff->Fit(polyFunctot,"Rt");
  TString xaxisText = "m_{" + schannel + "}";
  totgrEff->GetXaxis()->SetTitle(xaxisText);
  TString yaxisText = "Efficiency, " + sprocess + ", " + schannel;
  totgrEff->GetYaxis()->SetTitle(yaxisText);
  totgrEff->SetMinimum(0.0);
  totgrEff->SetMaximum(1.0);
  totgrEff->Draw("AP");
  polyFunctot->Draw("sames");
  ctot->Print(outname+".eps");
  ctot->Print(outname+".png"); // Does not work in batch?
  ctot->Print(outname+".pdf"); 
  ftot->cd();
  totgrEff->Write("TotalEfficiency");
  ftot->Close();

  cout << endl;
  cout << "------- Parameters for " << sprocess << " " << schannel << " sqrts=" << sqrts << endl;
  cout << "   a1 = " << polyFunctot->GetParameter(0) << endl;
  cout << "   a2 = " << polyFunctot->GetParameter(1) << endl;
  cout << "   a3 = " << polyFunctot->GetParameter(2) << endl;
  cout << "   a4 = " << polyFunctot->GetParameter(3) << endl;
  cout << "   b1 = " << polyFunctot->GetParameter(4) << endl;
  cout << "   b2 = " << polyFunctot->GetParameter(5) << endl;
  cout << "   b3 = " << polyFunctot->GetParameter(6) << endl;
  cout << "   g1 = " << polyFunctot->GetParameter(7) << endl;
  cout << "   g2 = " << polyFunctot->GetParameter(8) << endl;
  cout << "   g3 = " << polyFunctot->GetParameter(9) << endl;
  cout << "---------------------------" << endl << endl;

  string oftotprocess;
  if (process==1) oftotprocess="";
  if (process!=1) oftotprocess=sprocess;

  if (process==1) {
    oftot << endl;
    oftot << "## signal efficiency ##" << endl;
  }
  oftot << "signalEff " << oftotprocess << "a1  " << polyFunctot->GetParameter(0) << endl;
  oftot << "signalEff " << oftotprocess << "a2  " << polyFunctot->GetParameter(1) << endl;
  oftot << "signalEff " << oftotprocess << "a3  " << polyFunctot->GetParameter(2) << endl;
  oftot << "signalEff " << oftotprocess << "a4  " << polyFunctot->GetParameter(3) << endl;
  oftot << "signalEff " << oftotprocess << "b1  " << polyFunctot->GetParameter(4) << endl;
  oftot << "signalEff " << oftotprocess << "b2  " << polyFunctot->GetParameter(5) << endl;
  oftot << "signalEff " << oftotprocess << "b3  " << polyFunctot->GetParameter(6) << endl;
  oftot << "signalEff " << oftotprocess << "g1  " << polyFunctot->GetParameter(7) << endl;
  oftot << "signalEff " << oftotprocess << "g2  " << polyFunctot->GetParameter(8) << endl;
  oftot << "signalEff " << oftotprocess << "g3  " << polyFunctot->GetParameter(9) << endl;
  oftot << endl;
  oftot.close();
  
  cname = "eff" + sprocess + ssqrts + "_" + schannel + "_ratio";
  TCanvas *crat = new TCanvas(cname,cname);
  crat->SetGrid();

  outname = "sigFigs" + ssqrts +"/eff_" + sprocess + "_" + schannel + "_" + sjes + "_ratio";

  TF1 *ratiofit;
  if (process==1 || process==2) ratiofit = new TF1("ratiofit","([0]+[1]*x+[2]*x*x)",110.,xMax);
  if (process==3 || process==4 || process==5 ) ratiofit = new TF1("ratiofit","([0]+[1]*x)",110.,xMax);

  ratgrEff->Fit(ratiofit,"Rt");
  ratgrEff->GetXaxis()->SetTitle(xaxisText);
  TString yaxisratio = "Dijet ratio, " + sprocess + ", " + schannel;
  ratgrEff->GetYaxis()->SetTitle(yaxisratio);
  ratgrEff->SetMinimum(0.0);
  ratgrEff->SetMaximum(1.0);
  ratgrEff->Draw("AP");
  crat->Print(outname+".eps");
  crat->Print(outname+".png"); // Does not work in batch?
  crat->Print(outname+".pdf");
  fratio->cd();
  ratgrEff->Write("Ratio");
  fratio->Close();
  
  cout << endl;
  cout << "------- Parameters for " << sprocess << " " << schannel << " sqrts=" << sqrts << endl;
  cout << "   a1 = " << ratiofit->GetParameter(0) << endl;
  cout << "   a2 = " << ratiofit->GetParameter(1) << endl;
  if (process==1 || process==2) cout << "   a3 = " << ratiofit->GetParameter(2) << endl;
  cout << "---------------------------" << endl << endl;

  if (process==1) ofrat << "## signal efficiency ratios ##" << endl;
  ofrat << "signalEff tagged_" << sprocess << "_ratio " << ratiofit->GetParameter(0) << "+(" << ratiofit->GetParameter(1) << "*@0)";
  if (process==1 || process==2) ofrat << "+(" << ratiofit->GetParameter(2) << "*@0*@0)" << endl;
  else if (process==3 || process==4 ) ofrat << endl;
  else if (process==5) ofrat << endl << endl;
  ofrat.close();

  // deviations
  cout << "Deviations..." << endl;
  double maxResidual=0;
  for (int i = 0; i < nPoints; i++){
    double eval = polyFunctot->Eval(masses[i]);
    double residual = (eval - totefficiencyVal[i]);
    maxResidual = max(maxResidual,fabs(residual));
    if (verbose)    cout << "For mass, " << masses[i] << ": measured value is " << totefficiencyVal[i] << " and difference from function is " << residual <<endl;
  }
  cout << "Largest residual= " << maxResidual << endl;
}


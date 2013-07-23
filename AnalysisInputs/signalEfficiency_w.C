
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
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSystem.h"

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TROOT.h>
#include "Higgs/Higgs_CS_and_Width/include/HiggsCSandWidth.h"
#endif

#include "../CreateDatacards/include/tdrstyle.cc"

bool verbose = false;
bool useNewGGHPowheg = true;

using namespace std;
using namespace ROOT::Math;

//----------> SET INPUT VARIABLES in Config.h
#include "Config.h"
//<----------

enum Channel {fs4mu=1, fs4e=2, fs2e2mu=3};
enum Process {ggH=1, qqH=2, ZH=3, WH=4, ttH=5};

TFile* ftot,*fratio;

void signalEfficiency_w(int channel, double sqrts, int process, double JES, ofstream* txtYields=0);

struct Counts {
  Counts():
    totalCtr(0),
    numEventsRaw(0),    
    numEventsPowheg(0),
    numEventsPU(0),    
    numEventsDataMC(0),
    sumw2(0),
    eff_noweight(0),
    sumw_init2(0)       
  {}
  float totalCtr;         //  Actual efficiency
  int   numEventsRaw;     //  total raw events 
  float numEventsPowheg;  //  same, reweighted for hi-mass weight
  float numEventsPU;      //  same, reweighted for PU
  float numEventsDataMC;  //  same, reweighted for eff corr
  float sumw2;            //  square of weights, for error
  float eff_noweight;     //  efficiency, with no reweight
  float sumw_init2;       //  its error
  
  //Increase counters
  void incrCounters(float effW, float w_PUWeight, float w_powhegWeight, float w_dataMC){
    numEventsRaw++; 
    totalCtr+=effW;
    numEventsPowheg+=w_powhegWeight;
    numEventsPU+=w_powhegWeight*w_PUWeight;
    numEventsDataMC+=w_powhegWeight*w_PUWeight*w_dataMC;
    sumw2+=effW*effW;
    if (w_powhegWeight>0) {
      float w_initial = effW/(w_powhegWeight*w_PUWeight*w_dataMC); //initial weight (1/nevt_gen)
      eff_noweight += w_initial; 
      sumw_init2 += w_initial*w_initial;
    }
  }

};

  


// Run all final states and sqrts in one go
void signalEfficiency_w() {
  gSystem->Exec("mkdir -p sigFigs7TeV");
  gSystem->Exec("mkdir -p sigFigs8TeV");

  // // Not really needed
  //   gSystem->Load("libHiggsHiggs_CS_and_Width.so");
  //   gROOT->LoadMacro("$CMSSW_BASE/src/Higgs/Higgs_CS_and_Width/include/HiggsCSandWidth.h+");

  float JES = 0; // 1=JES up; -1=JES down

  //Create and open the text output file with yields per mass points
  ofstream fileOutYields;
  TString yieldsFileName = "yields.txt";
  fileOutYields.open(yieldsFileName);
  fileOutYields << "############### YIELDS ###############" << endl;
  fileOutYields << "Lumi 7 TeV: " << lumi7TeV << endl;
  fileOutYields << "Lumi 8 TeV: " << lumi8TeV << endl << endl;

  //ggH
  signalEfficiency_w(fs4mu,  7,ggH,JES,&fileOutYields);
  signalEfficiency_w(fs4e,   7,ggH,JES,&fileOutYields);
  signalEfficiency_w(fs2e2mu,7,ggH,JES,&fileOutYields);
  signalEfficiency_w(fs4mu,  8,ggH,JES,&fileOutYields);
  signalEfficiency_w(fs4e,   8,ggH,JES,&fileOutYields);
  signalEfficiency_w(fs2e2mu,8,ggH,JES,&fileOutYields);  

  //qqH
  signalEfficiency_w(fs4mu,  7,qqH,JES,&fileOutYields);
  signalEfficiency_w(fs4e,   7,qqH,JES,&fileOutYields);
  signalEfficiency_w(fs2e2mu,7,qqH,JES,&fileOutYields);
  signalEfficiency_w(fs4mu,  8,qqH,JES,&fileOutYields);
  signalEfficiency_w(fs4e,   8,qqH,JES,&fileOutYields);
  signalEfficiency_w(fs2e2mu,8,qqH,JES,&fileOutYields);

  //ZH
  signalEfficiency_w(fs4mu,  7,ZH,JES,&fileOutYields);
  signalEfficiency_w(fs4e,   7,ZH,JES,&fileOutYields);
  signalEfficiency_w(fs2e2mu,7,ZH,JES,&fileOutYields);
  signalEfficiency_w(fs4mu,  8,ZH,JES,&fileOutYields);
  signalEfficiency_w(fs4e,   8,ZH,JES,&fileOutYields);
  signalEfficiency_w(fs2e2mu,8,ZH,JES,&fileOutYields);

  //WH
  signalEfficiency_w(fs4mu,  7,WH,JES,&fileOutYields);
  signalEfficiency_w(fs4e,   7,WH,JES,&fileOutYields);
  signalEfficiency_w(fs2e2mu,7,WH,JES,&fileOutYields);
  signalEfficiency_w(fs4mu,  8,WH,JES,&fileOutYields);
  signalEfficiency_w(fs4e,   8,WH,JES,&fileOutYields);
  signalEfficiency_w(fs2e2mu,8,WH,JES,&fileOutYields);

  //ttH
  signalEfficiency_w(fs4mu,  7,ttH,JES,&fileOutYields);
  signalEfficiency_w(fs4e,   7,ttH,JES,&fileOutYields);
  signalEfficiency_w(fs2e2mu,7,ttH,JES,&fileOutYields);
  signalEfficiency_w(fs4mu,  8,ttH,JES,&fileOutYields);
  signalEfficiency_w(fs4e,   8,ttH,JES,&fileOutYields);
  signalEfficiency_w(fs2e2mu,8,ttH,JES,&fileOutYields);

  fileOutYields.close();
  
}


// The actual job
void signalEfficiency_w(int channel, double sqrts, int process, double JES, ofstream* txtYields) 
{
  TString schannel;
  if      (channel == 1) schannel = "4mu";
  else if (channel == 2) schannel = "4e";
  else if (channel == 3) schannel = "2e2mu";
  else cout << "Not a valid channel: " << channel << endl;

  TString schannelFilename = (schannel=="2e2mu"?"2mu2e":schannel); // adapt to naming convention of tree files

  TString sprocess;
  if      (process == ggH) sprocess = "ggH";
  else if (process == qqH) sprocess = "qqH";
  else if (process == ZH) sprocess = "ZH";
  else if (process == WH) sprocess = "WH";
  else if (process == ttH) sprocess = "ttH";
  else cout << "Not a valid channel: " << process << endl;  

  TString sjes;
  if      (JES==0.) sjes="";
  else if (JES>0.)  sjes="_up";
  else if (JES<0.)  sjes="_down";

  TString ssqrts = (long) sqrts + TString("TeV");


  // Print table with yields
  (*txtYields) << endl << endl 
	       << left << setw(7) << "*** Summary: " << sprocess << " process, sqrts = " << fixed << setprecision(0) <<sqrts << " TeV, channel = " << schannel << " ***" << endl << endl;
  (*txtYields) << left << setw(7) << "mH" << setw(13) << "XS*BR" << setw(13) 
	       << "Eff unt" << setw(13) << "Yield unt" << setw(13) << "Eff tag" 
	       << setw(13) << "Yield tag" << setw(13) << "Eff TOT" << setw(13) << "Yield TOT" 
	       << setw(13) << "Eff unt" << setw(13) << "Yield unt" << setw(13) << "Eff tag" 
	       << setw(13) << "Yield tag" << setw(13) << "Eff TOT" << setw(13) << "Yield TOT" 
	       << setw(13) << "n. raw" << setw(13) << "n. W pwg" << setw(13) << "n. W PU" 
	       << setw(13) << "n. W eff" << setw(13) 
	       << endl << left << setw(7) << " " << setw(13) << " " << setw(13)
	       << "(in MW)" << setw(13) << "(in MW)" << setw(13) << "(in MW)" 
	       << setw(13) << "(in MW)" << setw(13) << "(in MW)" << setw(13) << "(in MW)" 
	       << setw(13) << "(full)" << setw(13) << "(full)" << setw(13) << "(full)" 
	       << setw(13) << "(full)" << setw(13) << "(full)" << setw(13) << "(full)" 
	       << setw(13) << "(full U+T)" << setw(13) << "(full U+T)" << setw(13) << "(full U+T)" << setw(13) << "(full U+T)" 
	       << endl << endl;

  cout << "process = " << sprocess << " schannel = " << schannel << "  sqrts = " << sqrts << " JES = " << JES <<endl;

  TString totoutfile = "CardFragments/signalEfficiency_"  + ssqrts + "_" + schannel + sjes + ".txt";
  TString ratiooutfile = "CardFragments/signalEfficiency_" + ssqrts + "_" + schannel + sjes + "_ratio.txt";

  // Create card fragments using new powheg samples
  ofstream oftot;
  ofstream ofrat;

  if (process==ggH) {
    oftot.open(totoutfile,ios_base::out);
    ofrat.open(ratiooutfile,ios_base::out);
  } else {
    oftot.open(totoutfile,ios_base::out | ios_base::app);
    ofrat.open(ratiooutfile,ios_base::out | ios_base::app);
  }

  ftot = new TFile("sigFigs" + ssqrts +"/eff_" + sprocess + "_" + schannel + sjes + (useNewGGHPowheg ? ".root" : "_oldPwg.root"),"RECREATE");
  fratio = new TFile("sigFigs" + ssqrts +"/eff_" + sprocess + "_" + schannel + sjes + (useNewGGHPowheg ? "_ratio.root" : "_ratio_oldPwg.root"),"RECREATE");

  gSystem->AddIncludePath("-I$ROOFITSYS/include");
  setTDRStyle(false);
  gStyle->SetStatX(-0.5);

  int nPoints=0;
  int* masses=0;
  double* mHVal=0;

  // Pick the correct set of mass points, set subpath
  TString filepath;
  if (process==ggH){
    if (useNewGGHPowheg) {
      if (sqrts==7) {
	nPoints = nPoints7TeV_p15;
	masses  = masses7TeV_p15;
	mHVal   = mHVal7TeV_p15;
	filepath = filePath7TeV;
      } else if (sqrts==8) {
	nPoints = nPoints8TeV_p15;
	masses  = masses8TeV_p15;
	mHVal   = mHVal8TeV_p15;
	filepath = filePath8TeV;
      }
    } else { // OLD powheg samples
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
    }
  } else if (process==qqH) {
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
  } else if (process==ZH || process==WH || process==ttH) {
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
  double totefficiencyValInMW[arraySize];
  double totefficiencyErrInMW[arraySize];
  double dijetratioValInMW[arraySize];
  double dijetratioErrInMW[arraySize];


  // Define the object to compute XS and BRs
  HiggsCSandWidth *myCSW = new HiggsCSandWidth(gSystem->ExpandPathName("$CMSSW_BASE/src/Higgs/Higgs_CS_and_Width/txtFiles/"));
	
  TString infile;

  for (int i = 0; i < nPoints; i++){

    // Compute XS and BR
    double xsTimesBR = 0.;
    double BRH4e = myCSW->HiggsBR(12,masses[i]);
    double BRH2e2mu = myCSW->HiggsBR(13,masses[i]);
    double BRHZZ = myCSW->HiggsBR(11,masses[i]);
    double BR = BRHZZ;
    if (process==ggH || process==qqH) { 
      if (channel==fs4mu || channel==fs4e) BR = BRH4e;
      else BR = BRH2e2mu;
    }

    if (process==ggH)      xsTimesBR = BR*myCSW->HiggsCS(1,masses[i],sqrts);
    else if (process==qqH) xsTimesBR = BR*myCSW->HiggsCS(2,masses[i],sqrts);
    else if (process==ZH)  xsTimesBR = BR*myCSW->HiggsCS(3,masses[i],sqrts);
    else if (process==WH)  xsTimesBR = BR*myCSW->HiggsCS(4,masses[i],sqrts);
    else if (process==ttH) xsTimesBR = BR*myCSW->HiggsCS(5,masses[i],sqrts); 


    if (process==ggH) {
      if (useNewGGHPowheg){ 
	infile = filepath+ "/" + schannelFilename + "/HZZ4lTree_powheg15" + (masses[i]>200?"H":"jhuGenV3H") + (long)masses[i] + ".root";
      } else {
	infile = filepath+ "/" + schannelFilename + "/HZZ4lTree_H" + (long)masses[i] + ".root";
      }
    }
    
    else if (process==qqH) infile = filepath+ "/" + schannelFilename + "/HZZ4lTree_VBFH" + (long)masses[i] + ".root";
    else if (process==WH || process==ZH || process==ttH) infile = filepath+ "/" + schannelFilename + "/HZZ4lTree_" + sprocess + (long)masses[i] + ".root";    

    TFile *f = TFile::Open(infile) ;
    TTree *t1 = (TTree*) f->Get("SelectedTree");
    float MC_weight_norm, MC_weight_PUWeight, MC_weight_powhegWeight,  MC_weight_dataMC;
    float MC_weight_noxsec;
    float GenHPt;
    //int NJets;
    short genProcessId=0;
    //    short NJets30;
    vector<double> *JetPt=0;
    vector<double> *JetSigma=0;
    float ZZMass;
    t1->SetBranchAddress("MC_weight_norm",&MC_weight_norm); // For efficiency vs "proper" final state
    t1->SetBranchAddress("MC_weight_noxsec",&MC_weight_noxsec); // For efficiency vs all gen events
    t1->SetBranchAddress("MC_weight_powhegWeight",&MC_weight_powhegWeight);
    t1->SetBranchAddress("MC_weight_PUWeight",&MC_weight_PUWeight);
    t1->SetBranchAddress("MC_weight_dataMC",&MC_weight_dataMC);
    //t1->SetBranchAddress("NJets",&NJets);
    t1->SetBranchAddress("genProcessId",&genProcessId);
    t1->SetBranchAddress("JetPt",&JetPt);
    t1->SetBranchAddress("JetSigma",&JetSigma);
    //    t1->SetBranchAddress("NJets30",&NJets30);
    t1->SetBranchAddress("GenHPt",&GenHPt);
    t1->SetBranchAddress("ZZMass",&ZZMass);

    //Initialize counters for non-dijet events
    Counts* untagInMW = new Counts(); Counts* untagAll = new Counts();
    Counts* dijetInMW = new Counts(); Counts* dijetAll = new Counts();

    // Find window width // FIXME move to external function
    double valueWidth = myCSW->HiggsWidth(0,masses[i]);
    double windowVal = max(valueWidth,1.);
    
    double lowside = 100.;
    double highside = 1000.0;
      
    if (masses[i] >= 275){
      lowside = 180.0;
      highside = 650.0;
    }
    if (masses[i] >= 350){
      lowside = 200.0;
      highside = 900.0;
    }
    if (masses[i] >= 500){
      lowside = 250.0;
      highside = 1000.0;
    }
    if (masses[i] >= 700){
      lowside = 350.0;
      highside = 1400.0;
    }
    double low_M = max( (masses[i] - 20.*windowVal), lowside);
    double high_M = min((masses[i] + 15.*windowVal), highside);


    // // Load Higgs pT weights for old powheg
    //     TFile* fW; TH1D* h_HPtWeight; 
    //     TString fW_str = "./HPtWeights/weight_";
    //     fW_str += (long)masses[i];
    //     fW_str += (TString)".root";
    //     cout << fW_str << endl;
    //     if (process==ggH) {
    //       fW = TFile::Open(fW_str,"READ");
    //       h_HPtWeight = (TH1D*)fW->Get("h_weight");
    //     }


    for (int a = 0; a < t1->GetEntries(); a++){ 
      t1->GetEntry(a);
      // Skip VH events that do not belong to the right gen prod mechanism. This is no longer necessary with the proper VH samples
      if ((process==ZH && genProcessId!=24) || (process==WH && genProcessId!=26) || (process==ttH && (genProcessId!=121 && genProcessId!=122))) continue; 


      // We use the efficiency vs. generated events in the proper FS for ggH, VBF, and the efficiency vs all generated events for VH, ttH
      float effw = MC_weight_norm;
      if (process==ZH) {
	effw = MC_weight_noxsec*filter_eff_ZH_8TeV;
      }
      else if (process==WH){
	effw = MC_weight_noxsec*filter_eff_WH_8TeV;
      }
      else if (process==ttH){
	effw = MC_weight_noxsec*filter_eff_ttH_8TeV;
      }      

//       double HPtWeight = 1.;
//       if (process==ggH) HPtWeight = h_HPtWeight->GetBinContent(h_HPtWeight->FindBin(GenHPt));
//       //cout << "Higgs pT weight = " << HPtWeight << endl;
//       effw*=HPtWeight;
      

      int NJets=0;
      double jetptc=0;
      for (unsigned int j=0; j<JetPt->size();j++){
	if (JES==0.) jetptc=JetPt->at(j);
	else if (JES!=0.) jetptc=JetPt->at(j)*(1+JES*JetSigma->at(j));
	if (jetptc>30.) NJets++;
      }

      // Untagged
      if (NJets<2){
	untagAll->incrCounters(effw, MC_weight_PUWeight, MC_weight_powhegWeight, MC_weight_dataMC);
	if ( (ZZMass>low_M && ZZMass<high_M) ) untagInMW->incrCounters(effw, MC_weight_PUWeight, MC_weight_powhegWeight, MC_weight_dataMC);
      }
      else{ // Dijet
	dijetAll->incrCounters(effw, MC_weight_PUWeight, MC_weight_powhegWeight, MC_weight_dataMC);
	if ( (ZZMass>low_M && ZZMass<high_M) ) dijetInMW->incrCounters(effw, MC_weight_PUWeight, MC_weight_powhegWeight, MC_weight_dataMC);
      }
      
    }

    // FIXME: the 7TeV samples are assumed to have the ad-hoc correction factor for the mll>12 gen cut,
    // except for the 124,125,126 new samples. As this factor is accounted for in the x-section, we have to 
    // apply it here
    float m = masses[i];
    if (process==ggH && sqrts==7 && m>=123.9 &&  m<=126.1) {
      float mllCorr = 0.5 + 0.5*erf((m-80.85)/50.42);
      untagAll->totalCtr = untagAll->totalCtr/mllCorr;
      untagAll->eff_noweight=untagAll->eff_noweight/mllCorr;
      untagInMW->totalCtr = untagInMW->totalCtr/mllCorr;
      untagInMW->eff_noweight=untagInMW->eff_noweight/mllCorr;
      dijetAll->totalCtr = dijetAll->totalCtr/mllCorr;
      dijetAll->eff_noweight=dijetAll->eff_noweight/mllCorr;
      dijetInMW->totalCtr = dijetInMW->totalCtr/mllCorr;
      dijetInMW->eff_noweight=dijetInMW->eff_noweight/mllCorr;
    }

    if (verbose) {
      cout << " m = " << masses[i] 
	   << " :" <<endl;
      cout << "Selected non-dijet events (all) = " << untagAll->numEventsRaw 
	   << " Powheg Wtd= " << untagAll->numEventsPowheg
	   << " PU Wtd= " << untagAll->numEventsPU
	   << " Data/MC Wtd= " << untagAll->numEventsDataMC
	   << " Efficiency= " << untagAll->totalCtr
	   << endl;
      cout << "Selected non-dijet events (in mass window) = " << untagInMW->numEventsRaw 
	   << " Powheg Wtd= " << untagInMW->numEventsPowheg
	   << " PU Wtd= " << untagInMW->numEventsPU
	   << " Data/MC Wtd= " << untagInMW->numEventsDataMC
	   << " Efficiency= " << untagInMW->totalCtr
	   << endl;
      cout << "Selected dijet events (all) = " << dijetAll->numEventsRaw 
	   << " Powheg Wtd= " << dijetAll->numEventsPowheg
	   << " PU Wtd= " << dijetAll->numEventsPU
	   << " Data/MC Wtd= " << dijetAll->numEventsDataMC
	   << " Efficiency= " << dijetAll->totalCtr
	   << endl;
      cout << "Selected dijet events (in mass window) = " << dijetInMW->numEventsRaw 
	   << " Powheg Wtd= " << dijetInMW->numEventsPowheg
	   << " PU Wtd= " << dijetInMW->numEventsPU
	   << " Data/MC Wtd= " << dijetInMW->numEventsDataMC
	   << " Efficiency= " << dijetInMW->totalCtr
	   << endl;
    }

    // All events
    totefficiencyVal[i] = untagAll->totalCtr + dijetAll->totalCtr;
    cout << "All events:            " << sprocess << " " << m << " " << totefficiencyVal[i]<<endl;
    totefficiencyErr[i] = sqrt(untagAll->sumw2 + dijetAll->sumw2);
    dijetratioVal[i]=dijetAll->totalCtr/totefficiencyVal[i];
    dijetratioErr[i]=sqrt(pow(untagAll->totalCtr,2)*dijetAll->sumw2 + pow(dijetAll->totalCtr,2)*untagAll->sumw2)/pow(totefficiencyVal[i],2); // FIXME: misses 1 term 

    // Events inside the mass window
    totefficiencyValInMW[i] = untagInMW->totalCtr + dijetInMW->totalCtr;
    cout << "Events in mass window: " << sprocess << " " << m << " " << totefficiencyValInMW[i]<<endl;
    totefficiencyErrInMW[i] = sqrt(untagInMW->sumw2 + dijetInMW->sumw2);
    dijetratioValInMW[i]=dijetInMW->totalCtr/totefficiencyValInMW[i];
    dijetratioErrInMW[i]=sqrt(pow(untagInMW->totalCtr,2)*dijetInMW->sumw2 + pow(dijetInMW->totalCtr,2)*untagInMW->sumw2)/pow(totefficiencyValInMW[i],2);
    
    
    // Write yields to output file
    double lumi = -1.;
    sqrts == 7 ? lumi = lumi7TeV*1000 : lumi = lumi8TeV*1000;
    double yieldTot = xsTimesBR*lumi*totefficiencyVal[i];
    double yieldTag = xsTimesBR*lumi*dijetratioVal[i]*totefficiencyVal[i];
    double yieldUnt = xsTimesBR*lumi*untagAll->totalCtr;
    double yieldTotInMW = xsTimesBR*lumi*totefficiencyValInMW[i];
    double yieldTagInMW = xsTimesBR*lumi*dijetratioValInMW[i]*totefficiencyValInMW[i];
    double yieldUntInMW = xsTimesBR*lumi*untagInMW->totalCtr;

    
    (*txtYields) << left << setw(7) << fixed << setprecision(0) << masses[i] << setw(13) << fixed << setprecision(7) << xsTimesBR 
		 << setw(13) << fixed << setprecision(3) << untagInMW->totalCtr << setw(13) << yieldUntInMW << setw(13) << dijetratioValInMW[i]*totefficiencyValInMW[i] 
		 << setw(13) << yieldTagInMW << setw(13) << totefficiencyValInMW[i] << setw(13) << yieldTotInMW 
		 << setw(13) << fixed << setprecision(3) << untagAll->totalCtr << setw(13) << yieldUnt << setw(13) << dijetratioVal[i]*totefficiencyVal[i] 
		 << setw(13) << yieldTag << setw(13) << totefficiencyVal[i] << setw(13) << yieldTot 
		 << setw(13) << fixed << setprecision(0) << untagAll->numEventsRaw + dijetAll->numEventsRaw
		 << setw(13) << fixed << setprecision(2) << untagAll->numEventsRaw + dijetAll->numEventsPowheg
		 << setw(13) << untagAll->numEventsPU + dijetAll->numEventsPU
		 << setw(13) << untagAll->numEventsDataMC + dijetAll->numEventsDataMC
		 << endl;
     
  
    f->Close();
  }  

  (*txtYields) << endl << endl << endl;
	
  TGraphErrors* totgrEff;
  TGraphErrors* ratgrEff;

  if (process==ggH || process==qqH){
    totgrEff = new TGraphErrors( nPoints, mHVal, totefficiencyVal, 0, totefficiencyErr);
    ratgrEff = new TGraphErrors( nPoints, mHVal, dijetratioVal, 0, dijetratioErr);
  }
  else {
    totgrEff = new TGraphErrors( nPoints, mHVal, totefficiencyValInMW, 0, totefficiencyErrInMW);
    ratgrEff = new TGraphErrors( nPoints, mHVal, dijetratioValInMW, 0, dijetratioErrInMW);
  }
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

  if (channel==fs4mu && sqrts==7) {    
    polyFunctot->SetParLimits(7,0,0.035);
    polyFunctot->SetParLimits(8,160,210);
    polyFunctot->SetParLimits(9,30,50);
  }

  polyFunctot->SetLineColor(4);      
  TString cname = "eff" + sprocess + ssqrts + "_" + schannel;
  TCanvas *ctot = new TCanvas(cname,cname);
  ctot->SetGrid();

  TString outname = "sigFigs" + ssqrts +"/eff_" + sprocess + "_" + schannel + "_" + sjes;
  if (!useNewGGHPowheg) outname+="_oldPwg";

  totgrEff->Fit(polyFunctot,"Rt"); 
  TString xaxisText = "m_{" + schannel + "}";
  totgrEff->GetXaxis()->SetTitle(xaxisText);
  TString yaxisText = "Efficiency, " + sprocess + ", " + schannel;
  totgrEff->GetYaxis()->SetTitle(yaxisText);
  totgrEff->SetMinimum(0.0);
  totgrEff->SetMaximum(1.0);
  if (process>=3) totgrEff->SetMaximum(0.0035);
  totgrEff->Draw("AP");
  polyFunctot->Draw("sames");
  ctot->Print(outname+".eps");
  //ctot->Print(outname+".png"); // Does not work in batch?
  ctot->Print(outname+".pdf"); 
  //ctot->Print(outname+".root"); 
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


  // Create card fragments using new powheg samples
  string oftotprocess;
  if (process==ggH) oftotprocess="";
  else oftotprocess=sprocess;

  if (process==ggH) {
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
  if (!useNewGGHPowheg) outname+="_oldPwg";

  TF1 *ratiofit=0;
  if (process==ggH || process==qqH) ratiofit = new TF1("ratiofit","([0]+[1]*x+[2]*x*x)",110.,xMax);
  if (process==ZH || process==WH || process==ttH ) ratiofit = new TF1("ratiofit","([0]+[1]*x)",110.,xMax);

  ratgrEff->Fit(ratiofit,"Rt");
  ratgrEff->GetXaxis()->SetTitle(xaxisText);
  TString yaxisratio = "Dijet ratio, " + sprocess + ", " + schannel;
  ratgrEff->GetYaxis()->SetTitle(yaxisratio);
  ratgrEff->SetMinimum(0.0);
  ratgrEff->SetMaximum(1.0);
  ratgrEff->Draw("AP");
  crat->Print(outname+".eps");
  //crat->Print(outname+".png"); // Does not work in batch?
  crat->Print(outname+".pdf");
  //crat->Print(outname+".root");
  fratio->cd();
  ratgrEff->Write("Ratio");
  fratio->Close();
  
  cout << endl;
  cout << "------- Parameters for " << sprocess << " " << schannel << " sqrts=" << sqrts << endl;
  cout << "   a1 = " << ratiofit->GetParameter(0) << endl;
  cout << "   a2 = " << ratiofit->GetParameter(1) << endl;
  if (process==ggH || process==qqH) cout << "   a3 = " << ratiofit->GetParameter(2) << endl;
  cout << "---------------------------" << endl << endl;

  if (process==ggH) ofrat << "## signal efficiency ratios ##" << endl;
  ofrat << "signalEff tagged_" << sprocess << "_ratio " << ratiofit->GetParameter(0) << "+(" << ratiofit->GetParameter(1) << "*@0)";
  if (process==ggH || process==qqH) ofrat << "+(" << ratiofit->GetParameter(2) << "*@0*@0)" << endl;
  else if (process==ZH || process==WH ) ofrat << endl;
  else if (process==ttH) ofrat << endl << endl;
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


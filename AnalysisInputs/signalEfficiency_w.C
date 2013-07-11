
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

bool verbose = true;
bool useNewGGHPowheg = true;

using namespace std;
using namespace ROOT::Math;

//----------> SET INPUT VARIABLES in Config.h
#include "Config.h"
//<----------

TFile* ftot,*fratio;

void signalEfficiency_w(int channel, double sqrts, int process, double JES, ofstream* txtYields=0);

float WidthValue(float mHStarWidth);


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
  signalEfficiency_w(1,7,1,JES,&fileOutYields);
  signalEfficiency_w(2,7,1,JES,&fileOutYields);
  signalEfficiency_w(3,7,1,JES,&fileOutYields);
  if (!useNewGGHPowheg){
    //ggH (old powheg)
    signalEfficiency_w(1,8,1,JES,&fileOutYields);
    signalEfficiency_w(2,8,1,JES,&fileOutYields);
    signalEfficiency_w(3,8,1,JES,&fileOutYields);
  } else {
    //ggH (powheg15jhuGenV3)
    signalEfficiency_w(1,8,7,JES,&fileOutYields);
    signalEfficiency_w(2,8,7,JES,&fileOutYields);
    signalEfficiency_w(3,8,7,JES,&fileOutYields);
  }
  
//   //ggH (powheg15) // TEMPORARY for x-check
//   signalEfficiency_w(1,8,6,JES,&fileOutYields);
//   signalEfficiency_w(2,8,6,JES,&fileOutYields);
//   signalEfficiency_w(3,8,6,JES,&fileOutYields);


  //qqH
  signalEfficiency_w(1,7,2,JES,&fileOutYields);
  signalEfficiency_w(2,7,2,JES,&fileOutYields);
  signalEfficiency_w(3,7,2,JES,&fileOutYields);
  signalEfficiency_w(1,8,2,JES,&fileOutYields);
  signalEfficiency_w(2,8,2,JES,&fileOutYields);
  signalEfficiency_w(3,8,2,JES,&fileOutYields);
  //ZH
  signalEfficiency_w(1,7,3,JES,&fileOutYields);
  signalEfficiency_w(2,7,3,JES,&fileOutYields);
  signalEfficiency_w(3,7,3,JES,&fileOutYields);
  signalEfficiency_w(1,8,3,JES,&fileOutYields);
  signalEfficiency_w(2,8,3,JES,&fileOutYields);
  signalEfficiency_w(3,8,3,JES,&fileOutYields);
  //WH
  signalEfficiency_w(1,7,4,JES,&fileOutYields);
  signalEfficiency_w(2,7,4,JES,&fileOutYields);
  signalEfficiency_w(3,7,4,JES,&fileOutYields);
  signalEfficiency_w(1,8,4,JES,&fileOutYields);
  signalEfficiency_w(2,8,4,JES,&fileOutYields);
  signalEfficiency_w(3,8,4,JES,&fileOutYields);
  //ttH
  signalEfficiency_w(1,7,5,JES,&fileOutYields);
  signalEfficiency_w(2,7,5,JES,&fileOutYields);
  signalEfficiency_w(3,7,5,JES,&fileOutYields);
  signalEfficiency_w(1,8,5,JES,&fileOutYields);
  signalEfficiency_w(2,8,5,JES,&fileOutYields);
  signalEfficiency_w(3,8,5,JES,&fileOutYields);

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

  TString sprocess;
  if      (process == 1) sprocess = "ggH";
  else if (process == 2) sprocess = "qqH";
  else if (process == 3) sprocess = "ZH";
  else if (process == 4) sprocess = "WH";
  else if (process == 5) sprocess = "ttH";
  else if (process == 6) sprocess = "ggH_p15"; // TEMPORARY for x-check
  else if (process == 7) sprocess = "ggH_p15jhu"; // TEMPORARY for x-check
  else cout << "Not a valid channel: " << process << endl;

  TString sjes;
  if      (JES==0.) sjes="";
  else if (JES>0.)  sjes="_up";
  else if (JES<0.)  sjes="_down";

  TString ssqrts = (long) sqrts + TString("TeV");

  (*txtYields) << endl << endl 
	       << left << setw(10) << "*** Summary: " << sprocess << ", sqrts = " << fixed << setprecision(0) <<sqrts << " TeV ***" << endl << endl;
  (*txtYields) << left << setw(10) << "sqrts" << setw(13) << "Channel" << setw(16) << "Process" << setw(13) << "mH" 
	       << setw(13) << "Eff" << setw(13) << "XS*BR" << setw(13) << "Yield" << setw(13) 
	       << "unt-raw-TOT" << setw(13) << "unt-raw" << setw(13) << "unt-pwgW" 
	       << setw(13) << "unt-puW" << setw(13) << "unt-datamcW" << setw(13) 
	       << "tag-raw-TOT" << setw(13) << "tag-raw" << setw(13) << "tag-pwgW" 
	       << setw(13) << "tag-puW" << setw(13) << "tag-datamcW"
	       << endl << endl;

  cout << "process = " << sprocess << " schannel = " << schannel << "  sqrts = " << sqrts << " JES = " << JES <<endl;

  TString totoutfile = "CardFragments/signalEfficiency_"  + ssqrts + "_" + schannel + sjes + ".txt";
  TString ratiooutfile = "CardFragments/signalEfficiency_" + ssqrts + "_" + schannel + sjes + "_ratio.txt";

  // Create card fragments using new powheg samples
  ofstream oftot;
  if (process==1 || process==7) oftot.open(totoutfile,ios_base::out);
  if (process!=1 && process!=7) oftot.open(totoutfile,ios_base::out | ios_base::app);
  ofstream ofrat;
  if (process==1 || process==7) ofrat.open(ratiooutfile,ios_base::out);
  if (process!=1 && process!=7) ofrat.open(ratiooutfile,ios_base::out | ios_base::app);


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
  }  else if (process==6) {
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
  } else if (process==7) {
    if (sqrts==7) {
      cout << "No powheg15JhuGenV3 samples available at 7 TeV" << endl;
      exit(1);
    } else if (sqrts==8) {
      nPoints = nPoints8TeV_p15jhuGenV3;
      masses  = masses8TeV_p15jhuGenV3;
      mHVal   = mHVal8TeV_p15jhuGenV3;
      filepath = filePath8TeV;
    }
  }  

  // FIXME: need to skip ZH @ 7TeV for mH = 200 GeV
  if (process==3 && sqrts==7) nPoints = nPoints-1;


  float xMax = masses[nPoints-1]+10;


  const int arraySize=200;
  assert (arraySize>=nPoints);
  double totefficiencyVal[arraySize];
  double totefficiencyErr[arraySize];
  double dijetratioVal[arraySize];
  double dijetratioErr[arraySize];

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
    if (process==1 || process==2 || process>=6) {
      if (channel==1 || channel==2) BR = BRH4e;
      else BR = BRH2e2mu;
    }

    if (process==1)      xsTimesBR = BR*myCSW->HiggsCS(1,masses[i],sqrts);
    else if (process==2) xsTimesBR = BR*myCSW->HiggsCS(2,masses[i],sqrts);
    else if (process==3) xsTimesBR = BR*myCSW->HiggsCS(3,masses[i],sqrts);
    else if (process==4) xsTimesBR = BR*myCSW->HiggsCS(4,masses[i],sqrts);
    else if (process==5) xsTimesBR = BR*myCSW->HiggsCS(5,masses[i],sqrts);
    else if (process==6 || process==7) xsTimesBR = BR*myCSW->HiggsCS(1,masses[i],sqrts);


    if (process==1) infile = filepath+ "/" + (schannel=="2e2mu"?"2mu2e":schannel) + "/HZZ4lTree_H" + (long)masses[i] + ".root";
    else if (process==2) infile = filepath+ "/" + (schannel=="2e2mu"?"2mu2e":schannel) + "/HZZ4lTree_VBFH" + (long)masses[i] + ".root";
    else if (process==3 || process==4 || process==5) infile = filepath+ "/" + (schannel=="2e2mu"?"2mu2e":schannel) + "/HZZ4lTree_" + sprocess + (long)masses[i] + ".root";    
    else if (process==6) infile = filepath+ "/" + (schannel=="2e2mu"?"2mu2e":schannel) + "/HZZ4lTree_powheg15H" + (long)masses[i] + ".root";
    else if (process==7) infile = filepath+ "/" + (schannel=="2e2mu"?"2mu2e":schannel) + "/HZZ4lTree_powheg15" + (masses[i]>200?"H":"jhuGenV3H") + (long)masses[i] + ".root";
    TFile *f = TFile::Open(infile) ;
    TTree *t1 = (TTree*) f->Get("SelectedTree");
    float MC_weight_norm, MC_weight_PUWeight, MC_weight_powhegWeight,  MC_weight_dataMC;
    float MC_weight_noxsec;
    float GenHPt;
    //int NJets;
    short genProcessId=0;
    short NJets30;
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
    t1->SetBranchAddress("NJets30",&NJets30);
    t1->SetBranchAddress("GenHPt",&GenHPt);
    t1->SetBranchAddress("ZZMass",&ZZMass);

    //Initialize counters for non-dijet events
    int numndjEventsRaw = 0;
    int numndjEventsRawTot = 0;
    float numndjEventsPowheg =0;
    float numndjEventsPU =0;    
    float numndjEventsDataMC =0;
    float totalndjCtr=0;
    float ndjeff_noweight=0;
    float ndjsumw2=0;
    float ndjsumw_init2=0;
    //Initialize counters for dijet events
    int numdjEventsRaw = 0;
    int numdjEventsRawTot = 0;
    float numdjEventsPowheg =0;
    float numdjEventsPU =0;    
    float numdjEventsDataMC =0;
    float totaldjCtr=0;
    float djeff_noweight=0;
    float djsumw2=0;
    float djsumw_init2=0;


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
    //     if (process==1) {
    //       fW = TFile::Open(fW_str,"READ");
    //       h_HPtWeight = (TH1D*)fW->Get("h_weight");
    //     }


    for (int a = 0; a < t1->GetEntries(); a++){ 
      t1->GetEntry(a);
      if ((process==3 && genProcessId!=24) || (process==4 && genProcessId!=26) || (process==5 && (genProcessId!=121 && genProcessId!=122))) continue; 

      // We use the efficiency vs. generated events in the proper FS for ggH, VBF, and the efficiency vs all generated events for VH, ttH
      float effw = MC_weight_norm;
      if (process==3) {
	effw = MC_weight_noxsec*filter_eff_ZH_8TeV;
      }
      else if (process==4){
	effw = MC_weight_noxsec*filter_eff_WH_8TeV;
      }
      else if (process==5){
	effw = MC_weight_noxsec*filter_eff_ttH_8TeV;
      }      

//       double HPtWeight = 1.;
//       if (process==1) HPtWeight = h_HPtWeight->GetBinContent(h_HPtWeight->FindBin(GenHPt));
//       //cout << "Higgs pT weight = " << HPtWeight << endl;
//       effw*=HPtWeight;
      

      int NJets=0;
      double jetptc=0;
      for (unsigned int j=0; j<JetPt->size();j++){
	if (JES==0.) jetptc=JetPt->at(j);
	else if (JES!=0.) jetptc=JetPt->at(j)*(1+JES*JetSigma->at(j));
	if (jetptc>30.) NJets++;
      }

      if (NJets30!=NJets) cout << "ERROR: " << NJets30 << " " << NJets << endl;

      if(NJets30 < 2){
	// Increase counter of total events raw (counting also avents outside the mass window)
	++numndjEventsRawTot; 
	//Check if events falls inside the mass window (for VH)
	//cout << "lowM = " << low_M << ", highM = " << high_M << ", ZZMass = " << ZZMass << endl;
	if ( (process==3 || process==4) && (ZZMass<low_M || ZZMass>high_M) ) continue;
	//if ( (ZZMass<low_M || ZZMass>high_M) ) continue;
	//cout << "+++ lowM = " << low_M << ", highM = " << high_M << ", ZZMass = " << ZZMass << endl;
        totalndjCtr+=effw;      // Actual efficiency
	ndjsumw2 += effw*effw;  // square of weights, for error 
	// Numbers for x-check:
	++numndjEventsRaw;      // Number of raw events
	numndjEventsPowheg += MC_weight_powhegWeight; // Raw, hi-mass reweighted
	numndjEventsPU += MC_weight_powhegWeight*MC_weight_PUWeight; // Raw + hi-mass + PU rew
	numndjEventsDataMC += MC_weight_powhegWeight*MC_weight_PUWeight*MC_weight_dataMC; // Raw + hi-mass rew + PU +data/MC rew
	if (MC_weight_powhegWeight>0) {
	  float w_initial = effw/(MC_weight_powhegWeight*MC_weight_PUWeight*MC_weight_dataMC); //initial weight (1/nevt_gen)
	  ndjeff_noweight += w_initial; // efficiency, with no reweight
	  ndjsumw_init2 += w_initial*w_initial;
	}
      } else if(NJets > 1){
	// Increase counter of total events raw (counting also avents outside the mass window)
	++numdjEventsRawTot; 
	//Check if events falls inside the mass window (for VH)
	//cout << "lowM = " << low_M << ", highM = " << high_M << ", ZZMass = " << ZZMass << endl;
	if ( (process==3 || process==4) && (ZZMass<low_M || ZZMass>high_M) ) continue;
	//if ( (ZZMass<low_M || ZZMass>high_M) ) continue;
	//cout << "+++ lowM = " << low_M << ", highM = " << high_M << ", ZZMass = " << ZZMass << endl;
	totaldjCtr+=effw;      // Actual efficiency
	djsumw2 += effw*effw;  // square of weights, for error 
	// Numbers for x-check:
	++numdjEventsRaw;      // Number of raw events
	numdjEventsPowheg += MC_weight_powhegWeight; // Raw, hi-mass reweighted
	numdjEventsPU += MC_weight_powhegWeight*MC_weight_PUWeight; // Raw + hi-mass + PU rew
	numdjEventsDataMC += MC_weight_powhegWeight*MC_weight_PUWeight*MC_weight_dataMC; // Raw + hi-mass rew + PU +data/MC rew
	if (MC_weight_powhegWeight>0) {
	  float w_initial = effw/(MC_weight_powhegWeight*MC_weight_PUWeight*MC_weight_dataMC); //initial weight (1/nevt_gen)
	  djeff_noweight += w_initial; // efficiency, with no reweight
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
    cout << sprocess << " " << m << " " << totefficiencyVal[i]<<endl;
    totefficiencyErr[i] = sqrt(ndjsumw2 + djsumw2);
    dijetratioVal[i]=totaldjCtr/totefficiencyVal[i];
    dijetratioErr[i]=sqrt(pow(totalndjCtr,2)*djsumw2 + pow(totaldjCtr,2)*ndjsumw2)/pow(totefficiencyVal[i],2);
    
    // Write yields to output file
    double lumi = -1.;
    sqrts == 7 ? lumi = lumi7TeV*1000 : lumi = lumi8TeV*1000;
    double yield = xsTimesBR*lumi*totefficiencyVal[i];
    
    (*txtYields) << left << setw(10) << fixed << setprecision(0) << sqrts << setw(13) << schannel << setw(16) << sprocess << setw(13) << masses[i] 
		 << setw(13) << fixed << setprecision(4) << totefficiencyVal[i] << setw(13) << fixed << setprecision(5) 
		 << xsTimesBR << setw(13) << fixed << setprecision(5) << yield << setw(13) 
		 << fixed << setprecision(0) << setw(13) << numndjEventsRawTot << setw(13) << numndjEventsRaw << setw(13) 
		 << fixed << setprecision(2) << numndjEventsPowheg << setw(13) << numndjEventsPU << setw(13) 
		 << numndjEventsDataMC << setw(13) 
		 << fixed << setprecision(0) << setw(13) << numdjEventsRawTot << setw(13) << numdjEventsRaw << setw(13) 
		 << fixed << setprecision(2) << numdjEventsPowheg << setw(13) << numdjEventsPU << setw(13) 
		 << numdjEventsDataMC << setw(13) 
		 << endl;
    
  
    f->Close();
  }  

  (*txtYields) << endl << endl << endl;
	

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

  if (process!=6) totgrEff->Fit(polyFunctot,"Rt");
  TString xaxisText = "m_{" + schannel + "}";
  totgrEff->GetXaxis()->SetTitle(xaxisText);
  TString yaxisText = "Efficiency, " + sprocess + ", " + schannel;
  totgrEff->GetYaxis()->SetTitle(yaxisText);
  totgrEff->SetMinimum(0.0);
  totgrEff->SetMaximum(1.0);
  if (process>=3 && process!=6 && process!=7) totgrEff->SetMaximum(0.0035);
  totgrEff->Draw("AP");
  if (process!=6) polyFunctot->Draw("sames");
  ctot->Print(outname+".eps");
  //ctot->Print(outname+".png"); // Does not work in batch?
  ctot->Print(outname+".pdf"); 
  ftot->cd();
  totgrEff->Write("TotalEfficiency");
  ftot->Close();

  if (process!=6){
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
  if (process==1 || process==7) oftotprocess="";
  if (process!=1 && process!=7) oftotprocess=sprocess;

  if (process==1 || process==7) {
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
  }
  oftot.close();

  
  cname = "eff" + sprocess + ssqrts + "_" + schannel + "_ratio";
  TCanvas *crat = new TCanvas(cname,cname);
  crat->SetGrid();

  outname = "sigFigs" + ssqrts +"/eff_" + sprocess + "_" + schannel + "_" + sjes + "_ratio";

  TF1 *ratiofit=0;
  if (process==1 || process==2 || process==7) ratiofit = new TF1("ratiofit","([0]+[1]*x+[2]*x*x)",110.,xMax);
  if (process==3 || process==4 || process==5 ) ratiofit = new TF1("ratiofit","([0]+[1]*x)",110.,xMax);

  if (process!=6) ratgrEff->Fit(ratiofit,"Rt");
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
  
  if (process!=6){
  cout << endl;
  cout << "------- Parameters for " << sprocess << " " << schannel << " sqrts=" << sqrts << endl;
  cout << "   a1 = " << ratiofit->GetParameter(0) << endl;
  cout << "   a2 = " << ratiofit->GetParameter(1) << endl;
  if (process==1 || process==7 || process==2) cout << "   a3 = " << ratiofit->GetParameter(2) << endl;
  cout << "---------------------------" << endl << endl;

  if (process==1 || process==7) ofrat << "## signal efficiency ratios ##" << endl;
  ofrat << "signalEff tagged_" << (process==7?"ggH":sprocess) << "_ratio " << ratiofit->GetParameter(0) << "+(" << ratiofit->GetParameter(1) << "*@0)";
  if (process==1 || process==7 || process==2) ofrat << "+(" << ratiofit->GetParameter(2) << "*@0*@0)" << endl;
  else if (process==3 || process==4 ) ofrat << endl;
  else if (process==5) ofrat << endl << endl;
  }
  ofrat.close();

  // deviations
  cout << "Deviations..." << endl;
  double maxResidual=0;
  if (process!=6){
  for (int i = 0; i < nPoints; i++){
    double eval = polyFunctot->Eval(masses[i]);
    double residual = (eval - totefficiencyVal[i]);
    maxResidual = max(maxResidual,fabs(residual));
    if (verbose)    cout << "For mass, " << masses[i] << ": measured value is " << totefficiencyVal[i] << " and difference from function is " << residual <<endl;
  }
  }
  cout << "Largest residual= " << maxResidual << endl;
}


float WidthValue(float mHStarWidth)
{
  ostringstream MassString;
  MassString << mHStarWidth;
  
  ifstream widthFile("widthvalues.txt");
  
  string line;

  bool FindedMass = false;

  float Gamma_ggCal, Gamma_ZZCal, Gamma_TOTCal;
  
  while (getline(widthFile,line)) {
    if( line == "" || line[0] == '#' ) continue;
    
    if(line[0]== MassString.str()[0] && line[1]== MassString.str()[1] && line[2]== MassString.str()[2]){
      
      stringstream stringline;
      stringline << line;    
      string masschar;
      stringline>>masschar>>Gamma_ggCal>>Gamma_ZZCal>>Gamma_TOTCal;
      
      FindedMass = true;
    }
  }
  if(!FindedMass) abort();

  return Gamma_TOTCal;
}

/* 
 * Compute ZZ background rates and write them in a card fragment.
 * usage: 
 * -set all input variables in Config.h
 * -run with:
 * root -q -b ZZbackgroundRate.C+
 * This runs on 7 and 8 TeV and writes the output in a file (see stdout).
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
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include "TTree.h"
#include "TText.h"
#include "TROOT.h"
#include "TStyle.h"

//----------> SET INPUT VARIABLES in Config.h
#include "Config.h"
//<----------

using namespace std;
void compute(TString name, int sqrts, double lumi, bool VBFtag = false);
double sumWeights(TString name, double lumi, bool selVBF = false);


// Run both sqrts in one go
void ZZbackgroundRate() {

  compute(filePath7TeV, 7, lumi7TeV, true);
  compute(filePath8TeV, 8, lumi8TeV, true);

  compute(filePath7TeV, 7, lumi7TeV, false);
  compute(filePath8TeV, 8, lumi8TeV, false);

}


// The actual job
void compute(TString filePath, int sqrts, double lumi, bool VBFtag){

  double qqZZ[3];
  double ggZZ[3];
  
  cout<<"qqZZ 4mu"<<endl;  
  double qqZZ_4mu_ev4mu = sumWeights(filePath + "/4mu/HZZ4lTree_ZZTo4mu.root",lumi,VBFtag);
  double qqZZ_4mu_ev4e = sumWeights(filePath + "/4mu/HZZ4lTree_ZZTo4e.root",lumi,VBFtag);
  double qqZZ_4mu_ev2e2mu = sumWeights(filePath + "/4mu/HZZ4lTree_ZZTo2e2mu.root",lumi,VBFtag);
  double qqZZ_4mu_ev4tau = sumWeights(filePath + "/4mu/HZZ4lTree_ZZTo4tau.root",lumi,VBFtag);
  double qqZZ_4mu_ev2mu2tau = sumWeights(filePath + "/4mu/HZZ4lTree_ZZTo2mu2tau.root",lumi,VBFtag);
  double qqZZ_4mu_ev2e2tau = sumWeights(filePath + "/4mu/HZZ4lTree_ZZTo2e2tau.root",lumi,VBFtag);
  qqZZ[0]=qqZZ_4mu_ev4mu+qqZZ_4mu_ev4e+qqZZ_4mu_ev2e2mu+qqZZ_4mu_ev4tau+qqZZ_4mu_ev2mu2tau+qqZZ_4mu_ev2e2tau;
  cout<<qqZZ_4mu_ev4mu<<" "<<qqZZ_4mu_ev4e<<" "<<qqZZ_4mu_ev2e2mu<<" "<<qqZZ_4mu_ev4tau<<" "<<qqZZ_4mu_ev2mu2tau<<" "<<qqZZ_4mu_ev2e2tau
      <<" -> "<<qqZZ_4mu_ev4mu<<" + "
      <<qqZZ_4mu_ev4e+qqZZ_4mu_ev2e2mu+qqZZ_4mu_ev4tau+qqZZ_4mu_ev2mu2tau+qqZZ_4mu_ev2e2tau
      <<" -> "<< qqZZ[0]<< endl;
  

  cout<<"qqZZ 4e"<<endl;
  double qqZZ_4e_ev4e = sumWeights(filePath + "/4e/HZZ4lTree_ZZTo4e.root",lumi,VBFtag);
  double qqZZ_4e_ev4mu = sumWeights(filePath + "/4e/HZZ4lTree_ZZTo4mu.root",lumi,VBFtag);
  double qqZZ_4e_ev2e2mu = sumWeights(filePath + "/4e/HZZ4lTree_ZZTo2e2mu.root",lumi,VBFtag);
  double qqZZ_4e_ev4tau = sumWeights(filePath + "/4e/HZZ4lTree_ZZTo4tau.root",lumi,VBFtag);
  double qqZZ_4e_ev2mu2tau = sumWeights(filePath + "/4e/HZZ4lTree_ZZTo2mu2tau.root",lumi,VBFtag);
  double qqZZ_4e_ev2e2tau = sumWeights(filePath + "/4e/HZZ4lTree_ZZTo2e2tau.root",lumi,VBFtag);
  qqZZ[1]=qqZZ_4e_ev4e+qqZZ_4e_ev4mu+qqZZ_4e_ev2e2mu+qqZZ_4e_ev4tau+qqZZ_4e_ev2mu2tau+qqZZ_4e_ev2e2tau;
  cout<<qqZZ_4e_ev4e<<" "<<qqZZ_4e_ev4mu<<" "<<qqZZ_4e_ev2e2mu<<" "<<qqZZ_4e_ev4tau<<" "<<qqZZ_4e_ev2mu2tau<<" "<<qqZZ_4e_ev2e2tau
      <<" -> "<<qqZZ_4e_ev4e<<" + "
      <<qqZZ_4e_ev4mu+qqZZ_4e_ev2e2mu+qqZZ_4e_ev4tau+qqZZ_4e_ev2mu2tau+qqZZ_4e_ev2e2tau
      <<" -> "<< qqZZ[1] <<endl;

  cout<<"qqZZ 2e2mu"<<endl;
  double qqZZ_2e2mu_ev2e2mu = sumWeights(filePath + "/2mu2e/HZZ4lTree_ZZTo2e2mu.root",lumi,VBFtag);
  double qqZZ_2e2mu_ev4mu = sumWeights(filePath + "/2mu2e/HZZ4lTree_ZZTo4mu.root",lumi,VBFtag);
  double qqZZ_2e2mu_ev4e = sumWeights(filePath + "/2mu2e/HZZ4lTree_ZZTo4e.root",lumi,VBFtag);
  double qqZZ_2e2mu_ev4tau = sumWeights(filePath + "/2mu2e/HZZ4lTree_ZZTo4tau.root",lumi,VBFtag);
  double qqZZ_2e2mu_ev2mu2tau = sumWeights(filePath + "/2mu2e/HZZ4lTree_ZZTo2mu2tau.root",lumi,VBFtag);
  double qqZZ_2e2mu_ev2e2tau = sumWeights(filePath + "/2mu2e/HZZ4lTree_ZZTo2e2tau.root",lumi,VBFtag);
  qqZZ[2] = qqZZ_2e2mu_ev2e2mu+qqZZ_2e2mu_ev4e+qqZZ_2e2mu_ev4mu+qqZZ_2e2mu_ev4tau+qqZZ_2e2mu_ev2mu2tau+qqZZ_2e2mu_ev2e2tau; 
  cout<<qqZZ_2e2mu_ev2e2mu<<" "<<qqZZ_2e2mu_ev4e<<" "<<qqZZ_2e2mu_ev4mu<<" "<<qqZZ_2e2mu_ev4tau<<" "<<qqZZ_2e2mu_ev2mu2tau<<" "<<qqZZ_2e2mu_ev2e2tau
      <<" -> "<<qqZZ_2e2mu_ev2e2mu<<" + "
      <<qqZZ_2e2mu_ev4e+qqZZ_2e2mu_ev4mu+qqZZ_2e2mu_ev4tau+qqZZ_2e2mu_ev2mu2tau+qqZZ_2e2mu_ev2e2tau
      <<" -> "<<qqZZ[2]<<endl;
  
  cout<<"ggZZ 4mu"<<endl;
  double ggZZ_4mu_ev2l2l = sumWeights(filePath + "/4mu/HZZ4lTree_ggZZ2l2l.root",lumi,VBFtag);
  double ggZZ_4mu_ev4l = sumWeights(filePath + "/4mu/HZZ4lTree_ggZZ4l.root",lumi,VBFtag);
  ggZZ[0] = ggZZ_4mu_ev4l + ggZZ_4mu_ev2l2l;
  cout<<ggZZ_4mu_ev4l<<" + "<<ggZZ_4mu_ev2l2l<<" -> "<< ggZZ[0] <<endl;

  cout<<"ggZZ 4e"<<endl;
  double ggZZ_4e_ev2l2l = sumWeights(filePath + "/4e/HZZ4lTree_ggZZ2l2l.root",lumi,VBFtag);
  double ggZZ_4e_ev4l = sumWeights(filePath + "/4e/HZZ4lTree_ggZZ4l.root",lumi,VBFtag);
  ggZZ[1] = ggZZ_4e_ev4l + ggZZ_4e_ev2l2l;
  cout<<ggZZ_4e_ev4l<<" + "<<ggZZ_4e_ev2l2l<<" -> "<< ggZZ[1] <<endl;

  cout<<"ggZZ 2e2mu"<<endl;
  double ggZZ_2e2mu_ev2l2l = sumWeights(filePath + "/2mu2e/HZZ4lTree_ggZZ2l2l.root",lumi,VBFtag);
  double ggZZ_2e2mu_ev4l = sumWeights(filePath + "/2mu2e/HZZ4lTree_ggZZ4l.root",lumi,VBFtag);
  ggZZ[2] = ggZZ_2e2mu_ev4l + ggZZ_2e2mu_ev2l2l;
  cout<<ggZZ_2e2mu_ev4l<<" + "<<ggZZ_2e2mu_ev2l2l<<" -> "<< ggZZ[2] <<endl;

  TString schannel[3] = {"4mu","4e","2e2mu"};
  TString ssqrts = (long) sqrts + TString("TeV");
  for (int i=0; i<3; ++i) {
    TString outfile = "CardFragments/ZZRates_" + ssqrts + "_" + schannel[i] + "_" + Form("%d",int(VBFtag)) + ".txt";
    ofstream of(outfile,ios_base::out);
    of << "## rates --- format = chan N lumi ##" << endl
       << "## if lumi is blank, lumi for cards used ##" << endl;
    of << "rate qqZZ  " << qqZZ[i] << endl;
    of << "rate ggZZ  " << ggZZ[i] << endl;
    of.close();
    cout << "Output written to: " << outfile << endl;
  }
}


double sumWeights(TString name, double lumi, bool selVBF){
  TChain* tree = new TChain("SelectedTree");
  tree->Add((TString)(name));

  float MC_weight,ZZMass,mela,Z2Mass,Z1Mass;
  int NJets;
  //std::vector<double> *JetEta = new vector<double>;
  //std::vector<double> *JetPt = new vector<double>;
  //std::vector<double> *JetPhi = new vector<double>;
  //std::vector<double> *JetMass = new vector<double>;
  tree->SetBranchAddress("MC_weight",&MC_weight);  
  tree->SetBranchAddress("ZZMass",&ZZMass);  
  tree->SetBranchAddress("Z1Mass",&Z1Mass);  
  tree->SetBranchAddress("Z2Mass",&Z2Mass);  
  tree->SetBranchAddress("ZZLD",&mela);
  tree->SetBranchAddress("NJets",&NJets);
  //tree->SetBranchAddress("JetPt", &JetPt);
  //tree->SetBranchAddress("JetEta", &JetEta);
  //tree->SetBranchAddress("JetMass", &JetMass); 

  double totEvents=0;
  //Get Histos for LD
  for(int iEvt=0; iEvt<tree->GetEntries(); iEvt++){
    tree->GetEntry(iEvt);
    if( (selVBF == true && NJets > 1) || (selVBF == false && NJets < 2))
      {
	totEvents=totEvents+MC_weight;
      }
  }
  return totEvents*lumi;
}




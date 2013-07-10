/* 

plots MELA distribution for Zjets MC, data CR, qqZZ MC.  
compare dataCR and qqZZ, fit with line

to run:

root -l -n -b
.L ~/tdrstyle.C
setTDRStyle()
gStyle->SetOptStat(0);
.L ZXsystematics.C+
runSystematics()

*/


#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TF1.h"
#include "TTree.h"
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

//--------------------
// global variables

//TString discriminantName="ZZLD";
TString inputDir="root://lxcms02//data/Higgs/rootuplesOut/130702b/";
//TString inputDir="root://lxcms02//data/Higgs/rootuplesKD/310812/";

bool useSlopeError=true;

void setInputDir(TString name){

  inputDir=name;

}

/*void setDiscrimName(char* name){

  discriminantName=name;

  }*/

//=======================================================================
bool test_bit(int mask, unsigned int iBit){return (mask >>iBit) &1; }

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// function below plots qqZZ, data CR, MC CR shapes, calculates
// the ratio and fits the ratio with a line
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

pair<double,double> measureSystematics(double mzzLow=130, double mzzHigh=180){

  pair<double,double> results;
  results.first=-99;
  results.second=-99;

  // ----------
  // load trees 
  // ----------
  
  TChain* CRdata = new TChain("SelectedTree");
  CRdata->Add(inputDir+"/PRODFSR_8TeV/CR/HZZ4lTree_DoubleOr_CRZLLTree.root");
  CRdata->Add(inputDir+"/PRODFSR/CR/HZZ4lTree_Double*_CRZLLTree.root");


  TChain* CRmc = new TChain("SelectedTree");
  CRmc->Add(inputDir+"/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2*_CRZLLTree.root");
  CRmc->Add(inputDir+"/PRODFSR/CR/HZZ4lTree_DYJetsToLLTuneZ2*_CRZLLTree.root");


  TChain* qqZZ = new TChain("SelectedTree");
  qqZZ->Add(inputDir+"/PRODFSR/4mu/HZZ4lTree_ZZTo4mu.root");
  qqZZ->Add(inputDir+"/PRODFSR/4e/HZZ4lTree_ZZTo4e.root");
  qqZZ->Add(inputDir+"/PRODFSR/2mu2e/HZZ4lTree_ZZTo2e2mu.root");

  qqZZ->Add(inputDir+"/PRODFSR_8TeV/4mu/HZZ4lTree_ZZTo4mu.root");
  qqZZ->Add(inputDir+"/PRODFSR_8TeV/4e/HZZ4lTree_ZZTo4e.root");
  qqZZ->Add(inputDir+"/PRODFSR_8TeV/2mu2e/HZZ4lTree_ZZTo2e2mu.root");


  if( !CRmc || CRmc->GetEntries()<=0 ){
    cout << "problem loading CRmc files... " << endl;
    return results;
  }
  if( !CRdata || CRdata->GetEntries()<=0 ){
    cout << "problem loading CRdata files... " << endl;
    return results;
  }
  if( !qqZZ || qqZZ->GetEntries()<=0 ){
    cout << "problem loading qqZZ files... " << endl;
    return results;
  }

  // = = = = = = =
  // set branches 
  // = = = = = = =

  float mzz, w;
  float p0plus_VAJHU,bkg_VAMCFM;
  int CRflag;

  CRmc->SetBranchAddress("ZZMass",&mzz);
  //CRmc->SetBranchAddress(discriminantName,&D);
  CRmc->SetBranchAddress("p0plus_VAJHU",&p0plus_VAJHU);
  CRmc->SetBranchAddress("bkg_VAMCFM",&bkg_VAMCFM);
  CRmc->SetBranchAddress("MC_weight_noxsec",&w);
  CRmc->SetBranchAddress("CRflag",&CRflag);

  CRdata->SetBranchAddress("ZZMass",&mzz);
  //CRdata->SetBranchAddress(discriminantName,&D);
  CRdata->SetBranchAddress("p0plus_VAJHU",&p0plus_VAJHU);
  CRdata->SetBranchAddress("bkg_VAMCFM",&bkg_VAMCFM);
  CRdata->SetBranchAddress("CRflag",&CRflag);

  
  qqZZ->SetBranchAddress("ZZMass",&mzz);
  //qqZZ->SetBranchAddress(discriminantName,&D);
  qqZZ->SetBranchAddress("p0plus_VAJHU",&p0plus_VAJHU);
  qqZZ->SetBranchAddress("bkg_VAMCFM",&bkg_VAMCFM);
  qqZZ->SetBranchAddress("MC_weight_noxsec",&w);

  // = = = = = = = = = = = 
  // initialize histograms 
  // = = = = = = = = = = = 

  TH1F* h_CRmc = new TH1F("h_CRmc",";MELA;",30,0,1);
  h_CRmc->Sumw2();
  TH1F* h_CRdata = new TH1F("h_CRdata",";MELA;",30,0,1);
  h_CRdata->Sumw2();
  TH1F* h_qqZZ = new TH1F("h_qqZZ",";MELA;",30,0,1);
  h_qqZZ->Sumw2();

  TCanvas* can = new TCanvas("can","can",400,550);
  TPad* pad2 = new TPad("pad2","pad2",0.,0.,1.,0.3);
  TPad* pad1 = new TPad("pad1","pad1",0.,0.3,1.,1.);
  pad2->SetBottomMargin(0.22);
  
  pad1->Draw();
  pad2->Draw();

  // = = = = = =
  // fill histos
  // = = = = = = 


  for(int iEvt=0; iEvt<qqZZ->GetEntries(); iEvt++){
    
    qqZZ->GetEntry(iEvt);

    if(mzz>mzzLow&&mzz<mzzHigh){
      h_qqZZ->Fill(p0plus_VAJHU/(p0plus_VAJHU + bkg_VAMCFM),w);
    }
    
  }

  for(int iEvt=0; iEvt<CRmc->GetEntries(); iEvt++){
    
    CRmc->GetEntry(iEvt);

    if(mzz>mzzLow&&mzz<mzzHigh){
      if(test_bit(CRflag,5) || test_bit(CRflag,7) || test_bit(CRflag,9) || test_bit(CRflag,11)) 
	{
	  h_CRmc->Fill(p0plus_VAJHU/(p0plus_VAJHU + bkg_VAMCFM),w);
	}
    }
  }

  for(int iEvt=0; iEvt<CRdata->GetEntries(); iEvt++){
    
    CRdata->GetEntry(iEvt);

    if(mzz>mzzLow&&mzz<mzzHigh){
      if(test_bit(CRflag,5) || test_bit(CRflag,7) || test_bit(CRflag,9) || test_bit(CRflag,11)) 
	{
	  h_CRdata->Fill(p0plus_VAJHU/(p0plus_VAJHU + bkg_VAMCFM));
	}
    }
    
  }
  
  // = = = = = = 
  // draw histos
  // = = = = = = 

  h_qqZZ->Scale(h_CRdata->Integral()/h_qqZZ->Integral());
  h_qqZZ->SetLineColor(4);
  h_qqZZ->SetLineStyle(2);  
  h_qqZZ->SetLineWidth(2);

  h_CRmc->Scale(h_CRdata->Integral()/h_CRmc->Integral());
  h_CRmc->SetLineColor(2);
  h_CRmc->SetLineWidth(2);

  h_CRdata->SetMarkerStyle(8);

  pad1->cd();

  if(h_CRdata->GetMaximum()+h_CRdata->GetBinError(h_CRdata->GetMaximumBin()) > h_CRmc->GetMaximum()+h_CRmc->GetBinError(h_CRmc->GetMaximumBin())){
    h_CRdata->GetYaxis()->SetRangeUser(0,(h_CRdata->GetMaximum()+h_CRdata->GetBinError(h_CRdata->GetMaximumBin()))*1.5);
    h_CRdata->Draw("ep");
    h_CRmc->Draw("EhistSAME");
    h_qqZZ->Draw("EhistSAME");
  }else{
    h_CRmc->GetYaxis()->SetRangeUser(0,(h_CRmc->GetMaximum()+h_CRmc->GetBinError(h_CRmc->GetMaximumBin()))*1.5);
    h_CRmc->Draw("Ehist");
    h_CRdata->Draw("SAMEep");
    h_qqZZ->Draw("EhistSAME");
  }

  // ---------- LEGEND ---------------

  TLegend* leg = new TLegend(.5,.6,.90,.90);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h_CRdata,"Z+X (data)","lp");
  leg->AddEntry(h_CRmc,"Z+X (Madgraph MC)","l");
  leg->AddEntry(h_qqZZ,"qqZZ (Powheg MC)","l");
  
  leg->Draw();

  // ------------ RATIO PAD ---------------

  pad2->cd();
  
  TH1F* ratio_qqZZ   = new TH1F(*h_qqZZ);
  ratio_qqZZ->Divide(h_qqZZ);
  TH1F* ratio_CRmc   = new TH1F(*h_CRmc);
  ratio_CRmc->Divide(h_qqZZ);
  TH1F* ratio_CRdata = new TH1F(*h_CRdata);
  ratio_CRdata->Divide(h_qqZZ);

  ratio_CRdata->GetYaxis()->SetRangeUser(0.,2.);
  ratio_CRdata->Draw("p");

  // ------------- FIT RATIO -------------
  cout << "===================================================================" << endl;
  cout << "====================Fitting data CR " << mzzLow << "-" << mzzHigh << "========================" << endl;

  TF1* fline = new TF1("fline","[0]+[1]*x",mzzLow,mzzHigh);
  fline->SetLineWidth(2);
  fline->SetLineColor(kGreen+1);
  ratio_CRdata->Fit("fline","");
  //gStyle->SetOptFit(0);

  ratio_CRmc->Draw("EhistSAME");
  //ratio_qqZZ->Draw("EhistSAME");

  TF1* fline_plusError = new TF1("fline_plusError","[0]+[1]*x",mzzLow,mzzHigh);
  fline_plusError->SetLineWidth(2);
  fline_plusError->SetLineColor(kMagenta+1);
  fline_plusError->SetParameters(fline->GetParameter(0),fline->GetParameter(1)+
		       ( fline->GetParameter(1)/fabs(fline->GetParameter(1)) )*fline->GetParError(1));
  fline_plusError->Draw("SAME");
  // -------------------------------------

  results.first=fline->GetParameter(0);
  if(useSlopeError)
    results.second=fline->GetParameter(1)+
      ( fline->GetParameter(1)/fabs(fline->GetParameter(1)) )*fline->GetParError(1);
  else
    results.second=fline->GetParameter(1);
  //char temp[250];
  
  TString temp;

  //sprintf(temp,"Z+X_data_vs_MC_vs_qqZZ_%i-%i_%s.eps",(int)mzzLow,(int)mzzHigh,"VAKD");
  temp="Z+X_data_vs_MC_vs_qqZZ_";
  temp+=(int)mzzLow;
  temp+="-";
  temp+=(int)mzzHigh;
  temp+="_";
  temp+="VAKD";
  temp+=".eps";
  can->SaveAs(temp);
  //sprintf(temp,"Z+X_data_vs_MC_vs_qqZZ_%i-%i_%s.png",(int)mzzLow,(int)mzzHigh,"VAKD");
  temp="Z+X_data_vs_MC_vs_qqZZ_";
  temp+=(int)mzzLow;
  temp+="-";
  temp+=(int)mzzHigh;
  temp+="_";
  temp+="VAKD";
  temp+=".png";
  can->SaveAs(temp);

  delete h_CRmc;
  delete h_CRdata;
  delete h_qqZZ;
  
  delete pad1;
  delete pad2;
  delete can;
  
  delete CRdata;
  delete CRmc;
  delete qqZZ;

  return results;

}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// measure systematic effect in course binned mZZ windows
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

void runSystematics(){

  vector<double> low, high, slope, yInt;
  pair<double,double> fitResults;

  fitResults=measureSystematics(100,120);
  low.push_back(100);  high.push_back(120);
  slope.push_back(fitResults.second); yInt.push_back(fitResults.first);
  fitResults=measureSystematics(120,140);
  low.push_back(120);  high.push_back(140);
  slope.push_back(fitResults.second); yInt.push_back(fitResults.first);
  fitResults=measureSystematics(140,160);
  low.push_back(140);  high.push_back(160);
  slope.push_back(fitResults.second); yInt.push_back(fitResults.first);
  fitResults=measureSystematics(160,180);
  low.push_back(160);  high.push_back(180);
  slope.push_back(fitResults.second); yInt.push_back(fitResults.first);
  fitResults=measureSystematics(180,220);
  low.push_back(180);  high.push_back(220);
  slope.push_back(fitResults.second); yInt.push_back(fitResults.first);
  fitResults=measureSystematics(220,260);
  low.push_back(220);  high.push_back(260);
  slope.push_back(fitResults.second); yInt.push_back(fitResults.first);
  fitResults=measureSystematics(260,300);
  low.push_back(260);  high.push_back(300);
  slope.push_back(fitResults.second); yInt.push_back(fitResults.first);

  cout << "double low[" << low.size() << "]={";
  for(unsigned int i=0 ; i<low.size(); i++){
    if(i+1==low.size())
      cout << low[i] << "};" << endl;
    else
      cout << low[i] << ",";
  }

  cout << "double high[" << high.size() << "]={";
  for(unsigned int i=0 ; i<high.size(); i++){
    if(i+1==high.size())
      cout << high[i] << "};" << endl;
    else
      cout << high[i] << ",";
  }

  cout << "double slope[" << slope.size() << "]={";
  for(unsigned int i=0 ; i<slope.size(); i++){
    if(i+1==slope.size())
      cout << slope[i] << "};" << endl;
    else
      cout << slope[i] << ",";
  }

  cout << "double yInt[" << yInt.size() << "]={";
  for(unsigned int i=0 ; i<yInt.size(); i++){
    if(i+1==yInt.size())
      cout << yInt[i] << "};" << endl;
    else
      cout << yInt[i] << ",";
  }

}


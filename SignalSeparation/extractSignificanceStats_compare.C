#include <Riostream.h>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#include "TArrow.h"
#include "TGaxis.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TArrow.h"
#include "TColor.h"
#include "TROOT.h"
#include "Math/DistFunc.h"

void setTDRStyle();
void run(int);

TH1F* hPS_A;
TH1F* hSM_A;
TH1F* hObs_A;
double obsQ_A;
TH1F* hPS_B;
TH1F* hSM_B;
TH1F* hObs_B;
double obsQ_B;
TString inputA;
TString inputB;
bool unblind;

int extractSignificanceStats(TString legALT="0^{-}", TString nameALT="0m", TString input1="qmu_*", TString input2="/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/UFLresults/mu_float/sepExample_ggSpin0M_qmu.root"){

  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  setTDRStyle();

  inputA=input1;
  inputB=input2;
  unblind = true;

  cout<<endl<<"RUNNING FRAMEWROK A"<<endl;
  run(1);
  cout<<endl<<"RUNNING FRAMEWROK B"<<endl;
  run(2);
  cout<<endl<<"PLOTTING"<<endl;
  const float lumi7TeV=5.051;
  const float lumi8TeV=12.21;

  //Plotting
  gStyle->SetOptStat(0);
  TCanvas *c1=new TCanvas("c","c",500,500);
  c1->cd();
  hSM_A->Rebin(50);
  hPS_A->Rebin(50);
  hSM_A->SetXTitle(" -2 #times ln(L_{JP}/L_{0+})");
  hSM_A->SetYTitle("Pseudoexperiments");
  hPS_A->SetXTitle(" -2 #times ln(L_{JP}/L_{0+})");
  hPS_A->SetYTitle("Pseudoexperiments");
  hSM_A->SetLineColor(kRed+2);
  hSM_A->SetLineStyle(2);
  hSM_A->SetFillColor(798);
  hSM_A->SetLineWidth(2);
  hPS_A->SetFillColor(kAzure+7);
  hPS_A->SetLineColor(kBlue);
  hPS_A->SetLineWidth(1);
  hPS_A->SetFillStyle(3001);

  hObs_A->SetLineColor(kRed);
  hObs_A->SetLineWidth(2);

  TGraph *grObs=new TGraph();//dummy, just for the legend
  grObs->SetLineColor(kRed);
  grObs->SetLineWidth(1);
  
  hSM_A->GetXaxis()->SetRangeUser(-30.0,30.0);
  hSM_A->GetXaxis()->SetLabelFont(42);
  hSM_A->GetXaxis()->SetLabelOffset(0.007);
  hSM_A->GetXaxis()->SetLabelSize(0.045);
  hSM_A->GetXaxis()->SetTitleSize(0.05);
  hSM_A->GetXaxis()->SetTitleOffset(1.15);
  hSM_A->GetXaxis()->SetTitleFont(42);
  hSM_A->GetYaxis()->SetLabelFont(42);
  hSM_A->GetYaxis()->SetLabelOffset(0.007);
  hSM_A->GetYaxis()->SetLabelSize(0.045);
  hSM_A->GetYaxis()->SetTitleSize(0.05);
  hSM_A->GetYaxis()->SetTitleOffset(1.8);
  hSM_A->GetYaxis()->SetTitleFont(42); 
  //TGaxis::SetMaxDigits(2); 
  hSM_A->Scale(1./hSM_A->Integral("width"));
  hPS_A->Scale(1./hPS_A->Integral("width"));
  float maxhSM_A=hSM_A->GetBinContent(hSM_A->GetMaximumBin());
  float maxhPS_A=hPS_A->GetBinContent(hPS_A->GetMaximumBin());
  if(maxhPS_A>maxhSM_A){
    hSM_A->SetMaximum(maxhPS_A*1.3);
    hPS_A->SetMaximum(maxhPS_A*1.3);
  }
  else{
    hSM_A->SetMaximum(maxhSM_A*1.3);
    hPS_A->SetMaximum(maxhSM_A*1.3);
  }


  hSM_A->Draw();
  hPS_A->Draw("sames");

  TArrow *obsArrow=0;
  if(unblind)obsArrow=new TArrow(obsQ_A,hSM_A->GetMaximum()/5.0,obsQ_A,0.0,0.05,"|->");
  else obsArrow=new TArrow(0.0,hSM_A->GetMaximum()/5.0,0.0,0.0,0.05,"|->");
  obsArrow->SetLineColor(kRed);
  obsArrow->SetLineWidth(4.0);
  if(unblind)  obsArrow->Draw();

  //TLegend *leg = new TLegend(0.63,0.73,0.92,0.93);
  TLegend *leg = new TLegend(0.63,0.73,0.88,0.93);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);

  leg->AddEntry(hSM_A, "0^{+} framework A","f");
  leg->AddEntry(hPS_A, legALT+" framework A","f");
  if(unblind) leg->AddEntry(hObs_A,"Observed framework A","L");
  //  if(unblind) leg->AddEntry(hObs,"  Simulated data","L");

  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->SetLeftMargin(0.18);
  c1->SetRightMargin(0.05);
  c1->SetTopMargin(0.05);
  c1->SetBottomMargin(0.15);
  c1->SetFrameFillStyle(0);
  c1->SetFrameBorderMode(0);
  c1->SetFrameFillStyle(0);
  c1->SetFrameBorderMode(0);

  TPaveText *pt = new TPaveText(0.1577181,0.9562937,0.9580537,0.9947552,"brNDC");
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetFillStyle(0);
  pt->SetTextFont(42);
  pt->SetTextSize(0.03);
  TText *text = pt->AddText(0.01,0.5,"CMS preliminary");
  text = pt->AddText(0.3,0.6,"#sqrt{s} = 7 TeV, L = 5.1 fb^{-1}  #sqrt{s} = 8 TeV, L = 19.6 fb^{-1}");
  pt->Draw();   

  TArrow* UFLArrow=new TArrow(obsQ_B,hSM_A->GetMaximum()/5.0,obsQ_B,0.0,0.05,"|->");
  hSM_B->Rebin(50);
  hPS_B->Rebin(50);

  hPS_B->SetLineColor(kBlack);
  hPS_B->SetFillStyle(0);
  hSM_B->SetLineColor(kBlack);
  hSM_B->SetLineStyle(2);
  hSM_B->SetFillStyle(0);
  UFLArrow->SetLineColor(kGreen);
  UFLArrow->SetLineWidth(4);
  hObs_B->SetLineColor(kGreen);
  hObs_B->SetLineWidth(2);
  c1->cd();
  hPS_B->Scale(1./hPS_B->Integral("width"));
  hSM_B->Scale(1./hSM_B->Integral("width"));
  hPS_B->Draw("same");
  hSM_B->Draw("same");
  UFLArrow->Draw("same");

  TH1F* hObs2=new TH1F();
  hObs2->SetLineColor(kGreen);
  hObs2->SetLineWidth(3);
  leg->AddEntry(hSM_B, "0^{+} framework B","l");
  leg->AddEntry(hPS_B, legALT+" framework B","l");
  if(unblind) leg->AddEntry(hObs_B,"Observed framework B","L");
  leg->Draw();

  c1->SaveAs("sigsep_compare_"+nameALT+".eps");
  c1->SaveAs("sigsep_compare_"+nameALT+".png");
  c1->SaveAs("sigsep_compare_"+nameALT+".root");
  c1->SaveAs("sigsep_compare_"+nameALT+".C");

  return 0;
}//end main

void run(int framework){
 TChain* t = new TChain("q");
 if(framework==1){
   t->Add(inputA);
 }
 else if(framework==2){
   t->Add(inputB);
 }
 else{
   cout<<"wrong framework! Aborting..."<<endl;
   abort();
 }
  float q,m,w;
  int type;
  t->SetBranchAddress("q",&q);
  t->SetBranchAddress("mh",&m);
  t->SetBranchAddress("weight",&w);
  t->SetBranchAddress("type",&type);


  TH1F *hSM=new TH1F("hSM;S = -2 #times ln(L_{PS}/L_{SM});Number of Toys","",8000,-40,40);
  TH1F *hPS=new TH1F("hPS;S = -2 #times ln(L_{PS}/L_{SM});Number of Toys","",8000,-40,40);
  TH1F *hObs=new TH1F("hObserved","",8000,-40,40);

  std::vector<float> v_SM, v_PS,v_Obs;

  for(int i=0;i<t->GetEntries();i++){
    t->GetEntry(i);

    //if(i==0)cout<<"MASS in the TREE = "<<m<<endl<<endl;

    q*=2.0;
    if(type<0){ //SM hypothesis 
      hSM->Fill(-q);
      v_SM.push_back(-q);
    }
    else if(type>0){//ALT hypothesis
      hPS->Fill(-q);
      v_PS.push_back(-q);
    }
    else{
      hObs->Fill(-q);
      //cout<<"DATA -q ="<<(-q)<<endl;
      v_Obs.push_back(-q);
    }
  }//end loop on tree entries
  //cout<<"Finished to loop, sorting vectors "<<v_SM.size()<<" "<<v_PS.size()<<" "<<v_Obs.size()<<endl;
  sort(v_SM.begin(),v_SM.end());//sort in ascending order
  sort(v_PS.begin(),v_PS.end()); 
  sort(v_Obs.begin(),v_Obs.end());
  int ntoysSM= hSM->GetEntries();
  int ntoysPS= hPS->GetEntries();

  //we assume that SM is on the right and PS on the left of zero
  if(v_PS.at(0)>v_SM.at(ntoysSM-1)){
    cout<<"Swapped distributions !!! The alternative model shouldstay on the negative side of the significance."<<endl;
    cout<<"Please edit the code and change the sign of q when filling histos and vectors in the loop on tree entries"<<endl;
  }

  if((int(v_SM.size())!= ntoysSM)||(int(v_PS.size())!= ntoysPS)){
    cout<<"Mismatch in size of vectors and #entries of historgams ! v_SM.size()="<< v_SM.size() <<"  ntoysSM="<<ntoysSM<<endl;
  }

  float medianSM=v_SM.at(int(ntoysSM/2));
  float medianPS=v_PS.at(int(ntoysPS/2));
  cout<<"Toys generated "<<ntoysSM<<"\t"<<ntoysPS<<endl;
  //cout<<"Mean of SM/PS hypothesis: "<<hSM->GetMean()<<"\t"<<hPS->GetMean()<<endl;
  //cout<<"RMS  of SM/PS hypothesis: "<<hSM->GetRMS()<<"\t"<<hPS->GetRMS()<<endl;
  //cout<<"Median of SM/PS hypothesis: "<<medianSM<<"\t"<<medianPS<<endl;

  const float step=0.05;
  float coverage=0.0;
  float diff=10.0;
  float cut=v_PS.at(0)-step;
  float crosspoint=-99.0;
  int startSM=ntoysSM-1, startPS=0;
  //cout<<"Starting to loop with cut at "<<cut<<endl;

  while(cut<=v_SM.at(ntoysSM-1)+step){
    float cutSM=-1.0,cutPS=-1.0;
    for(int iSM=startSM;iSM>=0;iSM--){      
      if(v_SM.at(iSM)<cut){//gotcha
	cutSM=ntoysSM-iSM;
	break;
      }
    }

    for(int iPS=startPS;iPS<ntoysPS;iPS++){
      if(v_PS.at(iPS)>cut){//gotcha
	cutPS=iPS;
	break;
      }
    }

    if(cutSM>=0&&cutPS>=0){
      float fracSM=(ntoysSM-cutSM)/ntoysSM;
      float fracPS=(ntoysPS-cutPS)/ntoysPS;
      if(fabs(fracSM-fracPS)<diff){
	diff=fabs(fracSM-fracPS);
	coverage=fabs(fracSM+fracPS)/2.0;
	crosspoint=cut;
      }
    }
    cut+=step;
  }//end while loop
 
  //cout<<"Finished loop on vector elements, min diff is "<<diff<<", looped until cut_fin="<<cut<<endl;
  //cout<<"q value where SM and ALT distributions have same area on opposite sides: "<<crosspoint<<"  Coverage="<<coverage<<endl;
  float separation=2*ROOT::Math::normal_quantile_c(1.0 - coverage, 1.0);
  //cout<<"Separation from tail prob: "<<separation<<endl<<endl<<endl;
  

  float integralSM=hSM->Integral();
  float integralPS=hPS->Integral();
 
  float tailSM=hSM->Integral(1,hSM->FindBin(medianPS))/integralSM;
  float tailPS=hPS->Integral(hPS->FindBin(medianSM),hPS->GetNbinsX())/integralPS;

  diff=10.0;
  coverage=0.0;
  for(int i=1;i<hSM->GetNbinsX();i++){
    
    float fracSM=hSM->Integral(1,i) / integralSM;
    float fracPS=hPS->Integral(i,hPS->GetNbinsX()) / integralPS;
    if(fabs(fracSM-fracPS)<diff){
      diff=fabs(fracSM-fracPS);
      coverage=(fracSM+fracPS)/2.0;
    }

  }

  float sepH= 2*ROOT::Math::normal_quantile_c(1.0 - coverage, 1.0);
  cout<<"Separation from histograms = "<<sepH<<" with coverage "<<coverage<<endl;

  if(unblind){
    if(v_Obs.size()<1){
      cout<<"Ooops ! The size of the vector with the observed separation is not 1 but "<<v_Obs.size()<<" ! I am not going to plot the observed results."<<endl;
      unblind=false;
    }
    else{
      cout<<"q value in DATA:"<<v_Obs.at(0)<<endl;
      float obsTailSM=hSM->Integral(1,hSM->FindBin(v_Obs.at(0)))/integralSM;
      float obsTailPS=hPS->Integral(hPS->FindBin(v_Obs.at(0)),hPS->GetNbinsX())/integralPS;
      cout<<"P(SM < Obs): "<<obsTailSM<<"  ("<<ROOT::Math::normal_quantile_c(obsTailSM,1.0) <<" sigma)"<<endl;
      cout<<"P(PS > Obs): "<<obsTailPS<<"  ("<<ROOT::Math::normal_quantile_c(obsTailPS,1.0) <<" sigma)"<<endl;

      float obsCLsRatio = obsTailPS / (1.0 - obsTailSM);
      cout<<"CLs criterion P(PS > Obs) / P(SM > Obs) : "<<obsCLsRatio<<endl;//"  ("<<ROOT::Math::normal_quantile_c(obsCLsRatio,1.0) <<" sigma)"<<endl;
    }
  }//end if unblinding
  if(framework==1){
    hPS_A = (TH1F*) hPS->Clone("hPS_A");
    hSM_A = (TH1F*) hSM->Clone("hSM_A");
    obsQ_A= v_Obs.at(0);
    hObs_A= (TH1F*) hObs->Clone("hObs_A");
  }
  else if(framework==2){
    hPS_B = (TH1F*) hPS->Clone("hPS_B");
    hSM_B = (TH1F*) hSM->Clone("hSM_B");
    obsQ_B= v_Obs.at(0);
    hObs_B= (TH1F*) hObs->Clone("hObs_B");
  }
  else{
   cout<<"wrong framework! Aborting..."<<endl;
   abort();
  }
}

void plotAll(){
  //prefitmu
  cout<<"    #########  PREFIT MU, 0- ##########"<<endl;
  //extractSignificanceStats( "0^{-}", "prefitMu_0m", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/CMSSW_6_1_1/src/andrewFinalResults/cards_0-_8TeV_prefitMu/HCG/126/qmu_*.root", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/UFLFinalResults/mu_float/sepExample_ggSpin0M_2_qmu.root");
  //OLD UFL RESULTS
  extractSignificanceStats("0^{-}", "prefitMu_0m", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/CMSSW_6_1_1/src/andrewFinalResults/cards_0-_8TeV_prefitMu/HCG/126/qmu_*.root", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/UFLresults/mu_float/sepExample_ggSpin0M_qmu.root");

  cout<<"    #########  PREFIT MU, 0h+ ##########"<<endl;
  extractSignificanceStats( "0^{+}_{h}", "prefitMu_0hp", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/CMSSW_6_1_1/src/andrewFinalResults/cards_0h+_8TeV_prefitMu/HCG/126/qmu_*.root", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/UFLFinalResults/mu_float/sepExample_ggSpin0Ph_2_qmu.root");

  cout<<"    #########  PREFIT MU, 1+ ##########"<<endl;
  extractSignificanceStats( "1^{+}", "prefitMu_1p", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/CMSSW_6_1_1/src/andrewFinalResults/cards_1+_8TeV_prefitMu/HCG/126/qmu_*.root", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/UFLFinalResults/mu_float/sepExample_qqSpin1P_2_qmu.root");
  cout<<"    #########  PREFIT MU, 1- ##########"<<endl;
  extractSignificanceStats( "1^{-}", "prefitMu_1m", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/CMSSW_6_1_1/src/andrewFinalResults/cards_1-_8TeV_prefitMu/HCG/126/qmu_*.root", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/UFLFinalResults/mu_float/sepExample_qqSpin1M_2_qmu.root");
  
  cout<<"    #########  PREFIT MU, gg 2+m ##########"<<endl;
  //extractSignificanceStats( "2^{+}_{m}(gg)", "prefitMu_gg2mp", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/CMSSW_6_1_1/src/andrewFinalResults/cards_gg2m+_8TeV_prefitMu/HCG/126/qmu_*.root", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/UFLFinalResults/mu_float/sepExample_ggSpin2Pm_2_qmu.root"); //OLD UFL RESULTS
  extractSignificanceStats("2^{+}_{m}(gg)", "prefitMu_gg2mp", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/CMSSW_6_1_1/src/andrewFinalResults/cards_gg2m+_8TeV_prefitMu/HCG/126/qmu_*.root", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/UFLresults/mu_float/sepExample_ggSpin2Pm_qmu.root");

  cout<<"    #########  PREFIT MU, qq 2+m ##########"<<endl;
  extractSignificanceStats( "2^{+}_{m}(qq)", "prefitMu_qq2mp", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/CMSSW_6_1_1/src/andrewFinalResults/cards_qq2m+_8TeV_prefitMu/HCG/126/qmu_*.root", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/UFLFinalResults/mu_float/sepExample_qqSpin2Pm_2_qmu.root");

  //mu=1
  cout<<"    #########  MU=1, 0- ##########"<<endl;
  extractSignificanceStats( "0^{-}", "fixedMu_0m", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/CMSSW_6_1_1/src/andrewFinalResults/cards_0-_8TeV_fixedMu/HCG/126/qmu_*.root", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/UFLFinalResults/mu1/sepExample_ggSpin0M_qmu.root");
  cout<<"    #########  MU=1, 0h+ ##########"<<endl;
  extractSignificanceStats( "0^{+}_{h}", "fixedMu_0hp", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/CMSSW_6_1_1/src/andrewFinalResults/cards_0h+_8TeV_fixedMu/HCG/126/qmu_*.root", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/UFLFinalResults/mu1/sepExample_ggSpin0Ph_qmu.root");

  cout<<"    #########  MU=1, 1+ ##########"<<endl;
  extractSignificanceStats( "1^{+}", "fixedMu_1p", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/CMSSW_6_1_1/src/andrewFinalResults/cards_1+_8TeV_fixedMu/HCG/126/qmu_*.root", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/UFLFinalResults/mu1/sepExample_qqSpin1P_qmu.root");
  cout<<"    #########  MU=1, 1- ##########"<<endl;
  extractSignificanceStats( "1^{-}", "fixedMu_1m", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/CMSSW_6_1_1/src/andrewFinalResults/cards_1-_8TeV_fixedMu/HCG/126/qmu_*.root", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/UFLFinalResults/mu1/sepExample_qqSpin1M_qmu.root");

  cout<<"    #########  MU=1, gg 2+m ##########"<<endl;
  extractSignificanceStats( "2^{+}_{m}(gg)", "fixedMu_gg2mp", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/CMSSW_6_1_1/src/andrewFinalResults/cards_gg2m+_8TeV_fixedMu/HCG/126/qmu_*.root", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/UFLFinalResults/mu1/sepExample_ggSpin2Pm_qmu.root");
  cout<<"    #########  MU=1, qq 2+m ##########"<<endl;
  extractSignificanceStats( "2^{+}_{m}(qq)", "fixedMu_qq2mp", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/CMSSW_6_1_1/src/andrewFinalResults/cards_qq2m+_8TeV_fixedMu/HCG/126/qmu_*.root", "/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/UFLFinalResults/mu1/sepExample_qqSpin2Pm_qmu.root");
  

}

void setTDRStyle() {
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(750); //Height of canvas
  tdrStyle->SetCanvasDefW(1050); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  //  tdrStyle->SetErrorMarker(20);
  tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);

//For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.010);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.25);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.14);
  tdrStyle->SetPadRightMargin(0.04);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.005);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.25);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  //tdrStyle->SetNdivisions(505, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

// Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  Int_t n=35;
  Int_t *colors = new Int_t[n];
  for (Int_t i =0; i < n; i++) {
    //colors[i] = i+61;
    //colors[i] = 70-i;
    colors[i] = i+63;
    }
  /*colors[0]=61
  colors[1]=62;
  colors[2]=63;
  colors[3]=;
  colors[4]=;
  colors[5]=;
  colors[6]=;
  colors[7]=;
  colors[8]=;
  colors[9]=;
  colors[10]=;*/
  gStyle->SetPalette(n, colors);
  
  //gStyle->SetPalette(1, 0);
  tdrStyle->cd();
}
/*
void damnedPlot(){
TCanvas *c1 =  new TCanvas("paddy","paddy",700,700)
  c1->Divide(1,2,0,0);
 c1->cd(1).SetBottomMargin(0.001);
 c1->cd(1).SetTopMargin(0.01);
 c1->cd(1).SetRightMargin(0.01);
 c1.cd(2);
 c1.cd(2).SetTopMargin(0.001);
 c1.cd(2).SetRightMargin(0.01);

}
*/

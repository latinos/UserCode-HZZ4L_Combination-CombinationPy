#include <Riostream.h>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

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
#include "Math/DistFunc.h"

#include "tdrstyle.C"


double Lumi_7TeV = 5.1;
double Lumi_8TeV = 19.6;


int extractSignificanceStats(TString jobname=""){

  double Hist_Range_min=-35;
  double Hist_Range_max=35;

  setTDRStyle( true );

  char fileName[128];
  sprintf(fileName,"qmu.root");
  TFile *fq=new TFile(fileName,"READ");
  TTree *t=(TTree*)fq->Get("q");

  float q,m,w;
  int type;
  t->SetBranchAddress("q",&q);
  t->SetBranchAddress("mh",&m);
  t->SetBranchAddress("weight",&w);
  t->SetBranchAddress("type",&type);

  TH1F *hSM=new TH1F("hSM","",8000,Hist_Range_min,Hist_Range_max);
  TH1F *hPS=new TH1F("hPS","",8000,Hist_Range_min,Hist_Range_max);
  TH1F *hObs=new TH1F("hObserved","",1000,Hist_Range_min,Hist_Range_max);
  cout<<"Start to loop on tree in file "<<fileName<<endl;

  std::vector<float> v_SM, v_PS,v_Obs;

  for(int i=0;i<t->GetEntries();i++){
    t->GetEntry(i);
    if(i==0)cout<<"MASS in the TREE = "<<m<<endl<<endl;

    if(type<0){ //STD hypothesis
      hSM->Fill( -2*q );
      v_SM.push_back( -2*q );
    }
    else if(type>0){//ALT hypothesis (-> PS)
      hPS->Fill( -2*q );
      v_PS.push_back( -2*q );
    }
    else{
      hObs->Fill( -2*q );
      v_Obs.push_back( -2*q );
    }
  }//end loop on tree entries
  cout << "Actual data input: " << v_Obs[0] << endl;

  cout<<"Finished to loop, sorting vectors "<<v_SM.size()<<" "<<v_PS.size()<<" "<<v_Obs.size()<<endl;
  sort(v_SM.begin(),v_SM.end());//sort in ascending order
  sort(v_PS.begin(),v_PS.end()); 
  sort(v_Obs.begin(),v_Obs.end());
  int ntoysSM= hSM->GetEntries();
  int ntoysPS= hPS->GetEntries();

  //we assume that SM is on the right and PS on the left of zero
  if(v_PS.at(0)>v_SM.at(ntoysSM-1)){
    cout<<"Swapped distributions !!! The alternative model shouldstay on the negative side of the significance."<<endl;
    cout<<"Please edit the code and change the sign of q when filling histos and vectors in the loop on tree entries"<<endl;
    return 1;
  }
  float medianSM=v_SM.at(int(ntoysSM/2));
  float medianPS=v_PS.at(int(ntoysPS/2));
  cout<<"Toys generated "<<ntoysSM<<"\t"<<ntoysPS<<endl;
  cout<<"Mean of SM/PS hypothesis: "<<hSM->GetMean()<<"   /   "<<hPS->GetMean()<<endl;
  cout<<"RMS  of SM/PS hypothesis: "<<hSM->GetRMS()<<"   /   "<<hPS->GetRMS()<<endl;
  cout<<"Median of SM/PS hypothesis: "<<medianSM<<"   /   "<<medianPS<<endl;

  const float step=0.05;
  float coverage=0.0;
  float diff=10.0;
  float cut=v_PS.at(0)-step;
  float crosspoint=-99.0;
  int startSM=ntoysSM-1, startPS=0;

  // Possible scale location 1
//  hSM->Scale( 1/hSM->Integral() );
//  hPS->Scale( 1/hPS->Integral() );

  float integralSM=hSM->Integral();
  float integralPS=hPS->Integral();

  cout<<integralSM<<endl; 
  cout<<integralPS<<endl; 

  float tailSM=hSM->Integral(1,hSM->FindBin(medianPS))/integralSM;
  float tailPS=hPS->Integral(hPS->FindBin(medianSM),hPS->GetNbinsX())/integralPS;
  cout<<tailSM<<endl;
  cout<<tailPS<<endl;
  if(tailSM>0)cout<<"Prob( q < median(P) | S ) = "<<tailSM<<"  ("<<ROOT::Math::normal_quantile_c(tailSM,1.0) <<" sigma)"<<endl;
  if(tailPS>0)cout<<"Prob( q > median(S) | P ) = "<<tailPS<<"  ("<<ROOT::Math::normal_quantile_c(tailPS,1.0) <<" sigma)"<<endl;
  
  // Data point probability
  double Probability_Data_SM=hSM->Integral(0,hSM->FindBin(hObs->GetBinCenter( hObs->FindFirstBinAbove(0.99) ))-1)/integralSM;
  double Probability_Data_PS=hPS->Integral(0,hPS->FindBin(hObs->GetBinCenter( hObs->FindFirstBinAbove(0.99) ))-1)/integralPS;
  cout << "Data point at " << hObs->GetBinCenter( hObs->FindFirstBinAbove(0.99) ) << endl;

  if( Probability_Data_SM >= 1 || Probability_Data_SM <0 )
  {
    TF1 func( "f1","gaus",Hist_Range_min,Hist_Range_max );
    hSM->Fit( "f1", "QNO" );	//constant, mean, sigma
    Probability_Data_SM = ROOT::Math::normal_cdf( (hObs->GetBinCenter( hObs->FindFirstBinAbove(0.99) )-func.GetParameter(1))/func.GetParameter(2), 1, 0 );
  }
  if( Probability_Data_PS >= 1 || Probability_Data_PS <0 )
  {
    TF1 func( "f1","gaus",Hist_Range_min,Hist_Range_max );
    hPS->Fit( "f1", "QNO" );	//constant, mean, sigma
    Probability_Data_PS = ROOT::Math::normal_cdf( (hObs->GetBinCenter( hObs->FindFirstBinAbove(0.99) )-func.GetParameter(1))/func.GetParameter(2), 1, 0 );
  }

  cout << "Prob( q >= qObs | q~SM Higgs+Bkg ) = " << static_cast<double>(1-Probability_Data_SM) << "  (";
  if( Probability_Data_SM<1 ) cout << ROOT::Math::normal_quantile_c( 1.0 - Probability_Data_SM, 1.0 ) << " sigma)\n";
  else cout << "--- sigma)\n";
  cout << "Prob( q >= qObs | q~Alt. sig. + Bkg ) = " << static_cast<double>(1-Probability_Data_PS) << "  (";
  if( Probability_Data_PS<1 ) cout << ROOT::Math::normal_quantile_c( 1.0 - Probability_Data_PS, 1.0 ) << " sigma)\n";
  else cout << "--- sigma)\n";
  cout << "CLs = " << static_cast<double>((1-Probability_Data_PS)/(1-Probability_Data_SM)) << "  (";
  if( (1-Probability_Data_PS)/(1-Probability_Data_SM)>0 ) cout << ROOT::Math::normal_quantile_c( (1-Probability_Data_PS)/(1-Probability_Data_SM), 1.0 ) << " sigma)\n";
  else cout << "--- sigma)\n";
  
  
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

  // Possible scale location 2
  hSM->Scale( 1/hSM->Integral() );
  hPS->Scale( 1/hPS->Integral() );

//  hSM->Scale( hPS->Integral()/hSM->Integral() );

  //Fancy plot
  Int_t ci;   // for color index setting
  gStyle->SetOptStat(0);
  TCanvas *c1=new TCanvas("c1","c1",800,820);

   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetLeftMargin(0.18);
   c1->SetRightMargin(0.05);
   c1->SetTopMargin(0.05);
   c1->SetBottomMargin(0.15);
   c1->SetFrameFillStyle(0);
   c1->SetFrameBorderMode(0);
   c1->SetFrameFillStyle(0);
   c1->SetFrameBorderMode(0);

  c1->cd();
  hSM->Rebin(50);
  hPS->Rebin(50);

  if( (hSM->GetMaximum()) > (hPS->GetMaximum()) ) hSM->SetMaximum( (hSM->GetMaximum())*1.2  );
  else hSM->SetMaximum( (hPS->GetMaximum())*1.15 );

  hSM->SetXTitle("-2 #times ln(L_{JP}/L_{0+})");
  hSM->SetYTitle("Pseudoexperiments");
  //hPS->SetXTitle("Q = -2 #times ln(L_{SM+bkg}/L_{0^{-}+bkg})");
//   hSM->SetLineColor(kRed+1);
//   hSM->SetFillColor(kRed+1);
  hSM->SetLineWidth(2);
//   hSM->SetFillStyle(3003);	//3645
//   hPS->SetLineColor(kBlue-7);
  //hPS->SetFillColor(kBlue-7);
//   hPS->SetLineWidth(2);
  //hPS->SetFillStyle(3654);

  ci = TColor::GetColor("#ff0000");
  hObs->SetLineColor( ci );
  hObs->SetFillColor( ci );
  hObs->SetLineWidth( 2 );

  // PRL style
  ci = TColor::GetColor("#ffcc33");
  hSM->SetFillColor( ci );

  ci = TColor::GetColor("#990000");
  hSM->SetLineColor(ci);
  hSM->SetLineStyle(2);
  hSM->SetLineWidth(2);
  hSM->SetMarkerStyle(20);

  ci = TColor::GetColor("#0099ff");
  hPS->SetFillColor(ci);
  hPS->SetFillStyle(3001);

  ci = TColor::GetColor("#0000ff");
  hPS->SetLineColor(ci);
  hPS->SetLineStyle(0);
  hPS->SetMarkerStyle(20);

  hSM->Draw();
  hSM->GetXaxis()->SetTitleOffset( 1.4 );
  hSM->GetXaxis()->SetTitleSize( 0.04 );
  hSM->GetYaxis()->SetTitleOffset( 2.3 );
  hSM->GetYaxis()->SetTitleSize( 0.04 );
  hPS->Draw("sames");

  
  hObs->Scale( 200 );
//  hObs->Draw("sames");

  TLegend *leg = new TLegend(0.63,0.73,0.88,0.93,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hSM," 0^{+}_{m}","f");
  leg->AddEntry(hPS," 0^{+}_{h}","f");
  leg->AddEntry(hObs, " #frac{1}{2} #times CMS data","l");
  leg->Draw();


  TArrow Data_arrow( hObs->GetBinCenter( hObs->FindFirstBinAbove(0.99) ),0.2*hSM->GetMaximum(),hObs->GetBinCenter( hObs->FindFirstBinAbove(0.99) ),0.0, 0.05,"|->" );
  ci = TColor::GetColor("#ff0000");
//  Data_arrow.SetAngle(30);
  Data_arrow.SetLineWidth( 2 );
  Data_arrow.SetLineColor( ci );
  Data_arrow.SetFillColor( 1 );
  Data_arrow.SetFillStyle( 1001 );
  Data_arrow.Draw();


  TLatex *CP = new TLatex();
  CP->SetNDC(kTRUE);
  CP->SetTextSize(0.03);
  CP->SetTextAlign(31);
  CP->SetTextFont(42);
  CP->SetTextAlign(11);

  CP->DrawLatex( 0.145, 0.965, Form("CMS Preliminary #sqrt{s} = 7 TeV, L = %.1f fb^{-1}  #sqrt{s} = 8 TeV, L = %.1f fb^{-1}",
    Lumi_7TeV, Lumi_8TeV) );
//  CP->DrawLatex( 0.483, 0.925, Form("Cross-check results with alternative KD") );
//  CP->DrawLatex( 0.58, 0.925, Form("Alternative analysis with MEKD") );
  
  TLatex *Separation = new TLatex();
  Separation->SetNDC(kTRUE);
  Separation->SetTextSize(0.04);
  Separation->SetTextAlign(31);
  Separation->SetTextFont(42);
  Separation->SetTextAlign(11);

//   Separation->DrawLatex( 0.2, 0.86, Form("Expected separation power: %.2f#sigma", fabs(sepH)) );


//  TPaveText pt(0.16,0.93,0.8,0.99,"NDC");
//  pt.SetFillColor(0);
//  TString stmp; stmp.Form("Exp separation power = %.2f #sigma", fabs(sepH));
//  pt.AddText(stmp);
//  TPaveText pt2(0.55,0.90,0.99,0.94,"NDC");
//  pt2.SetFillColor(0);
//  pt2.AddText("#sqrt{s} = 7 TeV, L = 5.05 fb^{-1},  #sqrt{s} = 8 TeV, L = 12.21 fb^{-1}");
//  pt2.SetBorderSize(0);
//  pt.Draw();
//  pt2.Draw();

  c1->SaveAs("sigsep_"+jobname+".eps");
  c1->SaveAs("sigsep_"+jobname+".png");
  c1->SaveAs("sigsep_"+jobname+".root");

  return 0;
}//end main


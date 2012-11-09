#include "TH1F.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <vector>
#include "TCutG.h"
#include "TFile.h"
#include "TH2.h"
#include "TPad.h"

void getPvals(TFile *f, std::vector<double> &v_mh, std::vector<double> &v_obs);

// --------- Inputs ------- //
TString inputFile = "resultstemp/higgsCombineHZZ4L_PLP.root";
TString inputFile_1 = "results261012_1Dbranch05/higgsCombineHZZ4L_PLP.root";
TString inputFile_2 = "results261012_2Dbranch05/higgsCombineHZZ4L_PLP.root";//"results261012_2DnewZshape/higgsCombineHZZ4L_PLP.root";
TString inputFileExp ="resultstemp/higgsCombineHZZ4L_PLPE.root";
TString inputFileExp_1 = "results261012_1Dbranch05/higgsCombineHZZ4L_PLPE.root";
TString inputFileExp_2 = "results261012_2Dbranch05/higgsCombineHZZ4L_PLPE.root";//"results261012_2DnewZshape/higgsCombineHZZ4L_PLPE.root";
TString Graph1 = "2D Fit 7+8TeV";
TString Graph2 = "Observed (1D)";
TString Graph3 = "Observed";
const bool addObs = false;
const bool addObs_1 = true;
const bool addObs_2 = true;
const bool addExpected = false;
const bool addExpected_1 = true;
const bool addExpected_2 = true;
string method = "PLP";
Double_t xLow = 99.9;
Double_t xHigh = 1001.0;
Double_t yLow = 1e-7;
Double_t yHigh = 1.0;
TString xTitle = "m_{H} [GeV]";
TString yTitle = "local p-value";
const bool logy = true;
const bool logx = true;
const bool grid = false;
const bool gridOnTop = false;
const bool points = true;
const bool points_1 = true;
const bool points_2 = true;
const bool isTiny = false;
int canvasX = 700;
int canvasY = 700;
const bool _DEBUG_ = false;
string plotDir = "plots";
string dimension = "2D1D_no2l2tau";
// ----------------------- //


using namespace std;

void plot_Signif_3lines()
{

  gROOT->ProcessLine(".x tdrstyle.cc");
  gStyle->SetPadLeftMargin(0.15);//016
  gStyle->SetPadTopMargin(0.05);

  TFile *inFile;
  if(addObs) inFile = new TFile(inputFile,"READ");
  TFile *inFileExp;
  if(addExpected)inFileExp = new TFile(inputFileExp,"READ");

  if(!inFile && addObs){ cout << "Cannot find file " << inputFile << end; return;}
  if(addExpected && !inFileExp){cout << "Cannot find file " << inputFileExp << end; return;}


  TFile *inFile_1;
  if(addObs_1) inFile_1 = new TFile(inputFile_1,"READ");
  TFile *inFileExp_1;
  if(addExpected_1)inFileExp_1 = new TFile(inputFileExp_1,"READ");

  if(!inFile_1 && addObs_1){ cout << "Cannot find file " << inputFile_1 << end;return;}
  if(addExpected_1 && !inFileExp_1){cout <<"Cannot find file " << inputFileExp_1 << end; return;}


  TFile *inFile_2;
  if(addObs_2) inFile_2 = new TFile(inputFile_2,"READ");
  TFile *inFileExp_2;
  if(addExpected_2)inFileExp_2 = new TFile(inputFileExp_2,"READ");
  
  if(!inFile_2 && addObs_2){ cout << "Cannot find file " << inputFile_2 << end;return;}
  if(addExpected_2 && !inFileExp_2){cout <<"Cannot find file " << inputFileExp_2 << end; return;}
  
  
  // ------------------- Get Values -------------------- //
  
  vector<double> mH, Val_obs;
  vector<double> v_masses, v_obs;
  if(addObs){  
    getPvals(inFile,mH,Val_obs);
    for(unsigned int i = 0; i < mH.size(); i++)
      {
	v_masses.push_back( mH[i] );
	v_obs.push_back(Val_obs[i]);
	cout << "obs " <<  mH[i] << "  " << Val_obs[i] << "  ( "<< ROOT::Math::gaussian_quantile_c(Val_obs[i],1) << ")" << endl;
      }
  }
  
  
  vector<double> mH_1, Val_obs_1;
  vector<double> v_masses_1, v_obs_1;
  if(addObs_1){
    getPvals(inFile_1,mH_1,Val_obs_1);    
    for(unsigned int i = 0; i < mH_1.size(); i++)
      {
	v_masses_1.push_back( mH_1[i] );
	v_obs_1.push_back(Val_obs_1[i]);
	cout << "obs_1 " <<  mH_1[i] << "  " << Val_obs_1[i] << "  ( "<< ROOT::Math::gaussian_quantile_c(Val_obs_1[i],1) << ")" << endl;
      }
  }
  vector<double> mH_2, Val_obs_2;
  vector<double> v_masses_2, v_obs_2;
  if(addObs_2){
    getPvals(inFile_2,mH_2,Val_obs_2);
    
    for(unsigned int i = 0; i < mH_2.size(); i++)
      {
	v_masses_2.push_back( mH_2[i] );
	v_obs_2.push_back(Val_obs_2[i]);
	cout << "obs_2 " << mH_2[i] << "  " << Val_obs_2[i] << "  ( "<< ROOT::Math::gaussian_quantile_c(Val_obs_2[i],1) << ")" << endl;
      }
  }

  // ------------------ Get Expected -------------------- //
  vector<double> mH_exp, Val_exp;  
  vector<double> v_masses_exp, v_exp;
  if(addExpected)
    {
      getPvals(inFileExp,mH_exp,Val_exp);
      for(unsigned int i = 0; i < mH_exp.size(); i++)
      	{
	  v_masses_exp.push_back( mH_exp[i] );
	  v_exp.push_back(Val_exp[i]);
	  cout << "exp " << mH_exp[i] << "  " << Val_exp[i] << "  ( "<< ROOT::Math::gaussian_quantile_c(Val_exp[i],1) << ")" << endl;
	}
    }


  vector<double> mH_exp_1, Val_exp_1;  
  vector<double> v_masses_exp_1, v_exp_1;
  if(addExpected_1)
    {
      getPvals(inFileExp_1,mH_exp_1,Val_exp_1);
      cout << mH_exp_1.size() << endl;

      for(unsigned int i = 0; i < mH_exp_1.size(); i++)
      	{
	  v_masses_exp_1.push_back( mH_exp_1[i] );
	  v_exp_1.push_back(Val_exp_1[i]);
	  cout << "exp_1 " <<  mH_exp_1[i] << "  " << Val_exp_1[i] << "  ( "<< ROOT::Math::gaussian_quantile_c(Val_exp_1[i],1) << ")" << endl;
	}
    }
  

  vector<double> mH_exp_2, Val_exp_2;  
  vector<double> v_masses_exp_2, v_exp_2;
  if(addExpected_2)
    {
      getPvals(inFileExp_2,mH_exp_2,Val_exp_2);
      for(unsigned int i = 0; i < mH_exp_2.size(); i++)
      	{
	  v_masses_exp_2.push_back( mH_exp_2[i] );
	  v_exp_2.push_back(Val_exp_2[i]);
	  cout << "exp_2 " <<  mH_exp_2[i] << "  " << Val_exp_2[i] << "  ( "<< ROOT::Math::gaussian_quantile_c(Val_exp_2[i],1) << ")" << endl;
	}
    }


  int nMassEffExp=0;
  int nExcludedExp=0;
  const int sizeVExp = v_masses_exp.size();
  double a_masses_exp[sizeVExp], a_exp[sizeVExp];
  if(addExpected)
    {
      for(unsigned int n = 0; n < v_masses_exp.size(); n++)
	{
	  a_masses_exp[nMassEffExp] = v_masses_exp[n];
	  a_exp[nMassEffExp] = v_exp[n];
	  nMassEffExp++;
	}
      
    }

  int nMassEffExp_1=0;
  int nExcludedExp_1=0;
  const int sizeVExp_1 = v_masses_exp_1.size();
  double a_masses_exp_1[sizeVExp_1], a_exp_1[sizeVExp_1];
  if(addExpected_1)
    {
      for(unsigned int n = 0; n < v_masses_exp_1.size(); n++)
	{
	  a_masses_exp_1[nMassEffExp_1] = v_masses_exp_1[n];
	  a_exp_1[nMassEffExp_1] = v_exp_1[n];
	  nMassEffExp_1++;
	}
      
    }

  int nMassEffExp_2=0;
  int nExcludedExp_2=0;
  const int sizeVExp_2 = v_masses_exp_2.size();
  double a_masses_exp_2[sizeVExp_2], a_exp_2[sizeVExp_2];
  if(addExpected_2)
    {
      for(unsigned int n = 0; n < v_masses_exp_2.size(); n++)
	{
	  a_masses_exp_2[nMassEffExp_2] = v_masses_exp_2[n];
	  a_exp_2[nMassEffExp_2] = v_exp_2[n];
	  nMassEffExp_2++;
	}
      
    }


  // ------------------- Change Values to Arrays -------------------- //

  int nMassEff=0;
  int nExcluded=0;
  const int sizeV = v_masses.size();
  double a_masses[sizeV], a_obs[sizeV];
  for(unsigned int m = 0; m < v_masses.size(); m++)
    {
      a_masses[nMassEff] = v_masses[m];
      a_obs[nMassEff] = v_obs[m];
      nMassEff++;
      
    }


  int nMassEff_1=0;
  int nExcluded_1=0;
  const int sizeV_1 = v_masses_1.size();
  double a_masses_1[sizeV_1], a_obs_1[sizeV_1];
  for(unsigned int m = 0; m < v_masses_1.size(); m++)
    {
      a_masses_1[nMassEff_1] = v_masses_1[m];
      a_obs_1[nMassEff_1] = v_obs_1[m];
      nMassEff_1++;
      
    }


  int nMassEff_2=0;
  int nExcluded_2=0;
  const int sizeV_2 = v_masses_2.size();
  double a_masses_2[sizeV_2], a_obs_2[sizeV_2];
  for(unsigned int m = 0; m < v_masses_2.size(); m++)
    {
      a_masses_2[nMassEff_2] = v_masses_2[m];
      a_obs_2[nMassEff_2] = v_obs_2[m];
      nMassEff_2++;
      
    }

  cout << "Excluded " << nExcluded << " sick mass points!" << endl;

  // ------------------- Draw  -------------------- //

  TCanvas *c = new TCanvas("c","c",canvasX,canvasY);
  if(addObs){
    TGraph *grObs = new TGraph(nMassEff, a_masses, a_obs);
    grObs->SetLineWidth(2);
    grObs->SetLineColor(kRed);
    grObs->SetMarkerStyle(20);
    grObs->SetMarkerSize(0.7);
    grObs->SetMarkerColor(kRed);
  }
  if(addObs_1){
    TGraph *grObs_1 = new TGraph(nMassEff_1, a_masses_1, a_obs_1);
    grObs_1->SetLineWidth(2);
    grObs_1->SetLineColor(kBlue);
    grObs_1->SetMarkerStyle(20);
    grObs_1->SetMarkerSize(0.7);
    grObs_1->SetMarkerColor(kBlue);
  }
  if(addObs_2){
    TGraph *grObs_2 = new TGraph(nMassEff_2, a_masses_2, a_obs_2);
    grObs_2->SetLineWidth(2);
    grObs_2->SetLineColor(kBlack);
    grObs_2->SetMarkerStyle(20);
    grObs_2->SetMarkerSize(0.7);
    grObs_2->SetMarkerColor(kBlack);
  }
  if(addExpected)
    {
      TGraph *grExp = new TGraph(nMassEffExp, a_masses_exp, a_exp);
      grExp->SetLineWidth(2);
      grExp->SetLineColor(kRed);
      //grExp->SetLineStyle(2);
      grExp->SetLineStyle(7);
      grExp->SetName("Expected");
      grExp->SetTitle("Expected");
    }
  if(addExpected_1)
    {
      TGraph *grExp_1 = new TGraph(nMassEffExp_1, a_masses_exp_1, a_exp_1);
      grExp_1->SetLineWidth(2);
      grExp_1->SetLineColor(kBlue);
      grExp_1->SetLineStyle(7);
      grExp_1->SetName("Expected1");
      grExp_1->SetTitle("Expected1");
    }
  if(addExpected_2)
    {
      TGraph *grExp_2 = new TGraph(nMassEffExp_2, a_masses_exp_2, a_exp_2);
      grExp_2->SetLineWidth(2);
      grExp_2->SetLineColor(kBlack);
      grExp_2->SetLineStyle(7);
      grExp_2->SetName("Expected2");
      grExp_2->SetTitle("Expected2");
    }

  char outfileName[192];
  
  // --------------- Low Mass Zoom -------------- //

  double ptLow= 0.42, ptHigh = 0.69;

  TPaveText *pt = new TPaveText(0.1577181,0.9562937,0.9580537,0.9947552,"brNDC");
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetFillStyle(0);
  pt->SetTextFont(42);
  pt->SetTextSize(0.03);
  TText *text = pt->AddText(0.01,0.5,"CMS preliminary");
  text = pt->AddText(0.3,0.6,"#sqrt{s} = 7 TeV, L = 5.1 fb^{-1}  #sqrt{s} = 8 TeV, L = 12.2 fb^{-1}");
  pt->Draw();   
 //  TPaveText *pt = new TPaveText(ptLow,0.46,ptHigh,0.50,"NDC");
//   pt->SetFillColor(0);
//   pt->SetTextFont(42);
//   pt->AddText("CMS Preliminary");
//   TPaveText *pt4 = new TPaveText(ptLow,0.42,ptHigh,0.46,"NDC");
//   pt4->SetFillColor(0);
//   pt4->SetTextFont(42);
//   pt4->AddText("H #rightarrow ZZ #rightarrow 4L"); 

//   //TPaveText *pt2 = new TPaveText(0.69,0.94,0.98,0.99,"NDC");
//   TPaveText *pt2 = new TPaveText(ptLow,0.38,ptHigh,0.42,"NDC");
//   pt2->SetFillColor(0);
//   pt2->SetTextFont(42);
//   char lum[192];
//   sprintf(lum," #sqrt{s} = 7 TeV, L = %.2f fb^{-1}",5.051);
//   pt2->AddText(lum); 
//   //TPaveText *pt3 = new TPaveText(0.69,0.90,0.98,0.95,"NDC");
//   TPaveText *pt3 = new TPaveText(ptLow,0.34,ptHigh,0.38,"NDC");
//   pt3->SetFillColor(0);
//   pt3->SetTextFont(42);
//   char lum2[192];
//   sprintf(lum2," #sqrt{s} = 8 TeV, L = %.1f fb^{-1}",12.21);
//   pt3->AddText(lum2); 


  //TPaveText *oneSig = new TPaveText(0.94,0.85,0.97,0.89,"NDC");
  TPaveText *oneSig = new TPaveText(0.97,0.84,1,0.88,"NDC");
  oneSig->SetFillColor(0);
  oneSig->SetTextFont(42);
  //oneSig->SetTextColor(kRed);
  oneSig->AddText("1#sigma"); 

  //TPaveText *twoSig = new TPaveText(0.94,0.765,0.97,0.805,"NDC");
  TPaveText *twoSig = new TPaveText(0.97,0.74,1,0.78,"NDC");
  twoSig->SetFillColor(0);
  twoSig->SetTextFont(42);
  //twoSig->SetTextColor(kRed);
  twoSig->AddText("2#sigma"); 

  //TPaveText *threeSig = new TPaveText(0.94,0.63,0.97,0.67,"NDC");
  TPaveText *threeSig = new TPaveText(0.97,0.6,1,0.64,"NDC");
  threeSig->SetFillColor(0);
  threeSig->SetTextFont(42);
  //threeSig->SetTextColor(kRed);
  threeSig->AddText("3#sigma"); 

  //TPaveText *fourSig = new TPaveText(0.94,0.455,0.97,0.495,"NDC");
  TPaveText *fourSig = new TPaveText(0.97,0.405,1,0.445,"NDC");
  fourSig->SetFillColor(0);
  fourSig->SetTextFont(42);
  //fourSig->SetTextColor(kRed);
  fourSig->AddText("4#sigma"); 

  //TPaveText *fiveSig = new TPaveText(0.94,0.24,0.97,0.28,"NDC");
  TPaveText *fiveSig = new TPaveText(0.97,0.17,1,0.21,"NDC");
  fiveSig->SetFillColor(0);
  fiveSig->SetTextFont(42);
  //fiveSig->SetTextColor(kRed);
  fiveSig->AddText("5#sigma");

  TLegend * box2 = new TLegend(0.46,0.48,0.66,0.6);
  box2->SetFillColor(0);
  box2->SetTextFont(42);
  //box2->SetBorderSize(0);
  if(addObs)box2->AddEntry(grObs,Graph1,"l"); 
  if(addObs_1)box2->AddEntry(grObs_1,Graph2,"l"); 
  if(addObs_2)box2->AddEntry(grObs_2,Graph3,"l"); 
  if(addExpected_2)box2->AddEntry(grExp_2,"Expected","l");
  TLatex latex; latex.SetTextFont(42); latex.SetTextSize(0.045); latex.SetTextColor(2);



	
  if(grid) c->SetGrid();
   
  TH1F *hr = c->DrawFrame(105.0,yLow,180.0,yHigh);
  TLine *l1=new TLine();
  l1->SetLineStyle(1);
  l1->SetLineWidth(1);
  l1->SetLineColor(kRed);
  l1->DrawLine(105.0,ROOT::Math::normal_cdf_c(1, 1.0),180.0,ROOT::Math::normal_cdf_c(1, 1.0));
  TLine *l2=new TLine();
  l2->SetLineStyle(1);
  l2->SetLineWidth(1);
  l2->SetLineColor(kRed);
  l2->DrawLine(105.0,ROOT::Math::normal_cdf_c(2, 1.0),180.0,ROOT::Math::normal_cdf_c(2, 1.0));
  TLine *l3=new TLine();
  l3->SetLineStyle(1);
  l3->SetLineWidth(1);
  l3->SetLineColor(kRed);
  l3->DrawLine(105.0,ROOT::Math::normal_cdf_c(3, 1.0),180.0,ROOT::Math::normal_cdf_c(3, 1.0));
  TLine *l4=new TLine();
  l4->SetLineStyle(1);
  l4->SetLineWidth(1);
  l4->SetLineColor(kRed);
  l4->DrawLine(105.0,ROOT::Math::normal_cdf_c(4, 1.0),180.0,ROOT::Math::normal_cdf_c(4, 1.0));
  TLine *l5=new TLine();
  l5->SetLineStyle(1);
  l5->SetLineWidth(1);
  l5->SetLineColor(kRed);
  l5->DrawLine(105.0,ROOT::Math::normal_cdf_c(5, 1.0),180.0,ROOT::Math::normal_cdf_c(5, 1.0));


  if(addObs)grObs->Sort();
  if(addExpected)grExp->Sort();
  if(addObs_1)grObs_1->Sort();
  if(addExpected_1)grExp_1->Sort();
  if(addObs_2)grObs_2->Sort();
  if(addExpected_2)grExp_2->Sort();
  if(addObs)
    {
      if(points)grObs->Draw("LP");
      else grObs->Draw("L");
    }
  if(addObs_1)
    {
      if(points_1)grObs_1->Draw("LP");
      else grObs_1->Draw("L");
    }
  if(addObs_2)
    {
      if(points_2)grObs_2->Draw("LP");
      else grObs_2->Draw("L");
    }

  if(addExpected)grExp->Draw("L");
  if(addExpected_1)grExp_1->Draw("L");
  if(addExpected_2)grExp_2->Draw("L");
  //latex.DrawLatex(xHigh+(xHigh-xLow)*0.01, ROOT::Math::normal_cdf_c(1,1.0)*1.1,"1#sigma");

  pt->Draw("SAME");
//   pt2->Draw("SAME");
//   pt3->Draw("SAME");
//   pt4->Draw("SAME");


  hr->GetXaxis()->SetTitle(xTitle);
  hr->GetYaxis()->SetTitle(yTitle);
  hr->GetYaxis()->SetTitleOffset(1.2);		
  hr->GetXaxis()->SetLabelSize(0.04);
  hr->GetYaxis()->SetLabelSize(0.04);

  if(logy)gPad->SetLogy();
  c->Update();
  if(gridOnTop)gPad->RedrawAxis("g");

  oneSig->Draw("SAME");
  twoSig->Draw("SAME");
  threeSig->Draw("SAME");
  fourSig->Draw("SAME");
  fiveSig->Draw("SAME");

  box2->Draw();
  sprintf( outfileName,"%s/Pvals_%s_lowMass_%s.eps",plotDir.c_str(),method.c_str(),dimension.c_str() );
  c->SaveAs(outfileName);
  sprintf( outfileName,"%s/Pvals_%s_lowMass_%s.png",plotDir.c_str(),method.c_str(),dimension.c_str() );
  c->SaveAs(outfileName);
  sprintf( outfileName,"%s/Pvals_%s_lowMass_%s.root",plotDir.c_str(),method.c_str(),dimension.c_str() );
  c->SaveAs(outfileName);


  // --------------- Full Mass Range ---------------- //
	
  TLegend * box3 = new TLegend(0.6,0.18,0.85,0.31);
  box3->SetFillColor(0);
  box3->SetTextFont(42);
  //box3->SetBorderSize(0);
  if(addObs)box3->AddEntry(grObs,Graph1,"l"); 
  if(addObs_1)box3->AddEntry(grObs_1,Graph2,"l"); 
  if(addObs_2)box3->AddEntry(grObs_2,Graph3,"l"); 
  if(addExpected_2)box3->AddEntry(grExp_2,"Expected","l");

  TCanvas *cl = new TCanvas("cl","cl",canvasX,canvasY);
  cl->cd();
  if(grid) cl->SetGrid();
  
  TH1F *hrl = cl->DrawFrame(xLow,yLow,xHigh,yHigh);
  l1->DrawLine(xLow,ROOT::Math::normal_cdf_c(1, 1.0),xHigh,ROOT::Math::normal_cdf_c(1, 1.0));
  l2->DrawLine(xLow,ROOT::Math::normal_cdf_c(2, 1.0),xHigh,ROOT::Math::normal_cdf_c(2, 1.0));
  l3->DrawLine(xLow,ROOT::Math::normal_cdf_c(3, 1.0),xHigh,ROOT::Math::normal_cdf_c(3, 1.0));
  l4->DrawLine(xLow,ROOT::Math::normal_cdf_c(4, 1.0),xHigh,ROOT::Math::normal_cdf_c(4, 1.0));
  l5->DrawLine(xLow,ROOT::Math::normal_cdf_c(5, 1.0),xHigh,ROOT::Math::normal_cdf_c(5, 1.0));

  if(addObs)
    {
      if(points)grObs->Draw("LP");
      else grObs->Draw("L");
    }
  if(addObs_1)
    {
      if(points_1)grObs_1->Draw("LP");
      else grObs_1->Draw("L");
    }
  if(addObs_2)
    {
      if(points_2)grObs_2->Draw("LP");
      else grObs_2->Draw("L");
    }

  if(addExpected)grExp->Draw("L");
  if(addExpected_1)grExp_1->Draw("L");
  if(addExpected_2)grExp_2->Draw("L");
 
  pt->Draw("SAME");
//   pt2->Draw("SAME");
//   pt3->Draw("SAME");
//  pt4->Draw("SAME");

  hrl->GetXaxis()->SetTitle(xTitle);
  hrl->GetYaxis()->SetTitle(yTitle);
  hrl->GetYaxis()->SetTitleOffset(1.1);		

  if(logy)gPad->SetLogy();
  if(logx)
    {
      hrl->GetXaxis()->SetMoreLogLabels();
      hrl->GetXaxis()->SetNoExponent();
      gPad->SetLogx();
      TLine tick; tick.SetLineWidth(1); tick.SetLineColor(1);
      double dyh = yHigh * 0.08;
      double dyl = yLow * 0.08; //fabs(c1->PixeltoY(c1->VtoPixel(0.95)) - c1->PixeltoY(c1->VtoPixel(0.94)));
      if (gPad->GetLogy() && log(yHigh/yLow) > log(1e6)) { dyh *= 2; dyl *= 2; }
      if (gPad->GetLogy() == 0) { dyh = dyl = 0.01*(yHigh-yLow); }
      if (isTiny) { dyh *= 2; dyl *= 2; }
      for (int j = 100; j < 1000; j += 10)  {
	  if (j > 400 && j % 20 == 10) continue;
	  tick.DrawLine(j, yLow, j, yLow+(j % 100 == 0 ? 2*dyl : dyl));
	  tick.DrawLine(j, yHigh, j, yHigh-(j % 100 == 0 ? 2*dyh : dyh));
      }
    }
  cl->Update();
  hrl->GetXaxis()->SetRangeUser(xLow,xHigh);
  hrl->GetXaxis()->SetLabelSize(0.04);
  hrl->GetYaxis()->SetLabelSize(0.04);
  cl->Update();

  if(gridOnTop)gPad->RedrawAxis("g");
  box2->Draw("SAME");
  
  oneSig->Draw("SAME");
  twoSig->Draw("SAME");
  threeSig->Draw("SAME");
  fourSig->Draw("SAME");
  fiveSig->Draw("SAME");

  sprintf( outfileName,"%s/Pvals_%s_wholeMass_%s.eps",plotDir.c_str(),method.c_str(),dimension.c_str() );
  cl->SaveAs(outfileName);
  sprintf( outfileName,"%s/Pvals_%s_wholeMass_%s.png",plotDir.c_str(),method.c_str(),dimension.c_str() );
  cl->SaveAs(outfileName);
  sprintf( outfileName,"%s/Pvals_%s_wholeMass_%s.root",plotDir.c_str(),method.c_str(),dimension.c_str() );
  cl->SaveAs(outfileName);

	
}



void getPvals(TFile *f, std::vector<double> &v_mh,std::vector<double> &v_obs)
{

  TTree *tree =(TTree*)f->Get("limit");
  
  double mh,limit;
  float quant;
  tree->SetBranchAddress("mh",&mh);
  tree->SetBranchAddress("limit",&limit);
  tree->SetBranchAddress("quantileExpected",&quant);

  tree->BuildIndex("mh");
  TTreeIndex *index = (TTreeIndex*)tree->GetTreeIndex();
  
  for( int i = 0; i < index->GetN(); i++ )
    {
      Long64_t local = tree->LoadTree( index->GetIndex()[i] );
      tree->GetEntry(local);

      //double tmpMH = mh;
      //if( fabs(floor(tmpMH)-mh) == 0.5) continue;

      if(_DEBUG_)cout << "mH: " << mh << " limit: " << limit << " quantileExpected: " << quant << endl;  
      if(quant>-1.01&&quant<-0.99)
	{
	  v_obs.push_back(limit);
	  v_mh.push_back(mh);
	}
      else {cout<<"Error! Unknown Quantile =  " << quant << endl;}
    }
  
}

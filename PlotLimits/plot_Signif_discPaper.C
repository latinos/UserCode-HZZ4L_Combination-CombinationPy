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
TString inputFile = "results_2D_7TeV_CorScale/higgsCombineHZZ4L_PLP.root";
TString inputFile_1 = "results_2D_8TeV_CorScale/higgsCombineHZZ4L_PLP.root";
TString inputFile_2 = "results_2D_Combined_CorScale/higgsCombineHZZ4L_PLP.root";
TString inputFileExp ="results_2D_7TeV_CorScale/higgsCombineHZZ4L_PLPE.root";
TString inputFileExp_1 = "results_2D_8TeV_CorScale/higgsCombineHZZ4L_PLPE.root";
TString inputFileExp_2 = "results_2D_Combined_CorScale/higgsCombineHZZ4L_PLPE.root";
TString Graph1 = "7 TeV Observed";
TString Graph2 = "8 TeV Observed";
TString Graph3 = "Observed";
const bool addObs = true;
const bool addObs_1 = true;
const bool addObs_2 = true;
const bool addExpected = false;
const bool addExpected_1 = false;
const bool addExpected_2 = true;
string method = "PLP";
Double_t xLow = 99.9;
Double_t xHigh = 601.0;
Double_t yLow = 1e-5;
Double_t yHigh = 1.0;
TString xTitle = "m_{H} (GeV)";
TString yTitle = "local p-value";
const bool logy = true;
const bool logx = true;
const bool grid = false;
const bool gridOnTop = false;
const bool points = false;
const bool points_1 = false;
const bool points_2 = false;
const bool isTiny = false;
int canvasX = 900;
int canvasY = 700;
Double_t lumi = 5.051;
const bool _DEBUG_ = false;
string plotDir = "plots_final";
string dimension = "2D_forSMPaper";
// ----------------------- //


using namespace std;

void plot_Signif_discPaper()
{

  gROOT->ProcessLine(".x tdrstyle.cc");
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.06);
  gStyle->SetLabelOffset(0.012,"XY");

  TFile *inFile = new TFile(inputFile,"READ");
  TFile *inFileExp;
  if(addExpected)inFileExp = new TFile(inputFileExp,"READ");

  if(!inFile){ cout << "Cannot find file " << inputFile << end; return;}
  if(addExpected && !inFileExp){cout << "Cannot find file " << inputFileExp << end; return;}


  TFile *inFile_1 = new TFile(inputFile_1,"READ");
  TFile *inFileExp_1;
  if(addExpected_1)inFileExp_1 = new TFile(inputFileExp_1,"READ");

  if(!inFile_1){ cout << "Cannot find file " << inputFile_1 << end;return;}
  if(addExpected_1 && !inFileExp_1){cout <<"Cannot find file " << inputFileExp_1 << end; return;}


  TFile *inFile_2 = new TFile(inputFile_2,"READ");
  TFile *inFileExp_2;
  if(addExpected_2)inFileExp_2 = new TFile(inputFileExp_2,"READ");
  
  if(!inFile_2){ cout << "Cannot find file " << inputFile_2 << end;return;}
  if(addExpected_2 && !inFileExp_2){cout <<"Cannot find file " << inputFileExp_2 << end; return;}


  // ------------------- Get Values -------------------- //

  vector<double> mH, Val_obs;
  getPvals(inFile,mH,Val_obs);
  vector<double> v_masses, v_obs;

  for(unsigned int i = 0; i < mH.size(); i++)
    {
      v_masses.push_back( mH[i] );
      v_obs.push_back(Val_obs[i]);
      cout << "obs " <<  mH[i] << "  " << Val_obs[i] << endl;
    }


  vector<double> mH_1, Val_obs_1;
  getPvals(inFile_1,mH_1,Val_obs_1);
  vector<double> v_masses_1, v_obs_1;

  for(unsigned int i = 0; i < mH_1.size(); i++)
    {
      v_masses_1.push_back( mH_1[i] );
      v_obs_1.push_back(Val_obs_1[i]);
      cout << "obs_1 " <<  mH_1[i] << "  " << Val_obs_1[i] << endl;
    }


  vector<double> mH_2, Val_obs_2;
  getPvals(inFile_2,mH_2,Val_obs_2);
  vector<double> v_masses_2, v_obs_2;

  for(unsigned int i = 0; i < mH_2.size(); i++)
    {
      v_masses_2.push_back( mH_2[i] );
      v_obs_2.push_back(Val_obs_2[i]);
      cout << "obs_2 " << mH_2[i] << "  " << Val_obs_2[i] << endl;
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
	  cout << "exp " << mH_exp[i] << "  " << Val_exp[i] << endl;
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
	  cout << "exp_1 " <<  mH_exp_1[i] << "  " << Val_exp_1[i] << endl;
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
	  cout << "exp_2 " <<  mH_exp_2[i] << "  " << Val_exp_2[i] << endl;
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
  TGraph *grObs = new TGraph(nMassEff, a_masses, a_obs);
  grObs->SetLineWidth(3);
  grObs->SetLineColor(kRed);
  grObs->SetMarkerStyle(20);
  grObs->SetMarkerSize(1);
  grObs->SetMarkerColor(kRed);

  TGraph *grObs_1 = new TGraph(nMassEff_1, a_masses_1, a_obs_1);
  grObs_1->SetLineWidth(3);
  grObs_1->SetLineColor(kBlue);
  grObs_1->SetMarkerStyle(20);
  grObs_1->SetMarkerSize(1);
  grObs_1->SetMarkerColor(kBlue);

  TGraph *grObs_2 = new TGraph(nMassEff_2, a_masses_2, a_obs_2);
  grObs_2->SetLineWidth(3);
  grObs_2->SetLineColor(kBlack);
  grObs_2->SetMarkerStyle(20);
  grObs_2->SetMarkerSize(1.2);
  grObs_2->SetMarkerColor(kBlack);

  if(addExpected)
    {
      TGraph *grExp = new TGraph(nMassEffExp, a_masses_exp, a_exp);
      grExp->SetLineWidth(5);
      grExp->SetLineColor(kRed);
      grExp->SetLineStyle(7);
    }
  if(addExpected_1)
    {
      TGraph *grExp_1 = new TGraph(nMassEffExp_1, a_masses_exp_1, a_exp_1);
      grExp_1->SetLineWidth(5);
      grExp_1->SetLineColor(kBlue);
      grExp_1->SetLineStyle(7);
    }
  if(addExpected_2)
    {
      TGraph *grExp_2 = new TGraph(nMassEffExp_2, a_masses_exp_2, a_exp_2);
      grExp_2->SetLineWidth(5);
      grExp_2->SetLineColor(kBlack);
      grExp_2->SetLineStyle(7);
    }

  char outfileName[192];
  
  // --------------- Low Mass Zoom -------------- //

  double ptLow= 0.61, ptHigh = 0.89;


  TPaveText *pt = new TPaveText(0.16,0.94,0.23,1,"NDC");
  pt->SetFillColor(0);
  pt->SetTextFont(42);
  pt->AddText("CMS");
  TPaveText *pt2 = new TPaveText(0.45,0.94,0.71,1,"NDC");
  pt2->SetFillColor(0);
  pt2->SetTextFont(42);
  char lum[192];
  sprintf(lum," #sqrt{s} = 7 TeV, L = %.1f fb^{-1}",5.1);
  pt2->AddText(lum); 
  TPaveText *pt3 = new TPaveText(0.705,0.94,0.965,1,"NDC");
  pt3->SetFillColor(0);
  pt3->SetTextFont(42);
  char lum2[192];
  sprintf(lum2," #sqrt{s} = 8 TeV, L = %.1f fb^{-1}",5.3);
  pt3->AddText(lum2); 


  TPaveText *oneSig = new TPaveText(0.96,0.785,0.995,0.845,"NDC");
  oneSig->SetFillColor(0);
  oneSig->SetTextFont(42);
  //oneSig->SetTextColor(kRed);
  oneSig->AddText("1#sigma"); 

  TPaveText *twoSig = new TPaveText(0.96,0.65,1,0.71,"NDC");
  twoSig->SetFillColor(0);
  twoSig->SetTextFont(42);
  //twoSig->SetTextColor(kRed);
  twoSig->AddText("2#sigma"); 

  TPaveText *threeSig = new TPaveText(0.96,0.45,1,0.51,"NDC");
  threeSig->SetFillColor(0);
  threeSig->SetTextFont(42);
  //threeSig->SetTextColor(kRed);
  threeSig->AddText("3#sigma"); 

  TPaveText *fourSig = new TPaveText(0.96,0.19,1,0.25,"NDC");
  fourSig->SetFillColor(0);
  fourSig->SetTextFont(42);
  //fourSig->SetTextColor(kRed);
  fourSig->AddText("4#sigma"); 

  TLegend * box2 = new TLegend(0.55,0.22,0.92,0.46);
  box2->SetFillColor(0);
  box2->SetTextFont(42);
  box2->SetBorderSize(0);
  box2->SetTextSize(0.037);
  if(addObs_2)box2->AddEntry(grObs_2,Graph3,"l"); 
  if(addExpected_2)box2->AddEntry(grExp_2,"Exp. for SM Higgs Boson","l");
  if(addObs)box2->AddEntry(grObs,Graph1,"l"); 
  if(addObs_1)box2->AddEntry(grObs_1,Graph2,"l"); 


  TLatex latex; latex.SetTextFont(42); latex.SetTextSize(0.045); latex.SetTextColor(2);



	
  if(grid) c->SetGrid();
   
  TH1F *hr = c->DrawFrame(109.0,yLow,150.0,yHigh);
  TLine *l1=new TLine();
  l1->SetLineStyle(1);
  l1->SetLineWidth(1);
  l1->SetLineColor(kRed);
  l1->DrawLine(109.0,ROOT::Math::normal_cdf_c(1, 1.0),150.0,ROOT::Math::normal_cdf_c(1, 1.0));
  TLine *l2=new TLine();
  l2->SetLineStyle(1);
  l2->SetLineWidth(1);
  l2->SetLineColor(kRed);
  l2->DrawLine(109.0,ROOT::Math::normal_cdf_c(2, 1.0),150.0,ROOT::Math::normal_cdf_c(2, 1.0));
  TLine *l3=new TLine();
  l3->SetLineStyle(1);
  l3->SetLineWidth(1);
  l3->SetLineColor(kRed);
  l3->DrawLine(109.0,ROOT::Math::normal_cdf_c(3, 1.0),150.0,ROOT::Math::normal_cdf_c(3, 1.0));
  TLine *l4=new TLine();
  l4->SetLineStyle(1);
  l4->SetLineWidth(1);
  l4->SetLineColor(kRed);
  l4->DrawLine(109.0,ROOT::Math::normal_cdf_c(4, 1.0),150.0,ROOT::Math::normal_cdf_c(4, 1.0));


  grObs->Sort();
  if(addExpected)grExp->Sort();
  grObs_1->Sort();
  if(addExpected_1)grExp_1->Sort();
  grObs_2->Sort();
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
  pt2->Draw("SAME");
  pt3->Draw("SAME");
  //pt4->Draw("SAME");


  hr->GetXaxis()->SetTitle(xTitle);
  hr->GetYaxis()->SetTitle(yTitle);
  hr->GetYaxis()->SetTitleOffset(1.05);		
 

  if(logy)gPad->SetLogy();
  c->Update();
  if(gridOnTop)gPad->RedrawAxis("g");

  oneSig->Draw("SAME");
  twoSig->Draw("SAME");
  threeSig->Draw("SAME");
  fourSig->Draw("SAME");
 
  box2->Draw();
  sprintf( outfileName,"%s/Pvals_%s_lowMass_%s_7p8sep.eps",plotDir.c_str(),method.c_str(),dimension.c_str() );
  c->SaveAs(outfileName);
  sprintf( outfileName,"%s/Pvals_%s_lowMass_%s_7p8sep.png",plotDir.c_str(),method.c_str(),dimension.c_str() );
  c->SaveAs(outfileName);



  // --------------- Full Mass Range ---------------- //
	
  TLegend * box3 = new TLegend(0.6,0.18,0.85,0.31);
  box3->SetFillColor(0);
  box3->SetTextFont(42);
  //box3->SetBorderSize(0);
  if(addObs_2)box3->AddEntry(grObs_2,Graph3,"l"); 
  if(addExpected_2)box3->AddEntry(grExp_2,"SM Higgs Expected","l");
  if(addObs)box3->AddEntry(grObs,Graph1,"l"); 
  if(addObs_1)box3->AddEntry(grObs_1,Graph2,"l"); 

 

  TCanvas *cl = new TCanvas("cl","cl",canvasX,canvasY);
  cl->cd();
  if(grid) cl->SetGrid();
  
  TH1F *hrl = cl->DrawFrame(xLow,yLow,xHigh,yHigh);
  l1->DrawLine(xLow,ROOT::Math::normal_cdf_c(1, 1.0),xHigh,ROOT::Math::normal_cdf_c(1, 1.0));
  l2->DrawLine(xLow,ROOT::Math::normal_cdf_c(2, 1.0),xHigh,ROOT::Math::normal_cdf_c(2, 1.0));
  l3->DrawLine(xLow,ROOT::Math::normal_cdf_c(3, 1.0),xHigh,ROOT::Math::normal_cdf_c(3, 1.0));
  l4->DrawLine(xLow,ROOT::Math::normal_cdf_c(4, 1.0),xHigh,ROOT::Math::normal_cdf_c(4, 1.0));

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

  if(addExpected)grExp->Draw("C");
  if(addExpected_1)grExp_1->Draw("C");
  if(addExpected_2)grExp_2->Draw("C");
 
  pt->Draw("SAME");
  pt2->Draw("SAME");
  pt3->Draw("SAME");
  //pt4->Draw("SAME");

  hrl->GetXaxis()->SetTitle(xTitle);
  hrl->GetYaxis()->SetTitle(yTitle);
  hrl->GetYaxis()->SetTitleOffset(1.05);		

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
      for (int j = 100; j < 600; j += 10)  {
	  if (j > 400 && j % 20 == 10) continue;
	  tick.DrawLine(j, yLow, j, yLow+(j % 100 == 0 ? 2*dyl : dyl));
	  tick.DrawLine(j, yHigh, j, yHigh-(j % 100 == 0 ? 2*dyh : dyh));
      }
    }
  cl->Update();
  hrl->GetXaxis()->SetRangeUser(xLow,xHigh);
  cl->Update();

  if(gridOnTop)gPad->RedrawAxis("g");
  box2->Draw("SAME");
  
  oneSig->Draw("SAME");
  twoSig->Draw("SAME");
  threeSig->Draw("SAME");
  fourSig->Draw("SAME");
 
  sprintf( outfileName,"%s/Pvals_%s_wholeMass_%s_7p8sep.eps",plotDir.c_str(),method.c_str(),dimension.c_str() );
  cl->SaveAs(outfileName);
  sprintf( outfileName,"%s/Pvals_%s_wholeMass_%s_7p8sep.png",plotDir.c_str(),method.c_str(),dimension.c_str() );
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

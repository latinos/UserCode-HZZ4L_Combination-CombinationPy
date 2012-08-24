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
TString inputFile = "results_2D_Combined_CorScale/higgsCombineHZZ4L_PLP.root";
TString inputFileExp ="results_2D_Combined_CorScale/higgsCombineHZZ4L_PLPE.root";
const bool addExpected = true;
string method = "PLP";
Double_t xLow = 99.9;
Double_t xHigh = 601.0;
Double_t yLow = 1e-5;
Double_t yHigh = 1.0;
TString xTitle = "Higgs boson mass, GeV";
TString yTitle = "local #hat{p}-value";
const bool logy = true;
const bool logx = true;
const bool grid = true;
const bool gridOnTop = false;
const bool points = true;
const bool isTiny = false;
int canvasX = 900;
int canvasY = 700;
const bool _DEBUG_ = false;
string plotDir = "plots_TopUp";
string dimension = "2D";
int sqrts = 7;
Double_t lumi = 5.051;
// ----------------------- //


using namespace std;

void plot_Signif()
{

  gROOT->ProcessLine(".x tdrstyle.cc");
  gStyle->SetPadLeftMargin(0.16);

  TFile *inFile = new TFile(inputFile,"READ");
  TFile *inFileExp;
  if(addExpected)inFileExp = new TFile(inputFileExp,"READ");
  
  // ------------------- Get Values -------------------- //

  vector<double> mH, Val_obs;
  getPvals(inFile,mH,Val_obs);
  vector<double> v_masses, v_obs;
  
  for(unsigned int i = 0; i < mH.size(); i++)
    {
      v_masses.push_back( mH[i] );
      v_obs.push_back(Val_obs[i]);
      cout << mH[i] << "  " << Val_obs[i] << endl;

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

	  cout << mH_exp[i] << "  " << Val_exp[i] << endl;

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
  cout << "Excluded " << nExcluded << " sick mass points!" << endl;

  // ------------------- Draw  -------------------- //


  TCanvas *c = new TCanvas("c","c",canvasX,canvasY);
  TGraph *grObs = new TGraph(nMassEff, a_masses, a_obs);
  grObs->SetLineWidth(3);
  grObs->SetLineColor(kBlue);
  grObs->SetMarkerStyle(20);
  grObs->SetMarkerSize(1.2);
  grObs->SetMarkerColor(kBlue);

  if(addExpected)
    {
      TGraph *grExp = new TGraph(nMassEffExp, a_masses_exp, a_exp);
      grExp->SetLineWidth(1);
      grExp->SetLineColor(kBlack);
      grExp->SetLineStyle(2);
    }

  char outfileName[192];
  

  // --------------- Low Mass Zoom -------------- //
	
  TLegend * box2 = new TLegend(0.2,0.18,0.45,0.3);
  box2->SetFillColor(0);
  box2->SetTextFont(42);
  //box2->SetBorderSize(0);
   box2->AddEntry(grObs,"w/o m_{4l} uncertainties","l"); 
		
  TPaveText *pt = new TPaveText(0.16,0.95,0.44,0.99,"NDC");
  pt->SetFillColor(0);
  pt->SetTextFont(42);
  pt->AddText("CMS Preliminary 2012"); 
  TPaveText *pt2 = new TPaveText(0.69,0.95,0.98,0.99,"NDC");
  pt2->SetFillColor(0);
  pt2->SetTextFont(42);
  char lum[192];
  sprintf(lum," #sqrt{s} = %i TeV, L = %.2f fb^{-1}",sqrts,lumi);
  pt2->AddText(lum); 
	
  if(grid) c->SetGrid();
   
  TH1F *hr = c->DrawFrame(105.0,yLow,180.0,yHigh);

  if(points)grObs->Draw("LP");
  else grObs->Draw("L");
  if(addExpected)grExp->Draw("L");

  pt->Draw("SAME");
  pt2->Draw("SAME");
	
  hr->GetXaxis()->SetTitle(xTitle);
  hr->GetYaxis()->SetTitle(yTitle);
  hr->GetYaxis()->SetTitleOffset(1.2);		
  
  if(logy)gPad->SetLogy();
  c->Update();
  if(gridOnTop)gPad->RedrawAxis("g");
  TLine *l1=new TLine();
  l1->SetLineStyle(2);
  l1->SetLineWidth(3.0);
  l1->SetLineColor(kRed);
  l1->DrawLine(105.0,ROOT::Math::normal_cdf_c(1, 1.0),180.0,ROOT::Math::normal_cdf_c(1, 1.0));
  TLine *l2=new TLine();
  l2->SetLineStyle(2);
  l2->SetLineWidth(3.0);
  l2->SetLineColor(kRed);
  l2->DrawLine(105.0,ROOT::Math::normal_cdf_c(2, 1.0),180.0,ROOT::Math::normal_cdf_c(2, 1.0));
  TLine *l3=new TLine();
  l3->SetLineStyle(2);
  l3->SetLineWidth(3.0);
  l3->SetLineColor(kRed);
  l3->DrawLine(105.0,ROOT::Math::normal_cdf_c(3, 1.0),180.0,ROOT::Math::normal_cdf_c(3, 1.0));
  TLine *l4=new TLine();
  l4->SetLineStyle(2);
  l4->SetLineWidth(3.0);
  l4->SetLineColor(kRed);
  l4->DrawLine(105.0,ROOT::Math::normal_cdf_c(4, 1.0),180.0,ROOT::Math::normal_cdf_c(4, 1.0));
 
  //box2->Draw();
  sprintf( outfileName,"%s/Pvals_%s_lowMass_%s_%iTeV.eps",plotDir.c_str(),method.c_str(),dimension.c_str(),sqrts);
  //c->SaveAs(outfileName);
  sprintf( outfileName,"%s/Pvals_%s_lowMass_%s_%iTeV.png",plotDir.c_str(),method.c_str(),dimension.c_str(),sqrts);
  //c->SaveAs(outfileName);




  // --------------- Full Mass Range ---------------- //
	
  TLegend * box3 = new TLegend(0.2,0.18,0.45,0.3);
  box3->SetFillColor(0);
  box3->SetTextFont(42);
  //box3->SetBorderSize(0);
  box3->AddEntry(grObs,"w/o m_{4l} uncertainties","l"); 

  TCanvas *cl = new TCanvas("cl","cl",canvasX,canvasY);
  cl->cd();
  if(grid) cl->SetGrid();
  
  TH1F *hrl = cl->DrawFrame(xLow,yLow,xHigh,yHigh);
  if(points)grObs->Draw("LP");
  else grObs->Draw("L");
  if(addExpected)grExp->Draw("L");

  pt->Draw("SAME");
  pt2->Draw("SAME");
	
  hrl->GetXaxis()->SetTitle(xTitle);
  hrl->GetYaxis()->SetTitle(yTitle);
  hrl->GetYaxis()->SetTitleOffset(1.2);		

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
  //TLine *l1=new TLine();
  //l1->SetLineStyle(3);
  //l1->SetLineWidth(2.0);
  //l1->SetLineColor(kRed);
  l1->DrawLine(xLow,ROOT::Math::normal_cdf_c(1, 1.0),xHigh,ROOT::Math::normal_cdf_c(1, 1.0));
  //TLine *l2=new TLine();
  //l2->SetLineStyle(3);
  //l2->SetLineWidth(2.0);
  //l2->SetLineColor(kRed);
  l2->DrawLine(xLow,ROOT::Math::normal_cdf_c(2, 1.0),xHigh,ROOT::Math::normal_cdf_c(2, 1.0));
  //TLine *l3=new TLine();
  //l3->SetLineStyle(3);
  //l3->SetLineWidth(2.0);
  //l3->SetLineColor(kRed);
  l3->DrawLine(xLow,ROOT::Math::normal_cdf_c(3, 1.0),xHigh,ROOT::Math::normal_cdf_c(3, 1.0));
  //TLine *l4=new TLine();
  //l4->SetLineStyle(3);
  //l4->SetLineWidth(2.0);
  //l4->SetLineColor(kRed);
  l4->DrawLine(xLow,ROOT::Math::normal_cdf_c(4, 1.0),xHigh,ROOT::Math::normal_cdf_c(4, 1.0));
  if(gridOnTop)gPad->RedrawAxis("g");
  //box3->Draw();
  sprintf( outfileName,"%s/Pvals_%s_wholeMass_%s_%iTeV.eps",plotDir.c_str(),method.c_str(),dimension.c_str(),sqrts);
  //cl->SaveAs(outfileName);
  sprintf( outfileName,"%s/Pvals_%s_wholeMass_%s_%iTeV.png",plotDir.c_str(),method.c_str(),dimension.c_str(),sqrts);
  //cl->SaveAs(outfileName);


	
}



void getPvals(TFile *f, std::vector<double> &v_mh,std::vector<double> &v_obs)
{

  TTree *tree =(TTree*)f->Get("limit");
  
  double mh,limit;
  float quant;
  tree->SetBranchAddress("mh",&mh);
  tree->SetBranchAddress("limit",&limit);
  tree->SetBranchAddress("quantileExpected",&quant);
  
  for(int i=0;i<tree->GetEntries();i++)
    {
      tree->GetEntry(i);

      double tmpMH = mh;
      if( fabs(floor(tmpMH)-mh) == 0.5) continue;

      if(_DEBUG_)cout << "mH: " << mh << " limit: " << limit << " quantileExpected: " << quant << endl;  
      if(quant>-1.01&&quant<-0.99)
	{
	  v_obs.push_back(limit);
	  v_mh.push_back(mh);
	}
      else {cout<<"Error! Unknown Quantile =  " << quant << endl;}
    }
  
}

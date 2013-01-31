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

void getLimits(TFile *f, std::vector<double> &v_mh,std::vector<double> &v_mean,std::vector<double> &v_68l,std::vector<double> &v_68h,std::vector<double> &v_95l,std::vector<double> &v_95h,std::vector<double> &v_obs);

// --------- Inputs ------- //
TString inputFile = "../../plot/higgsCombineHZZ4L_ASCLS_HCP_MorLumi.root";
const bool addObsLimit = false;
const bool isXSxBR = false;
const bool _DEBUG_ = false;
string method = "FREQ";
Double_t xLow = 99.9;
Double_t xHigh = 1000.0;
Double_t yLow = 0.1;
Double_t yHigh = 25.0;
TString xTitle = "m_{H} [GeV]";
TString yTitle = "95% CL limit on #sigma/#sigma_{SM}";
const bool logy = true;
const bool logx = true;
const bool grid = true;
const bool gridOnTop = true;
const bool points = true;
const bool isTiny = false;
int canvasX = 900;
int canvasY = 700;
//double sqrts = 8.0;
//Double_t lumi = 2.968;
double sqrts = 8.0;
Double_t lumi = 5.261;
std::string plotDir = "plots_test";
std::string append = "2D";
TString LimitType = "2D";
bool printLimits = true;
// ----------------------- //


using namespace std;

void plot_Exclusion()
{

  gROOT->ProcessLine(".x tdrstyle.cc");
  gStyle->SetPadLeftMargin(0.16);

  if(isXSxBR)
    {
      gSystem->Load("../CreateDatacards/include/HiggsCSandWidth_cc");
      gSystem->Load("../CreateDatacards/include/HiggsCSandWidthSM4_cc");
      HiggsCSandWidth *myCSW = new HiggsCSandWidth();
      HiggsCSandWidthSM4 *myCSWSM4 = new HiggsCSandWidthSM4();
    }

  TFile *inFile = new TFile(inputFile,"READ");
  
  // ------------------- Get Values -------------------- //

  vector<double> mH, Val_obs, Val_mean, Val_68h, Val_68l, Val_95h, Val_95l;
  getLimits(inFile,mH,Val_mean,Val_68l,Val_68h,Val_95l,Val_95h,Val_obs);
  vector<double> v_masses, v_means, v_lo68, v_hi68, v_lo95, v_hi95, v_obs;
  vector<double> expExclusion,obsExclusion;
  for(unsigned int i = 0; i < mH.size(); i++)
    {
      v_masses.push_back( mH[i] );
      v_means.push_back( Val_mean[i] );
      v_lo68.push_back( min( Val_68l[i], Val_68h[i]) );
      v_hi68.push_back( max( Val_68h[i], Val_68l[i]) );
      v_lo95.push_back( min( Val_95l[i], Val_95h[i]) );
      v_hi95.push_back( max( Val_95h[i], Val_95l[i]) );
      v_obs.push_back(Val_obs[i]);
      if(Val_mean[i] < 1.0) expExclusion.push_back(mH[i]);
      if(Val_obs[i] < 1.0 && addObsLimit) obsExclusion.push_back(mH[i]);
      if(printLimits)cout << mH[i] << "  " << Val_mean[i] << "  "  << Val_obs[i] << endl;
    }

  // ------------------- For XSxBR --------------------- //
  
  double mHxs;
  double CStot,  CStot_p,  CStot_m;
  double BRHZZ4L, XSBRHZZ4L;
  const int sizeVXSBR = 490;
  
  double a_mhxs[sizeVXSBR], a_xsbr_err_h[sizeVXSBR], a_xsbr_err_l[sizeVXSBR];
  double a_xsbr[sizeVXSBR];
  
  if(isXSxBR)
    {
      int kk = 0;
      for( int k = 110; k < 600; k++)
	{ 
	  mHxs = k;
	  CStot = myCSW->HiggsCS(0,mHxs,sqrts,true);
	  CStot_p = myCSW->HiggsCSErrPlus(1,mHxs,sqrts);
	  CStot_m = myCSW->HiggsCSErrMinus(1,mHxs,sqrts);
	  
	  BRHZZ4L = myCSW->HiggsBR(14,mHxs,true);
	  XSBRHZZ4L = BRHZZ4L*CStot;
	  
	  a_mhxs[kk] = mHxs;
	  a_xsbr_err_h[kk] = XSBRHZZ4L*CStot_p*1000;
	  a_xsbr_err_l[kk] = -1*XSBRHZZ4L*CStot_m*1000;
	  a_xsbr[kk] = XSBRHZZ4L*1000;
	  kk++;
	}
      TGraph *gr_XSBR=new TGraphAsymmErrors(kk,a_mhxs,a_xsbr);
      TGraphAsymmErrors *gr_XSBR68=new TGraphAsymmErrors(kk,a_mhxs,a_xsbr,0,0,a_xsbr_err_l,a_xsbr_err_h);
    }

  // ------------------- Change Values to Arrays -------------------- //

  int nMassEff=0;
  int nExcluded=0;
  const int sizeV = v_masses.size();
  double a_masses[sizeV], a_means[sizeV], a_lo68[sizeV], a_hi68[sizeV], a_lo95[sizeV], a_hi95[sizeV], a_obs[sizeV];
  for(unsigned int m = 0; m < v_masses.size(); m++)
    {
      if(v_hi68.at(m)>=v_hi95.at(m) || v_lo68.at(m)<=v_lo95.at(m))
	{
	  cout << "Point at M = " << v_masses.at(m) << " excluded" << endl;
	  nExcluded++;
	  continue;
	}
      
      a_masses[nMassEff] = v_masses[m];
      a_means[nMassEff] = v_means[m];
      a_lo68[nMassEff] = v_lo68[m];
      a_hi68[nMassEff] = v_hi68[m];
      a_lo95[nMassEff] = v_lo95[m];
      a_hi95[nMassEff] = v_hi95[m];
      a_obs[nMassEff] = v_obs[m];
      if(isXSxBR)
	{
	  double xs = myCSW->HiggsCS(0,a_masses[nMassEff],sqrts,true);
	  double br = myCSW->HiggsBR(14,a_masses[nMassEff],true);
	  double XSxBR_factor = xs*br*1000;
	  cout << XSxBR_factor << endl;

	  a_obs[nMassEff]*=XSxBR_factor;
	  a_means[nMassEff]*=XSxBR_factor;
	  a_hi68[nMassEff]*=XSxBR_factor;
	  a_lo68[nMassEff]*=XSxBR_factor;
	  a_hi95[nMassEff]*=XSxBR_factor;
	  a_lo95[nMassEff]*=XSxBR_factor;
	}
      nMassEff++;
      
    }
  cout << "Excluded " << nExcluded << " sick mass points!" << endl;

  // --------------- Excluded Regions --------------- //

  for(int p = 0; p < expExclusion.size(); p++)
    {
      cout << "Expected Exclusion: " <<  expExclusion[p] << endl;
    }

  if(addObsLimit)
    {
      for(int q = 0; q < obsExclusion.size(); q++)
	{
	  cout << "Observed Exclusion: " <<  obsExclusion[q] << endl;
	}
    }

  // ------------------- Draw  -------------------- //




  TCanvas *c = new TCanvas("c","c",canvasX,canvasY);
  TGraph *gr = new TGraph(nMassEff, a_masses, a_means);
  TGraph* grshade_68 = new TGraph(2*nMassEff);
  TGraph* grshade_95 = new TGraph(2*nMassEff);
  TGraph *grObs = new TGraph(nMassEff, a_masses, a_obs);
  grObs->SetLineWidth(3);
  grObs->SetLineColor(kBlack);
  grObs->SetMarkerStyle(21);
  grObs->SetMarkerSize(0.8);
 

  gr->SetLineStyle(2);
  gr->SetLineWidth(3);
  gr->SetLineColor(kBlue);
  grshade_68->SetFillColor(kGreen);
  grshade_95->SetFillColor(kYellow);		
  grshade_68->SetLineStyle(2);
  grshade_95->SetLineStyle(2);
  grshade_68->SetLineWidth(3);
  grshade_95->SetLineWidth(3);
  grshade_68->SetLineColor(kBlue);
  grshade_95->SetLineColor(kBlue);
  
  for (int a = 0; a < nMassEff; a++)
    {
      grshade_68->SetPoint(a,a_masses[a],a_lo68[a]);
      grshade_68->SetPoint(sizeV+a,a_masses[nMassEff-a-1],a_hi68[nMassEff-a-1]);
      grshade_95->SetPoint(a,a_masses[a],a_lo95[a]);
      grshade_95->SetPoint(nMassEff+a,a_masses[nMassEff-a-1],a_hi95[nMassEff-a-1]);
    }
	
	
  char outfileName[192];

  TF1* oneLine = new TF1("oneLine","1",xLow,xHigh);
  oneLine->SetLineColor(kRed);
  oneLine->SetLineWidth(3);

  // --------------- Low Mass Zoom -------------- //
	
  //  TLegend * box2 = new TLegend(0.35,0.7,0.75,0.90);
  TLegend * box2 = new TLegend(0.52,0.7,0.92,0.9);
  box2->SetFillColor(0);
  //box2->SetBorderSize(0);
  if (addObsLimit){ box2->AddEntry(grObs,"Observed Asym. CLs","l"); }
  //box2->AddEntry(gr,LimitType+" Fit Expected Asym. CLs","l");
  box2->AddEntry(gr,"Expected Asym. CLs","l");
  //box2->AddEntry(grshade_68,"68% expectation","f");
  //box2->AddEntry(grshade_95,"95% expectation","f");
  box2->AddEntry(grshade_68,"Expected #pm 1#sigma","lf");
  box2->AddEntry(grshade_95,"Expected #pm 2#sigma","lf");
  //box2->AddEntry(oneLine,"#sigma / #sigma_{SM}","l");

  double ptLow= 0.26, ptHigh = 0.51;
		
  TPaveText *pt = new TPaveText(ptLow,0.84,ptHigh,0.88,"NDC");
  pt->SetFillColor(0);
  pt->SetTextFont(42);
  pt->AddText("CMS Preliminary");
  TPaveText *pt4 = new TPaveText(ptLow,0.8,ptHigh,0.84,"NDC");
  pt4->SetFillColor(0);
  pt4->SetTextFont(42);
  pt4->AddText("H #rightarrow ZZ #rightarrow 4L"); 

  //TPaveText *pt2 = new TPaveText(0.69,0.94,0.98,0.99,"NDC");
  TPaveText *pt2 = new TPaveText(ptLow,0.76,ptHigh,0.8,"NDC");
  pt2->SetFillColor(0);
  pt2->SetTextFont(42);
  char lum[192];
  sprintf(lum," #sqrt{s} = 8 TeV, L = %.2f fb^{-1}",5.261);
  //sprintf(lum," #sqrt{s} = 7 TeV, L = %.2f fb^{-1}",5.051);
  pt2->AddText(lum); 
  //TPaveText *pt3 = new TPaveText(0.69,0.90,0.98,0.95,"NDC");
  TPaveText *pt3 = new TPaveText(ptLow,0.29,ptHigh,0.33,"NDC");
  pt3->SetFillColor(0);
  pt3->SetTextFont(42);
  char lum2[192];
  sprintf(lum2," #sqrt{s} = 8 TeV, L = %.2f fb^{-1}",5.261);
  pt3->AddText(lum2); 
	
  if(grid) c->SetGrid();
  if(logy)gPad->SetLogy();
   
  TH1F *hr = c->DrawFrame(105.0,yLow,180.0,yHigh);
	
  
  grshade_95->Draw("f");
  grshade_68->Draw("f");
  gr->Draw("C");
  if(isXSxBR)
    {
      gr_XSBR68->SetFillColor(kRed);
      gr_XSBR68->SetFillStyle(3013);
      gr_XSBR68->Draw("3same");
      gr_XSBR->SetLineColor(kRed);
      gr_XSBR->Draw("C");
    }
  if(addObsLimit)
    {
      if(points)grObs->Draw("CP");
      else grObs->Draw("C");
    }

  if(!isXSxBR)oneLine->Draw("LSAME");


  hr->GetXaxis()->SetTitle(xTitle);
  hr->GetYaxis()->SetTitle(yTitle);
  hr->GetYaxis()->SetTitleOffset(1.2);		
  
  
  c->Update();
  if(gridOnTop)gPad->RedrawAxis("g");
  box2->Draw();
  pt->Draw("SAME");
  pt2->Draw("SAME");
  pt4->Draw("SAME");


  sprintf( outfileName,"%s/UpperLimit_%s_%.0fTeV_lowMass_%s.eps",plotDir.c_str(),method.c_str(),sqrts, append.c_str() );
  //c->SaveAs(outfileName);
  sprintf( outfileName,"%s/UpperLimit_%s_%.0fTeV_lowMass_%s.png",plotDir.c_str(),method.c_str(),sqrts, append.c_str() );	
  //c->SaveAs(outfileName);


  // --------------- Full Mass Range ---------------- //
	
  //TLegend * box3 = new TLegend(0.4,0.72,0.8,0.92);
  TLegend * box3 = new TLegend(0.52,0.7,0.92,0.9);
  box3->SetFillColor(0);
  //box3->SetBorderSize(0);
  if (addObsLimit){ box3->AddEntry(grObs,"Observed Asym. CLs","l"); }
  //box3->AddEntry(gr,LimitType+" Fit Expected Asym. CLs","l");
  box3->AddEntry(gr,"Expected Asym. CLs","l");
  box3->AddEntry(grshade_68,"Expected #pm 1#sigma","lf");
  box3->AddEntry(grshade_95,"Expected #pm 2#sigma","lf");
  //box3->AddEntry(oneLine,"#sigma / #sigma_{SM}","l");
  

  TCanvas *cl = new TCanvas("cl","cl",canvasX,canvasY);
  cl->cd();
  if(grid) cl->SetGrid();
  
  TH1F *hrl = cl->DrawFrame(xLow,yLow,xHigh,yHigh);
  grshade_95->Draw("f");
  grshade_68->Draw("f");
  gr->Draw("C");
  if(isXSxBR)
    {
      gr_XSBR68->SetFillColor(kRed);
      gr_XSBR68->SetFillStyle(3013);
      gr_XSBR68->Draw("3same");
      gr_XSBR->SetLineColor(kRed);
      gr_XSBR->Draw("C");
    }
  if(addObsLimit)
    {
      if(points)grObs->Draw("CP");
      else grObs->Draw("C");
    }

 
  if(!isXSxBR)oneLine->Draw("CSAME");
 

  hrl->GetXaxis()->SetTitle(xTitle);
  hrl->GetYaxis()->SetTitle(yTitle);
  hrl->GetYaxis()->SetTitleOffset(1.2);		
  //hrl->GetXaxis()->SetXmax(xHigh);
  

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
  box3->Draw();

  pt->Draw("SAME");
  pt2->Draw("SAME");
  pt4->Draw("SAME");

  sprintf( outfileName,"%s/UpperLimit_%s_%.0fTeV_wholeMass_%s.eps",plotDir.c_str(),method.c_str(), sqrts, append.c_str() );
  //cl->SaveAs(outfileName);
  sprintf( outfileName,"%s/UpperLimit_%s_%.0fTeV_wholeMass_%s.png",plotDir.c_str(),method.c_str(), sqrts, append.c_str() );	
  //cl->SaveAs(outfileName);	

  // ---------------- Root File --------------- //

  //TFile outf("testUL_single.root","RECREATE");
  //outf.cd();
  //gr->Write();
  //outf.Write();

	
}



void getLimits(TFile *f, std::vector<double> &v_mh,std::vector<double> &v_mean,std::vector<double> &v_68l,std::vector<double> &v_68h,std::vector<double> &v_95l,std::vector<double> &v_95h,std::vector<double> &v_obs)
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
      if(_DEBUG_)cout << "mH: " << mh << " limit: " << limit << " quantileExpected: " << quant << endl;  
      if(quant>-1.01&&quant<-0.99)
	{
	  v_obs.push_back(limit);
	  v_mh.push_back(mh);
	}
      else if(quant>0.024 && quant<0.026) v_95l.push_back(limit);
      else if(quant>0.15  && quant<0.17 ) v_68l.push_back(limit);
      else if(quant>0.49  && quant<0.51 ) v_mean.push_back(limit);
      else if(quant>0.83  && quant<0.85 ) v_68h.push_back(limit);
      else if(quant>0.974 && quant<0.976) v_95h.push_back(limit);
      else {cout<<"Error! Unknown Quantile =  " << quant << endl;}
      
    }
  
}

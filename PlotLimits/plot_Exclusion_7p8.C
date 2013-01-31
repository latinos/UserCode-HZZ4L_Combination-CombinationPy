#include "TGraphAsymmErrors.h"
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

void getLimits(TFile *f, std::vector<double> &v_mh,std::vector<double> &v_mean,std::vector<double> &v_68l,std::vector<double> &v_68h,std::vector<double> &v_95l,std::vector<double> &v_95h,std::vector<double> &v_obs,std::vector<double> &v_obs95);

TGraphAsymmErrors *slidingWindowAverage(TGraphAsymmErrors *input, int slidingWindow);
TGraph *slidingWindowAverage2(TGraph *input, int slidingWindow);
TGraphAsymmErrors * removeGlitches(TGraphAsymmErrors *out);
TGraph * removeGlitches2(TGraph *out);


// --------- Inputs ------- //
TString inputFile = "../CreateDatacards/CMSSW_5_2_5/src/Moriond_low_v1/results/higgsCombineHZZ4L_ASCLS.root";
TString inputFile2 = "../CreateDatacards/CMSSW_5_2_5/src/Moriond_high_v1/results/higgsCombineHZZ4L_ASCLS.root";
TString inputFileHCP = "../CreateDatacards/CMSSW_5_2_5/src/HCP_MorLumi/results/higgsCombineHZZ4L_ASCLS.root";
const bool addObsLimit = false;
const bool isXSxBR = false;
const bool _DEBUG_ = false;
string method = "FREQ";
Double_t xLow = 99.9;
Double_t xHigh = 1001.0;
Double_t yLow = 0.1;
Double_t yHigh = 20.0;
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
//Double_t lumi = 1.616;
double sqrts = 7.0;
Double_t lumi = 5.051;
std::string plotDir = "plots";
std::string append = "2D_no2l2tau";
TString LimitType = "3D";
// ----------------------- //


using namespace std;

void plot_Exclusion_7p8()
{

  gROOT->ProcessLine(".x tdrstyle.cc");
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadTopMargin(0.05);

  if(isXSxBR)
    {
      gSystem->Load("../CreateDatacards/include/HiggsCSandWidth_cc");
      gSystem->Load("../CreateDatacards/include/HiggsCSandWidthSM4_cc");
      HiggsCSandWidth *myCSW = new HiggsCSandWidth();
      HiggsCSandWidthSM4 *myCSWSM4 = new HiggsCSandWidthSM4();
    }

  TFile *inFile = new TFile(inputFile,"READ");
  TFile *inFile2 = new TFile(inputFile2, "READ");
  TFile *inFileHCP = new TFile(inputFileHCP,"READ");
  
  // ------------------- Get Values -------------------- //

  vector<double> mH, Val_obs, Val_mean, Val_68h, Val_68l, Val_95h, Val_95l, Val_obs95;
  getLimits(inFile,mH,Val_mean,Val_68l,Val_68h,Val_95l,Val_95h,Val_obs,Val_obs95);
  getLimits(inFile2,mH,Val_mean,Val_68l,Val_68h,Val_95l,Val_95h,Val_obs,Val_obs95);
  vector<double> v_masses, v_means, v_lo68, v_hi68, v_lo95, v_hi95, v_obs;
  vector<double> expExclusion,obsExclusion,expExcl95,obsExcl95;
  for(unsigned int i = 1; i < mH.size(); i++)
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
      if( max( Val_95h[i], Val_95l[i])< 1.0)expExcl95.push_back(mH[i]);
      if(Val_obs95[i]<1.0)obsExcl95.push_back(mH[i]);
    }

  vector<double> mH2, Val_obs2, Val_mean2, Val_68h2, Val_68l2, Val_95h2, Val_95l2, Val_obs952;
  getLimits(inFileHCP,mH2,Val_mean2,Val_68l2,Val_68h2,Val_95l2,Val_95h2,Val_obs2,Val_obs952);
  //getLimits(inFile2,mH,Val_mean,Val_68l,Val_68h,Val_95l,Val_95h,Val_obs,Val_obs95);
  vector<double> v_masses2, v_means2, v_lo682, v_hi682, v_lo952, v_hi952, v_obs2;
  vector<double> expExclusion2,obsExclusion2,expExcl952,obsExcl952;
  for(unsigned int i = 1; i < mH2.size(); i++)
    {
      v_masses2.push_back( mH2[i] );
      v_means2.push_back( Val_mean2[i] );
      v_lo682.push_back( min( Val_68l2[i], Val_68h2[i]) );
      v_hi682.push_back( max( Val_68h2[i], Val_68l2[i]) );
      v_lo952.push_back( min( Val_95l2[i], Val_95h2[i]) );
      v_hi952.push_back( max( Val_95h2[i], Val_95l2[i]) );
      v_obs2.push_back(Val_obs2[i]);
      if(Val_mean2[i] < 1.0) expExclusion2.push_back(mH2[i]);
      if(Val_obs2[i] < 1.0 && addObsLimit) obsExclusion2.push_back(mH2[i]);
      if( max( Val_95h2[i], Val_95l[i])< 1.0)expExcl952.push_back(mH2[i]);
      if(Val_obs952[i]<1.0)obsExcl952.push_back(mH2[i]);
    }
  //this is because point at 1TeV is read first->TBC
  /*v_masses.push_back( mH[0] );
  v_means.push_back( Val_mean[0] );
  v_lo68.push_back( min( Val_68l[0], Val_68h[0]) );
  v_hi68.push_back( max( Val_68h[0], Val_68l[0]) );
  v_lo95.push_back( min( Val_95l[0], Val_95h[0]) );
  v_hi95.push_back( max( Val_95h[0], Val_95l[0]) );
  v_obs.push_back(Val_obs[0]);
  if(Val_mean[0] < 1.0) expExclusion.push_back(mH[0]);
  if(Val_obs[0] < 1.0 && addObsLimit) obsExclusion.push_back(mH[0]);
  if(v_hi95[0] < 1.0)expExcl95.push_back(mH[0]);
  if(Val_obs95[0]<1.0)obsExcl95.push_back(mH[0]);
  */
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
  double a_masses[sizeV], a_means[sizeV], a_lo68[sizeV], a_hi68[sizeV], a_lo95[sizeV], a_hi95[sizeV], a_obs[sizeV],a_zero[sizeV];
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
      a_zero[nMassEff] = 0;
      //cout << v_masses[m] << "  " << v_means[m] << endl;
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

  int nMassEff2=0;
  int nExcluded2=0;
  const int sizeV2 = v_masses2.size();
  double a_masses2[sizeV2], a_means2[sizeV2], a_lo682[sizeV2], a_hi682[sizeV2], a_lo952[sizeV2], a_hi952[sizeV2], a_obs2[sizeV2],a_zero2[sizeV2];
  for(unsigned int m = 0; m < v_masses2.size(); m++)
    {
      if(v_hi682.at(m)>=v_hi952.at(m) || v_lo682.at(m)<=v_lo952.at(m))
	{
	  cout << "Point at M = " << v_masses2.at(m) << " excluded" << endl;
	  nExcluded2++;
	  continue;
	}
      
      a_masses2[nMassEff2] = v_masses2[m];
      a_means2[nMassEff2] = v_means2[m];
      a_lo682[nMassEff2] = v_lo682[m];
      a_hi682[nMassEff2] = v_hi682[m];
      a_lo952[nMassEff2] = v_lo952[m];
      a_hi952[nMassEff2] = v_hi952[m];
      a_obs2[nMassEff2] = v_obs2[m];
      a_zero2[nMassEff2] = 0;

      nMassEff2++;
      
    }
  cout << "Excluded2 " << nExcluded2 << " sick mass points!" << endl;

  // --------------- Excluded Regions --------------- //
  cout<<"*******EXCLUSION*********"<<endl;
  for(int p = 0; p < expExclusion.size(); p++)
    {
      cout << "Expected Exclusion: " <<  expExclusion[p] << endl;
    }

  for(int s=0;s<expExcl95.size();s++)cout<<"Expected Excluxion 95% "<< expExcl95[s]<<endl;

  if(addObsLimit)
    {
      for(int q = 0; q < obsExclusion.size(); q++)
	{
	  cout << "Observed Exclusion: " <<  obsExclusion[q] << endl;
	}
      for(int t=0;t<obsExcl95.size();t++)cout<<"Observed Excluxion 95% "<< obsExcl95[t]<<endl;
    }

  cout<<endl;
  cout<<"*******LIMITS*********"<<endl;

  for(int r=0;r<nMassEff;r++){
    cout<<"Expected Limit "<<a_masses[r]<<": "<<a_means[r]<<endl;
  }

 if(addObsLimit)
    {
      for(int q2 = 0; q2 < Val_obs.size(); q2++)
	{
	  cout << "Observed limit " <<a_masses[q2]<<": "<< Val_obs[q2] << endl;
	}
      for(int t2=0;t2<Val_obs95.size();t2++)cout<<"Observed limit 95% " <<a_masses[t2]<<": "<< Val_obs95[t2] << endl;
    }
  
  // ------------------- Draw  -------------------- //



  TCanvas *c = new TCanvas("c","c",canvasX,canvasY);
  TGraphAsymmErrors* gr = new TGraphAsymmErrors(nMassEff, a_masses, a_means);
  TGraphAsymmErrors* grshade_68 = new TGraphAsymmErrors(2*nMassEff);
  TGraphAsymmErrors* grshade_95 = new TGraphAsymmErrors(2*nMassEff);
  gr->SetName("Expected");gr->SetTitle("Expected");
  grshade_68->SetName("68");grshade_68->SetTitle("68");
  grshade_95->SetName("95");grshade_95->SetTitle("95");
  TGraph *grObs = new TGraph(nMassEff, a_masses, a_obs);
  grObs->SetLineWidth(3);
  grObs->SetLineColor(kBlack);
  grObs->SetMarkerStyle(21);
  grObs->SetMarkerSize(0.8);
 
  //cout<<"nMasses: "<<nMassEff<<endl;
  for (int a = 0; a < nMassEff; a++)
    {
      grshade_68->SetPoint(a,a_masses[a],a_lo68[a]);
      grshade_68->SetPoint(sizeV+a,a_masses[nMassEff-a-1],a_hi68[nMassEff-a-1]);
      //cout<<"setting point "<<sizeV+a<<endl;
      grshade_95->SetPoint(a,a_masses[a],a_lo95[a]);
      grshade_95->SetPoint(nMassEff+a,a_masses[nMassEff-a-1],a_hi95[nMassEff-a-1]);
    }

  TGraphAsymmErrors* gr2 = new TGraphAsymmErrors(nMassEff2, a_masses2, a_means2);
  TGraphAsymmErrors* grshade_682 = new TGraphAsymmErrors(2*nMassEff2);
  TGraphAsymmErrors* grshade_952 = new TGraphAsymmErrors(2*nMassEff2);
  gr2->SetName("Expected2");gr2->SetTitle("Expected2");
  grshade_682->SetName("68_2");grshade_682->SetTitle("68_2");
  grshade_952->SetName("95_2");grshade_952->SetTitle("95_2");
  TGraph *grObs2 = new TGraph(nMassEff2, a_masses2, a_obs2);
  grObs2->SetLineWidth(3);
  grObs2->SetLineColor(kBlack);
  grObs2->SetMarkerStyle(21);
  grObs2->SetMarkerSize(0.8);
 
  //cout<<"nMasses: "<<nMassEff<<endl;
  for (int a = 0; a < nMassEff2; a++)
    {
      grshade_682->SetPoint(a,a_masses2[a],a_lo682[a]);
      grshade_682->SetPoint(sizeV2+a,a_masses2[nMassEff2-a-1],a_hi682[nMassEff2-a-1]);
      cout<<"setting point "<<sizeV2+a<<endl;
      grshade_952->SetPoint(a,a_masses2[a],a_lo952[a]);
      grshade_952->SetPoint(nMassEff2+a,a_masses2[nMassEff2-a-1],a_hi952[nMassEff2-a-1]);
    }
  //for(int a = 0; a < nMassEff; a++)cout<<"mass: "<<a_masses[a]<<endl;
  //TGraphAsymmErrors* grshade_68m = removeGlitches(grshade_68n);
  //TGraphAsymmErrors* grshade_95m = removeGlitches(grshade_95n);
  //TGraphAsymmErrors* gr = removeGlitches(grn);

  //TGraphAsymmErrors* grshade_68 = slidingWindowAverage(grshade_68m,2);
  //TGraphAsymmErrors* grshade_95 = slidingWindowAverage(grshade_95m,2);
  //TGraph* gr = slidingWindowAverage2(grm,2);
  //gr->Sort();
  //grshade_68->Sort();
  //grshade_95->Sort();
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

  gr2->SetLineStyle(2);
  gr2->SetLineWidth(3);
  gr2->SetLineColor(kBlack);
  grshade_682->SetFillColor(6);
  grshade_952->SetFillColor(7);		
  grshade_682->SetLineStyle(2);
  grshade_952->SetLineStyle(2);
  grshade_682->SetLineWidth(3);
  grshade_952->SetLineWidth(3);
  grshade_682->SetLineColor(6);
  grshade_682->SetFillStyle(3013);
  grshade_952->SetLineColor(7);
	
	
  char outfileName[192];

  TF1* oneLine = new TF1("oneLine","1",xLow,xHigh);
  oneLine->SetLineColor(kRed);
  oneLine->SetLineWidth(3);

  // --------------- Low Mass Zoom -------------- //
	
  TLegend * box2 = new TLegend(0.56,0.7,0.85,0.9);
  box2->SetFillColor(0);
  //box2->SetBorderSize(0);
  if (addObsLimit){ box2->AddEntry(grObs,"Observed","l"); }
  //box2->AddEntry(gr,LimitType+" Fit Expected Asym. CLs","l");
  box2->AddEntry(gr,"Expected Moriond","l");
  box2->AddEntry(gr2,"Expected HCP at MoriondLumi","l");
  //box2->AddEntry(grshade_68,"68% expectation","f");
  //box2->AddEntry(grshade_95,"95% expectation","f");
  box2->AddEntry(grshade_68,"Expected #pm 1#sigma","lf");
  box2->AddEntry(grshade_95,"Expected #pm 2#sigma","lf");
  box2->AddEntry(grshade_682,"Expected #pm 1#sigma HCP at MoriondLumi","lf");
  //box2->AddEntry(oneLine,"#sigma / #sigma_{SM}","l");


  double ptLow= 0.28, ptHigh = 0.54;
		
  TPaveText *pt = new TPaveText(ptLow,0.84,ptHigh,0.88,"NDC");
  pt->SetFillColor(0);
  pt->SetTextFont(42);
  pt->AddText("CMS Preliminary");
  TPaveText *pt4 = new TPaveText(ptLow,0.8,ptHigh,0.84,"NDC");
  pt4->SetFillColor(0);
  pt4->SetTextFont(42);
  pt4->AddText("H #rightarrow ZZ #rightarrow 4L Moriond strategy"); 

  //TPaveText *pt2 = new TPaveText(0.69,0.94,0.98,0.99,"NDC");
  TPaveText *pt2 = new TPaveText(ptLow,0.76,ptHigh,0.8,"NDC");
  pt2->SetFillColor(0);
  pt2->SetTextFont(42);
  char lum[192];
  sprintf(lum," #sqrt{s} = 7 TeV, L = %.2f fb^{-1}",5.051);
  pt2->AddText(lum); 
  //TPaveText *pt3 = new TPaveText(0.69,0.90,0.98,0.95,"NDC");
  TPaveText *pt3 = new TPaveText(ptLow,0.72,ptHigh,0.76,"NDC");
  pt3->SetFillColor(0);
  pt3->SetTextFont(42);
  char lum2[192];
  sprintf(lum2," #sqrt{s} = 8 TeV, L = %.2f fb^{-1}",19.62);
  pt3->AddText(lum2); 

  if(grid) c->SetGrid();
   
  TH1F *hr = c->DrawFrame(105.0,yLow,180.0,yHigh);
	
  gr->Sort();
  //grshade_68->Sort();
  //grshade_95->Sort();
  
  grshade_95->Draw("f");
  grshade_68->Draw("f");
  grshade_682->Draw("f");
  gr->Draw("C");
  gr2->Draw("C");
  if(isXSxBR)
    {
      gr_XSBR68->Sort();
      gr_XSBR68->SetFillColor(kRed);
      gr_XSBR68->SetFillStyle(3013);
      gr_XSBR68->Draw("3same");
      gr_XSBR->SetLineColor(kRed);
      gr_XSBR->Draw("C");
    }
  if(addObsLimit)
    {
      grObs->Sort();
      if(points)grObs->Draw("CP");
      else grObs->Draw("C");
    }

  if(!isXSxBR)oneLine->Draw("LSAME");


  hr->GetXaxis()->SetTitle(xTitle);
  hr->GetYaxis()->SetTitle(yTitle);
  hr->GetYaxis()->SetTitleOffset(1.2);		
  if(logy)gPad->SetLogy();
  
  
  c->Update();
  if(gridOnTop)gPad->RedrawAxis("g");
  pt->Draw("SAME");
  pt2->Draw("SAME");
  pt3->Draw("SAME");
  pt4->Draw("SAME");

  box2->Draw();
  sprintf( outfileName,"%s/UpperLimit_%s_7p8TeV_lowMass_%s.eps",plotDir.c_str(),method.c_str(), append.c_str() );
  c->SaveAs(outfileName);
  sprintf( outfileName,"%s/UpperLimit_%s_7p8TeV_lowMass_%s.png",plotDir.c_str(),method.c_str(), append.c_str() );	
  c->SaveAs(outfileName);
  sprintf( outfileName,"%s/UpperLimit_%s_7p8TeV_lowMass_%s.root",plotDir.c_str(),method.c_str(), append.c_str() );
  c->SaveAs(outfileName);

  // --------------- Full Mass Range ---------------- //
	
  //TLegend * box3 = new TLegend(0.64,0.67,0.94,0.87);
  TLegend * box3 = new TLegend(0.56,0.7,0.86,0.9);
    box3->SetFillColor(0);
  //box3->SetBorderSize(0);
  if (addObsLimit){ box3->AddEntry(grObs,"Observed","l"); }
  //box3->AddEntry(gr,LimitType+" Fit Expected Asym. CLs","l");
  box3->AddEntry(gr,"Expected","l");
  box3->AddEntry(gr2,"Expected HCP at MoriondLumi","l");
  box3->AddEntry(grshade_68,"Expected #pm 1#sigma","lf");
  box3->AddEntry(grshade_95,"Expected #pm 2#sigma","lf");
  box3->AddEntry(grshade_682,"Expected #pm 1#sigma HCP at MoriondLumi","lf");
  //box3->AddEntry(oneLine,"#sigma / #sigma_{SM}","l");
  

  TCanvas *cl = new TCanvas("cl","cl",canvasX,canvasY);
  cl->cd();
  if(grid) cl->SetGrid();
  
  TH1F *hrl = cl->DrawFrame(xLow,yLow,xHigh,yHigh-5);

  grshade_95->Draw("f");
  grshade_68->Draw("f");

  grshade_682->Draw("f");
  gr->Draw("C");
  gr2->Draw("C");

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
      double dyh = (yHigh-5) * 0.08;
      double dyl = yLow * 0.08; //fabs(c1->PixeltoY(c1->VtoPixel(0.95)) - c1->PixeltoY(c1->VtoPixel(0.94)));
      if (gPad->GetLogy() && log((yHigh-5)/yLow) > log(1e6)) { dyh *= 2; dyl *= 2; }
      if (gPad->GetLogy() == 0) { dyh = dyl = 0.01*((yHigh-5)-yLow); }
      if (isTiny) { dyh *= 2; dyl *= 2; }
      for (int j = 100; j < 600; j += 10)  {
	  if (j > 400 && j % 20 == 10) continue;
	  tick.DrawLine(j, yLow, j, yLow+(j % 100 == 0 ? 2*dyl : dyl));
	  tick.DrawLine(j, (yHigh-5), j, yHigh-5-(j % 100 == 0 ? 2*dyh : dyh));
      }
    }
  cl->Update();
  hrl->GetXaxis()->SetRangeUser(xLow,xHigh);
  cl->Update();

  if(gridOnTop)gPad->RedrawAxis("g");

  pt->Draw("SAME");
  pt2->Draw("SAME");
  pt3->Draw("SAME");
  pt4->Draw("SAME");
  
  box3->Draw();
  sprintf( outfileName,"%s/UpperLimit_%s_7p8TeV_wholeMass_%s.eps",plotDir.c_str(),method.c_str(), append.c_str() );
  cl->SaveAs(outfileName);
  sprintf( outfileName,"%s/UpperLimit_%s_7p8TeV_wholeMass_%s.png",plotDir.c_str(),method.c_str(), append.c_str() );	
  cl->SaveAs(outfileName);
  sprintf( outfileName,"%s/UpperLimit_%s_7p8TeV_wholeMass_%s.root",plotDir.c_str(),method.c_str(), append.c_str() );
  cl->SaveAs(outfileName);	

  // ---------------- Root File --------------- //

  TFile outf("testUL_single.root","RECREATE");
  outf.cd();
  gr->Write();
  outf.Write();

	
}



void getLimits(TFile *f, std::vector<double> &v_mh,std::vector<double> &v_mean,std::vector<double> &v_68l,std::vector<double> &v_68h,std::vector<double> &v_95l,std::vector<double> &v_95h,std::vector<double> &v_obs,std::vector<double> &v_obs95)
{

  TTree *tree =(TTree*)f->Get("limit");
  
  double mh,limit,errObs;
  float quant;
  tree->SetBranchAddress("mh",&mh);
  tree->SetBranchAddress("limit",&limit);
  tree->SetBranchAddress("quantileExpected",&quant);
  tree->SetBranchAddress("limitErr",&errObs);

  int nentry=tree->GetEntries();

  for(int i=0;i<nentry;i++)
    {
      tree->GetEntry(i);
      if(_DEBUG_)cout << "mH: " << mh << " limit: " << limit << " quantileExpected: " << quant << endl;  
      if(quant>-1.01&&quant<-0.99)
	{
	  v_obs.push_back(limit);
	  v_mh.push_back(mh);
	  v_obs95.push_back(limit+2*errObs/limit);
	}
      else if(quant>0.024 && quant<0.026) v_95l.push_back(limit);
      else if(quant>0.15  && quant<0.17 ) v_68l.push_back(limit);
      else if(quant>0.49  && quant<0.51 ) v_mean.push_back(limit);
      else if(quant>0.83  && quant<0.85 ) v_68h.push_back(limit);
      else if(quant>0.974 && quant<0.976) v_95h.push_back(limit);
      else {cout<<"Error! Unknown Quantile =  " << quant << endl;}
      
    }
  
}






TGraphAsymmErrors *slidingWindowAverage(TGraphAsymmErrors *input, int slidingWindow) {

  //cout << 11 << endl;

    TGraphAsymmErrors *out  = (TGraphAsymmErrors *) input->Clone();
    //bool isLogPlot = TString(input->GetName()).Contains("smcls");

    //cout << 22 << endl;

    bool isLogPlot = true;
    for (int i = 0, n = input->GetN(); i < n; ++i) {
        /*if (i >= 1 && i < n-1 && (input->GetX()[i] < input->GetX()[i-1] || input->GetX()[i] > input->GetX()[i+1])) {
            out->GetX()[i] = 0.5*(input->GetX()[i-1] + input->GetX()[i+1]);
            out->GetY()[i] = 0.5*(input->GetY()[i-1] + input->GetY()[i+1]);
            out->GetEYlow()[i] = 0.5*(input->GetEYlow()[i-1] + input->GetEYlow()[i+1]);
            out->GetEYhigh()[i] = 0.5*(input->GetEYhigh()[i-1] + input->GetEYhigh()[i+1]);
            continue;
        }*/

      //cout << 33 << endl;

        double y0 = input->GetY()[i];
        double sum = 0, sumhi = 0, sumlo = 0, sumw = 0;
        for (int j = i-slidingWindow; j <= i+slidingWindow; ++j) 
	  {

            if (j < 0 || j >= n) continue;
            double y = input->GetY()[j], w = 1.0; // /(1.0 + abs(i-j));
	    //cout << j << "  " << y0 << "  " << y << "  " << slidingWindow << endl;
            if (isLogPlot) 
	      {
		if (fabs(log(y0/y)) > log(2)) continue;    
	      }
	    else {
	      if (fabs(y0-y) > 0.1*y0) continue;    
            }
            if (isLogPlot) 
	      {
                sum   += w*log(y);
                sumlo += w*log(input->GetEYlow()[j]);
                sumhi += w*log(input->GetEYhigh()[j]);
	      } 
	    else {
	      sum   += w*y;
	      sumlo += w*input->GetEYlow()[j];
	      sumhi += w*input->GetEYhigh()[j];
            }
            sumw  += w;
	  }
        if (sumw == 0) continue;
        if (isLogPlot) 
	  {
	    out->GetY()[i] = exp(sum/sumw);
	    out->GetEYlow()[i] = exp(sumlo/sumw);
	    out->GetEYhigh()[i] = exp(sumhi/sumw);
	  } 
	else {
	  out->GetY()[i] = sum/sumw;
	  out->GetEYlow()[i] = sumlo/sumw;
	  out->GetEYhigh()[i] = sumhi/sumw;
        }
    }
    return out;
}

TGraph *slidingWindowAverage2(TGraph *input, int slidingWindow) {

  cout << 11 << endl;

    TGraph *out  = (TGraph*) input->Clone();
    //bool isLogPlot = TString(input->GetName()).Contains("smcls");

    cout << 22 << endl;

    bool isLogPlot = true;
    for (int i = 0, n = input->GetN(); i < n; ++i) {
        /*if (i >= 1 && i < n-1 && (input->GetX()[i] < input->GetX()[i-1] || input->GetX()[i] > input->GetX()[i+1])) {
            out->GetX()[i] = 0.5*(input->GetX()[i-1] + input->GetX()[i+1]);
            out->GetY()[i] = 0.5*(input->GetY()[i-1] + input->GetY()[i+1]);
            out->GetEYlow()[i] = 0.5*(input->GetEYlow()[i-1] + input->GetEYlow()[i+1]);
            out->GetEYhigh()[i] = 0.5*(input->GetEYhigh()[i-1] + input->GetEYhigh()[i+1]);
            continue;
        }*/

      cout << 33 << endl;

        double y0 = input->GetY()[i];
        double sum = 0, sumhi = 0, sumlo = 0, sumw = 0;
        for (int j = i-slidingWindow; j <= i+slidingWindow; ++j) 
	  {

            if (j < 0 || j >= n) continue;
            double y = input->GetY()[j], w = 1.0; // /(1.0 + abs(i-j));
	    cout << j << "  " << y0 << "  " << y << "  " << slidingWindow << endl;
            if (isLogPlot) 
	      {
		if (fabs(log(y0/y)) > log(2)) continue;    
	      }
	    else {
	      if (fabs(y0-y) > 0.1*y0) continue;    
            }
            if (isLogPlot) 
	      {
                sum   += w*log(y);
                sumlo += w*log(input->GetEYlow()[j]);
                sumhi += w*log(input->GetEYhigh()[j]);
	      } 
	    else {
	      sum   += w*y;
	      sumlo += w*input->GetEYlow()[j];
	      sumhi += w*input->GetEYhigh()[j];
            }
            sumw  += w;
	  }
        if (sumw == 0) continue;
        if (isLogPlot) 
	  {
	    out->GetY()[i] = exp(sum/sumw);
	    out->GetEYlow()[i] = exp(sumlo/sumw);
	    out->GetEYhigh()[i] = exp(sumhi/sumw);
	  } 
	else {
	  out->GetY()[i] = sum/sumw;
	  out->GetEYlow()[i] = sumlo/sumw;
	  out->GetEYhigh()[i] = sumhi/sumw;
        }
    }
    return out;
}





TGraphAsymmErrors * removeGlitches(TGraphAsymmErrors *out) {
  do {
    int bad = -1; bool hasgood = false;
    for (int i = 0; i < out->GetN(); ++i) {
      if (out->GetEYlow()[i] == 0 && out->GetEYhigh()[i] == 0) {
	bad = i;
      } else {
	hasgood = true;
      }
    }
    if (!hasgood) return out;
    if (bad == -1) return out;
    out->RemovePoint(bad);
  } while (1);
}

TGraph * removeGlitches2(TGraph *out) {
  do {
    int bad = -1; bool hasgood = false;
    for (int i = 0; i < out->GetN(); ++i) {
      if (out->GetEYlow()[i] == 0 && out->GetEYhigh()[i] == 0) {
	bad = i;
      } else {
	hasgood = true;
      }
    }
    if (!hasgood) return out;
    if (bad == -1) return out;
    out->RemovePoint(bad);
  } while (1);
}

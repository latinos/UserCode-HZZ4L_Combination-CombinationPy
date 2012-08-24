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
#include "TSpline.h"

void getLimits(TFile *f, std::vector<double> &v_mh,std::vector<double> &v_mean,std::vector<double> &v_68l,std::vector<double> &v_68h,std::vector<double> &v_95l,std::vector<double> &v_95h,std::vector<double> &v_obs);

// --------- Inputs ------- //
TString inputFile = "higgsCombine_ASCLS_1D_FSR_7TeV.root";
TString inputFile2D = "higgsCombine_ASCLS_1D_FSR_7TeV.root";
TString inputFilePRL = "higgsCombine_ASCLS_PRL_5.021_06.06.root";
const bool addObsLimit = false;
const bool isXSxBR = false;
const bool _DEBUG_ = false;
string method = "ASCLS";
Double_t xLow = 99.9;
Double_t xHigh = 601.0;
Double_t yLow = 0.1;
Double_t yHigh = 15.0;
TString xTitle = "Higgs boson mass, GeV";
TString yTitle = "95% CL limit on #sigma/#sigma_{SM}";
const bool logy = true;
const bool logx = true;
const bool grid = true;
const bool gridOnTop = true;
const bool points = true;
const bool isTiny = false;
int canvasX = 900;
int canvasY = 700;
double sqrts = 7.0;
Double_t lumi = 5.051;
std::string plotDir = "Limits_App";
TString TplotDir = plotDir;
// ----------------------- //


using namespace std;

void plot_Exclusion_PRL1d2d()
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
  TFile *inFile2D = new TFile(inputFile2D,"READ");
  TFile *inFilePRL = new TFile(inputFilePRL,"READ");

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
    }



  vector<double> mH2, Val_obs2, Val_mean2, Val_68h2, Val_68l2, Val_95h2, Val_95l2;
  getLimits(inFile2D,mH2,Val_mean2,Val_68l2,Val_68h2,Val_95l2,Val_95h2,Val_obs2);
  vector<double> v_masses2, v_means2, v_lo682, v_hi682, v_lo952, v_hi952, v_obs2;
   for(unsigned int i = 0; i < mH2.size(); i++)
    {
      v_masses2.push_back( mH2[i] );
      v_means2.push_back( Val_mean2[i] );
      // v_obs2.push_back(Val_obs2[i]);
    }


  vector<double> mH3, Val_obs3, Val_mean3, Val_68h3, Val_68l3, Val_95h3, Val_95l3;
  getLimits(inFilePRL,mH3,Val_mean3,Val_68l3,Val_68h3,Val_95l3,Val_95h3,Val_obs3);
  vector<double> v_masses3, v_means3, v_lo683, v_hi683, v_lo953, v_hi953, v_obs3;
   for(unsigned int i = 0; i < mH3.size(); i++)
    {
      v_masses3.push_back( mH3[i] );
      v_means3.push_back( Val_mean3[i] );
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
	  cout << "Observed Exclusion: " <<  obsExlusion[q] << endl;
	}
    }

  // ------------------- Draw  -------------------- //


  int nMassEff2=0;
  const int sizeV2 = v_masses2.size();
  double a_masses2[sizeV], a_means2[sizeV];
  for(unsigned int m = 0; m < v_masses2.size(); m++)
    {
      a_masses2[nMassEff2] = v_masses2[m];
      a_means2[nMassEff2] = v_means2[m];
      nMassEff2++;
    }

  int nMassEff3=0;
  const int sizeV3 = v_masses3.size();
  double a_masses3[sizeV], a_means3[sizeV];
  for(unsigned int m = 0; m < v_masses3.size(); m++)
    {
      a_masses3[nMassEff3] = v_masses3[m];
      a_means3[nMassEff3] = v_means3[m];
      nMassEff3++;
    }


  TCanvas *c = new TCanvas("c","c",canvasX,canvasY);
  c->cd();
  TPad *pad2 = new TPad("pad2", "foo", 0, 0,   1, 0.3);
  TPad *pad1 = new TPad("pad1", "foo", 0, 0.3, 1, 1);
  pad2->SetBottomMargin(0.22);
  //pad1->Draw();
  //pad2->Draw();
  pad2->SetGrid();

  TGraph *gr = new TGraph(nMassEff, a_masses, a_means);
  TGraph *gr2D = new TGraph(nMassEff2, a_masses2, a_means2);
  TGraph *grPRL = new TGraph(nMassEff3, a_masses3, a_means3);
  TGraph* grshade_68 = new TGraph(2*nMassEff);
  TGraph* grshade_95 = new TGraph(2*nMassEff);
  TGraph *grObs = new TGraph(nMassEff, a_masses, a_obs);
  grObs->SetLineWidth(3);
  grObs->SetLineColor(kBlack);
  grObs->SetMarkerStyle(21);
  grObs->SetMarkerSize(0.8);


  int nRatio = 0, nRatioPRL = 0;
  int nRatio2 = 0, nRatioPRL2 = 0;
  const int nME2 = nMassEff2;
  const int nME3 = nMassEff3;
  double ratio2D1D[nME2],ratio1DPRL[nME3];
  double ratio2D1D2[nME2],ratio1DPRL2[nME3];
  for(int qq = 0; qq < nME2-4; qq++)
    {
      ratio2D1D[nRatio] = 100*(grPRL->Eval(a_masses2[qq],0,"S")-gr2D->Eval(a_masses2[qq],0,"S"))/grPRL->Eval(a_masses2[qq],0,"S");
      nRatio++;
    }
  for(int rr = 0; rr < nME3-4; rr++)
    {
      ratio1DPRL[nRatioPRL] = 100*(grPRL->Eval(a_masses2[rr],0,"S")-gr->Eval(a_masses2[rr],0,"S"))/grPRL->Eval(a_masses2[rr],0,"S");
      nRatioPRL++;
    }

 for(int qq = 0; qq < nME2-4; qq++)
    {
      ratio2D1D2[nRatio2] = gr2D->Eval(a_masses2[qq],0,"S")/grPRL->Eval(a_masses2[qq],0,"S");
      nRatio2++;
    }
  for(int rr = 0; rr < nME3-4; rr++)
    {
      ratio1DPRL2[nRatioPRL2] = gr->Eval(a_masses2[rr],0,"S")/grPRL->Eval(a_masses2[rr],0,"S");
      nRatioPRL2++;
    }

  TGraph *grRatio = new TGraph(nRatio, a_masses2, ratio2D1D);
  TGraph *grRatioPRL = new TGraph(nRatioPRL, a_masses2, ratio1DPRL);

  TGraph *grRatio2 = new TGraph(nRatio2, a_masses2, ratio2D1D2);
  TGraph *grRatioPRL2 = new TGraph(nRatioPRL2, a_masses2, ratio1DPRL2);

  gr->SetLineStyle(2);
  gr->SetLineWidth(3);
  gr->SetLineColor(kBlue);
  gr2D->SetLineStyle(1);
  gr2D->SetLineWidth(3);
  gr2D->SetLineColor(kBlack);
  grPRL->SetLineStyle(1);
  grPRL->SetLineWidth(3);
  grPRL->SetLineColor(kRed);
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
	
  TLegend * box2 = new TLegend(0.35,0.7,0.75,0.90);
  box2->SetFillColor(0);
  //box2->SetBorderSize(0);
  if (addObsLimit){ box2->AddEntry(grObs,"Observed Asym. CLs","l"); }
  box2->AddEntry(grPRL,"PRL Expected Asym. CLs","l");
  box2->AddEntry(gr,"1D Fit Expected Asym. CLs","l");
  box2->AddEntry(gr2D,"2D Fit Expected Asym. CLs","l");
  //box2->AddEntry(grshade_68,"68% expectation","f");
  //box2->AddEntry(grshade_95,"95% expectation","f");
  box2->AddEntry(grshade_68,"Expected #pm 1#sigma","lf");
  box2->AddEntry(grshade_95,"Expected #pm 2#sigma","lf");
  //box2->AddEntry(oneLine,"#sigma / #sigma_{SM}","l");

		
  TPaveText *pt = new TPaveText(0.16,0.95,0.44,0.99,"NDC");
  pt->SetFillColor(0);
  pt->SetTextFont(42);
  pt->AddText("CMS Preliminary 2012"); 
  TPaveText *pt2 = new TPaveText(0.69,0.95,0.98,0.99,"NDC");
  pt2->SetFillColor(0);
  pt2->SetTextFont(42);
  char lum[192];
  sprintf(lum," #sqrt{s} = %.0f TeV, L = %.2f fb^{-1}",sqrts,lumi);
  pt2->AddText(lum); 
	
  if(grid) pad1->SetGrid();

  pad1->cd();
  TH1F *hr = pad1->DrawFrame(105.0,yLow,180.0,yHigh);
	
  
  grshade_95->Draw("f");
  grshade_68->Draw("f");
  gr->Draw("L");
  gr2D->Draw("L");
  grPRL->Draw("L");
  if(isXSxBR)
    {
      gr_XSBR68->SetFillColor(kRed);
      gr_XSBR68->SetFillStyle(3013);
      gr_XSBR68->Draw("3same");
      gr_XSBR->SetLineColor(kRed);
      gr_XSBR->Draw("L");
    }
  if(addObsLimit)
    {
      if(points)grObs->Draw("LP");
      else grObs->Draw("L");
    }

  if(!isXSxBR)oneLine->Draw("LSAME");

  pt->Draw("SAME");
  pt2->Draw("SAME");
	
  hr->GetXaxis()->SetTitle(xTitle);
  hr->GetYaxis()->SetTitle(yTitle);
  hr->GetYaxis()->SetTitleOffset(1.2);		

  pad2->cd();
  TH1F *hr2 = pad2->DrawFrame(105.0,yLow,180.0,yHigh);
  hr2->GetYaxis()->SetTitleOffset(0.06);

  ///DRAW RATIOs

  
  c->Update();
  if(gridOnTop)gPad->RedrawAxis("g");
  box2->Draw();
  sprintf( outfileName,"%s/UpperLimit_%s_lowMass_PRL.eps",plotDir.c_str(),method.c_str() );
  c->SaveAs(outfileName);
  //sprintf( outfileName,"%s/UpperLimit_%s_lowMass_PRL.png",plotDir.c_str(),method.c_str() );
  //c->SaveAs(outfileName);
  //sprintf( outfileName,"plots/UpperLimit_%s_lowMass_2D.C",method.c_str() );	
  //c->SaveAs(outfileName);


  // --------------- Full Mass Range ---------------- //
	
  TLegend * box3 = new TLegend(0.4,0.72,0.8,0.92);
  box3->SetFillColor(0);
  //box3->SetBorderSize(0);
  if (addObsLimit){ box3->AddEntry(grObs,"Observed Asym. CLs","l"); }
  box3->AddEntry(grPRL,"PRL Expected Asym. CLs","l");
  box3->AddEntry(gr,"1D Fit Expected Asym. CLs","l");
  box3->AddEntry(gr2D,"2D Fit Expected Asym. CLs","l");
  box3->AddEntry(grshade_68,"Expected #pm 1#sigma","lf");
  box3->AddEntry(grshade_95,"Expected #pm 2#sigma","lf");
  //box3->AddEntry(oneLine,"#sigma / #sigma_{SM}","l");
  

  TCanvas *cl = new TCanvas("cl","cl",canvasX,canvasY);
  cl->cd();


  TPad *padl2 = new TPad("padl2", "foo", 0, 0,   1, 0.3);
  TPad *padl1 = new TPad("padl1", "foo", 0, 0.3, 1, 1);

  padl2->SetBottomMargin(0.22);

  //padl1->Draw();
  //padl2->Draw();

  if(grid) padl1->SetGrid();
  padl1->cd();
  TH1F *hrl = padl1->DrawFrame(xLow,yLow,xHigh,yHigh);

  grshade_95->Draw("f");
  grshade_68->Draw("f");
  gr->Draw("L");
  gr2D->Draw("L");
  grPRL->Draw("L");
  if(isXSxBR)
    {
      gr_XSBR68->SetFillColor(kRed);
      gr_XSBR68->SetFillStyle(3013);
      gr_XSBR68->Draw("3same");
      gr_XSBR->SetLineColor(kRed);
      gr_XSBR->Draw("L");
    }
  if(addObsLimit)
    {
      if(points)grObs->Draw("LP");
      else grObs->Draw("L");
    }

 
  if(!isXSxBR)oneLine->Draw("LSAME");
 
  pt->Draw("SAME");
  pt2->Draw("SAME");
	
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
  padl1->Update();
  hrl->GetXaxis()->SetRangeUser(xLow,xHigh);
  padl1->Update();


  padl2->cd();
  TH1F *hrl2 = padl2->DrawFrame(xLow,yLow,xHigh,yHigh);
  hrl2->GetYaxis()->SetTitleOffset(0.6);


  ///DRAW RATIOs


  if(gridOnTop)gPad->RedrawAxis("g");
  box3->Draw();
  sprintf( outfileName,"%s/UpperLimit_%s_wholeMass_PRLCompare.eps",plotDir.c_str(),method.c_str() );
  cl->SaveAs(outfileName);
  //sprintf( outfileName,"%s/UpperLimit_%s_wholeMass_PRLCompare.png",plotDir.c_str(),method.c_str() );
  //cl->SaveAs(outfileName);
  //sprintf( outfileName,"plots/UpperLimit_%s_wholeMass_2D_05.30.C",method.c_str() );	
  //cl->SaveAs(outfileName);	

  // ---------------- Root File --------------- //


  TLegend *legRatio = new TLegend(0.8,0.8,0.99,0.99);
  legRatio->SetFillColor(kWhite);
  legRatio->AddEntry(grRatio,"2D to PRL","L");
  legRatio->AddEntry(grRatioPRL,"1D to PRL","L");



  gStyle->SetOptTitle(0);
  TCanvas *cRatio = new TCanvas("cRatio","cRatio",canvasX,canvasY);
  cRatio->cd();
  TH1F *hrRatio = cRatio->DrawFrame(99.0,0.0,601.0,50);
  grRatio->SetLineWidth(3);
  grRatio->SetLineColor(kRed);
  grRatioPRL->SetLineWidth(3);
  grRatioPRL->SetLineColor(kBlue);
  grRatio->Draw("C");
  grRatioPRL->Draw("C");
  legRatio->Draw("SAME");
  cRatio->Update();
  cRatio->Modified();
  hrRatio->GetYaxis()->SetTitle("Percent Improvement");
  hrRatio->GetXaxis()->SetTitle(xTitle);
  cRatio->Update();
  cRatio->SaveAs(TplotDir+"/PRL_Improvement_06.06.eps");
  cRatio->SaveAs(TplotDir+"/PRL_Improvement_06.06.png");


  TCanvas *cRatio2 = new TCanvas("cRatio2","cRatio",canvasX,canvasY);
  cRatio2->cd();
  TH1F *hrRatio2 = cRatio->DrawFrame(99.0,0.0,601.0,2.0);
  grRatio2->SetLineWidth(3);
  grRatio2->SetLineColor(kRed);
  grRatioPRL2->SetLineWidth(3);
  grRatioPRL2->SetLineColor(kBlue);
  grRatio2->Draw("C");
  grRatioPRL2->Draw("C");
  legRatio->Draw("SAME");
  cRatio2->Update();
  cRatio2->Modified();
  hrRatio2->GetYaxis()->SetTitle("Ratio");
  hrRatio2->GetXaxis()->SetTitle(xTitle);
  cRatio2->Update();
  cRatio2->SaveAs(TplotDir+"/Ratio_06.06.eps");
  cRatio2->SaveAs(TplotDir+"/Ratio_06.06.png");
 


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

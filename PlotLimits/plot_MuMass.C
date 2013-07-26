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
TString inputFile = "results/higgsCombineHZZ4L_ML.root";//"results261012_1DoldZshape/higgsCombineHZZ4L_ASCLS.root";
const bool addObsLimit = true;
const bool isXSxBR = false;
const bool _DEBUG_ = false;
string method = "";
Double_t xLow = 99.9;
Double_t xHigh = 1001.0;
Double_t yLow = -0.4;
Double_t yHigh = 1.8;
TString xTitle = "m_{H} [GeV]";
TString yTitle = "Best fit #sigma/#sigma_{SM}";
const bool logy = false;
const bool logx = true;
const bool grid = true;
const bool gridOnTop = true;
const bool points = false;
const bool isTiny = false;
int canvasX = 700;
int canvasY = 700;
//double sqrts = 8.0;
//Double_t lumi = 1.616;
double sqrts = 7.0;
Double_t lumi = 5.1;
std::string plotDir = "plots";
std::string append = "3D_no2l2tau";
// ----------------------- //


using namespace std;

void plot_MuMass()
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
  
  // ------------------- Get Values -------------------- //

  vector<double> mH, Val_obs, Val_mean, Val_68h, Val_68l, Val_95h, Val_95l, Val_obs95;
  getLimits(inFile,mH,Val_mean,Val_68l,Val_68h,Val_95l,Val_95h,Val_obs,Val_obs95);

  //return;

  vector<double> v_masses, v_means, v_lo68, v_hi68, v_lo95, v_hi95, v_obs;
  vector<double> expExclusion,obsExclusion,expExcl95,obsExcl95;
  printf("%p %p\n",Val_95h,Val_95l);
  for(unsigned int i = 0; i < mH.size(); i++)
    {
      v_masses.push_back( mH[i] );
      v_means.push_back( Val_mean[i] );
      v_lo68.push_back( min( Val_68l[i], Val_68h[i]) );
      v_hi68.push_back( max( Val_68h[i], Val_68l[i]) );
      v_obs.push_back(Val_obs[i]);
      if(Val_mean[i] < 1.0) expExclusion.push_back(mH[i]);
      if(Val_obs[i] < 1.0 && addObsLimit) obsExclusion.push_back(mH[i]);
    }

  // ------------------- Change Values to Arrays -------------------- //

  int nMassEff=0;
  int nExcluded=0;
  const int sizeV = v_masses.size();
  double a_masses[sizeV], a_means[sizeV], a_lo68[sizeV], a_hi68[sizeV], a_obs[sizeV],a_zero[sizeV];
  for(unsigned int m = 0; m < v_masses.size(); m++)
    {
      //if(v_hi68.at(m)>=v_hi95.at(m) || v_lo68.at(m)<=v_lo95.at(m))
      //	{
      //	  cout << "Point at M = " << v_masses.at(m) << " excluded" << endl;
      //	  nExcluded++;
      //	  continue;
      //	}
      
      a_masses[nMassEff] = v_masses[m];
      a_means[nMassEff] = v_means[m];
      a_lo68[nMassEff] = v_lo68[m];
      a_hi68[nMassEff] = v_hi68[m];
      a_obs[nMassEff] = v_obs[m];
      a_zero[nMassEff] = 0;
      //cout << v_masses[m] << "  " << v_means[m] << endl;
      nMassEff++;
      
    }
  cout << "Excluded " << nExcluded << " sick mass points!" << endl;

  // ------------------- Draw  -------------------- //



  TCanvas *c = new TCanvas("c","c",canvasX,canvasY);
  TGraphAsymmErrors* gr = new TGraphAsymmErrors(nMassEff, a_masses, a_means);
 //  TGraphAsymmErrors* grshade_68 = new TGraphAsymmErrors(2*nMassEff);
//   TGraphAsymmErrors* grshade_95 = new TGraphAsymmErrors(2*nMassEff);
  TGraphAsymmErrors* grshade_68 = new TGraphAsymmErrors(nMassEff);
  TGraphAsymmErrors* grshade_95 = new TGraphAsymmErrors(nMassEff);
  gr->SetName("Expected");gr->SetTitle("Expected");
  grshade_68->SetName("68");grshade_68->SetTitle("68");
  grshade_95->SetName("95");grshade_95->SetTitle("95");
  TGraph *grObstemp = new TGraph(nMassEff, a_masses, a_obs);
  grObstemp->Sort();
  double *xtemp=grObstemp->GetX();
  double *ytemp=grObstemp->GetY();
  TGraph *grObs = new TGraph(0);
  for(int igo=0;igo<grObstemp->GetN();igo++){
    if(igo>0)
      if(xtemp[igo]-xtemp[igo-1]<0.5){
	float delta = TMath::Abs(xtemp[igo]-(int)xtemp[igo]);
	bool okd = false;
	if(TMath::Abs(delta-0.5)<0.01)okd=true;
	if(delta<0.01)okd=true;
	if(!okd)continue;
      }
    grObs->Set(grObs->GetN()+1);
    grObs->SetPoint(grObs->GetN()-1,xtemp[igo],ytemp[igo]);
  }
  grObs->SetLineWidth(3);
  grObs->SetLineColor(kBlack);
  grObs->SetMarkerStyle(20);
  grObs->SetMarkerSize(0.7);
 
  for (int a = 0; a < nMassEff; a++)
    {
      grshade_68->SetPoint(a,a_masses[a],a_means[a]);
      grshade_68->SetPointError(a,0,0,a_means[a]-a_lo68[a],-a_means[a]+a_hi68[a]);
    }
  //TGraphAsymmErrors* grshade_68m = removeGlitches(grshade_68n);
  //TGraphAsymmErrors* grshade_95m = removeGlitches(grshade_95n);
  //TGraphAsymmErrors* gr = removeGlitches(grn);
  
  //TGraphAsymmErrors* grshade_68 = slidingWindowAverage(grshade_68m,2);
  //TGraphAsymmErrors* grshade_95 = slidingWindowAverage(grshade_95m,2);
  //TGraph* gr = slidingWindowAverage2(grm,2);
  gr->Sort();
  grshade_68->Sort();
  gr->SetLineStyle(2);
  gr->SetLineWidth(3);
  gr->SetLineColor(kBlue);
  grshade_68->SetFillColor(kGreen);
  grshade_68->SetLineStyle(2);
  grshade_68->SetLineWidth(3);
  grshade_68->SetLineColor(kGreen);
	
  char outfileName[192];

  TF1* oneLine = new TF1("oneLine","1",xLow,xHigh);
  oneLine->SetLineColor(kRed);
  oneLine->SetLineWidth(3);

  // --------------- Low Mass Zoom -------------- //
	
  TLegend * box2 = new TLegend(0.7,0.85,0.93,0.92);
  box2->SetFillColor(0);
  //box2->SetBorderSize(0);
  //if (addObsLimit){ box2->AddEntry(grObs,"Observed","l"); }
  //box2->AddEntry(gr,LimitType+" Fit Expected Asym. CLs","l");
  //box2->AddEntry(gr,"Expected","l");
  box2->AddEntry(grshade_68,"68% CL band","lf");
  //box2->AddEntry(oneLine,"#sigma / #sigma_{SM}","l");


  double ptLow= 0.28, ptHigh = 0.54;

  TPaveText *pt = new TPaveText(0.1577181,0.9562937,0.9580537,0.9947552,"brNDC");
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetFillStyle(0);
  pt->SetTextFont(42);
  pt->SetTextSize(0.03);
  TText *text = pt->AddText(0.01,0.5,"CMS preliminary");
  text = pt->AddText(0.3,0.6,"#sqrt{s} = 7 TeV, L = 5.1 fb^{-1}  #sqrt{s} = 8 TeV, L = 19.8 fb^{-1}");
		

  if(grid) c->SetGrid();
   
  TH1F *hr = c->DrawFrame(110.0,yLow,145.0,yHigh);
	
  gr->Sort();
 
  grshade_68->Draw("e3");
  //gr->Draw("L");

  if(addObsLimit)
    {
      grObs->Sort();
      if(points)grObs->Draw("CP");
      else grObs->Draw("C");
    }

  if(!isXSxBR)oneLine->Draw("LSAME");


  hr->GetXaxis()->SetTitle(xTitle);
  hr->GetYaxis()->SetTitle(yTitle);
  hr->GetXaxis()->SetTitleSize(0.05);
  hr->GetYaxis()->SetTitleSize(0.05);
  hr->GetXaxis()->SetLabelSize(0.04);
  hr->GetYaxis()->SetLabelSize(0.04);
  hr->GetYaxis()->SetTitleOffset(1.2);		
  if(logy)gPad->SetLogy();
  
  
  c->Update();
  if(gridOnTop)gPad->RedrawAxis("g");
  pt->Draw("SAME");
 //  pt2->Draw("SAME");
//   pt3->Draw("SAME");
//   pt4->Draw("SAME");

  box2->Draw();
  sprintf( outfileName,"%s/MuScan%s_lowMass_%s.eps",plotDir.c_str(),method.c_str(), append.c_str() );
  c->SaveAs(outfileName);
  sprintf( outfileName,"%s/MuScan%s_lowMass_%s.png",plotDir.c_str(),method.c_str(), append.c_str() );	
  c->SaveAs(outfileName);
  sprintf( outfileName,"%s/MuScan%s_lowMass_%s.root",plotDir.c_str(),method.c_str(), append.c_str() );
  c->SaveAs(outfileName);

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
      if(mh==113||mh==118||mh==119.5||mh==128.5||mh==132||mh==134||mh==141||mh==142||mh==143)continue;
      if(quant>-1.01&&quant<-0.99)
	{
	  v_obs.push_back(limit);
	  v_mh.push_back(mh);
	  v_obs95.push_back(limit+2*errObs/limit);
	}
      //else if(quant>0.024 && quant<0.026) {v_95l.push_back(limit);printf("95l\n");}
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

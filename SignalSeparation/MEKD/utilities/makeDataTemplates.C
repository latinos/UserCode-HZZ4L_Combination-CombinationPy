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
#include "TH2F.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "Math/DistFunc.h"

int makeDataTemplates(TString dataLoc = "", TString outLoc = ""){


  if (dataLoc == "") dataLoc = "CMSdata";
  if (outLoc == "") outLoc = "dataTemplates";

  gROOT->ProcessLine(".!mkdir -p "+outLoc);

  const int N = 6, M = 1;
  
  TString outName[N] = {"Hist_Scatter_Data_2e2m_7TeV.dat.root","Hist_Scatter_Data_4m_7TeV.dat.root","Hist_Scatter_Data_4e_7TeV.dat.root",
			"Hist_Scatter_Data_2e2m_8TeV.dat.root","Hist_Scatter_Data_4m_8TeV.dat.root","Hist_Scatter_Data_4e_8TeV.dat.root"};
  TString fileName[N] = {"hzz2e2mu_5.051.root","hzz4mu_5.051.root","hzz4e_5.051.root","hzz2e2mu_19.63.root","hzz4mu_19.63.root","hzz4e_19.63.root"};
  TString treeName = "data_obs";
  TString histName[M] = {"hist"};
  TString yName[M] = {"CMS_zz4l_pseudoKD"};
  TString xName[M] = {"CMS_zz4l_smd"};

  TFile *f = new TFile("templates2D/Dsignal_ALT_2e2mu.root","READ");
  TH2F *hist = (TH2F*)f->Get("h_superDpsD");
  TH2F *clonedHist = hist->Clone();
  clearHist(clonedHist);

  for(int i = 0; i < N; i++)
    {
      TFile *outFile = new TFile(outLoc+"/"+outName[i],"RECREATE");
      
      for( int j = 0; j < M; j++)
	{
	  
	  makeTemplate(dataLoc+"/"+fileName[i],treeName,xName[j],yName[j],clonedHist);
	  outFile->cd();
	  clonedHist->SetNameTitle(histName[j],histName[j]);
	  clonedHist->Write();
	  clearHist(clonedHist);
	}

      outFile->Close();
    }

  return 0;

}






void makeTemplate(TString fileName, TString treeName, TString xName, TString yName, TH2F* &hist){

  double massLow = 106;
  double massHigh = 141;

  double xKD, yKD, mass4l;

  TFile *f = new TFile(fileName,"READ");
  TTree *tree = (TTree*)f->Get(treeName);

  for (int i = 0; i < tree->GetEntries(); i++)
    {

      tree->GetEntry(i);

      tree->SetBranchAddress("CMS_zz4l_mass",&mass4l);
      tree->SetBranchAddress(xName,&xKD);
      tree->SetBranchAddress(yName,&yKD);


      if (mass4l < massLow || mass4l > massHigh) continue;

      hist->Fill(xKD,yKD);
	
    }

  f->Close();
  

  
}


void clearHist(TH2F* &hist){

  for(int i = 0; i < hist->GetNbinsX(); i++)
    {
      for(int j = 0; j < hist->GetNbinsY(); j++)
	{
	  hist->SetBinContent(i,j,0);
	}
    }
  
}

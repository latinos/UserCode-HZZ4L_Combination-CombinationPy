// Directly produce MC templates from samples: to be used with new MC

// C++ includes
#include <iostream>
#include <string>
#include <stdio.h>
#include <sstream>
// #include <algorithm>

// ROOT includes
#include <TROOT.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TProfile.h>
#include <TChain.h>

const TString destDir = "./selRootFiles/";
TString whichVar = "PTRestricted";
static const int nsamp = 8;
const TString dataFileNames[nsamp] = {"gg","vbf","wh","zh","tth","zz","zx","ggzz"};
TString systSources[nsamp][5];

using namespace std;

void adjustHistogram(TH2F* hist) {

  // fill with larger bins if not enough events
  hist->Smooth(1,"k5b"); 

  // fill empty bins  
  double floor = ((hist->Integral())/(hist->GetNbinsX()*hist->GetNbinsY()))*0.001;
  for(int i = 1; i <= hist->GetNbinsX(); i++) {
    for(int j = 1; j <= hist->GetNbinsY(); j++) {
      double orig = hist->GetBinContent(i,j);
      if (orig < 0.) orig = 0.;
      hist->SetBinContent(i,j,(orig+floor));
    }
  } 

  // normalize slices
  double norm;
  TH1F* tempProj;
 
  for(int i=1; i<=hist->GetNbinsX(); i++){
    
    tempProj = (TH1F*)hist->ProjectionY("tempProj",i,i);
    norm = tempProj->Integral();
 
    if (norm>0) { // Avoid introducing NaNs in the histogram
      for(int j=1; j<=hist->GetNbinsY(); j++){
	hist->SetBinContent(i,j,hist->GetBinContent(i,j)/norm);
      }
    }

  }

  return;
}

void makeErrorTemplatesPtSyst(TString channel="2e2mu") { 
    
  systSources[0][0] = "Resummation";
  systSources[0][1] = "TopMass";
  systSources[0][2] = "Mela";
  
  systSources[1][0] = "PDF-VBF";
  systSources[1][1] = "scale-VBF";
  systSources[1][2] = "Mela";
  
  systSources[5][0] = "SingleZ";
  systSources[5][1] = "PDF-ZZ";
  systSources[5][2] = "scale-ZZ";
  systSources[5][3] = "Mela"; 
 
  systSources[2][0] = "NLOLO_WH";
  systSources[3][0] = "NLOLO_ZH";

  TString aLongString;
  char fileName[200];
  char fileName2[200];

  //  for (Int_t lhc=8; lhc<9; lhc++) {
  for (Int_t lhc=7; lhc<9; lhc++) {
    for (Int_t k=0; k<nsamp; k++) {

      TString lhcs = "7";
      if (lhc==8) lhcs="8";
 
      aLongString = destDir + "/" + whichVar + "_" + dataFileNames[k] + "_TEMPL_" + channel + "_" + lhcs + "TeV.root";  
      TFile* ftemp = new TFile(aLongString,"UPDATE");
      
      TString whichVar2 = whichVar;
      whichVar2.ToLower();
      sprintf(fileName,"%sH_Default",whichVar2.Data());
      cout << aLongString << " " << fileName << endl;
      TH2F* baseHist = (TH2F*)ftemp->Get(fileName);
      cout << "Histo found OK " << ftemp << " " << baseHist << endl;

      TH2F* upHist = (TH2F*)baseHist->Clone();
      TH2F* downHist = (TH2F*)baseHist->Clone();

      //MC stats
      for (int i=1; i<=baseHist->GetNbinsX(); i++) { 
	for (int j=1; j<=baseHist->GetNbinsY(); j++) {
	  upHist->SetBinContent(i,j,(baseHist->GetBinError(i,j))*(baseHist->GetBinError(i,j)));
	  downHist->SetBinContent(i,j,(baseHist->GetBinError(i,j))*(baseHist->GetBinError(i,j)));
	}
      }
      cout << "Errors filled OK" << endl;
      cout << dataFileNames[k] << " entries " << baseHist->GetBinContent(100,10) << " "  << upHist->GetBinContent(100,10) << " " << downHist->GetBinContent(100,10) << endl ;

      //Other systs
      for (int ss = 0; ss < 5; ss++) {
	if (systSources[k][ss] != "") {

	  TH2F* thisHistUp;
	  TH2F* thisHistDown;

	  if (systSources[k][ss] == "Resummation") {
	    sprintf(fileName,"%sH_ResummationUp",whichVar2.Data());
            sprintf(fileName2,"%sH_ResummationDown",whichVar2.Data());
	  } else if (systSources[k][ss] == "Mela") {
            sprintf(fileName2,"%sH_Mela00-03",whichVar2.Data());
            sprintf(fileName,"%sH_Mela06-10",whichVar2.Data());
          } else {
	    sprintf(fileName,"%sH_%s",whichVar2.Data(),systSources[k][ss].Data());
	    sprintf(fileName2,"%sH_%s",whichVar2.Data(),systSources[k][ss].Data());  
	  }
	  thisHistUp = (TH2F*)((TH2F*)ftemp->Get(fileName))->Clone();
          cout << "Histo up found OK" << endl;
          thisHistDown = (TH2F*)((TH2F*)ftemp->Get(fileName2))->Clone();
          cout << "Histo down found OK" << endl;

	  for (int i=1; i<=baseHist->GetNbinsX(); i++) { 
	    for (int j=1; j<=baseHist->GetNbinsY(); j++) {
	      float upval = upHist->GetBinContent(i,j) + pow(thisHistUp->GetBinContent(i,j)-baseHist->GetBinContent(i,j),2);
	      float downval = upHist->GetBinContent(i,j) + pow(thisHistDown->GetBinContent(i,j)-baseHist->GetBinContent(i,j),2);
	      upHist->SetBinContent(i,j,upval);
	      downHist->SetBinContent(i,j,downval);
	    }
	  }

          if (k == 6) cout << "zx entries after syst " << systSources[k][ss] << " " << upHist->GetBinContent(100,10) << " " << downHist->GetBinContent(100,10) << endl;
	}
      }

      // Subtract or add
      for (int i=1; i<=baseHist->GetNbinsX(); i++) { 
	for (int j=1; j<=baseHist->GetNbinsY(); j++) {
	  float upval = baseHist->GetBinContent(i,j) + sqrt(upHist->GetBinContent(i,j));
          float downval = baseHist->GetBinContent(i,j) - sqrt(downHist->GetBinContent(i,j));
	  upHist->SetBinContent(i,j,upval);
	  downHist->SetBinContent(i,j,downval);
	}
      }

      if (k == 6 && lhc == 8) cout << "zx entries after sqrt " << baseHist->GetBinContent(100,10) << " " << upHist->GetBinContent(100,10) << " " << downHist->GetBinContent(100,10) << endl;

      ftemp->cd();
      upHist->SetName(whichVar2 + "H_TotalUp");
      adjustHistogram(upHist);
      upHist->Write();
      downHist->SetName(whichVar2 + "H_TotalDown");
      adjustHistogram(downHist); 
      downHist->Write();
    }      
  }
  return;
}


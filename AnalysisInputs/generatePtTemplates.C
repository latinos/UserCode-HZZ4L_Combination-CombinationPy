/* 
 * Create 2D (mass, LD) templates. Script imported from: http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/JHU/MELA/scripts/generateTemplates.C?revision=1.1.2.1&view=markup&pathrev=post_unblinding
  * usage: 
 * -set input paths variables in Config.h
 * -run with:
 * > root -q -b 
 * > .L RooModifTsallis.cc+ 
 * > .x generatePtTemplates.C+
 * 2D templates are written to "destDir"
 *
 */

#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TString.h"
#include <sstream>
#include <vector>

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TROOT.h>
//#include "ZZMatrixElement/MELA/interface/Mela.h"
#endif

#include "RooRealVar.h"
#include "RooModifTsallis.h"
//----------> SET INPUT VARIABLES in Config.h
#include "Config.h"
//<----------

//---
int useSqrts=0;              //0=use 7+8TeV; 1=use 7TeV only, 2 use 8TeV only
TString melaName = "ZZLD"; // name of MELA branch to be used. Possibilities are ZZLD,ZZLD_analBkg,ZZLD_postICHEP,ZZLD_PtY,pseudoMelaLD, spinTwoMinimalMelaLD 

bool extendToHighMass = true; // Include signal samples above 600 GeV

float highMzz=(extendToHighMass?1600:800);
float mBinSize=2.;

const TString destDir = "../CreateDatacards/templates2D/";
static const int nsamp = 8;
const TString dataFileNames[nsamp] = {"gg","vbf","wh","zh","tth","zz","zx","ggzz"};
TString systSources[nsamp][5];

//=======================================================================

TH2F* fillTemplate(TString channel="4mu", int sampleIndex=0,int LHCsqrts= 7,bool isLowMass=true, TString systName = "Default", bool down = false){  
 
  TH2F* bkgHist;
  if(!isLowMass)
    bkgHist = new TH2F("bkgHisto","bkgHisto",int((highMzz-100.)/mBinSize+0.5),100,highMzz,500,0,4);
  else
    bkgHist = new TH2F("bkgHisto","bkgHisto",int(35/mBinSize+0.5),106,141,500,0,4);

  bkgHist->Sumw2();

  if (sampleIndex < 0 || sampleIndex > nsamp) {
    cout << "Samples are numbered from 0 to " << nsamp-1 << endl;
    return bkgHist;
  }

  // fill histogram
  RooRealVar* ptoverm = new RooRealVar("ptoverm","p_{T}/M^{4l}",0.,4.,"GeV/c");

  RooRealVar mup("mup","emme", 1.,-1000000.,1000000.);
  RooRealVar nup("nup","enne", 0.93, -1000000.,1000000.);
  RooRealVar n2up("n2up","enne2", 0.75, -1000000.,1000000.);
  RooRealVar bbup("bbup","bibi",0.02, -1000000.,1000000.);
  RooRealVar Tup("Tup","tti",0.2,-1000000.,1000000.);
  RooRealVar bb2up("bb2up","bibi2",0.02, -1000000.,1000000.);
  RooRealVar fexpup("fexpup","f_exp",0.02, -1000000.,1000000.);
 
  RooModifTsallis* rtup = new RooModifTsallis("rtup","rtup",*ptoverm,
					    mup,nup,n2up,bbup,bb2up,Tup,fexpup);
  RooRealVar m("m","emme", 1.,-1000000.,1000000.);
  RooRealVar n("n","enne", 0.93, -1000000.,1000000.);
  RooRealVar n2("n2","enne2", 0.75, -1000000.,1000000.);
  RooRealVar bb("bb","bibi",0.02, -1000000.,1000000.);
  RooRealVar T("T","tti",0.2,-1000000.,1000000.);
  RooRealVar bb2("bb2","bibi2",0.02, -1000000.,1000000.);
  RooRealVar fexp("fexp","f_exp",0.02, -1000000.,1000000.);
 
  RooModifTsallis* rt = new RooModifTsallis("rt","rt",*ptoverm,
					    m,n,n2,bb,bb2,T,fexp);

  // default
  char fileName[200];
  int nXbins=bkgHist->GetNbinsX();
  int nYbins=bkgHist->GetNbinsY();
    
  sprintf(fileName,"tsallisPars/paramsPTOverMCJLST_%s_%dTeV_Default.txt",dataFileNames[sampleIndex].Data(),LHCsqrts);
  if (sampleIndex == 4) sprintf(fileName,"tsallisPars/paramsPTOverMCJLST_%s_%dTeV_Default.txt",dataFileNames[6].Data(),LHCsqrts);
  ifstream theFile(fileName);
  (RooArgSet(mup,nup,n2up,bbup,bb2up,fexpup,Tup)).readFromStream(theFile,false);   
  if (systName == "Default") {
    
    for (Int_t i=1; i<=nYbins; i++) {      
      ptoverm->setVal(bkgHist->GetYaxis()->GetBinCenter(i));
      // if (sampleIndex == 0 && LHCsqrts == 8) cout << i << " " << ptoverm->getVal() << " " << rt->getVal() << endl; 
      for (Int_t j=1; j<=nXbins; j++) {
	bkgHist->SetBinContent(j,i,rtup->getVal());
      }
    } 

  } else if (systName == "Mela") {
    
    if (down) {
      sprintf(fileName,"tsallisPars/paramsPTOverMCJLST_%s_8TeV_Mela00-03.txt",dataFileNames[sampleIndex].Data());
    } else {
      sprintf(fileName,"tsallisPars/paramsPTOverMCJLST_%s_8TeV_Mela06-10.txt",dataFileNames[sampleIndex].Data());
    }

    ifstream theFile2(fileName);
    (RooArgSet(m,n,n2,bb,bb2,fexp,T)).readFromStream(theFile2,false);  

    for (Int_t i=1; i<=nYbins; i++) {      
      ptoverm->setVal(bkgHist->GetYaxis()->GetBinCenter(i));
      // if (sampleIndex == 0 && LHCsqrts == 8) cout << i << " " << ptoverm->getVal() << " " << rt->getVal() << endl; 
      for (Int_t j=1; j<=nXbins; j++) {
	bkgHist->SetBinContent(j,i,rt->getVal());
      }
    }
 
  } else {
    
    sprintf(fileName,"tsallisPars/paramsPTOverMCJLST_%s_%dTeV_%s.txt",dataFileNames[sampleIndex].Data(),LHCsqrts,systName.Data());
    ifstream theFile2(fileName);
    (RooArgSet(m,n,n2,bb,bb2,fexp,T)).readFromStream(theFile2,false);     

    if (down) {
      m.setVal(fabs(2*mup.getVal() - m.getVal()));
      n.setVal(fabs(2*nup.getVal() - n.getVal()));
      n2.setVal(fabs(2*n2up.getVal() - n2.getVal()));
      bb.setVal(fabs(2*bbup.getVal() - bb.getVal()));
      T.setVal(fabs(2*Tup.getVal() - T.getVal()));
      bb2.setVal(fabs(2*bb2up.getVal() - bb2.getVal()));
      fexp.setVal(fabs(2*fexpup.getVal() - fexp.getVal()));
    }
    for (Int_t i=1; i<=nYbins; i++) {      
      ptoverm->setVal(bkgHist->GetYaxis()->GetBinCenter(i));
      // if (sampleIndex == 0 && LHCsqrts == 8) cout << i << " " << ptoverm->getVal() << " " << rt->getVal() << endl; 
      for (Int_t j=1; j<=nXbins; j++) {
	bkgHist->SetBinContent(j,i,rt->getVal());
      }
    } 
  }
  // normalize slices

  double norm;
  TH1F* tempProj;
  
  for(int i=1; i<=nXbins; i++){
    
    tempProj = (TH1F*)bkgHist->ProjectionY("tempProj",i,i);
    norm = tempProj->Integral();
    // if (sampleIndex == 0 && LHCsqrts == 8) cout << i << " " << norm << endl;

    if (norm>0) { // Avoid introducing NaNs in the histogram
      for(int j=1; j<=nYbins; j++){
	bkgHist->SetBinContent(i,j,bkgHist->GetBinContent(i,j)/norm);
      }
    }

  }
  
  return bkgHist;
  
}

//=======================================================================

void makeTemplate(TString channel="4mu",bool isLowMass=true){

  char fileName[200];

  for (Int_t lhc=7; lhc<9; lhc++) {
    for (Int_t k=0; k<nsamp; k++) {
     
      TString lhcs = "7";
      if (lhc==8) lhcs="8";

      TFile* ftemp = new TFile(destDir + "PtoverM_mZZ_" + dataFileNames[k] + "_" + channel + "_"+ lhcs + "TeV.root","RECREATE");

      cout << "*** Now filling: " << dataFileNames[k] << ", Default template" << endl;
      TH2F* h_mzzD = fillTemplate(channel,k,lhc,isLowMass);
     
      ftemp->cd();
  
      h_mzzD->Write("h_Ptmzz_mzz");

      for (int ss = 0; ss < 5; ss++) {
	if (systSources[k][ss] != "") {

	  cout << "*** Now filling: " << dataFileNames[k] << ", " << systSources[k][ss] << " templates" << endl;

	  TH2F* h_mzzDup = fillTemplate(channel,k,lhc,isLowMass,systSources[k][ss],false);
	  sprintf(fileName,"h_Ptmzz_mzz_%s_up",systSources[k][ss].Data());
          h_mzzDup->Write(fileName);

	  TH2F* h_mzzDdown = fillTemplate(channel,k,lhc,isLowMass,systSources[k][ss],true);
	  sprintf(fileName,"h_Ptmzz_mzz_%s_down",systSources[k][ss].Data());
          h_mzzDdown->Write(fileName);
	}
      }

      ftemp->Close();
    }
  }

}

//=======================================================================

void generatePtTemplates(){

  bool isLowMass = false;

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

  makeTemplate("4mu",isLowMass);
  makeTemplate("4e",isLowMass);
  makeTemplate("2e2mu",isLowMass);

}




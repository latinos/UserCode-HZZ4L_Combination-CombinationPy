//Fit the PT curves for sig and BKG
//
//By: Roberto Covarelli (U of Rochester)

#include <math.h>
#include <string>
#include <vector>
#include <iostream>
#include <stdio.h>

#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TPad.h"
#include "TLatex.h"

#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooBinning.h"
#include "RooHistPdf.h"
#include "RooFitResult.h"
#include "RooRandom.h"

// #include "PDFs/RooTsallis.h"
// #include "PDFs/RooTsallis2.h"
#include "PDFs/RooModifTsallis.h" 

using namespace RooFit;

RooDataHist* withSmartBinning(TH1F* source, RooRealVar* var, float min, float max, float minWeight, float sizeSig, int type) {
  
  // DEFINE BINNING
  // type:
  // 0 = background
  // 1 = HRes gg signal
  // 2 = POWHEG gg signal
  // 2 = VBF signal
  float lowlimit;
  std::vector<float> sizes;
  std::vector<float> highlimits;  
  if (type == 0) {
    lowlimit = min;
    sizes.push_back(1.);
    highlimits.push_back(max);
  } else if (type == 1 || type == 2) {
    lowlimit = min;
    sizes.push_back(sizeSig);        
    // sizes.push_back(8.);        
    // sizes.push_back(4.);
    // sizes.push_back(8.);
    // highlimits.push_back(150.); 
    // highlimits.push_back(200.);
    // highlimits.push_back(320.);
    highlimits.push_back(max);
  } else {
    lowlimit = min;
    sizes.push_back(sizeSig);
    highlimits.push_back(max);
  }

  // ROOREALVAR BINNING
  var->setMin(lowlimit);
  var->setMax(highlimits.back());
  RooBinning rb(lowlimit,highlimits.back());
  float leftlimit = lowlimit;
  for (unsigned int ii = 0; ii < highlimits.size(); ii++) {
    float rightlimit = highlimits.at(ii);
    int howManyBins = (rightlimit-leftlimit)/sizes.at(ii);
    rb.addUniform(howManyBins,leftlimit,rightlimit);
    leftlimit = rightlimit;
  }
  var->setBinning(rb);

  // FILL ROODATAHIST
  RooDataHist* result = new RooDataHist("result","A dataset",RooArgList(*var));
  float nWeight = source->GetEntries();
  if (type == 2) nWeight = 0.00215;
  float binWidth = source->GetBinWidth(1);

  float thisWeight = 0.;
  float thisWeightsum = 0.;
  int j = 0;
  int whichInterval = 0;
  for (Int_t i=1; i<=source->GetNbinsX(); i++) {

    float thispt = source->GetXaxis()->GetBinCenter(i);
    if (thispt > lowlimit && thispt < highlimits.back()) {
      bool lastbin = false;
      if (source->GetXaxis()->GetBinCenter(i+1) >= highlimits.back()) lastbin = true;
      if (thispt > highlimits.at(whichInterval)) whichInterval++;
      var->setVal(thispt);
      thisWeight += source->GetBinContent(i)*nWeight;
      thisWeightsum += source->GetBinError(i)*nWeight*source->GetBinError(i)*nWeight;
      j++;
      cout << "Entries in TH1 bin " << i << " = " << source->GetBinContent(i)*nWeight << " +/- " << source->GetBinError(i)*nWeight << endl;   
      if (binWidth*j >= sizes.at(whichInterval) || lastbin) {
        // cout << "Entro" << endl; 
	if (thisWeight <= 0) thisWeight = minWeight; 
        if (fabs(thisWeight/sizes.at(whichInterval) - 0.0829038) < 0.00001) thisWeight /= 10.;
	if (sqrt(thisWeightsum) > thisWeight) thisWeightsum = thisWeight*thisWeight*0.9;
	cout << " " << var->getVal() << " " << thisWeight << " " << thisWeightsum << " " << sizes.at(whichInterval) << endl; 
	result->add(RooArgSet(*var),
		    thisWeight/sizes.at(whichInterval),
		    thisWeightsum/(sizes.at(whichInterval)*sizes.at(whichInterval)));
	thisWeight = 0.;
	thisWeightsum = 0.;
	j = 0;
      }
    }
  } 

  return result;
}

void fitPtOverMCJLST(int mass = 125, int LHCsqrts = 7, int whichtype = 1, 
		     bool correctErrors = false, /* string changeParName = "", */
		     bool showErrorPDFs = false, string systString = "Default")

// whichtype
// 0 - gg Signal
// 1 - VBF Signal
// 2 - ZZ
// 3 - ZX
// 4 - ggZZ
// 5 - WH
// 6 - ZH
// 7 - ttH

{

  string changeParName = "";
  if (systString == "Default") changeParName = "up";

  string nameSample[8] = {"gg","vbf","zz","zx","ggzz","wh","zh","tth"};
  float maxType[8] = {2.4,3.2,1.6,1.6,1.6,3.2,3.2,3.2};   
  float rebinType[8] = {1,2,1,1,4,10,10,40};
  
  for (int t = 0; t < 8; t++) {
    if (mass > 150) maxType[t] = int(117.90*maxType[t]/sqrt(mass-10.91))/10.;
  }

  char fileToOpen[200];
  sprintf(fileToOpen,"selRootFiles/PToverM_%s%d_SEL_%dTeV.root",nameSample[whichtype].c_str(),mass,LHCsqrts);
  // if (whichtype == 3) sprintf(fileToOpen,"PTOVERM_%s_SEL_allTeV.root",nameSample[whichtype].c_str());

  RooRealVar* ptoverm = new RooRealVar("ptoverm","p_{T}/M^{4l}",0.,maxType[whichtype],"GeV/c");
 
  TFile input(fileToOpen);
  // if (systString == "Mass" || systString == "Mela") {
  //  sprintf(fileToOpen,"ptovermH_%sUp",systString.c_str());
  // } else {
  sprintf(fileToOpen,"ptovermH_%s",systString.c_str());
  //}
  TH1F* ptovermH = (TH1F*)input.Get(fileToOpen);
  
  if (rebinType[whichtype] > 1) ptovermH->Rebin(rebinType[whichtype]);
  if (maxType[whichtype] < ptovermH->GetBinLowEdge(ptovermH->GetNbinsX() + 1) - ptovermH->GetBinWidth(1)) {
    int theBin = ptovermH->FindBin(maxType[whichtype]);
    ptovermH->GetXaxis()->SetRange(1,theBin-1);
  }

  gROOT->ProcessLine(".L mytdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");
  
  // cout << endl << "Signal " << endl;   
  ptoverm->setBins(ptovermH->GetNbinsX());

  RooDataHist* rdh = new RooDataHist("rdh","Some dataset",RooArgList(*ptoverm),Import(*ptovermH,kFALSE));
 
  // fit definitions
  // RooWorkspace *ws = new RooWorkspace("ws");

  RooRealVar m("m","emme", 1.,0.01, 40.);
  RooRealVar n("n","enne", 0.93, 0.05, 15.);
  RooRealVar n2("n2","enne2", 0.75, 0.5, 15.);
  RooRealVar bb("bb","bibi",0.02, 0.000005, 20.0);
  RooRealVar T("T","tti",0.2,0.00000005,1.);
  RooRealVar bb2("bb2","bibi2",0.02, 0.0005, 10.0);
  RooRealVar fexp("fexp","f_exp",0.02, 0.0, 1.0);
  if (whichtype == 0) {
    if (LHCsqrts == 8) {
      m.setVal(3.319);   // m.setConstant(kTRUE);
      n.setVal(0.7606);    if (systString != "Default" || mass != 125) n.setConstant(kTRUE); 
      n2.setVal(0.8061);   n2.setConstant(kTRUE);
      bb.setVal(3.728);   // bb.setConstant(kTRUE);
      T.setVal(0.00333);   // T.setConstant(kTRUE);
      bb2.setVal(1.7172);    // bb2.setConstant(kTRUE);
      fexp.setVal(0.002144);   if (systString != "Default" || mass != 125) fexp.setConstant(kTRUE);
    } else {
      m.setVal(0.061);   // m.setConstant(kTRUE);
      n.setVal(1.6141);   if (systString == "Resummation" || systString == "TopMass")  n.setConstant(kTRUE);
      n2.setVal(1.3294);   n2.setConstant(kTRUE);
      bb.setVal(4.2761);   // bb.setConstant(kTRUE);
      T.setVal(0.0361);   // T.setConstant(kTRUE);
      bb2.setVal(1.6643);   bb2.setConstant(kTRUE);
      fexp.setVal(0.0004);   // fexp.setConstant(kTRUE);
    }
  } else if (whichtype == 1) {
    m.setVal(1.006);   // m.setConstant(kTRUE);
    n.setVal(10.939);   n.setConstant(kTRUE);
    n2.setVal(1.1448);   n2.setConstant(kTRUE);
    bb.setVal(0.02048);   bb.setConstant(kTRUE);
    T.setVal(0.16115);   if (systString.find("Mela") != string::npos) T.setConstant(kTRUE); // T.setConstant(kTRUE);
    bb2.setVal(1.0024);   bb2.setConstant(kTRUE);
    fexp.setVal(0.005);   fexp.setConstant(kTRUE);
    if (mass > 300) {
      fexp.setVal(0.0);   fexp.setConstant(kFALSE);
    }
    if (mass > 500) {
      bb2.setVal(5.0);  //  bb2.setConstant(kFALSE);
    }
    if (mass > 500) {
      bb.setVal(15.0);  //  bb.setConstant(kFALSE);
    }
  } else if (whichtype == 2) {
    if (LHCsqrts == 8) {
      m.setVal(1.0476);   // m.setConstant(kTRUE);    
      bb.setVal(3.3088);  // if (mass != 140) bb.setConstant(kTRUE);
      n2.setVal(0.7146);   n2.setConstant(kTRUE);
      n.setVal(0.9518);      n.setConstant(kTRUE);
      bb2.setVal(100000.);  bb2.setConstant(kTRUE);
      T.setVal(0.0021889);    if (systString.find("Mela") != string::npos || mass != 140) T.setConstant(kTRUE);
      fexp.setVal(0.0);    fexp.setConstant(kTRUE);
    } else {
      m.setVal(1.028);   // m.setConstant(kTRUE);    
      bb.setVal(2.91); // bb.setConstant(kTRUE);
      n2.setVal(0.7146);  n2.setConstant(kTRUE);
      n.setVal(0.9518);     n.setConstant(kTRUE);
      bb2.setVal(100000.);  bb2.setConstant(kTRUE);
      T.setVal(0.002248);   if (systString.find("Mela") != string::npos) T.setConstant(kTRUE);
      fexp.setVal(0.0);    fexp.setConstant(kTRUE);
    }
  } else if (whichtype == 3) {
    m.setVal(1.411);   // m.setConstant(kTRUE);
    n.setVal(3.4523);     n.setConstant(kTRUE);
    n2.setVal(0.6910);    n2.setConstant(kTRUE);
    bb.setVal(0.00039);  // bb.setConstant(kTRUE);
    T.setVal(0.118);   // T.setConstant(kTRUE);
    bb2.setVal(0.0224);   bb2.setConstant(kTRUE);
    fexp.setVal(0.0);   fexp.setConstant(kTRUE);
  } else if (whichtype == 4) {
    m.setVal(1.411);   // m.setConstant(kTRUE);
    n.setVal(5.756);    // n.setConstant(kTRUE);
    n2.setVal(0.8738);   // n2.setConstant(kTRUE);
    bb.setVal(0.00039);  // bb.setConstant(kTRUE);
    T.setVal(0.118);   // T.setConstant(kTRUE);
    bb2.setVal(0.0224);   bb2.setConstant(kTRUE);
    fexp.setVal(0.0);   fexp.setConstant(kTRUE);
  } else if (whichtype == 5 && LHCsqrts == 8) {
    m.setVal(1.006);   // m.setConstant(kTRUE);
    n.setVal(10.939);    n.setConstant(kTRUE);
    n2.setVal(1.1448);   n2.setConstant(kTRUE);
    bb.setVal(3.897);   bb.setConstant(kTRUE);
    T.setVal(0.1009);  // T.setConstant(kTRUE);
    bb2.setVal(1.0224);   bb2.setConstant(kTRUE);
    fexp.setVal(0.01);   fexp.setConstant(kTRUE);
  } else {
    // cout << "Entro qui" << endl;
    m.setVal(1.006);   // m.setConstant(kTRUE);
    n.setVal(10.939);    n.setConstant(kTRUE);
    n2.setVal(1.1448);   n2.setConstant(kTRUE);
    bb.setVal(0.0129);  bb.setConstant(kTRUE);
    T.setVal(0.1009);    // T.setConstant(kTRUE);
    bb2.setVal(1.0224);   bb2.setConstant(kTRUE);
    fexp.setVal(0.01);   fexp.setConstant(kTRUE);
  }
 
  
  RooModifTsallis* rt3 = new RooModifTsallis("rt3","rt3",*ptoverm,m,n,n2,bb,bb2,T,fexp);
  // ws->import(*rt3);

  // fit
  RooFitResult* fit = rt3->fitTo(*rdh,Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(1));  

  float mVal = m.getVal();   
  float nVal = n.getVal();
  float n2Val = n2.getVal();
  float bbVal = bb.getVal();
  float bb2Val = bb2.getVal();
  float fexpVal = fexp.getVal();
  float TVal = T.getVal();

  if (correctErrors) {
    // Tsallis errors not reliable, use toy MC

    TH1F* mHist = new TH1F("mHist","m",21,-0.5*mVal,0.5*mVal);
    TH1F* nHist = new TH1F("nHist","n",21,-0.2*nVal,0.2*nVal);
    TH1F* n2Hist = new TH1F("n2Hist","n2",21,-0.2*n2Val,0.2*n2Val);
    TH1F* bbHist = new TH1F("bbHist","bb",21,-0.2*bbVal,0.2*bbVal);
    TH1F* bb2Hist = new TH1F("bb2Hist","bb2",21,-0.2*bb2Val,0.2*bb2Val);
    TH1F* fexpHist = new TH1F("fexpHist","fexp",21,-0.2*fexpVal-0.000001,0.2*fexpVal+0.000001);
    TH1F* THist = new TH1F("THist","T",21,-0.5*TVal,0.5*TVal);
    mHist->GetXaxis()->SetTitle("m-m_{gen}");
    nHist->GetXaxis()->SetTitle("n-n_{gen}");
    n2Hist->GetXaxis()->SetTitle("n2-n2_{gen}");
    bbHist->GetXaxis()->SetTitle("bb-bb_{gen}");
    bb2Hist->GetXaxis()->SetTitle("bb2-bb2_{gen}");
    THist->GetXaxis()->SetTitle("T-T_{gen}");
    fexpHist->GetXaxis()->SetTitle("fexp-fexp_{gen}");

    for (unsigned int iToy = 0; iToy < 200; iToy++) {

      cout << endl << "####" << endl;
      cout << "Generating toy experiment n. " << iToy+1 << endl;

      m.setVal(mVal);
      n.setVal(nVal);
      n2.setVal(n2Val);
      bb.setVal(bbVal);
      bb2.setVal(bb2Val);
      fexp.setVal(fexpVal);
      T.setVal(TVal);

      TDatime *now = new TDatime();
      Int_t seed = now->GetDate() + now->GetTime();
      cout << "RooFit Generation Seed = " << seed+iToy << endl;
      RooRandom::randomGenerator()->SetSeed(seed+iToy);
      cout << "####" << endl << endl;

      RooDataSet *dataToy = rt3->generate(RooArgSet(*ptoverm),ptovermH->GetEntries());
      RooDataHist *dataToyH = new RooDataHist("dataToyH","toy",RooArgSet(*ptoverm),*dataToy);
      
      rt3->fitTo(*dataToyH,Minos(0),SumW2Error(kTRUE),NumCPU(1));  
  
      if (fit->floatParsFinal().find("m")) mHist->Fill(m.getVal()-mVal);
      if (fit->floatParsFinal().find("n")) nHist->Fill(n.getVal()-nVal);
      if (fit->floatParsFinal().find("n2")) n2Hist->Fill(n2.getVal()-n2Val);
      if (fit->floatParsFinal().find("bb")) bbHist->Fill(bb.getVal()-bbVal);
      if (fit->floatParsFinal().find("bb2")) bb2Hist->Fill(bb2.getVal()-bb2Val);
      if (fit->floatParsFinal().find("fexp")) fexpHist->Fill(fexp.getVal()-fexpVal);
      if (fit->floatParsFinal().find("T")) THist->Fill(T.getVal()-TVal);
    }

    TCanvas cant("cant","Test canvas",5.,5.,900.,500.);
    cant.Divide(4,2);
    cant.cd(1);   mHist->Draw();
    cant.cd(2);   nHist->Draw();
    cant.cd(3);   n2Hist->Draw();
    cant.cd(4);   bbHist->Draw();
    cant.cd(5);   bb2Hist->Draw();
    cant.cd(6);   fexpHist->Draw();
    cant.cd(7);   THist->Draw();
    // cant.SaveAs("figs/testToys.pdf");
    cant.SaveAs("newfigs/testToys.pdf");

    if (fit->floatParsFinal().find("m")) m.setError(mHist->GetRMS());
    if (fit->floatParsFinal().find("n")) n.setError(nHist->GetRMS());
    if (fit->floatParsFinal().find("n2")) n2.setError(n2Hist->GetRMS());
    if (fit->floatParsFinal().find("bb")) bb.setError(bbHist->GetRMS());
    if (fit->floatParsFinal().find("bb2")) bb2.setError(bb2Hist->GetRMS());
    if (fit->floatParsFinal().find("fexp")) fexp.setError(fexpHist->GetRMS());
    if (fit->floatParsFinal().find("T")) T.setError(THist->GetRMS());
  }

  m.setVal(mVal);
  n.setVal(nVal);
  n2.setVal(n2Val);
  bb.setVal(bbVal);
  bb2.setVal(bb2Val);
  fexp.setVal(fexpVal);
  T.setVal(TVal);

  char fileToSave[200];
  // if (changeParName != "") 
  //  sprintf(fileToSave,"text/paramsPTOverMCJLST_%s_%dTeV_%s_%s.txt",nameSample[whichtype].c_str(),LHCsqrts,systString.c_str(),changeParName.c_str()); 
  // else 
  sprintf(fileToSave,"text/paramsPTOverMCJLST_%s%d_%dTeV_%s.txt",nameSample[whichtype].c_str(),mass,LHCsqrts,systString.c_str());
  ofstream os1(fileToSave);
  if (changeParName != "") {
    sprintf(fileToSave,"m%s",changeParName.c_str());  m.SetName(fileToSave);
    sprintf(fileToSave,"n%s",changeParName.c_str());  n.SetName(fileToSave);
    sprintf(fileToSave,"n2%s",changeParName.c_str());  n2.SetName(fileToSave);
    sprintf(fileToSave,"bb%s",changeParName.c_str());  bb.SetName(fileToSave);
    sprintf(fileToSave,"bb2%s",changeParName.c_str());  bb2.SetName(fileToSave);
    sprintf(fileToSave,"fexp%s",changeParName.c_str());  fexp.SetName(fileToSave);
    sprintf(fileToSave,"T%s",changeParName.c_str());  T.SetName(fileToSave);
  }
  (RooArgSet(m,n,n2,bb,bb2,fexp,T)).writeToStream(os1,false);
  os1.close();

  RooRealVar mup("mup","emme", 1.,0.01, 30.);
  RooRealVar nup("nup","enne", 0.93, 0.5, 15.);
  RooRealVar n2up("n2up","enne2", 0.75, 0.5, 15.);
  RooRealVar bbup("bbup","bibi",0.02, 0.00005, 20.0);
  RooRealVar Tup("Tup","tti",0.2,0.00000005,1.);
  RooRealVar bb2up("bb2up","bibi2",0.02, 0.0005, 10.0);
  RooRealVar fexpup("fexpup","f_exp",0.02, 0.0, 1.0);
 
  RooModifTsallis* rt3up = new RooModifTsallis("rt3up","rt3up",*ptoverm,mup,nup,n2up,bbup,bb2up,Tup,fexpup);
  // ws->import(*rt3up);
 
  RooRealVar mdown("mdown","emme", 1.,0.01, 30.);
  RooRealVar ndown("ndown","enne", 0.93, 0.5, 15.);
  RooRealVar n2down("n2down","enne2", 0.75, 0.5, 15.);
  RooRealVar bbdown("bbdown","bibi",0.02, 0.00005, 20.0);
  RooRealVar Tdown("Tdown","tti",0.2,0.00000005,1.);
  RooRealVar bb2down("bb2down","bibi2",0.02, 0.0005, 10.0);
  RooRealVar fexpdown("fexpdown","f_exp",0.02, 0.0, 1.0);

  RooModifTsallis* rt3down = new RooModifTsallis("rt3down","rt3down",*ptoverm,mdown,ndown,n2down,bbdown,bb2down,Tdown,fexpdown);
  // ws->import(*rt3down);

  RooPlot *frame = ptoverm->frame();

  char reducestr[300];
  sprintf(reducestr,"ptoverm > %f && ptoverm < %f",ptoverm->getMin(),ptoverm->getMax());
  
  rdh->plotOn(frame,DataError(RooAbsData::SumW2),Cut(reducestr));
  static RooHist *hpull;
  float chi2 = 0.;

  if (changeParName == "") {
    sprintf(fileToSave,"text/paramsPTOverMCJLST_%s%d_%dTeV_Default.txt",nameSample[whichtype].c_str(),mass,LHCsqrts);
    ifstream is1(fileToSave);
    (RooArgSet(mup,nup,n2up,bbup,bb2up,fexpup,Tup)).readFromStream(is1,false);

    mdown.setVal(fabs(3*mup.getVal() - 2*m.getVal()));
    ndown.setVal(fabs(3*nup.getVal() - 2*n.getVal()));
    n2down.setVal(fabs(3*n2up.getVal() - 2*n2.getVal()));
    bbdown.setVal(fabs(3*bbup.getVal() - 2*bb.getVal()));
    Tdown.setVal(fabs(3*Tup.getVal() - 2*T.getVal()));
    bb2down.setVal(fabs(3*bb2up.getVal() - 2*bb2.getVal()));
    fexpdown.setVal(fabs(3*fexpup.getVal() - 2*fexp.getVal()));

    if (showErrorPDFs) {
      rt3->plotOn(frame,LineColor(kRed),LineStyle(kDashed),Normalization(rdh->sumEntries(),RooAbsReal::NumEvent));
      hpull = frame->pullHist();
      rt3up->plotOn(frame,LineColor(kBlue),Normalization(rdh->sumEntries(),RooAbsReal::NumEvent));
      if (systString.find("Mela") == string::npos) rt3down->plotOn(frame,LineColor(kRed),LineStyle(kDashed),Normalization(rdh->sumEntries(),RooAbsReal::NumEvent));
    } else {
      rt3->plotOn(frame,LineColor(kBlue),Normalization(rdh->sumEntries(),RooAbsReal::NumEvent));
      hpull = frame->pullHist();
    }
  } else {
    mup.setVal(m.getVal() + m.getError());   cout << "mup = " << mup.getVal() << endl;
    nup.setVal(n.getVal() + n.getError());
    n2up.setVal(n2.getVal() + n2.getError());
    bbup.setVal(bb.getVal() + bb.getError());
    Tup.setVal(T.getVal() + T.getError());
    bb2up.setVal(bb2.getVal() + bb2.getError());
    fexpup.setVal(fexp.getVal() + fexp.getError());

    mdown.setVal(m.getVal() - m.getError());  cout << "mdown = " << mdown.getVal() << endl;
    ndown.setVal(n.getVal() - n.getError());
    n2down.setVal(n2.getVal() - n2.getError());
    bbdown.setVal(bb.getVal() - bb.getError());
    Tdown.setVal(T.getVal() - T.getError());
    bb2down.setVal(bb2.getVal() - bb2.getError());
    fexpdown.setVal(fexp.getVal() - fexp.getError());

    rt3->plotOn(frame,LineColor(kBlue),Normalization(rdh->sumEntries(),RooAbsReal::NumEvent));
    hpull = frame->pullHist();
    if (showErrorPDFs) {
      rt3up->plotOn(frame,LineColor(kRed),LineStyle(kDashed),Normalization(rdh->sumEntries(),RooAbsReal::NumEvent));
      rt3down->plotOn(frame,LineColor(kRed),LineStyle(kDashed),Normalization(rdh->sumEntries(),RooAbsReal::NumEvent));
    }
  }

  double *ypulls = hpull->GetY();
  unsigned int nBins = rdh->numEntries();
  unsigned int nFullBins = 0;
  for (unsigned int i = 0; i < nBins; i++) {
    cout << "Pull of bin " << i << " = " << ypulls[i] << endl;
    if (fabs(ypulls[i]) < 5.0) chi2 += ypulls[i]*ypulls[i]; 
    cout << "Partial chi2 = " << chi2 << endl;
    if (fabs(ypulls[i]) > 0.0001 && fabs(ypulls[i]) < 5.0) nFullBins++;
  }
  for (unsigned int i = 0; i < nBins; i++) {
    if (fabs(ypulls[i]) < 0.0001) ypulls[i] = 999.; 
    hpull->SetPointError(i,0.,0.,0.,0.);
  } 
  int nFitPar = fit->floatParsFinal().getSize() - 1;

  TCanvas can("can","The canvas",5.,5.,500.,900.); 
  can.Divide(1,3);

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(22);
  t->SetTextSize(0.06);

  can.cd(1);
  gPad->SetBottomMargin(0.0);
  frame->Draw();
  // gPad->SetLogy(); 
  // Htest->Draw();
  sprintf(fileToSave,"%s %d GeV at %d TeV",nameSample[whichtype].c_str(),mass,LHCsqrts);
  t->DrawLatex(0.6,0.8,fileToSave); 

  can.cd(2);
  gPad->SetLogy(); 
  gPad->SetTopMargin(0.0);
  frame->Draw();
 
  RooPlot* pull = ptoverm->frame(Title("Pull Distribution")) ;
  pull->GetYaxis()->SetTitle("Pull");
  /* pull->SetLabelSize(0.08,"XYZ");
  pull->SetTitleSize(0.08,"XYZ");
  pull->SetTitleOffset(0.6,"Y");
  pull->SetTitleOffset(1.0,"X"); */
  pull->addPlotable(hpull,"P") ; 
  pull->SetMinimum(-6.); 
  pull->SetMaximum(6.); 

  can.cd(3);
  gPad->SetGridy();
  pull->Draw();
  sprintf(fileToSave,"#chi^{2}/n_{DoF} = %4.1f/%d",chi2,nFullBins - nFitPar);
  if (chi2 < 1000.) t->DrawLatex(0.80,0.86,fileToSave);

  // sprintf(fileToSave,"figs/fitPTOverMCJLST_%s%d_%dTeV_%s.pdf",nameSample[whichtype].c_str(),mass,LHCsqrts,systString.c_str());
  sprintf(fileToSave,"newfigs/fitPTOverMCJLST_%s%d_%dTeV_%s.pdf",nameSample[whichtype].c_str(),mass,LHCsqrts,systString.c_str());
  can.SaveAs(fileToSave);

}




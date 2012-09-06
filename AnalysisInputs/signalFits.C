/* 
 * Compute shape parameters parametrization as a function of the invariant mass  for signals and write them in a card fragment.
 * usage: 
 * -set all input variables in Config.h
 * -run with:
 * root -q -b signalFits.C
 * This runs on the 3 final states for 7 and 8 TeV and writes the output in a file (see stdout).
 *
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

/*
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TSystem.h"

#include "RooGlobalFunc.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooFFTConvPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooCBShape.h"
#include "RooPlot.h"
*/

#include "HZZ4LRooPdfs.cc"

//----------> SET INPUT VARIABLES in Config.h
#include "Config.h"
//<----------

using namespace RooFit ;
using namespace std;

//Declaration
void signalFits(int channel, int sqrts);
float WidthValue(float mHStarWidth);

void signalFits()
{
  gSystem->Exec("mkdir -p sigFigs7TeV");
  gSystem->Exec("mkdir -p sigFigs8TeV");

  signalFits(1,7);
  signalFits(2,7);
  signalFits(3,7);

  signalFits(1,8);
  signalFits(2,8);
  signalFits(3,8);

  return;
}

//The actual job
void signalFits(int channel, int sqrts)
{
  string schannel;
  if (channel == 1) schannel = "4mu";
  if (channel == 2) schannel = "4e";
  if (channel == 3) schannel = "2mu2e";
  cout << "Final state = " << schannel << " and sqrt(s) = " << sqrts << endl;

  //Pick the correct mass points and paths
  TString filePath;
  int nPoints;
  double *masses;

  if(sqrts == 7){
    filePath = filePath7TeV;
    nPoints = nPoints7TeV;
    masses = mHVal7TeV;
  }
  else if(sqrts == 8){
    filePath = filePath8TeV;
    nPoints = nPoints8TeV;
    masses = mHVal8TeV;
  }
  else abort();

  filePath.Append(schannel);

  //Prepare to store all the shape parameters for the mass points	
  const int arraySize=200;
  assert(arraySize >= nPoints);

  Double_t a_meanCB[arraySize];  Double_t a_meanCB_err[arraySize]; 
  Double_t a_sigmaCB[arraySize]; Double_t a_sigmaCB_err[arraySize];
  Double_t a_alphaCB[arraySize]; Double_t a_alphaCB_err[arraySize];
  Double_t a_nCB[arraySize];     Double_t a_nCB_err[arraySize];

  Double_t a_fitCovQual[arraySize];
  Double_t a_fitEDM[arraySize];
  Double_t a_fitStatus[arraySize];

  char outfile[192];
  sprintf(outfile,"sigFigs%iTeV",sqrts);
                
  //Loop over the mass points
  for (int i = 0; i < nPoints; i++){
		
    //Open input file with shapes and retrieve the tree
    char tmp_finalInPath[200];
    sprintf(tmp_finalInPath,"/HZZ4lTree_H%i.root",masses[i]);
    string finalInPath = filePath + tmp_finalInPath;

    TFile *f = new TFile(finalInPath.c_str()) ;
    TTree *tree= (TTree*) f->Get("SelectedTree");
    if(tree==NULL){
      cout << "Impossible to retrieve the tree for mass point " << masses[i] <<" GeV " << endl;
      abort();
    }

    double valueWidth = WidthValue(masses[i]);
    double windowVal = max(valueWidth,1.);
    double lowside = 100.;
    if(masses[i] > 300) lowside = 200.;
    double low_M = max( (masses[i] - 15.*windowVal), lowside) ;
    double high_M = min( (masses[i] + 10.*windowVal), 1400.) ;
    cout << "lowM = " << low_M << ", highM = " << high_M << endl;

    //Set the observable and get the RooDataSomething
    RooRealVar ZZMass("ZZMass","ZZMass",low_M,high_M);
    RooDataSet *set2 = new RooDataSet("data","data", tree, RooArgSet(ZZMass));
    RooDataHist *set = (RooDataHist*)set2->binnedClone("datahist","datahist");

    //Theoretical signal model  
    RooRealVar MHStar("MHStar","MHStar",masses[i],0.,2000.);
    MHStar.setConstant(true);
    RooRealVar Gamma_TOT("Gamma_TOT","Gamma_TOT",valueWidth,0.,500.);
    Gamma_TOT.setConstant(true);
    
    //RooGenericPdf SignalTheor("SignalTheor","1./( pow(pow(@1,2)-pow(@0,2),2) + (pow(@0,4)/pow(@1,2))*pow(@2,2) )",RooArgSet(ZZMass,MHStar,Gamma_TOT));
    RooGenericPdf SignalTheor("SignalTheor","1./( pow(pow(@0,2)-pow(@1,2),2) + (pow(@0,4)/pow(@1,2))*pow(@2,2) )",RooArgSet(ZZMass,MHStar,Gamma_TOT));

    //Experimental resolution
    RooRealVar meanCB("meanCB","meanCB",0.,-20.,20.);
    RooRealVar sigmaCB("sigmaCB","sigmaCB",1.,0.01,30.);
    RooRealVar alphaCB("alphaCB","alphaCB",3.,-10.,10.);
    RooRealVar nCB("nCB","nCB",2.,-10.,10.);

    //Initialize to decent values
    float m = masses[i];
    if(channel == 1) sigmaCB.setVal(11.282-0.213437*m+0.0015906*m*m-5.18846e-06*m*m*m+8.05552e-09*m*m*m*m -4.69101e-12*m*m*m*m*m);
    else if(channel == 2) sigmaCB.setVal(3.58777+-0.0252106*m+0.000288074*m*m+-8.11435e-07*m*m*m+7.9996e-10*m*m*m*m);
    else if(channel == 3) sigmaCB.setVal(7.42629+-0.100902*m+0.000660553*m*m+-1.52583e-06*m*m*m+1.2399e-09*m*m*m*m);
    else abort();

    RooCBShape massRes("massRes","crystal ball",ZZMass,meanCB,sigmaCB,alphaCB,nCB);

    //Convolute theoretical shape and resolution
    RooFFTConvPdf sigPDF("sigPDF","sigPDF",ZZMass, SignalTheor,massRes);
    sigPDF.setBufferFraction(0.2);

    //Fit the shape
    RooFitResult *fitRes = sigPDF.fitTo(*set,Save(1));

    a_fitEDM[i] = fitRes->edm();
    a_fitCovQual[i] = fitRes->covQual();
    a_fitStatus[i] = fitRes->status();

    a_meanCB[i]  = meanCB.getVal();
    a_sigmaCB[i] = sigmaCB.getVal();
    a_alphaCB[i]  = alphaCB.getVal();
    a_nCB[i]     = nCB.getVal();

    a_meanCB_err[i]  = meanCB.getError();
    a_sigmaCB_err[i] = sigmaCB.getError();
    a_alphaCB_err[i]  = alphaCB.getError();
    a_nCB_err[i]     = nCB.getError();

    //Plot in the figures directory
    RooPlot *xplot = ZZMass.frame();
    set->plotOn(xplot);
    sigPDF.plotOn(xplot);

    TCanvas canv;
    canv.cd();
    xplot->Draw();

    string tmp_plotFileTitle;
    tmp_plotFileTitle.insert(0,outfile);
    tmp_plotFileTitle += "/fitMass_H";
    char tmp2_plotFileTitle[200];
    sprintf(tmp2_plotFileTitle,"%i_%iTeV_",masses[i],sqrts);
    string plotFileTitle = tmp_plotFileTitle + tmp2_plotFileTitle + schannel + ".gif";

    canv.SaveAs(plotFileTitle.c_str());

  }

  Double_t x_err[arraySize];

  TGraphErrors* gr_meanCB  = new TGraphErrors(nPoints, masses, a_meanCB, x_err, a_meanCB_err);
  TGraphErrors* gr_sigmaCB = new TGraphErrors(nPoints, masses, a_sigmaCB, x_err, a_sigmaCB_err);
  TGraphErrors* gr_alphaCB = new TGraphErrors(nPoints, masses, a_alphaCB, x_err, a_alphaCB_err);
  TGraphErrors* gr_nCB     = new TGraphErrors(nPoints, masses, a_nCB, x_err, a_nCB_err);

  gr_meanCB->Fit("pol4");
  gr_sigmaCB->Fit("pol4");
  gr_alphaCB->Fit("pol4");
  gr_nCB->Fit("pol4");

  TF1 *fit_meanCB  = gr_meanCB->GetListOfFunctions()->First();
  TF1 *fit_sigmaCB = gr_sigmaCB->GetListOfFunctions()->First();
  TF1 *fit_alphaCB = gr_alphaCB->GetListOfFunctions()->First();
  TF1 *fit_nCB     = gr_nCB->GetListOfFunctions()->First();

  gr_meanCB->GetXaxis()->SetTitle("Mean value of the CB function");
  gr_sigmaCB->GetXaxis()->SetTitle("Sigma of the CB function");
  gr_alphaCB->GetXaxis()->SetTitle("Alpha parameter of the CB function");
  gr_nCB->GetXaxis()->SetTitle("n parameter of the CB function");

  gr_meanCB->SetTitle("");
  gr_sigmaCB->SetTitle("");
  gr_alphaCB->SetTitle("");
  gr_nCB->SetTitle("");

  TCanvas canv2;
  canv2.Divide(2,2);

  canv2.cd(1); gr_meanCB->Draw("A*");  fit_meanCB->Draw("SAME");
  canv2.cd(2); gr_sigmaCB->Draw("A*"); fit_sigmaCB->Draw("SAME");
  canv2.cd(3); gr_alphaCB->Draw("A*"); fit_alphaCB->Draw("SAME");
  canv2.cd(4); gr_nCB->Draw("A*");     fit_nCB->Draw("SAME");

  string tmp_paramPlotFileTitle;
  tmp_paramPlotFileTitle.insert(0,outfile);
  tmp_paramPlotFileTitle += "/fitParam_";
  char tmp2_paramPlotFileTitle[200];
  sprintf(tmp2_paramPlotFileTitle,"%iTeV_",sqrts);
  string paramPlotFileTitle = tmp_paramPlotFileTitle + tmp2_paramPlotFileTitle + schannel + ".gif";

  canv2.SaveAs(paramPlotFileTitle.c_str());

  char tmp_outCardName[200];
  sprintf(tmp_outCardName,"CardFragments/signalFunctions_%iTeV_",sqrts);
  string outCardName = tmp_outCardName + schannel + ".txt";
  ofstream ofsCard(outCardName.c_str(),fstream::out);
  ofsCard << "## signal functions --- no spaces! ##" << endl;
  ofsCard << "signalShape n_CB " << fit_nCB->GetParameter(0) << "+(" << fit_nCB->GetParameter(1) << "*@0)+(" << fit_nCB->GetParameter(2) << "*@0*@0)+(" << fit_nCB->GetParameter(3) << "*@0*@0*@0)+(" << fit_nCB->GetParameter(4) << "*@0*@0*@0*@0)" << endl;
  ofsCard << "signalShape alpha_CB " << fit_alphaCB->GetParameter(0) << "+(" << fit_alphaCB->GetParameter(1) << "*@0)+(" << fit_alphaCB->GetParameter(2) << "*@0*@0)+(" << fit_alphaCB->GetParameter(3) << "*@0*@0*@0)+(" << fit_alphaCB->GetParameter(4) << "*@0*@0*@0*@0)" <<endl;
  ofsCard << "signalShape mean_CB " << fit_meanCB->GetParameter(0) << "+(" << fit_meanCB->GetParameter(1) << "*@0)+(" << fit_meanCB->GetParameter(2) << "*@0*@0)+(" << fit_meanCB->GetParameter(3) << "*@0*@0*@0)+(" <<  fit_meanCB->GetParameter(4) << "*@0*@0*@0*@0)" <<endl;
  ofsCard << "signalShape sigma_CB " << fit_sigmaCB->GetParameter(0) << "+(" << fit_sigmaCB->GetParameter(1) << "*@0)+(" << fit_sigmaCB->GetParameter(2) << "*@0*@0)+(" << fit_sigmaCB->GetParameter(3) << "*@0*@0*@0)+(" <<  fit_sigmaCB->GetParameter(4) << "*@0*@0*@0*@0)" <<endl;
  ofsCard << endl;
  return;
}

float WidthValue(float mHStarWidth)
{
  ostringstream MassString;
  MassString << mHStarWidth;
  
  ifstream widthFile("widthvalues.txt");
  
  string line;

  bool FindedMass = false;

  float Gamma_ggCal, Gamma_ZZCal, Gamma_TOTCal;
  
  while (getline(widthFile,line)) {
    if( line == "" || line[0] == '#' ) continue;
    
    if(line[0]== MassString.str()[0] && line[1]== MassString.str()[1] && line[2]== MassString.str()[2]){
      
      stringstream stringline;
      stringline << line;    
      string masschar;
      stringline>>masschar>>Gamma_ggCal>>Gamma_ZZCal>>Gamma_TOTCal;
      
      FindedMass = true;
    }
  }
  if(!FindedMass) abort();

  return Gamma_TOTCal;
}

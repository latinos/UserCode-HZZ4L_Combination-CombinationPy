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

//----------> SET INPUT VARIABLES in Config.h
#include "Config.h"
//<----------

//Parameters to choose which samples to examine
//HCP option is useggH=true, useVBF,useVH = false and usedijet,usenondijet=true
bool debug = false;
bool useggH = true;
bool useVBF = false;
bool useVH = false;
bool usedijet = true;
bool usenondijet = true;

TFile* massfits;

using namespace RooFit;
using namespace std;

//Declaration
void signalFits(int channel, int sqrts);
float WidthValue(float mHStarWidth);
float highparameter(TF1 *generatedfit, int gamma);

void signalFits(){
  gSystem->Exec("mkdir -p sigFigs7TeV");
  gSystem->Exec("mkdir -p sigFigs8TeV");

  gSystem->Load("../CreateDatacards/CMSSW_5_2_5/lib/slc5_amd64_gcc462/libHiggsAnalysisCombinedLimit.so");

  if (!usedijet && !usenondijet){
    cout << "Neither dijet tagging category was chosen. Please choose one or both." << endl;
    abort();
  }

  if (!useggH && !useVBF && !useVH){
    cout << "Please choose at least one production method: ggH, VBF, VH."<<endl;
    abort();
  }

  signalFits(1,7);
  signalFits(2,7);
  signalFits(3,7);
  signalFits(1,8);
  signalFits(2,8);
  signalFits(3,8);

  return;
}

//The actual job
void signalFits(int channel, int sqrts){
  string schannel;
  if (channel == 1) schannel = "4mu";
  if (channel == 2) schannel = "4e";
  if (channel == 3) schannel = "2e2mu";
  cout << "Final state = " << schannel << " and sqrt(s) = " << sqrts << endl;

  //Pick the correct mass points and paths
  TString filePath;
  int nPointsggH,nPointsVBF,nPointsVH;
  double* massesggH,*massesVBF,*massesVH;

  if (sqrts==7) {
    nPointsggH = nPoints7TeV;
    massesggH  = mHVal7TeV;
    nPointsVBF = nVBFPoints7TeV;
    massesVBF  = mHVBFVal7TeV;
    nPointsVH = nVHPoints7TeV;
    massesVH  = mHVHVal7TeV;
    filePath = filePath7TeV;  
  } else if (sqrts==8) {
    nPointsggH = nPoints8TeV;
    massesggH  = mHVal8TeV;
    nPointsVBF = nVBFPoints8TeV;
    massesVBF  = mHVBFVal8TeV;
    nPointsVH = nVHPoints8TeV;
    massesVH  = mHVHVal8TeV;
    filePath = filePath8TeV;
  }
  else abort();

  int nPoints;
  double masses[200];
  for (int i=0;i<200;i++){
    masses[i]=-1;
  }
  if (useggH){
    for (int i=0; i<nPointsggH; i++){
      masses[i]=massesggH[i];
    }
    nPoints=nPointsggH;
  }else if(!useggH && useVBF){
    for (int i=0; i<nPointsVBF; i++){
      masses[i]=massesVBF[i];
    }
    nPoints=nPointsVBF;
  }else if(!useggH && !useVBF && useVH){
    for (int i=0; i<nPointsVH; i++){
      masses[i]=massesVH[i];
    }
      nPoints=nPointsVH;
  }


  bool flag=0;
  //If any mass points are in VBF/VH but not in ggH, this will include them
  if (useggH && useVBF){
    for (int i=0; i<nPointsVBF; i++){
      flag=1;
      for (int j=0; j<nPoints; j++){
	if (massesVBF[i]==masses[j]) flag=0;
      }
      if (flag==1){
	masses[nPoints]=massesVBF[i];
	nPoints++;
      }
    }
  }
  if ((useggH || useVBF) && useVH){
    for (int i=0; i<nPointsVH; i++){
      flag=1;
      for (int j=0; j<nPoints; j++){
	if (massesVH[i]==masses[j]) flag=0;
      }
      if (flag==1){
	masses[nPoints]=massesVH[i];
	nPoints++;
      }
    }
  }

  filePath.Append(schannel=="2e2mu"?"2mu2e":schannel);

  //Prepare to store all the shape parameters for the mass points	
  const int arraySize=200;
  assert(arraySize >= nPoints);

  Double_t a_meanCB[arraySize];    Double_t a_meanCB_err[arraySize]; 
  Double_t a_sigmaCB[arraySize];   Double_t a_sigmaCB_err[arraySize];
  Double_t a_alphaCB_1[arraySize]; Double_t a_alphaCB_1_err[arraySize];
  Double_t a_nCB_1[arraySize];     Double_t a_nCB_1_err[arraySize];
  Double_t a_Gamma[arraySize];     Double_t a_Gamma_err[arraySize];
  Double_t a_alphaCB_2[arraySize]; Double_t a_alphaCB_2_err[arraySize];
  Double_t a_nCB_2[arraySize];     Double_t a_nCB_2_err[arraySize];

  Double_t a_fitCovQual[arraySize];
  Double_t a_fitEDM[arraySize];
  Double_t a_fitStatus[arraySize];

  char outfile[192];
  sprintf(outfile,"sigFigs%iTeV",sqrts);
                
  //Loop over the mass points
  for (int i = 0; i < nPoints; i++){

    bool flagVBF==0;
    bool flagVH==0;

    if (debug && masses[i]!=126.) continue;
		
    //Open input file with shapes and retrieve the tree
    char tmp_finalInPathggH[200],tmp_finalInPathVBF[200],tmp_finalInPathZH[200],tmp_finalInPathWH[200],tmp_finalInPathttH[200];
    sprintf(tmp_finalInPathggH,"/HZZ4lTree_H%i.root",masses[i]);
    sprintf(tmp_finalInPathVBF,"/HZZ4lTree_VBFH%i.root",masses[i]);
    sprintf(tmp_finalInPathZH,"/HZZ4lTree_ZH%i.root",masses[i]);
    sprintf(tmp_finalInPathWH,"/HZZ4lTree_WH%i.root",masses[i]);
    sprintf(tmp_finalInPathttH,"/HZZ4lTree_ttH%i.root",masses[i]);
    string finalInPathggH = filePath + tmp_finalInPathggH;
    string finalInPathVBF = filePath + tmp_finalInPathVBF;
    string finalInPathZH = filePath + tmp_finalInPathZH;
    string finalInPathWH = filePath + tmp_finalInPathWH;
    string finalInPathttH = filePath + tmp_finalInPathttH;

    TChain *f = new TChain("SelectedTree");
    if (useggH){
      if (i<nPointsggH) f->Add(finalInPathggH.c_str());
      else{
	cout<<"No ggH sample at this mass point."<<endl;
      }
    }
    else if (!useggH && useVBF){
      if (i<nPointsVBF) f->Add(finalInPathVBF.c_str());
      else{
	cout<<"No qqH sample at this mass point."<<endl;
      }
    }
    else if (!useggH && !useVBF && useVH){
      if (i<nPointsVH){
	f->Add(finalInPathZH.c_str());
	f->Add(finalInPathWH.c_str());
	f->Add(finalInPathttH.c_str());
      }
      else{
	cout<<"No VH samples at this mass point."<<endl;
      }
    }    
    if (useggH && useVBF){
      for (int j=0;j<nPointsVBF;j++){
	if (massesVBF[j]=masses[i]) flagVBF=1;
      }
      if (flagVBF==1) f->Add(finalInPathVBF.c_str());
      if (flagVBF==0) cout<<"No qqH sample at this mass point."<<endl;
    }
    if ((useggH || useVBF) && useVH){
      for (int j=0;j<nPointsVH;j++){
	if (massesVH[j]=masses[i]) flagVH=1;
      }
      if (flagVH==1){
	f->Add(finalInPathZH.c_str());
	f->Add(finalInPathWH.c_str());
	f->Add(finalInPathttH.c_str());
      }
      if (flagVH==0) cout<<"No VH samples at this mass point."<<endl;
    }

    double valueWidth = WidthValue(masses[i]);
    double windowVal = max(valueWidth,1.);
    double lowside = 100.;
    if(masses[i] > 300) lowside = 200.;
    double low_M = max( (masses[i] - 15.*windowVal), lowside) ;
    double high_M = min( (masses[i] + 10.*windowVal), 1400.);

    if(masses[i] > 399.){
      low_M = max( (masses[i] - 2.*windowVal), 250.) ;
      high_M = min( (masses[i] + 2.*windowVal), 1600.);
    }

    cout << "lowM = " << low_M << ", highM = " << high_M << endl;

    //Set the observable and get the RooDataSomething
    RooRealVar ZZMass("ZZMass","ZZMass",low_M,high_M);
    RooRealVar MC_weight("MC_weight","MC_weight",0.,10.);
    RooRealVar NJets("NJets","NJets",0.,100.);
    RooRealVar genProcessId("genProcessId","genProcessId",0.,150.);

    if(channel == 2) ZZMass.setBins(50);
    if(channel == 3) ZZMass.setBins(50);

    RooDataSet* set;

    if (!flagVH){
      if (usedijet && !usenondijet){
	set = new RooDataSet("data","data", f, RooArgSet(ZZMass,MC_weight,NJets), "NJets>1", "MC_weight");
      }
      else if (usenondijet && !usedijet){
	set = new RooDataSet("data","data", f, RooArgSet(ZZMass,MC_weight,NJets), "NJets<2", "MC_weight");
      } else{
	set = new RooDataSet("data","data", f, RooArgSet(ZZMass,MC_weight), "", "MC_weight");
      }
    }
    if (flagVH){
      if (usedijet && !usenondijet){
	set = new RooDataSet("data","data", f, RooArgSet(ZZMass,MC_weight,NJets,genProcessId), "NJets>1 && (genProcessId==24 || genProcessId==26 || genProcessId==121 || genProcessId==122 || genProcessId== 10011 )", "MC_weight");
      }
      else if (usenondijet && !usedijet){
	set = new RooDataSet("data","data", f, RooArgSet(ZZMass,MC_weight,NJets,genProcessId), "NJets<2 && (genProcessId==24 || genProcessId==26 || genProcessId==121 || genProcessId==122 || genProcessId== 10011)", "MC_weight");
      } else{
	set = new RooDataSet("data","data", f, RooArgSet(ZZMass,MC_weight,NJets,genProcessId), "(genProcessId==24 || genProcessId==26 || genProcessId==121 || genProcessId==122 || genProcessId== 10011)", "MC_weight");
      }
    }

    RooDataHist *hist = (RooDataHist*)set->binnedClone("datahist","datahist"); 

    //Theoretical signal model  
    RooRealVar MHStar("MHStar","MHStar",masses[i],0.,2000.);
    MHStar.setConstant(true);
    RooRealVar Gamma_TOT("Gamma_TOT","Gamma_TOT",valueWidth,0.,700.);
    if(masses[i] < 399.) Gamma_TOT.setConstant(true);
    RooRealVar one("one","one",1.0);
    one.setConstant(kTRUE);

    RooGenericPdf SignalTheor("SignalTheor","(@0)/( pow(pow(@0,2)-pow(@1,2),2) + pow(@0,2)*pow(@2,2) )",RooArgSet(ZZMass,MHStar,Gamma_TOT));
    RooRelBWUFParam SignalTheorLM("signalTheorLM","signalTheorLM",ZZMass,MHStar,one);

    //Experimental resolution
    RooRealVar meanCB("meanCB","meanCB",0.,-25.,25.);
    RooRealVar sigmaCB("sigmaCB","sigmaCB",1.,0.,5.);
    RooRealVar sigmaCB_high("sigmaCB_high","sigmaCB_high",5.,0.,150.);
    RooRealVar alphaCB_1("alphaCB_1","alphaCB_1",2.,0.,50.);
    RooRealVar nCB_1("nCB_1","nCB_1",2.,0.,12.);
    RooRealVar alphaCB_2("alphaCB_2","alphaCB_2",2.,0.,50.);
    RooRealVar nCB_2("nCB_2","nCB_2",20.);
    nCB_2.setConstant(kTRUE);

    //Initialize to decent values
    float m = masses[i];

    /*if(channel == 1){
      sigmaCB_R.setVal(11.282-0.213437*m+0.0015906*m*m-5.18846e-06*m*m*m+8.05552e-09*m*m*m*m -4.69101e-12*m*m*m*m*m);
      sigmaCB_L.setVal(11.282-0.213437*m+0.0015906*m*m-5.18846e-06*m*m*m+8.05552e-09*m*m*m*m -4.69101e-12*m*m*m*m*m);
    }
    else if(channel == 2) sigmaCB_R.setVal(3.58777+-0.0252106*m+0.000288074*m*m+-8.11435e-07*m*m*m+7.9996e-10*m*m*m*m);
    else if(channel == 3) sigmaCB_R.setVal(7.42629+-0.100902*m+0.000660553*m*m+-1.52583e-06*m*m*m+1.2399e-09*m*m*m*m);*/

    sigmaCB.setVal(-4.56178+0.123209*m-0.00107193*m*m+4.5413e-06*m*m*m-8.19429e-09*m*m*m*m+4.75955e-12*m*m*m*m*m);
    sigmaCB_high.setVal(151.967-0.939938*m+0.00173551*m*m-8.26677e-07*m*m*m);

    //else abort();

    RooDoubleCB massRes("massRes","Double Crystal Ball",ZZMass,meanCB,sigmaCB,alphaCB_1,nCB_1,alphaCB_2,nCB_2);
    RooDoubleCB massResH("massResH","DCB Highmass",ZZMass,meanCB,sigmaCB_high,alphaCB_1,nCB_1,alphaCB_2,nCB_2);

    //Convolute theoretical shape and resolution
    RooFFTConvPdf *sigPDF;
    if(masses[i] < 399.) sigPDF = new RooFFTConvPdf("sigPDF","sigPDF",ZZMass,SignalTheorLM,massRes);
    else sigPDF = new RooFFTConvPdf("sigPDF","sigPDF",ZZMass,SignalTheor,massResH);
    sigPDF->setBufferFraction(0.2);

    RooPlot *xplot = ZZMass.frame();
    TCanvas *canv = new TCanvas("canv","canv",1200,800);

    string tmp_plotFileTitle;
    tmp_plotFileTitle.insert(0,outfile);
    tmp_plotFileTitle += "/fitMass_";
    char tmp2_plotFileTitle[200];
    sprintf(tmp2_plotFileTitle,"%i_%iTeV_",masses[i],sqrts);
    string plotFileTitle = tmp_plotFileTitle + tmp2_plotFileTitle + schannel;
    TString rootTitle = tmp_plotFileTitle + tmp2_plotFileTitle + schannel;
    if (useggH){
      plotFileTitle+="_ggH";
      rootTitle+="_ggH";
    }
    if (useVBF){
      plotFileTitle+="_VBF";
      rootTitle+="_VBF";
    }
    if (useVH){
      plotFileTitle+="_VH";
      rootTitle+="_VH";
    }
    if (usedijet && !usenondijet){
      plotFileTitle += "_1";
      rootTitle += "_1";
    }
    if (usenondijet && !usedijet){
      plotFileTitle += "_0";
      rootTitle += "_0";
    }

    double mass,mean,sigma,a1,n1,a2,n2,gamma;
    TH1F* parameters;

    //Fit the shape
    massfits = new TFile(rootTitle + ".root","RECREATE");
    
    RooFitResult *fitRes = sigPDF->fitTo(*hist,Save(1), SumW2Error(kTRUE));

    a_fitEDM[i] = fitRes->edm();
    a_fitCovQual[i] = fitRes->covQual();
    a_fitStatus[i] = fitRes->status();

    mass = masses[i];
    mean = meanCB.getVal();
    if (mass > 399.) sigma = sigmaCB_high.getVal();
    else sigma = sigmaCB.getVal();
    a1=alphaCB_1.getVal();
    n1=nCB_1.getVal();
    a2=alphaCB_2.getVal();
    n2=nCB_2.getVal();
    gamma=Gamma_TOT.getVal();

    a_meanCB[i]  = mean;
    a_sigmaCB[i] = sigma;
    a_alphaCB_1[i]  = a1;
    a_nCB_1[i]     = n1;
    a_alphaCB_2[i]  = a2;
    a_nCB_2[i]     = n2;
    a_Gamma[i] = gamma;
    
    a_meanCB_err[i]  = meanCB.getError();
    if (masses[i] > 399.) a_sigmaCB_err[i] = sigmaCB_high.getError();
    else a_sigmaCB_err[i] = sigmaCB.getError();
    a_alphaCB_1_err[i]  = alphaCB_1.getError();
    a_nCB_1_err[i]     = nCB_1.getError();
    a_alphaCB_2_err[i]  = alphaCB_2.getError();
    a_nCB_2_err[i]     = 0;
    if(masses[i] > 399.) a_Gamma_err[i] = Gamma_TOT.getError();
    else a_Gamma_err[i] = 0.;

    //Plot in the figures directory
    hist->plotOn(xplot);
    sigPDF->plotOn(xplot);
    canv->cd();
    xplot->Draw();

    string plotgif = plotFileTitle + ".gif";
    canv->SaveAs(plotgif.c_str());
    massfits->cd();
    set->Write("MassData");
    parameters = new TH1F("","",8,0,8);
    parameters->Fill(0,mass);
    parameters->SetBinError(1,0);
    parameters->Fill(1,mean);
    parameters->SetBinError(2,a_meanCB_err[i]);
    parameters->Fill(2,sigma);
    parameters->SetBinError(3,a_sigmaCB_err[i]);
    parameters->Fill(3,a1);
    parameters->SetBinError(4,a_alphaCB_1_err[i]);
    parameters->Fill(4,n1);
    parameters->SetBinError(5,a_nCB_1_err[i]);
    parameters->Fill(5,a2);
    parameters->SetBinError(6,a_alphaCB_2_err[i]);
    parameters->Fill(6,n2);
    parameters->SetBinError(7,a_nCB_2_err[i]);
    parameters->Fill(7,gamma);
    parameters->SetBinError(8,a_Gamma_err[i]);
    parameters->Write("Parameters");
    massfits->Close();
    
  }

  if (debug) return;

  TGraph* gr_meanCB  = new TGraph(nPoints, masses, a_meanCB);
  TGraph* gr_sigmaCB = new TGraph(nPoints, masses, a_sigmaCB);
  TGraph* gr_alphaCB_1 = new TGraph(nPoints, masses, a_alphaCB_1);
  TGraph* gr_nCB_1     = new TGraph(nPoints, masses, a_nCB_1);
  TGraph* gr_alphaCB_2 = new TGraph(nPoints, masses, a_alphaCB_2);
  TGraph* gr_nCB_2     = new TGraph(nPoints, masses, a_nCB_2);
  TGraph* gr_Gamma   = new TGraph(nPoints, masses, a_Gamma);

  TF1 *paramfit = new TF1("paramfit","(x<400)*([0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x)+(x>=400)*(([0]+([1]-[6])*400+([2]-[7])*pow(400,2)+([3]-[8])*pow(400,3)+[4]*pow(400,4)+[5]*pow(400,5))+[6]*x+[7]*x*x+[8]*x*x*x)",115,1000);
  TF1 *gammafit = new TF1("gammafit","(x<400)*([0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x)+(x>=400)*(([0]+([6]-[2])*pow(400,2)+2*([7]-[3])*pow(400,3)-3*[4]*pow(400,4)-4*[5]*pow(400,5))+([1]+2*([2]-[6])*400+3*([3]-[7])*pow(400,2)+4*[4]*pow(400,3)+5*[5]*pow(400,4))*x+[6]*x*x+[7]*x*x*x)",115,1000);

  gr_meanCB->Fit("paramfit");
  gr_sigmaCB->Fit("paramfit");
  gr_alphaCB_1->Fit("paramfit");
  gr_nCB_1->Fit("paramfit");
  gr_alphaCB_2->Fit("paramfit");
  gr_nCB_2->Fit("pol0");
  gr_Gamma->Fit("gammafit");

  TF1 *fit_meanCB  = gr_meanCB->GetListOfFunctions()->First();
  TF1 *fit_sigmaCB = gr_sigmaCB->GetListOfFunctions()->First();
  TF1 *fit_alphaCB_1 = gr_alphaCB_1->GetListOfFunctions()->First();
  TF1 *fit_nCB_1     = gr_nCB_1->GetListOfFunctions()->First();
  TF1 *fit_alphaCB_2 = gr_alphaCB_2->GetListOfFunctions()->First();
  TF1 *fit_nCB_2     = gr_nCB_2->GetListOfFunctions()->First();
  TF1 *fit_Gamma   = gr_Gamma->GetListOfFunctions()->First();

  gr_meanCB->GetXaxis()->SetTitle("Mean value of the DCB function");
  gr_sigmaCB->GetXaxis()->SetTitle("Sigma of the DCB function");
  gr_alphaCB_1->GetXaxis()->SetTitle("Alpha parameter of the R leg of DCB function");
  gr_nCB_1->GetXaxis()->SetTitle("n parameter of the R leg of DCB function");
  gr_alphaCB_2->GetXaxis()->SetTitle("Alpha parameter of the L leg of DCB function");
  gr_nCB_2->GetXaxis()->SetTitle("n parameter of the L leg of DCB function");
  gr_Gamma->GetXaxis()->SetTitle("#Gamma of the BW function");

  gr_meanCB->SetTitle("");
  gr_sigmaCB->SetTitle("");
  gr_alphaCB_1->SetTitle("");
  gr_nCB_1->SetTitle("");
  gr_alphaCB_2->SetTitle("");
  gr_nCB_2->SetTitle("");
  gr_Gamma->SetTitle("");

  TCanvas *canv2 = new TCanvas("canv2","canv2",1600,800);
  canv2->Divide(4,2);

  canv2->cd(1); gr_meanCB->Draw("A*");  fit_meanCB->Draw("SAME");
  canv2->cd(2); gr_sigmaCB->Draw("A*"); fit_sigmaCB->Draw("SAME");
  canv2->cd(3); gr_alphaCB_1->Draw("A*"); fit_alphaCB_1->Draw("SAME");
  canv2->cd(4); gr_nCB_1->Draw("A*");     fit_nCB_1->Draw("SAME");
  canv2->cd(5); gr_alphaCB_2->Draw("A*"); fit_alphaCB_2->Draw("SAME");
  canv2->cd(6); gr_nCB_2->Draw("A*");     fit_nCB_2->Draw("SAME");
  canv2->cd(7); gr_Gamma->Draw("A*");   fit_Gamma->Draw("SAME");

  string tmp_paramPlotFileTitle;
  tmp_paramPlotFileTitle.insert(0,outfile);
  tmp_paramPlotFileTitle += "/fitParam_";
  char tmp2_paramPlotFileTitle[200];
  sprintf(tmp2_paramPlotFileTitle,"%iTeV_",sqrts);
  string paramPlotFileTitle = tmp_paramPlotFileTitle + tmp2_paramPlotFileTitle + schannel;
  if (useggH) paramPlotFileTitle+="_ggH";
  if (useVBF) paramPlotFileTitle+="_VBF";
  if (useVH) paramPlotFileTitle+="_VH";
  string paramgif=paramPlotFileTitle;
  if (usedijet && !usenondijet){
    paramgif+="_0.gif";
  }else if (usenondijet && !usedijet){
    paramgif+="_1.gif";
  }else if (usenondijet && usedijet){
    paramgif+="deriv5.gif";
  }
  canv2->SaveAs(paramgif.c_str());

  char tmp_outCardName[200];
  sprintf(tmp_outCardName,"%iTeV_",sqrts);
  string prependName = "CardFragments/signalFunctions_";
  string appendName = ".txt";
  string outCardName =  prependName + tmp_outCardName + schannel + appendName;

  float highn1,higha1,higha2,highmean,highsigma,highgamma1,highgamma2;
  highn1=highparameter(fit_nCB_1,0);
  higha1=highparameter(fit_alphaCB_1,0);
  higha2=highparameter(fit_alphaCB_2,0);
  highmean=highparameter(fit_meanCB,0);
  highsigma=highparameter(fit_sigmaCB,0);  
  highgamma1=highparameter(fit_Gamma,1);
  highgamma2=highparameter(fit_Gamma,2);

  ofstream ofsCard;
  if (usedijet && usenondijet){
    ofsCard.open(outCardName.c_str(),fstream::out);
    ofsCard << "## signal functions --- no spaces! ##" << endl;
    ofsCard << "signalShape n_CB " << fit_nCB_1->GetParameter(0) << "+(" << fit_nCB_1->GetParameter(1) << "*@0)+(" << fit_nCB_1->GetParameter(2) << "*@0*@0)+(" << fit_nCB_1->GetParameter(3) << "*@0*@0*@0)+(" << fit_nCB_1->GetParameter(4) << "*@0*@0*@0*@0)+(" << fit_nCB_1->GetParameter(5) << "*@0*@0*@0*@0*@0)" << endl;
    ofsCard << "signalShape alpha_CB " << fit_alphaCB_1->GetParameter(0) << "+(" << fit_alphaCB_1->GetParameter(1) << "*@0)+(" << fit_alphaCB_1->GetParameter(2) << "*@0*@0)+(" << fit_alphaCB_1->GetParameter(3) << "*@0*@0*@0)+(" << fit_alphaCB_1->GetParameter(4) << "*@0*@0*@0*@0)+(" << fit_alphaCB_1->GetParameter(5) << "*@0*@0*@0*@0*@0)" << endl;
    ofsCard << "signalShape n2_CB " << fit_nCB_2->GetParameter(0) << endl;
    ofsCard << "signalShape alpha2_CB " << fit_alphaCB_2->GetParameter(0) << "+(" << fit_alphaCB_2->GetParameter(1) << "*@0)+(" << fit_alphaCB_2->GetParameter(2) << "*@0*@0)+(" << fit_alphaCB_2->GetParameter(3) << "*@0*@0*@0)+(" << fit_alphaCB_2->GetParameter(4) << "*@0*@0*@0*@0)+(" << fit_alphaCB_2->GetParameter(5) << "*@0*@0*@0*@0*@0)" << endl;
    ofsCard << "signalShape mean_CB " << fit_meanCB->GetParameter(0) << "+(" << fit_meanCB->GetParameter(1) << "*@0)+(" << fit_meanCB->GetParameter(2) << "*@0*@0)+(" << fit_meanCB->GetParameter(3) << "*@0*@0*@0)+(" << fit_meanCB->GetParameter(4) << "*@0*@0*@0*@0)+(" << fit_meanCB->GetParameter(5) << "*@0*@0*@0*@0*@0)" << endl;
    ofsCard << "signalShape sigma_CB " << fit_sigmaCB->GetParameter(0) << "+(" << fit_sigmaCB->GetParameter(1) << "*@0)+(" << fit_sigmaCB->GetParameter(2) << "*@0*@0)+(" << fit_sigmaCB->GetParameter(3) << "*@0*@0*@0)+(" << fit_sigmaCB->GetParameter(4) << "*@0*@0*@0*@0)+(" << fit_sigmaCB->GetParameter(5) << "*@0*@0*@0*@0*@0)" << endl;
    //ofsCard << "signalShape gamma_BW " << fit_Gamma->GetParameter(0) << "+(" << fit_Gamma->GetParameter(1) << "*@0)+(" << fit_Gamma->GetParameter(2) << "*@0*@0)+(" << fit_Gamma->GetParameter(3) << "*@0*@0*@0)+(" << fit_Gamma->GetParameter(4) << "*@0*@0*@0*@0)+(" << fit_Gamma->GetParameter(5) << "*@0*@0*@0*@0*@0)" << endl;
    ofsCard << "HighMasssignalShape n_CB " << highn1 << "+(" << fit_nCB_1->GetParameter(6) << "*@0)+(" << fit_nCB_1->GetParameter(7) << "*@0*@0)+(" << fit_nCB_1->GetParameter(8) << "*@0*@0*@0)" <<endl;
    ofsCard << "HighMasssignalShape alpha_CB " << higha1 << "+(" << fit_alphaCB_1->GetParameter(6) << "*@0)+(" << fit_alphaCB_1->GetParameter(7) << "*@0*@0)+(" << fit_alphaCB_1->GetParameter(8) << "*@0*@0*@0)" << endl;
    ofsCard << "HighMasssignalShape n2_CB " << fit_nCB_2->GetParameter(0) << endl;
    ofsCard << "HighMasssignalShape alpha2_CB " << higha2 << "+(" << fit_alphaCB_2->GetParameter(6) << "*@0)+(" << fit_alphaCB_2->GetParameter(7) << "*@0*@0)+(" << fit_alphaCB_2->GetParameter(8) << "*@0*@0*@0)" << endl;
    ofsCard << "HighMasssignalShape mean_CB " << highmean << "+(" << fit_meanCB->GetParameter(6) << "*@0)+(" << fit_meanCB->GetParameter(7) << "*@0*@0)+(" << fit_meanCB->GetParameter(8) << "*@0*@0*@0)" << endl;
    ofsCard << "HighMasssignalShape sigma_CB " << highsigma << "+(" << fit_sigmaCB->GetParameter(6) << "*@0)+(" << fit_sigmaCB->GetParameter(7) << "*@0*@0)+(" << fit_sigmaCB->GetParameter(8) << "*@0*@0*@0)" << endl;
    ofsCard << "HighMasssignalShape gamma_BW " << highgamma1 << "+(" << highgamma2 << "*@0)+(" << fit_Gamma->GetParameter(6) << "*@0*@0)+(" << fit_Gamma->GetParameter(7) << "*@0*@0*@0)" << endl;
    ofsCard << endl;
  }

  return;
}

float WidthValue(float mHStarWidth){
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

float highparameter(TF1 *generatedfit, int gamma){
  float highzero;

  if (gamma==0){
    highzero=generatedfit->GetParameter(0);
    highzero+=(generatedfit->GetParameter(1)-generatedfit->GetParameter(6))*400;
    highzero+=(generatedfit->GetParameter(2)-generatedfit->GetParameter(7))*pow(400,2);
    highzero+=(generatedfit->GetParameter(3)-generatedfit->GetParameter(8))*pow(400,3);
    highzero+=generatedfit->GetParameter(4)*pow(400,4);
    highzero+=generatedfit->GetParameter(5)*pow(400,5);
  } else if (gamma==1){
    highzero=generatedfit->GetParameter(0);
    highzero+=(generatedfit->GetParameter(6)-generatedfit->GetParameter(2))*pow(400,2);
    highzero+=2*(generatedfit->GetParameter(7)-generatedfit->GetParameter(3))*pow(400,3);
    highzero-=3*(generatedfit->GetParameter(4))*pow(400,4);
    highzero-=4*(generatedfit->GetParameter(5))*pow(400,5);
  } else if (gamma==2){
    highzero=generatedfit->GetParameter(1);
    highzero+=2*(generatedfit->GetParameter(2)-generatedfit->GetParameter(6))*400;
    highzero+=3*(generatedfit->GetParameter(3)-generatedfit->GetParameter(7))*pow(400,2);
    highzero+=4*generatedfit->GetParameter(4)*pow(400,3);
    highzero+=5*generatedfit->GetParameter(5)*pow(400,4);
  }
  return highzero;
}

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

using namespace RooFit ;
using namespace std;

//Declaration
void signalFits(int channel, int sqrts, int icomp);
float WidthValue(float mHStarWidth);
bool writeFits = false;
bool writeParam = true;

void compareSignalFits()
{
  gSystem->Exec("mkdir -p compFigs7TeV");
  gSystem->Exec("mkdir -p compFigs8TeV130613");

  gSystem->Load("../CreateDatacards/CMSSW_6_1_1/lib/slc5_amd64_gcc472/libHiggsAnalysisCombinedLimit.so");

  for(int ich=1;ich<4;ich++)
    for(int ien=8;ien<9;ien++)
      for(int icomp=0;icomp<3;icomp++)
	signalFits(ien,ich,icomp);
  //signalFits(8,1,1);
}

//The actual job
void signalFits(int sqrts, int channel, int icomp)
{
  TCanvas canv;
  string schannel;
  if (channel == 1) schannel = "4mu";
  if (channel == 2) schannel = "4e";
  if (channel == 3) schannel = "2e2mu";

  string scomp="H";
  if (icomp == 0) scomp = "ZH";
  if (icomp == 1) scomp = "WH";
  if (icomp == 2) scomp = "ttH";
  if (icomp == 3) scomp = "VBFH";

  cout << "Final state = " << schannel << " and sqrt(s) = " << sqrts << " "<<scomp<<endl;

  //Pick the correct mass points and paths
  TString filePath;
  int nPoints;
  double *masses;
  TString filePathH;

  if(sqrts == 7){
    filePath = filePath7TeV;
    nPoints = nPoints7TeV;
    masses = mHVal7TeV;
  }
  else if(sqrts == 8){
    filePath =filePath8TeV;//"root://lxcms02//data/Higgs/rootuplesOut/130613/PRODFSR_8TeV/";
    filePathH=filePath8TeV;//"root://lxcms02//data/Higgs/rootuplesOut/130613/PRODFSR_8TeV/";
    nPoints = nPoints8TeV;
    masses = mHVal8TeV;
  }
  else abort();

  filePathH.Append(schannel=="2e2mu"?"2mu2e":schannel);
  filePath.Append(schannel=="2e2mu"?"2mu2e":schannel);

  //Prepare to store all the shape parameters for the mass points	
  const int arraySize=200;
  assert(arraySize >= nPoints);

  Double_t a_meanCB[arraySize];  Double_t a_meanCB_err[arraySize]; 
  Double_t a_sigmaCB[arraySize]; Double_t a_sigmaCB_err[arraySize];
  Double_t a_alphaCB[arraySize]; Double_t a_alphaCB_err[arraySize];
  Double_t a_nCB[arraySize];     Double_t a_nCB_err[arraySize];
  Double_t a_Gamma[arraySize];   Double_t a_Gamma_err[arraySize];

  Double_t a_fitCovQual[arraySize];
  Double_t a_fitEDM[arraySize];
  Double_t a_fitStatus[arraySize];

  char outfile[192];
  sprintf(outfile,"compFigs%iTeV130613",sqrts);

  //Loop over the mass points
  for (int i = 0; i < nPoints; i++){
    //if(masses[i]<299||masses[i]>301)continue;          
    //Open input file with shapes and retrieve the tree
    char tmp_finalInPath[200];
    sprintf(tmp_finalInPath,"/HZZ4lTree_%s%i.root",scomp.c_str(),masses[i]);
    string finalInPath = filePath + tmp_finalInPath;
    cout<<finalInPath.c_str()<<endl;
    TFile *f = TFile::Open(finalInPath.c_str()) ;
    TTree *tree= (TTree*) f->Get("SelectedTree");
    if(tree==NULL){
      cout << "Impossible to retrieve the tree for mass point " << masses[i] <<" GeV " << endl;
      abort();
    }

    sprintf(tmp_finalInPath,"/HZZ4lTree_H%i.root",masses[i]);
    finalInPath = filePathH + tmp_finalInPath;
    cout<<finalInPath.c_str()<<endl;
    TFile *fH = TFile::Open(finalInPath.c_str()) ;
    if(!fH)continue;
    TTree *treeH= (TTree*) fH->Get("SelectedTree");
    if(treeH==NULL){
      cout << "Impossible to retrieve the tree for mass point " << masses[i] <<" GeV " << endl;
      abort();
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

    if(channel == 2) ZZMass.setBins(50);
    if(channel == 3) ZZMass.setBins(50);

    RooDataSet *set2 = new RooDataSet("data","data", tree, RooArgSet(ZZMass,MC_weight), "", "MC_weight");
    RooDataHist *set = (RooDataHist*)set2->binnedClone("datahist","datahist");

    RooRealVar ZZMassH("ZZMass","ZZMassH",low_M,high_M);
    RooRealVar MC_weightH("MC_weight","MC_weightH",0.,10.);

    if(channel == 2) ZZMassH.setBins(50);
    if(channel == 3) ZZMassH.setBins(50);

    RooDataSet *set2H = new RooDataSet("data ggH","dataH", treeH, RooArgSet(ZZMassH,MC_weightH), "", "MC_weightH");
    RooDataHist *setH = (RooDataHist*)set2H->binnedClone("datahistH","datahistH");

    double sumset = (double)set->sumEntries();
    double sumsetH = (double)setH->sumEntries();

    cout<<"setH "<<setH->sumEntries()<<" set "<<set->sumEntries()<<endl;
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
    RooRealVar meanCB("meanCB","meanCB",0.,-20.,20.);
    RooRealVar sigmaCB("sigmaCB","sigmaCB",1.,0.01,100.);
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
    RooFFTConvPdf *sigPDF;
    if(masses[i] < 399.) sigPDF = new RooFFTConvPdf("sigPDF","sigPDF",ZZMass, SignalTheorLM,massRes);
    else sigPDF = new RooFFTConvPdf("sigPDF","sigPDF",ZZMass, SignalTheor,massRes);
    sigPDF->setBufferFraction(0.2);

    double pm=0,ps=1,pa=3,pn=2,ga=1;

    if(!writeFits){
      TString finname;finname.Form("foutFit%d%d%d.root",sqrts,channel,icomp);
      TFile *finFits = TFile::Open(finname.Data());
      TF1 * fitm = (TF1*)finFits->Get("mean");
      TF1 * fits = (TF1*)finFits->Get("sigma");
      TF1 * fita = (TF1*)finFits->Get("alpha");
      TF1 * fitn = (TF1*)finFits->Get("n");
      TF1 * fitg = (TF1*)finFits->Get("gamma");

      pm=fitm->Eval(masses[i]);
      pa=fita->Eval(masses[i]);
      ps=fits->Eval(masses[i]);
      pn=fitn->Eval(masses[i]);
      ga=fitg->Eval(masses[i]);
    }
    RooRealVar meanCB_p("meanCB_p","meanCB_p",pm);
    RooRealVar sigmaCB_p("sigmaCB_p","sigmaCB_p",ps);
    RooRealVar alphaCB_p("alphaCB_p","alphaCB_p",pa);
    RooRealVar nCB_p("nCB_p","nCB_p",pn);
    RooRealVar pgamma("pgamma","pgamma",ga);

    RooCBShape paramCB("mparamCB","crystal ball param",ZZMass,meanCB_p,sigmaCB_p,alphaCB_p,nCB_p);
    RooRelBWUFParam paramBW("paramBW","paramBW",ZZMass,MHStar,pgamma);
    RooFFTConvPdf *paramPDF =  new RooFFTConvPdf("fitparamPDF","fitparamPDF",ZZMass, paramBW,paramCB);

    //Fit the shape
    RooFitResult *fitRes = sigPDF->fitTo(*set,Save(1), SumW2Error(kTRUE));

    a_fitEDM[i] = fitRes->edm();
    a_fitCovQual[i] = fitRes->covQual();
    a_fitStatus[i] = fitRes->status();

    a_meanCB[i]  = meanCB.getVal();
    a_sigmaCB[i] = sigmaCB.getVal();
    a_alphaCB[i]  = alphaCB.getVal();
    a_nCB[i]     = nCB.getVal();
    a_Gamma[i] = Gamma_TOT.getVal();

    a_meanCB_err[i]  = meanCB.getError();
    a_sigmaCB_err[i] = sigmaCB.getError();
    a_alphaCB_err[i]  = alphaCB.getError();
    a_nCB_err[i]     = nCB.getError();
    if(masses[i] > 399.) a_Gamma_err[i] = Gamma_TOT.getError();
    else a_Gamma_err[i] = 0.;

    //Plot in the figures directory
    RooPlot *xplot = ZZMass.frame();
    double ymax = xplot->GetYaxis()->GetXmax();
    double scale = sumset/sumsetH;
    set->plotOn(xplot);
    sigPDF->plotOn(xplot);
    if(!writeFits && writeParam)paramPDF->plotOn(xplot,LineColor(kRed+1),LineStyle(2));
    //setH->plotOn(xplot,MarkerStyle(24),Rescale(scale));
    //xplot->GetYaxis()->SetRange(0.,ymax);
    RooPlot *xplotH = ZZMassH.frame();
    setH->plotOn(xplotH,MarkerStyle(24),Rescale(scale));
    //sigPDF->plotOn(xplotH);

    //canv.Divide(2,1);
    canv.cd();
    xplot->Draw();
    //xplot->GetYaxis()->SetRange(0.,ymax);
    //canv.cd(2);
    xplotH->Draw("SAME");

    TLegend *leg1 = new TLegend(0.15,0.65,0.35,0.85);
    leg1->SetFillColor(0);
    leg1->SetLineColor(0);
    leg1->SetFillStyle(0);
    TLegendEntry *edata = leg1->AddEntry("data","data","lpe");
    edata->SetMarkerStyle(20);
    TLegendEntry *edataggh = leg1->AddEntry("data ggH","data ggH","lpe");
    edataggh->SetMarkerStyle(24);
    edataggh->SetFillColor(0);
    TLegendEntry *esigpdf = leg1->AddEntry("sigPDF","PDF fit","l");
    esigpdf->SetLineColor(kBlue);
    if(writeParam){
      TLegendEntry *eparpdf = leg1->AddEntry("paramPDF","Parametric PDF","l");
      eparpdf->SetLineColor(kRed+1);
      eparpdf->SetLineStyle(2);
    }
    leg1->Draw("SAME");

    string tmp_plotFileTitle;
    tmp_plotFileTitle.insert(0,outfile);
    tmp_plotFileTitle += "/fitMass_";
    char tmp2_plotFileTitle[200];
    sprintf(tmp2_plotFileTitle,"%s%i_%iTeV_",scomp.c_str(),masses[i],sqrts);
    string plotFileTitle = tmp_plotFileTitle + tmp2_plotFileTitle + schannel + ".png";

    canv.SaveAs(plotFileTitle.c_str());
  }

  if(writeFits){
    Double_t x_err[arraySize];
    
    TGraph* gr_meanCB  = new TGraph(nPoints, masses, a_meanCB);
    TGraph* gr_sigmaCB = new TGraph(nPoints, masses, a_sigmaCB);
    TGraph* gr_alphaCB = new TGraph(nPoints, masses, a_alphaCB);
    TGraph* gr_nCB     = new TGraph(nPoints, masses, a_nCB);
    TGraph* gr_Gamma   = new TGraph(nPoints, masses, a_Gamma);

    gr_meanCB->Fit("pol3");
    gr_sigmaCB->Fit("pol3");
    gr_alphaCB->Fit("pol3");
    gr_nCB->Fit("pol3");
    gr_Gamma->Fit("pol3");
    
    TF1 *fit_meanCB  = gr_meanCB->GetListOfFunctions()->First();
    TF1 *fit_sigmaCB = gr_sigmaCB->GetListOfFunctions()->First();
    TF1 *fit_alphaCB = gr_alphaCB->GetListOfFunctions()->First();
    TF1 *fit_nCB     = gr_nCB->GetListOfFunctions()->First();
    TF1 *fit_Gamma   = gr_Gamma->GetListOfFunctions()->First();

    TString foutname =   "foutFit";
    foutname.Append(Form("%d",sqrts) ) ;
    foutname.Append(Form("%d",channel))  ;
    foutname.Append(Form("%d.root",icomp))  ;

    printf("creating file %s\n",foutname.Data());
    TFile *foutFit = new TFile(foutname.Data(),"RECREATE");
    foutFit->cd();

    fit_meanCB->SetName("mean");
    fit_sigmaCB->SetName("sigma");
    fit_alphaCB->SetName("alpha");
    fit_nCB->SetName("n");
    fit_Gamma->SetName("gamma");

    fit_meanCB->Write();
    fit_sigmaCB->Write();
    fit_alphaCB->Write();
    fit_nCB->Write();
    fit_Gamma->Write();

    foutFit->Close();
  }

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

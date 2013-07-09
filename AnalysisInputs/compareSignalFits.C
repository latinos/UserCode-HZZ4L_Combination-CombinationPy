#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include "WWAnalysis/TreeModifiers/macro/SignalInterpolationStrings.h"
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
void signalFits(int channel, int sqrts, int icomp,int njets =0,TCanvas *canv=0x0, bool isLongRange = false);
float WidthValue(float mHStarWidth);
bool writeFits = false;
bool plotParam = false;
bool plotSingleFit = false;
bool plotMoriondPdf= false;
bool plotPdfFromFile = true;
bool plotOverFlowBins = false;
bool saveRoots = true;

float debugMass = 0;

void compareSignalFits()
{
  gSystem->Exec("mkdir -p compFigs7TeV130702b");//if change folders, change as well in line 121
  gSystem->Exec("mkdir -p compFigs8TeV130702b");

  gSystem->Load("../CreateDatacards/CMSSW_6_1_1/lib/slc5_amd64_gcc472/libHiggsAnalysisCombinedLimit.so");

  TCanvas *canv = new TCanvas();

  for(int ich=1;ich<4;ich++){
    for(int ien=7;ien<9;ien++){
      for(int icomp=0;icomp<6;icomp++){
	signalFits(ien,ich,icomp,0,canv);
	signalFits(ien,ich,icomp,2,canv);
	signalFits(ien,ich,icomp,-1,canv);//no jet tagging
	signalFits(ien,ich,icomp,0,canv,1);//jet tagging, long mass range
	signalFits(ien,ich,icomp,2,canv,1);//jet tagging, long mass range
	
      }
    }
  }
  //SignalFits(8,1,3,-1,canv);//for tests
}

//The actual job
void signalFits(int sqrts, int channel, int icomp, int njets, TCanvas *canv, bool isLongRange)
{

 
  if(!canv)canv = new TCanvas();
  string schannel;
  if (channel == 1) schannel = "4mu";
  if (channel == 2) schannel = "4e";
  if (channel == 3) schannel = "2e2mu";

  string scomp="H";
  if (icomp == 0) scomp = "ZH";
  if (icomp == 1) scomp = "WH";
  if (icomp == 2) scomp = "ttH";
  if (icomp == 3) scomp = "VBFH";
  if (icomp == 4) scomp = "powheg15H";
  if (icomp == 5) scomp = "powheg15jhuGenV3H";

  if(icomp==5 && sqrts ==7)return;
  cout << "Final state = " << schannel << " and sqrt(s) = " << sqrts << " "<<scomp<<endl;

  //Pick the correct mass points and paths
  TString filePath;
  int nPoints;
  double *masses;
  TString filePathH;

  if(sqrts == 7){
    filePath =filePath7TeV;
    filePathH=filePath7TeV;
    nPoints = nPoints7TeV;
    masses = mHVal7TeV;
  }
  else if(sqrts == 8){
    filePath =filePath8TeV;
    filePathH=filePath8TeV;
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
  Double_t a_alpha2CB[arraySize]; Double_t a_alpha2CB_err[arraySize];
  Double_t a_n2CB[arraySize];     Double_t a_n2CB_err[arraySize];
  Double_t a_Gamma[arraySize];   Double_t a_Gamma_err[arraySize];

  Double_t a_fitCovQual[arraySize];
  Double_t a_fitEDM[arraySize];
  Double_t a_fitStatus[arraySize];

  char outfile[192];
  sprintf(outfile,"compFigs%iTeV130702b",sqrts);

  //Loop over the mass points
  for (int i = 0; i < nPoints; i++){
    if(debugMass>100 &&(masses[i]< debugMass-0.5||masses[i]>debugMass+0.5))continue;  //For tests
    if( masses[i]>399 &&( icomp<3||plotParam)) continue;  
    if( masses[i]>200 && icomp==5) continue;  
    //if(masses[i]<125 || masses[i]>126)continue;
    //Open input file with shapes and retrieve the tree
    char tmp_finalInPath[200];
    sprintf(tmp_finalInPath,"/HZZ4lTree_%s%i.root",scomp.c_str(),masses[i]);
    string finalInPath = filePath + tmp_finalInPath;
    cout<<finalInPath.c_str()<<endl;
    TFile *f = TFile::Open(finalInPath.c_str()) ;
    if(!f)continue;
    TTree *tree= (TTree*) f->Get("SelectedTree");
    if(tree==NULL){
      cout << "Impossible to retrieve the tree for mass point " << masses[i] <<" GeV " << endl;
      continue;
    }


    //Open ggH file and retrieve the tree
    sprintf(tmp_finalInPath,"/HZZ4lTree_H%i.root",masses[i]);
    finalInPath = filePathH + tmp_finalInPath;
    cout<<finalInPath.c_str()<<endl;
    TFile *fH = TFile::Open(finalInPath.c_str()) ;
    if(!fH)continue;
    TTree *treeH= (TTree*) fH->Get("SelectedTree");
    if(treeH==NULL){
      cout << "Impossible to retrieve the tree for mass point " << masses[i] <<" GeV " << endl;
      continue;
    }


    //Mass Range
    double valueWidth = WidthValue(masses[i]);
    double windowVal = max(valueWidth,1.);

    double lowside = 100.;
    double highside = 1000.0;
        
      if (masses[i] >= 275){
	lowside = 180.0;
	highside = 650.0;
      }
    if (masses[i] >= 350){
      lowside = 200.0;
      highside = 900.0;
    }
    if (masses[i] >= 500){
      lowside = 250.0;
      highside = 1000.0;
    }
    if (masses[i] >= 700){
      lowside = 350.0;
      highside = 1400.0;
    }
    double low_M = max( (masses[i] - 20.*windowVal), lowside);
    double high_M = min((masses[i] + 15.*windowVal), highside);

    if(isLongRange){
      low_M = 70;
      high_M = 600;
    }

    //cout << "lowM = " << low_M << ", highM = " << high_M << endl;

    //Set the observable and get the RooDataSomething
    RooRealVar ZZMass("ZZMass","ZZMass",low_M,high_M);
    RooRealVar MC_weight("MC_weight","MC_weight",0.,10.);
    RooRealVar NJets30("NJets30","NJets30",0,100);

    ZZMass.setBins(100);
    if(channel == 2) ZZMass.setBins(60);
    if(channel == 3) ZZMass.setBins(60);
 
    if(icomp==0)ZZMass.setBins(100);
    if(icomp==0 && channel == 2)ZZMass.setBins(60);
    if(icomp==0 && channel == 3)ZZMass.setBins(60);

    if(isLongRange) ZZMass.setBins(200);

    RooDataSet *set2;
    RooArgSet ntupleVarSet(ZZMass,NJets30,MC_weight);
    if(njets==-1) set2 = new RooDataSet("data","data", tree, RooArgSet(ZZMass,MC_weight), "", "MC_weight");
    else {  
      set2 = new RooDataSet("data","data",ntupleVarSet,WeightVar("MC_weight"));
      
      Float_t myMC,myMass;
      Short_t myNJets;
      int nentries = tree->GetEntries();
      
      tree->SetBranchAddress("ZZMass",&myMass);
      tree->SetBranchAddress("MC_weight",&myMC);
      tree->SetBranchAddress("NJets30",&myNJets);
      
      for(int itr =0;itr<nentries;itr++) {
	tree->GetEntry(itr);
	if(njets==0 && myNJets>1) continue;
	if(njets==2 && myNJets<2) continue;
	if(!plotOverFlowBins && (myMass>high_M || myMass<low_M))continue;
	ntupleVarSet.setRealValue("ZZMass",myMass);
	ntupleVarSet.setRealValue("MC_weight",myMC);
	ntupleVarSet.setRealValue("NJets30",(double)myNJets);
	
	set2->add(ntupleVarSet, myMC);
      }
    }

    RooDataHist *set = (RooDataHist*)set2->binnedClone("datahist","datahist");
    
    RooRealVar ZZMassH("ZZMass","ZZMassH",low_M,high_M);
    RooRealVar MC_weightH("MC_weight","MC_weightH",0.,10.);
    
    ZZMassH.setBins(100);
    if(channel == 2) ZZMassH.setBins(60);
    if(channel == 3) ZZMassH.setBins(60);
  
    if(icomp==0)ZZMassH.setBins(100);
    if(icomp==0 && channel == 2)ZZMassH.setBins(60);
    if(icomp==0 && channel == 3)ZZMassH.setBins(60);
    
    if(isLongRange) ZZMassH.setBins(200);
    
    RooDataSet *set2H = new RooDataSet("data ggH","dataH", treeH, RooArgSet(ZZMassH,MC_weightH), "", "MC_weight");
    RooDataHist *setH = (RooDataHist*)set2H->binnedClone("datahistH","datahistH");
    
    double sumset = (double)set->sumEntries();
    double sumsetH = (double)setH->sumEntries();
    
    RooDataSet setM;
    RooRealVar *CMS_zz4l_mass;
    double sumsetM =1;
    
    //cout<<"setH "<<setH->sumEntries()<<" set "<<set->sumEntries()<<endl;
    
    //Theoretical signal model  
    RooRealVar MHStar("MHStar","MHStar",masses[i],0.,2000.);
    MHStar.setConstant(true);
    RooRealVar Gamma_TOT("Gamma_TOT","Gamma_TOT",valueWidth,0.,700.);
    if(masses[i] < 399.) Gamma_TOT.setConstant(true);
    RooRealVar one("one","one",1.0);
    //one.setConstant(kTRUE);
    
    RooGenericPdf SignalTheor("SignalTheor","(@0)/( pow(pow(@0,2)-pow(@1,2),2) + pow(@0,2)*pow(@2,2) )",RooArgSet(ZZMass,MHStar,Gamma_TOT));
    RooRelBWUFParam SignalTheorLM("signalTheorLM","signalTheorLM",ZZMass,MHStar,one);
    
    //Experimental resolution
    RooRealVar meanCB("meanCB","meanCB",0.,-20.,20.);
    RooRealVar sigmaCB("sigmaCB","sigmaCB",1.,0.01,100.);
    RooRealVar alphaCB("alphaCB","alphaCB",3.,-10.,10.);
    RooRealVar nCB("nCB","nCB",2.,-10.,10.);
    RooRealVar alpha2CB("alpha2CB","alpha2CB",3.,-10.,10.);
    RooRealVar n2CB("n2CB","n2CB",2.,-10.,10.);
  
    //Initialize to decent values
    float m = masses[i];
    if(channel == 1) sigmaCB.setVal(11.282-0.213437*m+0.0015906*m*m-5.18846e-06*m*m*m+8.05552e-09*m*m*m*m -4.69101e-12*m*m*m*m*m);
    else if(channel == 2) sigmaCB.setVal(3.58777+-0.0252106*m+0.000288074*m*m+-8.11435e-07*m*m*m+7.9996e-10*m*m*m*m);
    else if(channel == 3) sigmaCB.setVal(7.42629+-0.100902*m+0.000660553*m*m+-1.52583e-06*m*m*m+1.2399e-09*m*m*m*m);
    else abort();

    //RooCBShape massRes("massRes","crystal ball",ZZMass,meanCB,sigmaCB,alphaCB,nCB);
    RooDoubleCB massRes("massRes","crystal ball",ZZMass,meanCB,sigmaCB,alphaCB,nCB,alpha2CB,n2CB);

    //Convolute theoretical shape and resolution
    RooFFTConvPdf *sigPDF;
    if(masses[i] < 399.) sigPDF = new RooFFTConvPdf("sigPDF","sigPDF",ZZMass, SignalTheorLM,massRes);
    else sigPDF = new RooFFTConvPdf("sigPDF","sigPDF",ZZMass, SignalTheor,massRes);
    sigPDF->setBufferFraction(0.2);

    //take the pdf from the parametrization
    RooFFTConvPdf *paramPDF;
    double pm=0,ps=1,pa=3,pn=2,ga=1,pa2=3,pn2=2;
    if(!writeFits && plotParam){
      TString finname;finname.Form("foutFit%d%d%d.root",sqrts,channel,icomp);
      TFile *finFits = TFile::Open(finname.Data());
      if(finFits){
	TF1 * fitm = (TF1*)finFits->Get("mean");
	TF1 * fits = (TF1*)finFits->Get("sigma");
	TF1 * fita = (TF1*)finFits->Get("alpha");
	TF1 * fitn = (TF1*)finFits->Get("n");
	TF1 * fita2 = (TF1*)finFits->Get("alpha2");
	TF1 * fitn2 = (TF1*)finFits->Get("n2");
	TF1 * fitg = (TF1*)finFits->Get("gamma");
	
	pm=fitm->Eval(masses[i]);
	pa=fita->Eval(masses[i]);
	pa2=fita2->Eval(masses[i]);
	ps=fits->Eval(masses[i]);
	pn=fitn->Eval(masses[i]);
	pn2=fitn2->Eval(masses[i]);
	ga=fitg->Eval(masses[i]);
      }
    

      RooRealVar meanCB_p("meanCB_p","meanCB_p",pm);
      RooRealVar sigmaCB_p("sigmaCB_p","sigmaCB_p",ps);
      RooRealVar alphaCB_p("alphaCB_p","alphaCB_p",pa);
      RooRealVar nCB_p("nCB_p","nCB_p",pn);
      RooRealVar alpha2CB_p("alpha2CB_p","alpha2CB_p",pa2);
      RooRealVar n2CB_p("n2CB_p","n2CB_p",pn2);
      RooRealVar pgamma("pgamma","pgamma",ga);

      //RooCBShape paramCB("mparamCB","crystal ball param",ZZMass,meanCB_p,sigmaCB_p,alphaCB_p,nCB_p);
      RooDoubleCB paramCB("paramCB","crystal ball param",ZZMass,meanCB_p,sigmaCB_p,alphaCB_p,nCB_p,alpha2CB_p,n2CB_p);
      RooRelBWUFParam paramBW("paramBW","paramBW",ZZMass,MHStar,pgamma);
      RooFFTConvPdf *paramPDF =  new RooFFTConvPdf("fitparamPDF","fitparamPDF",ZZMass, paramBW,paramCB);
    }
    //take the pdf from the (Moriond) datacards
    //STILL BUG FIXING THIS PART!!!!!!!!!!!!!!!
    RooDoubleCB MoriondCB;
    if(plotMoriondPdf){
      int mass=masses[i];
      TString cardName;
      cardName.Form("/afs/cern.ch/work/g/gortona/MoriondWS/hzz4l/%i/hzz4l_%sS_%iTeV_0.input.root",mass,schannel.c_str(),sqrts);
      TFile *fff = TFile::Open(cardName.Data());
      if(fff){
	RooWorkspace w = (RooWorkspace)fff->Get("w");
	//w.Print();
	//abort();
	MoriondCB = (RooDoubleCB)w.pdf("signalCB_ggH_SM");
	CMS_zz4l_mass = (RooRealVar*)w.var("CMS_zz4l_mass");
	//setM = (RooDataSet)w.data("CMS_zz4l_mass");
	//abort();
      }
    }

    //take the pdf from committed RooFormulaVar (parameters from Emanuele are here: UserCode/Mangano/WWAnalysis/TreeModifiers/macro/SignalInterpolationStrings.h
    RooRealVar zero("zero","zero",0.0);
    RooRealVar HMass("HMass","HMass",masses[i]);

    bool ien=1;
    if(sqrts==8)ien=0;
    bool doFFT =1;//default low mas =1

    TString sigstring = getSignalCBMeanString(masses[i],channel-1,ien,doFFT);
    RooFormulaVar meanCB_f("meanCB_f"    ,sigstring.Data(),RooArgList(HMass,zero));
    sigstring = getSignalCBSigmaString(masses[i],channel-1,ien);
    RooFormulaVar sigmaCB_f("sigmaCB_f"  ,sigstring.Data(),RooArgList(HMass,zero));
    sigstring = getSignalCBAlphaLString(masses[i],channel-1,ien);
    RooFormulaVar alphaCB_f("alphaCB_f"  ,sigstring.Data(),RooArgList(HMass));
    sigstring = getSignalCBNLString(masses[i],channel-1,ien);
    RooFormulaVar nCB_f("nCB_f"          ,sigstring.Data(),RooArgList(HMass));
    sigstring = getSignalCBAlphaRString(masses[i],channel-1,ien);
    RooFormulaVar alpha2CB_f("alpha2CB_f",sigstring.Data(),RooArgList(HMass));
    sigstring = getSignalCBNRString(masses[i],channel-1,ien);
    RooFormulaVar n2CB_f("n2CB_f"        ,sigstring.Data(),RooArgList(HMass));
    sigstring = getSignalBWGammaString(masses[i],channel-1,ien);
    RooFormulaVar fgamma("fgamma"        ,sigstring.Data(),RooArgList(HMass,zero));

    RooDoubleCB formulaCB("formulaCB","crystal ball from Formulas",ZZMass,meanCB_f,sigmaCB_f,alphaCB_f,nCB_f,alpha2CB_f,n2CB_f);
    RooFFTConvPdf *formulaPDF;// = new RooFFTConvPdf("fitFormulaPDF","fitFormulaPDF",ZZMass, formulaBW,formulaCB) ;
     
    if(masses[i]>399){
      RooRelBWHighMass formulaBW("formulaBW","formulaBW",ZZMass,HMass,fgamma);
      formulaPDF = new RooFFTConvPdf("fitFormulaPDF","fitFormulaPDF",ZZMass, formulaBW,formulaCB);
      //fbw = (RooRelBWUFParam)formulaBW;
    }else{
      RooRelBWUFParam formulaBWUF("formulaBW","formulaBW",ZZMass,HMass,one);
      formulaPDF = new RooFFTConvPdf("fitFormulaPDF","fitFormulaPDF",ZZMass, formulaBWUF,formulaCB);
      //fbw = (RooRelBWUFParam)formulaBW;
    }
    
    //    RooFFTConvPdf *formulaPDF = new RooFFTConvPdf("fitFormulaPDF","fitFormulaPDF",ZZMass, formulaBW,formulaCB);

    //Fit the shape
    RooFitResult *fitRes = sigPDF->fitTo(*set,Save(1), SumW2Error(kTRUE));

    a_fitEDM[i] = fitRes->edm();
    a_fitCovQual[i] = fitRes->covQual();
    a_fitStatus[i] = fitRes->status();

    a_meanCB[i]  = meanCB.getVal();
    a_sigmaCB[i] = sigmaCB.getVal();
    a_alphaCB[i]  = alphaCB.getVal();
    a_nCB[i]     = nCB.getVal();
    a_alpha2CB[i]  = alpha2CB.getVal();
    a_n2CB[i]     = n2CB.getVal();
    a_Gamma[i] = Gamma_TOT.getVal();

    a_meanCB_err[i]  = meanCB.getError();
    a_sigmaCB_err[i] = sigmaCB.getError();
    a_alphaCB_err[i]  = alphaCB.getError();
    a_nCB_err[i]     = nCB.getError();
    a_alpha2CB_err[i]  = alpha2CB.getError();
    a_n2CB_err[i]     = n2CB.getError();
    if(masses[i] > 399.) a_Gamma_err[i] = Gamma_TOT.getError();
    else a_Gamma_err[i] = 0.;

    //Plot in the figures directory
    RooPlot *xplot = ZZMass.frame();
    double scale = sumset/sumsetH;
    set->plotOn(xplot);
    if(plotSingleFit) sigPDF->plotOn(xplot);
    if(!writeFits && plotParam) paramPDF->plotOn(xplot,LineColor(kRed+1),LineStyle(2));
    if(plotPdfFromFile) {
      //formulaCB.plotOn(xplot,LineColor(kYellow+2),LineStyle(2));
      formulaPDF->plotOn(xplot,LineColor(kGreen+2),LineStyle(2));
    }
    RooPlot *xplotH = ZZMassH.frame();
    setH->plotOn(xplotH,MarkerStyle(24),Rescale(scale));

    canv->cd();
    xplot->Draw();
    xplotH->Draw("SAME");

    if(plotMoriondPdf){
      cout<<"fino a qui 0"<<endl;
      RooPlot *mplot = CMS_zz4l_mass->frame();
      double mscale = sumset/sumsetM;
      cout<<"fino a qui"<<endl;
      MoriondCB.plotOn(mplot,LineColor(kOrange+7),LineStyle(4)/*,Rescale(mscale)*/);
      mplot->Draw("SAME");
    }
  

    TLegend *leg1 = new TLegend(0.15,0.65,0.35,0.85);
    leg1->SetFillColor(0);
    leg1->SetLineColor(0);
    leg1->SetFillStyle(0);
    TString dataname = scomp.c_str();
    TString jetname;
    if(njets==2)jetname.Form(" %i+ jets",njets);
    if(njets==0)jetname.Form(" 0/1 jets");
    dataname.Append(jetname.Data());
    TLegendEntry *edata = leg1->AddEntry(dataname.Data(),dataname.Data(),"lpe");
    edata->SetMarkerStyle(20);
    TLegendEntry *edataggh = leg1->AddEntry("ggH","ggH","lpe");
    edataggh->SetMarkerStyle(24);
    edataggh->SetFillColor(0);
    if(plotSingleFit){
      TLegendEntry *esigpdf = leg1->AddEntry("sigPDF","PDF fit","l");
      esigpdf->SetLineColor(kBlue);
    }
    if(plotParam){
      TLegendEntry *eparpdf = leg1->AddEntry("paramPDF","Parametric PDF","l");
      eparpdf->SetLineColor(kRed+1);
      eparpdf->SetLineStyle(2);
    }
    if(plotMoriondPdf){
      TLegendEntry *emorpdf = leg1->AddEntry("Moriond PDF","Moriond PDF","l");
      emorpdf->SetLineColor(kGreen+2);
      emorpdf->SetLineStyle(4);
    }
    if(plotPdfFromFile){
      TLegendEntry *eforpdf = leg1->AddEntry("Emanuele Shape","Emanuele Shape","l");
      eforpdf->SetLineColor(kGreen+2);
      eforpdf->SetLineStyle(2);
    }
    leg1->Draw("SAME");

    string tmp_plotFileTitle;
    tmp_plotFileTitle.insert(0,outfile);
    tmp_plotFileTitle += "/fitMass_";
    char tmp2_plotFileTitle[200];
    sprintf(tmp2_plotFileTitle,"%s%i_%iTeV_",scomp.c_str(),masses[i],sqrts);
    string njetstr = "";
    if(njets==0)njetstr+="0jets_";
    if(njets==2)njetstr+="2plusjets_";
    if(isLongRange)njetstr+="longRange_";
    string plotFileTitle = tmp_plotFileTitle + tmp2_plotFileTitle + njetstr + schannel + ".png";
    canv->SaveAs(plotFileTitle.c_str());
    if(saveRoots){
      plotFileTitle = tmp_plotFileTitle + tmp2_plotFileTitle + njetstr + schannel + ".root";
      canv->SaveAs(plotFileTitle.c_str());
    }
    //delete newTree;
    //fo->Close();
    //delete fo;
    delete sigPDF;
    delete paramPDF;
    delete formulaPDF;
    delete set2;
    delete set2H;
    delete leg1;
  }

  if(writeFits){
    Double_t x_err[arraySize];
    
    TGraph* gr_meanCB  = new TGraph(nPoints, masses, a_meanCB);
    TGraph* gr_sigmaCB = new TGraph(nPoints, masses, a_sigmaCB);
    TGraph* gr_alphaCB = new TGraph(nPoints, masses, a_alphaCB);
    TGraph* gr_nCB     = new TGraph(nPoints, masses, a_nCB);
    TGraph* gr_alpha2CB = new TGraph(nPoints, masses, a_alpha2CB);
    TGraph* gr_n2CB     = new TGraph(nPoints, masses, a_n2CB);
    TGraph* gr_Gamma   = new TGraph(nPoints, masses, a_Gamma);

    gr_meanCB->Fit("pol3");
    gr_sigmaCB->Fit("pol3");
    gr_alphaCB->Fit("pol3");
    gr_nCB->Fit("pol3");
    gr_alpha2CB->Fit("pol3");
    gr_n2CB->Fit("pol3");
    gr_Gamma->Fit("pol3");
    
    TF1 *fit_meanCB  = gr_meanCB->GetListOfFunctions()->First();
    TF1 *fit_sigmaCB = gr_sigmaCB->GetListOfFunctions()->First();
    TF1 *fit_alphaCB = gr_alphaCB->GetListOfFunctions()->First();
    TF1 *fit_nCB     = gr_nCB->GetListOfFunctions()->First();
    TF1 *fit_alpha2CB = gr_alpha2CB->GetListOfFunctions()->First();
    TF1 *fit_n2CB     = gr_n2CB->GetListOfFunctions()->First();
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
    fit_alpha2CB->SetName("alpha2");
    fit_n2CB->SetName("n2");
    fit_Gamma->SetName("gamma");

    fit_meanCB->Write();
    fit_sigmaCB->Write();
    fit_alphaCB->Write();
    fit_nCB->Write();
    fit_alpha2CB->Write();
    fit_n2CB->Write();
    fit_Gamma->Write();

    foutFit->Close();
  }
  printf("%s\n",sigstring.Data());
  printf("%f\n",fgamma.getVal());
  printf("%f\n",HMass.getVal());
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

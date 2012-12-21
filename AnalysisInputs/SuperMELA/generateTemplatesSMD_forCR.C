/* 
 * Create 2D (mass, LD) templates. Script imported from: http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/JHU/MELA/scripts/generateTemplates.C?revision=1.1.2.1&view=markup&pathrev=post_unblinding
 * Requires ZZMatrixElement/MELA to have been checked out and compiled.
 * usage: 
 * -set input paths variables in Config.h
 * -run with:
 * root -q -b ../loadMELA.C generateTemplatesSMD_forCR.C+
 * 2D templates are written to "destDir"
 *
 */

#include "TFile.h"
#include "TStyle.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TString.h"
#include <sstream>
#include <Riostream.h>
#include <vector>
#include <utility>

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TROOT.h>
#include "ZZMatrixElement/MELA/interface/Mela.h"
#endif


//----------> SET INPUT VARIABLES in Config.h
#include "ConfigSMD.h"
//<----------

//--- Flags to control re-computation of KD
bool recompute_ = false;       // Recompute KD (instead of taking it from the tree); when true, the following flags apply:
bool usePowhegTemplate=false;  // false use analytic bg
bool withPt_ = false;          // Include pT in KD
bool withY_  = false;          //    "    Y  "  "
int sqrts    = 7;              // sqrts, used only for withPt_/withY_
Mela* myMELA; //used if recompute is true

////////////////////////////////////
//--- Really important params --- //
const int mH=125;
const float mzzCutLow=105;
const float mzzCutHigh=140;
const int useSqrts=2;              //0=use 7+8TeV; 1=use 7TeV only, 2 use 8TeV only
TString melaName = "pseudoLD"; // name of KD branch to be used.
const TString destDir = "../../CreateDatacards/templates2D_smd_8TeV_20121105_NOIntRew/"; //it must already exist !
bool makePSTemplate = true;
bool makeAltSignal = true;
const float melaCut=-1.0; //if negative, it is deactivated
const bool applyInterferenceRew=false;
string fInterferenceName="./1DinterferenceReweight.root";
TString melaCutName = melaName; //distinction btw LD and mela needed because 
                               //we might want both psMELA (for 2D template) and MELA (for cut)
//-----


bool extendToHighMass = false; // Include signal samples above 600 GeV
float highMzz=(extendToHighMass?1000:800);
float mBinSize=2.;



//-- binning of 2D template
const int nbinsX=21;
float binsX[nbinsX+1]={0.000, 0.030, 0.060, 0.100, 0.200, 0.300, 0.400, 0.500, 0.550, 0.600, 
		       0.633, 0.666, 0.700, 0.733, 0.766, 0.800, 0.833, 0.866, 0.900, 0.933,
		       0.966, 1.000};
const int nbinsYps=25;
float binsYps[nbinsYps+1]={0.000, 0.100, 0.150, 0.200, 0.233, 0.266, 0.300, 0.333, 0.366, 0.400, 
		       0.433, 0.466, 0.500, 0.533, 0.566, 0.600, 0.633, 0.666, 0.700, 0.733, 
		       0.766, 0.800, 0.850, 0.900, 0.950, 1.000};

const int nbinsYgrav=29;
float binsYgrav[nbinsYgrav+1]={0.000, 0.100, 0.150, 0.175 , 0.200, 0.225, 0.250, 0.275, 0.300, 0.325, 
		       0.350, 0.375, 0.400, 0.425 , 0.450, 0.475, 0.500, 0.525, 0.575, 0.600, 
		       0.633, 0.666, 0.700, 0.733 , 0.766, 0.800, 0.850, 0.900, 0.950, 1.000};

/*
const int nbinsY=31;
float binsY[nbinsY+1]={0.000, 0.100, 0.150, 0.200, 0.233, 0.266, 0.300, 0.320, 0.340, 0.360,
		       0.380, 0.400, 0.420, 0.440, 0.460, 0.480, 0.500, 0.533, 0.566, 0.600,
		       0.633, 0.666, 0.700, 0.733, 0.766, 0.800, 0.833, 0.866, 0.900, 0.933,
		       0.966, 1.000};
*/
//--
pair<TH2F*,TH2F*> reweightForCRunc(TH2F* temp);
void buildChainSingleMass(TChain* bkgMC, TString channel, int sampleIndex=0, int mh=125) ;
TH2F* reweightForInterference(TH2F* temp);
double calcInterfRew(TH1 *h,double KD );
void makePlot1D( TH1 *h ,TString label );
void makePlot2D( TH2 *h ,TString label );
TH2F* fillTemplate(TString channel="4mu", int sampleIndex=0,TString superMelaName="superLD",TString templateName="bkgHisto",  bool smooth=false);
TH1F *fillKDhisto(TString channel="4mu", int sampleIndex=0,float mzzLow=0.0,float mzzHigh=99999.0,  bool smooth=false);
void makeTemplate(TString channel="4mu");
void storeLDDistribution();
void generateTemplatesSMD_forCR() ;

//Following are NOT USED:
void buildChain(TChain* bkgMC, TString channel, int sampleIndex=0);
TH2F* mergeTemplates(TH2F* lowTemp, TH2F* highTemp);

//=======================================================================

pair<TH2F*,TH2F*> reweightForCRunc(TH2F* temp){

  cout << "reweightForCRunc" << endl;

  TH2F* tempUp = new TH2F(*temp);
  TH2F* tempDn = new TH2F(*temp);

  pair<TH2F*,TH2F*> histoPair(0,0);

  // ---------------------
  // functions for scaling
  // ---------------------
  
  double oldTempValue=0;
  double newTempValue=0;
  int point=-1;


  //  int numPtmp=0;
  const int numPtmpPS=5;
  const int numPtmp=8;

  const int numPoints=makePSTemplate? numPtmpPS : numPtmp;
  double low[numPoints] ;
  double high[numPoints] ;

  // ================ binning for pseudoMELA ==============================
  double lowBinsPS[numPtmpPS]   ={100.,        120.,        140.,         160.,     180.  };
  double highBinsPS[numPtmpPS]  ={120.,        140.,        160.,         180.,     1002. };
  // =======================================================================

  // ================ binning for MELA ==============================
  double lowBins[numPtmp]   ={100.,        120.,        140.,         160.,     180.,     220.,     260.,     300. }; 
  double highBins[numPtmp]  ={120.,        140.,        160.,         180.,     220.,     260.,     300.,     1002.};
  // ======================================================================


  // ================ systematics for pseudoMELA ==========================
  double slopePS_syst[numPtmpPS] ={-3.32705e-01, -1.90814e-01, -9.77189e-01, -3.81680e-01, 0.0 };
  double yIntrPS_syst[numPtmpPS] ={ 9.05727e-01, 9.95995e-01,  1.40367e+00,  1.12690,      1.0 }; 
  //  ==================================================================

  // ================ systematics for MELA ==========================
  double slope_syst[numPtmp] ={4.71836e-01, 1.17671e-01, -3.81680e-01, -1.20481, -1.21944, -2.06928, -1.35337, 0.0 };
  double yIntr_syst[numPtmp] ={6.83860e-01, 9.38454e-01, 1.12690,      1.24502,  1.72764,  2.11050,  1.52771,  1.0 }; 
  //==================================================================

  double slope[numPoints], yIntr[numPoints];

  if(makePSTemplate){
    for(int ib=0;ib<numPoints;ib++){
      low[ib]=lowBinsPS[ib];
      high[ib]=highBinsPS[ib];
      slope[ib]=slopePS_syst[ib];
      yIntr[ib]=yIntrPS_syst[ib];
    }
  }
  else{
    for(int ib=0;ib<numPoints;ib++){
      low[ib]=lowBins[ib];
      high[ib]=highBins[ib];
      slope[ib]=slope_syst[ib];
      yIntr[ib]=yIntr_syst[ib];
    }
  }


  for(int i=1; i<=temp->GetNbinsX(); i++){
    point = -1;

    // choose correct scale factor
    for(int p=0; p<numPoints; p++){
      //float m=(i*2.+101.); // NA: This is the center of bin i+1 and not of bin i... why?
      float m=temp->GetBinCenter(i+1); 
      if( m>=low[p] && m<high[p] ){
	point = p;
      }
    }
    if(point == -1){
      cout << "ERROR: could not find correct scale factor"<< endl;
      return histoPair;
    }

    for(int j=1; j<=temp->GetNbinsY(); j++){

      oldTempValue = temp->GetBinContent(i,j);
      newTempValue = oldTempValue*(slope[point]*(double)j/30.+yIntr[point]);
      tempUp->SetBinContent(i,j,newTempValue);
      newTempValue = oldTempValue*(-slope[point]*(double)j/30.+2.-yIntr[point]);
      tempDn->SetBinContent(i,j,newTempValue);

    }// end loop over Y bins

    // -------------- normalize mZZ slice ----------------

    double norm_up=(tempUp->ProjectionY("temp",i,i))->Integral();
    double norm_dn=(tempDn->ProjectionY("temp",i,i))->Integral();


    for(int j=1; j<=temp->GetNbinsY(); j++){
      
      tempUp->SetBinContent(i,j,tempUp->GetBinContent(i,j)/norm_up);
      tempDn->SetBinContent(i,j,tempDn->GetBinContent(i,j)/norm_dn);

    }

    // ---------------------------------------------------

  }// end loop over X bins

  histoPair.first  = tempUp;
  histoPair.second = tempDn;

  return histoPair;

}

//=======================================================================
TH2F* reweightForInterference(TH2F* temp){

  cout << "reweightForInterference" << endl;

  // for interference reweighting of MELA
  TF1* reweightFunc =0;

  if(makePSTemplate){// or makeAltSignal
  // ===================================================
  // for interference reweighting of pseudo-MELA
    reweightFunc = new TF1("reweightFunc","([0]+[1]*(x-110) )*0.5*(1 + TMath::Erf([2]*(x -[3]) ))*0.5*(1 + TMath::Erf([4]*([5]-x) ))  ",100,200);
    
    reweightFunc->SetParameter(0,-5.66409e-01);
    reweightFunc->SetParameter(1, 1.22591e-02);
    reweightFunc->SetParameter(2, 1.64942e+00);
    reweightFunc->SetParameter(3, 1.10080e+02);
    reweightFunc->SetParameter(4, 2.10905e+00);
    reweightFunc->SetParameter(5, 1.78529e+02);
    //  ==================================================== 
  }
  else{
    reweightFunc =new TF1("reweightFunc","gaus",100,1000);
    
    reweightFunc->SetParameter(0,0.354258);
    reweightFunc->SetParameter(1,114.909);
    reweightFunc->SetParameter(2,17.1512);
  }

  TH2F* newTemp = new TH2F(*temp);
  
  // ---------------------
  // functions for scaling
  // ---------------------
  
  double oldTempValue=0;
  double newTempValue=0;

  double slope;

  for(int i=1; i<=temp->GetNbinsX(); i++){

    // choose correct scale factor

    // for reweighting MELA
    if(makePSTemplate || i>8){
      slope=reweightFunc->Eval((double)((i-1)*2+101));
    }
    else      slope=.354;
     

    /* ==============================================
    // for reweighting pseudo-MELA
    slope = reweightFunc->Eval((double)((i-1)*2+101));
    ============================================== */

    for(int j=1; j<=temp->GetNbinsY(); j++){
      
      oldTempValue = temp->GetBinContent(i,j);
      newTempValue = oldTempValue*(1+slope*((double)j/30.-.5));
      newTemp->SetBinContent(i,j,newTempValue);

    }// end loop over Y bins

    // -------------- normalize mZZ slice ----------------

    double norm=(newTemp->ProjectionY("temp",i,i))->Integral();

    for(int j=1; j<=temp->GetNbinsY(); j++){
      
      newTemp->SetBinContent(i,j,newTemp->GetBinContent(i,j)/norm);

    }

    // ---------------------------------------------------

  }// end loop over X bins

  return newTemp;

}

double calcInterfRew(TH1 *h,double KD ){
  return h->GetBinContent(h->FindBin(KD));
}

void makePlot1D( TH1 *h ,TString label ){

  gStyle->SetOptStat(1);
  
  TCanvas *c1D=new TCanvas("c1d",("CANVAS "+label).Data());
  c1D->cd();
  h->SetXTitle("superMELA");
  h->SetYTitle("Norm. to unity");
  h->GetXaxis()->SetLabelSize(0.035);
  h->GetYaxis()->SetLabelSize(0.035);
  h->GetYaxis()->SetTitleOffset(1.15);
  h->SetLineWidth(2);
  h->Draw("HIST");
  c1D->SaveAs(("can_template1D_SuperMELA_"+label+".root").Data());
  delete c1D;
}


void makePlot2D( TH2 *h ,TString label ){

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  TCanvas *c2D=new TCanvas("c2d",("CANVAS "+label).Data());
  c2D->cd();
  h->SetXTitle("superMELA");
  h->SetYTitle("pseudoMELA");
  h->GetXaxis()->SetLabelSize(0.035);
  h->GetYaxis()->SetLabelSize(0.035);
  h->GetYaxis()->SetTitleOffset(1.15);
  h->Draw("col");
  c2D->SaveAs(("can_template_SMDvsPSD_"+label+".root").Data());
  delete c2D;
}

//=======================================================================

void buildChainSingleMass(TChain* bkgMC, TString channel, int sampleIndex, int mh) {

  //  TString sample[4]={"H*","ZZTo*","ggZZ*","H*Pse"};
  //  TString sampleName[4]={"signal","qqZZ","ggZZ","signal_PS"};

  if(useSqrts!=1 &&useSqrts!=2){
    cout<<"Error ! cannot build templates for superMELA mixing 7 and 8 TeV samples. Sqrt(s) set is "<<useSqrts<<endl;
    return;
  }

  TString chPath =channel;// (channel=="2e2mu"?"2mu2e":channel); // Adapt to different naming convention...

  char mch[32];
  sprintf(mch,"%d",int(mh));
  string strM=mch;
  //An error is issued on missing files; if a single file is missing in one set it can be safely ignored.

  string suffix="_withSMD_doubleCBonly.root";

  if(sampleIndex==0){
    //7TeV
    if(useSqrts==1)bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H"+strM+suffix);
    else     bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H"+strM+suffix);
    
  } 
  else if (sampleIndex==1){
     if(useSqrts==1){
      cout<<"Readign in 7 TeV for bkgd"<<endl;
    //7TeV
    bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ZZTo2e2mu"+suffix);
    bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ZZTo2e2tau"+suffix);
    bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ZZTo2mu2tau"+suffix);
    bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ZZTo4e"+suffix);
    bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ZZTo4mu"+suffix);
    bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ZZTo4tau"+suffix);
    }
    else{
      //8TeV
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_ZZTo2e2mu"+suffix);
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_ZZTo2e2tau"+suffix);
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_ZZTo2mu2tau"+suffix);
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_ZZTo4e"+suffix);
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_ZZTo4mu"+suffix);
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_ZZTo4tau"+suffix);
    }

  } 
  else if (sampleIndex==2){
     if(useSqrts==1){
    //7TeV
    bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ggZZ2l2l"+suffix);
    bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ggZZ4l"+suffix);
    }
    else{    //8TeV
    bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_ggZZ2l2l"+suffix);
    bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_ggZZ4l"+suffix);
    }
    
  } 
  else if(sampleIndex==3){ //this is for alternative signal samples

     if(useSqrts==1){   //7TeV
      cout<<"Readign in 7 TeV for Alt signal"<<endl;
      bkgMC->Add(filePath7TeVPS + "/" + chPath +"/HZZ4lTree_jhuPseH"+strM+""+suffix);
    }
    else{   //8TeV
      bkgMC->Add(filePath8TeVPS + "/" + chPath +"/HZZ4lTree_jhuPseH"+strM+""+suffix);
    }
  }
  else if(sampleIndex==4){ //this is for another alternative signal samples

     if(useSqrts==1){   //7TeV
      cout<<"Readign in 7 TeV for Alt signal (2)"<<endl;
      bkgMC->Add(filePath7TeVPS + "/" + chPath +"/HZZ4lTree_graviH"+strM+""+suffix);
    }
    else{   //8TeV
      bkgMC->Add(filePath8TeVPS + "/" + chPath +"/HZZ4lTree_graviH"+strM+""+suffix);
    }
  }
  else if (sampleIndex==5){
     if(useSqrts==1){
      cout<<"Readign in 7 TeV for Z+X data control regions"<<endl;
    //7TeV
    bkgMC->Add(filePath7TeV + "/CR/HZZ4lTree_DoubleEle"+suffix);
    bkgMC->Add(filePath7TeV + "/CR/HZZ4lTree_DoubleMu"+suffix);
    bkgMC->Add(filePath7TeV + "/CR/HZZ4lTree_DoubleOr"+suffix);
    }
    else{
cout<<"Readign in 8 TeV for Z+X data control regions"<<endl;
      //8TeV
    bkgMC->Add(filePath8TeV + "/CR/HZZ4lTree_DoubleEle"+suffix);
    bkgMC->Add(filePath8TeV + "/CR/HZZ4lTree_DoubleMu"+suffix);
    bkgMC->Add(filePath8TeV + "/CR/HZZ4lTree_DoubleOr"+suffix);
    cout<<"Finished to add files to chain: #entries is "<<bkgMC->GetEntries()<<endl;
    }

  }//end if sampleindex==5
    
}


//=======================================================================

TH2F* fillTemplate(TString channel, int sampleIndex,TString superMelaName,TString templateName,  bool smooth){

  TChain* bkgMC = new TChain("SelectedTree");//"data_obs" "SelectedTree"


  //  buildChain(bkgMC, channel, sampleIndex);
  bool use8TeV=false;
  buildChainSingleMass(bkgMC, channel, sampleIndex,mH);
  cout << "Chain for " << channel << " " << sampleIndex << "  "<<(useSqrts==2 ? "8TeV" : "7TeV") << " ===>" << bkgMC->GetEntries() << endl;
  // bkgMC->ls();

  float mzz,KD,KD_cut,w=0;
  double sKD;
  float m1=0, m2=0, costheta1=0, costheta2=0, costhetastar=0, phi=0, phi1=0;
  float pt4l=0, Y4l=0;
  float psig=0, pbkg=0;
  
  
  
  bool cutSameVar=false;
  bkgMC->SetBranchAddress(melaName.Data(),&KD);
  if (melaName!=melaCutName) {
    bkgMC->SetBranchAddress(melaCutName.Data(),&KD_cut);
  }
 else {
   //    KD_cut=KD;
   cutSameVar=true;
   // bkgMC->SetBranchAddress(melaName.Data(),&KD_cut);
  }
  

  bkgMC->SetBranchAddress("ZZMass",&mzz);
  bkgMC->SetBranchAddress("MC_weight_noxsec",&w);
  //  bkgMC->SetBranchAddress("superLD",&sKD);
  bkgMC->SetBranchAddress(superMelaName,&sKD);

  if (recompute_) {
    bkgMC->SetBranchAddress("Z1Mass",&m1);
    bkgMC->SetBranchAddress("Z2Mass",&m2);
    bkgMC->SetBranchAddress("helcosthetaZ1",&costheta1);
    bkgMC->SetBranchAddress("helcosthetaZ2",&costheta2);
    bkgMC->SetBranchAddress("costhetastar",&costhetastar);
    bkgMC->SetBranchAddress("helphi",&phi);
    bkgMC->SetBranchAddress("phistarZ1",&phi1);
    bkgMC->SetBranchAddress("ZZPt",&pt4l);
    bkgMC->SetBranchAddress("ZZRapidity",&Y4l);
  }
  

  bool usePseudoLDbin=melaName.Contains("pseudo");
  const int nbinsY=(usePseudoLDbin? nbinsYps : nbinsYgrav);
  float binsY[nbinsY+1];
  for(int ib=0;ib<=nbinsY;ib++){
    if(usePseudoLDbin) binsY[ib]=binsYps[ib];
    else binsY[ib]=binsYgrav[ib];
  }
  TH2F* bkgHist = new TH2F(templateName,templateName,nbinsX,binsX,nbinsY,binsY);
  // TH2F* bkgHist = new TH2F(templateName,templateName,50,0.0,1.0, 25,0.0,1.0);
  // const int nBinsFine=100;
  // float xfine[nBinsFine+1];
  // for(int i=0;i<=nBinsFine;i++)xfine[i]=i*(1.0/nBinsFine);
  //  TH2F* bkgHist = new TH2F(templateName,templateName,nBinsFine,xfine,nbinsY,binsY);



  bkgHist->Sumw2();


  //if asked to applyinterf reweighting, load the TH1
  TH1F *hInterfRewX=0,*hInterfRewY=0;
  TFile *fInterfRew=0;
  if(applyInterferenceRew){
    if((sampleIndex==0||sampleIndex==3||sampleIndex==4)
       && (channel=="4mu" || channel=="4e")){

      string curPath=gDirectory->GetPath();
      //  cout<<"Current Path is "<<curPath.c_str()<<endl;
      fInterfRew=new TFile( fInterferenceName.c_str());
      string hIntNameX="";

      if(sampleIndex==0)hIntNameX+="scalar_";
      else if(sampleIndex==3)hIntNameX+="pseudoscalar_";
      else if(sampleIndex==4)hIntNameX+="graviton_";
      else{
	cout<<"Error from fillTemplate: unrecognized sample -> "<<sampleIndex<<endl;
	hIntNameX+="unknownSample_";
      }

      string hIntNameY=hIntNameX;
      hIntNameX+="superMELA";
      if(melaName.Contains("pseudo"))hIntNameY+="pseudoMELA";
      else  if(melaName.Contains("gravi"))hIntNameY+="graviMELA";
      else{
	cout<<"Error from fillTemplate: unrecognized variable name -> "<<melaName.Data()<<endl;
	hIntNameY+="unknownVar";
      }
      hInterfRewX=(TH1F*)fInterfRew->Get(hIntNameX.c_str());
      hInterfRewY=(TH1F*)fInterfRew->Get(hIntNameY.c_str());
      gDirectory->cd(curPath.c_str());
      //      cout<<"Current Path is "<<gDirectory->GetPath()<<endl;
    }
  }



  // fill histogram
  cout<<"Input Chain has "<<bkgMC->GetEntries()<<"  entries"<<endl;
  for(int i=0; i<bkgMC->GetEntries(); i++){

    bkgMC->GetEntry(i);

    // if(i%5000==0) cout << "event: " << i << "/" << bkgMC->GetEntries() << endl;
    if(cutSameVar)KD_cut=KD;

    bool cutPassed= (melaCut>0.0) ? (KD_cut>melaCut) : true;
    if( cutPassed &&sKD>=0.0 && mzz>mzzCutLow&&mzz<mzzCutHigh){
      

      if(recompute_){

	/*
	cout << "===========================" << endl;
	cout << "mzz: " << mzz << endl;
	cout << "m1: " << m1 << endl;
	cout << "m2: " << m2 << endl;
	cout << "costheta1: " << costheta1 << endl;
	cout << "costheta2: " << costheta2 << endl;
	cout << "costhetastar: " << costhetastar << endl;
	cout << "phi: " << phi << endl;
	cout << "phi1: " << phi1 << endl;
	cout << "pt4l: " << pt4l << endl;
	cout << "Y4l: " << Y4l << endl;
	*/

	myMELA->computeKD(mzz,m1,m2,
			  costhetastar,
			  costheta1,
			  costheta2,
			  phi,
			  phi1,
			  KD,psig,pbkg,
			  withPt_,pt4l,
			  withY_, Y4l);
	
	//cout << "LD: " << LD << endl;
	
      }


      if(hInterfRewX!=0&&hInterfRewY!=0){
	w*=calcInterfRew(hInterfRewX,sKD);//reweight in the supermela direction
	w*=calcInterfRew(hInterfRewY,KD);//reweight in the sig sep KD direction
      }
      bkgHist->Fill(sKD,KD,w);
      //   bkgHist->Fill(mzz,KD,w);

      bkgHist->Fill(sKD,KD,w);
      //   bkgHist->Fill(mzz,KD,w);

    }

  }

  int nXbins=bkgHist->GetNbinsX();
  int nYbins=bkgHist->GetNbinsY();
    

  cout<<"Before smoothing, TH2 named "<<bkgHist->GetName()<<" has "<<bkgHist->GetEntries()<<"  entries and integral="<<bkgHist->Integral()<<endl;

// smooth 


    // bkgHist->Smooth();
  for(int i=1; i<=bkgHist->GetNbinsX(); i++){
    for(int j=1; j<=bkgHist->GetNbinsY(); j++){
      if(bkgHist->GetBinContent(i,j)<0.000001)
	bkgHist->SetBinContent(i,j,0.000001);
    }// for(int j=1; j<=nYbins; j++){
  }// for(int i=1; i<=nXbins; i++){

  if(smooth)  bkgHist->Smooth(1,"k5b");
 


  // normalize TH2
  double totArea=bkgHist->Integral();
  bkgHist->Scale(1.0/totArea);
  cout<<"TH2 named "<<bkgHist->GetName()<<" has "<<bkgHist->GetEntries()<<"  entries and integral="<<bkgHist->Integral()<<endl;
  /*
  //normalize in slices
  double norm;
  TH1F* tempProj;
  
  for(int i=1; i<=nXbins; i++){
    
    tempProj = (TH1F*) bkgHist->ProjectionY("tempProj",i,i);
    norm=tempProj->Integral();

    if (norm>0) { // Avoid introducing NaNs in the histogram
      for(int j=1; j<=nYbins; j++){
	bkgHist->SetBinContent(i,j, bkgHist->GetBinContent(i,j)/norm   );
      }
    }

  }
  */


  cout<<"Finishing fillTemplate for sample "<<sampleIndex<<endl;
  return bkgHist;
  
}


//===================================================
TH1F *fillKDhisto(TString channel, int sampleIndex,float mzzLow,float mzzHigh,  bool smooth){


  TChain* bkgMC = new TChain("SelectedTree");//"data_obs" "SelectedTree"
  buildChainSingleMass(bkgMC, channel, sampleIndex,mH);
  float mzz=-1.0,mela=-444.0,w=0;
  double LD=-999.;
  char yVarName[32];
  bkgMC->SetBranchAddress("ZZMass",&mzz);
  bkgMC->SetBranchAddress("ZZLD",&mela);
  bkgMC->SetBranchAddress("superLD",&LD);
  bkgMC->SetBranchAddress("MC_weight_noxsec",&w);


  char hTitle[128];
  sprintf(hTitle,"Distribution of superMELA KD with M_{4l} in [%d, %d]",int(mzzLow),int(mzzHigh));
  TH1F* outHist=new TH1F("finHisto",hTitle,200,0.0,1.0);
  outHist->Sumw2();
 // fill histogram
 // cout<<"Looping on tree netries (fillKDHisto) "<<bkgMC->GetEntries()<<endl;
  for(int i=0; i<bkgMC->GetEntries(); i++){
    bkgMC->GetEntry(i);
    //   if(i%200==0)cout<<"Entry "<<i<<"  mzz="<<mzz<<"  weight="<<w<<"  LD="<<LD<<endl;
    bool cutPassed= (melaCut>0.0) ? (mela>melaCut) : true;
    if( cutPassed &&mzz>mzzLow&&mzz<mzzHigh){
      outHist->Fill(LD,w);
    }
  }
  if(smooth) outHist->Smooth(1,"R"); 
  //normalize to unity as this is supposed to be a pdf
  outHist->Scale(1.0/outHist->Integral());

  return outHist;
}

//=======================================================================

void makeTemplate(TString channel){

  
  // =======================================
  // Z+X CR template
  TFile* fqqZZ = new TFile(destDir + "Dbackground_ZJetsCR_AllChans.root","RECREATE");
  pair<TH2F*,TH2F*> histoPair;
  TH2F* oldTemp;
  TH2F* low,*high,*h_mzzD;
  TH2F *h_mzzD_syst1Up,*h_mzzD_syst1Down,*h_mzzD_syst2Up,*h_mzzD_syst2Down;
  TH1F *h_D;
  
  low = fillTemplate(channel,5,"superLD","zjetsHist",true);
  //  high = fillTemplate(channel,1,false);
  h_mzzD =   (TH2F*)low->Clone("h_mzzD");// mergeTemplates(low,high);
  h_D = fillKDhisto(channel,5,mzzCutLow,mzzCutHigh,true);//last two are cuts on mZZ
 
  // ---------- apply interference reweighting --------
  
  oldTemp = new TH2F(*h_mzzD);
  oldTemp->SetName("oldTemp");

  //systematic on MELA shape for Z+X currently not available, to be estimated
  //  cout << "apply systematics for zjets control region" << endl;
  //  histoPair = reweightForCRunc(h_mzzD);
  histoPair.first=h_mzzD;
  histoPair.second=h_mzzD;
  // --------------------------------------------------
  string strSqrt="7TeV";
  if(useSqrts==2)strSqrt="8TeV";
  makePlot1D( h_D,strSqrt+"_ZJetsCR" );
  makePlot2D( h_mzzD,strSqrt+"_ZJetsCR" );

  fqqZZ->cd();
  h_mzzD->Write("h_superDpsD");
  oldTemp->Write("oldTemp");
  h_D->Write("h_superD");

  
  fqqZZ->Close();


}

//=======================================================================

void storeLDDistribution(){

  makeTemplate("4mu");

}


void generateTemplatesSMD_forCR() {
  
  storeLDDistribution();
}



// NOT USED
void buildChain(TChain* bkgMC, TString channel, int sampleIndex) {

  //  TString sample[4]={"H*","ZZTo*","ggZZ*","H*Pse"};
  //  TString sampleName[4]={"signal","qqZZ","ggZZ","signal_PS"};

  TString chPath = (channel=="2e2mu"?"2mu2e":channel); // Adapt to different naming convention...
  int mySqrts=useSqrts;

  if(recompute_ && (withPt_||withY_) &&useSqrts==0){
    cout<<"ERROR: use either 7 or 8 TeV to run Pt/Y templates"<<endl;
    if(sqrts==7)mySqrts=1;
    else mySqrts=2;
    cout<<"Setting data sets to "<<sqrts<<"TeV. Please check this is correct!"<<endl;
  }
  if(sampleIndex==0){
    //7TeV
    if(mySqrts<2){
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H120.root");
      //    bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H125.root"); // Skip: Special sample with a different composition
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H130.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H140.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H150.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H160.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H170.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H180.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H190.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H200.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H210.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H220.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H250.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H275.root"); // Missing in 240612
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H300.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H325.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H350.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H400.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H425.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H450.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H475.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H525.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H550.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H575.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H600.root");
      if (extendToHighMass) {
	bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H650.root");
	bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H700.root");
	bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H750.root");
	bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H800.root");  // Sample to be debugged, unexpected efficiency in 2e2mu
	bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H850.root");  // Missing in 310812
	bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H900.root");
	bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H950.root");
	bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H1000.root"); // Missing in 310812
      }
    }
    if(mySqrts%2==0){
      //8TeV
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H115.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H116.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H117.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H118.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H119.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H120.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H121.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H122.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H123.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H124.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H125.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H126.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H127.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H128.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H129.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H130.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H145.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H150.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H180.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H200.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H250.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H300.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H325.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H350.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H400.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H450.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H500.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H550.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H600.root");
      if (extendToHighMass) {
	bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H650.root");
	bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H700.root");
	bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H750.root");
	bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H800.root");
	bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H850.root");
	bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H900.root");
	bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H950.root");
	bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H1000.root");
      }
    }
    
  } else if (sampleIndex==1){
    //7TeV
    if(mySqrts<2){
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ZZTo2e2mu.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ZZTo2e2tau.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ZZTo2mu2tau.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ZZTo4e.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ZZTo4mu.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ZZTo4tau.root");
    }
    //8TeV
    if(mySqrts%2==0){
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_ZZTo2e2mu.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_ZZTo2e2tau.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_ZZTo2mu2tau.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_ZZTo4e.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_ZZTo4mu.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_ZZTo4tau.root");
    }
  } else if (sampleIndex==2){
    //7TeV
    if(useSqrts<2){
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ggZZ2l2l.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ggZZ4l.root");
    }
    //8TeV
    if(mySqrts%2==0){
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_ggZZ2l2l.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_ggZZ4l.root");
    }
    
  } else if(sampleIndex==3){ //this is for alternative signal samples
    abort(); // Standard location of these files still being arranged.
    //       sprintf(temp,"CJLSTtree_Jun25_2012/JHUsignal/HZZ%sTree_%s.root",channel,sample[sampleIndex].c_str());
    //       bkgMC->Add(temp);
  }
    
}


TH2F* mergeTemplates(TH2F* lowTemp, TH2F* highTemp){

  int nYbins=lowTemp->GetNbinsY();
  if (highTemp->GetNbinsY()!=nYbins) {
    cout << "ERROR: mergeTemplates: incorrect binning " << endl;
    abort();
  }

  TH2F* h_mzzD = new TH2F("h_mzzD","h_mzzD",int((highMzz-100.)/mBinSize +0.5),100,highMzz,nYbins,0,1);

  // copy lowmass into h_mzzD
  for(int i=1; i<=lowTemp->GetNbinsX(); ++i){
    for(int j=1; j<=nYbins; ++j){
      h_mzzD->SetBinContent(i,j, lowTemp->GetBinContent(i,j)  );
    }// end loop over D
  }// end loop over mZZ

  // copy high mass into h_mzzD
  for(int i=1; i<=highTemp->GetNbinsX(); ++i){
    for(int j=1; j<=nYbins; ++j){
      h_mzzD->SetBinContent(i+lowTemp->GetNbinsX(),j, highTemp->GetBinContent(i,j)  );
    }// end loop over D
  }// end loop over mZZ

  return h_mzzD;

}

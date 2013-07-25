/* 
 * Create 2D (mass, LD) templates. Script imported from: http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/JHU/MELA/scripts/generateTemplates.C?revision=1.1.2.1&view=markup&pathrev=post_unblinding
 * Requires ZZMatrixElement/MELA to have been checked out and compiled.
 * usage: 
 * -set input paths variables in Config.h
 * -run with:
 * root -q -b loadMELA.C generateTemplates.C+
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


//----------> SET INPUT VARIABLES in Config.h
#include "Config.h"
//<----------

//--- Flags to control re-computation of KD
bool recompute_ = false;       // Recompute KD (instead of taking it from the tree); when true, the following flags apply:
bool usePowhegTemplate=false;  // false use analytic bg
bool withPt_ = false;          // Include pT in KD
bool withY_  = false;          //    "    Y  "  "
int sqrts    = 7;              // sqrts, used only for withPt_/withY_

//---
int useSqrts=0;              //0=use 7+8TeV; 1=use 7TeV only, 2 use 8TeV only
TString SigName = "p0plus_VAJHU"; // name of MELA branch to be used. Possibilities are ZZLD,ZZLD_analBkg,ZZLD_postICHEP,ZZLD_PtY,pseudoMelaLD, spinTwoMinimalMelaLD 
TString BkgName = "bkg_VAMCFM";
TString melaName = "VAKD";

bool makePSTemplate = false;
bool makeAltSignal = false;
const float melaCut=-1.0;
bool extendToHighMass = true; // Include signal samples above 600 GeV

float highMzz=(extendToHighMass?1600:800);
float mBinSize=2.;

//Mela* myMELA;

const TString destDir = "../CreateDatacards/templates2D/";

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

  // ================ binning for MELA ====================================
  double lowBins[numPtmp]={100.,120.,140.,160.,180.,220.,260.,300.};
  double highBins[numPtmp]={120.,140.,160.,180.,220.,260.,300.,1602.};
  // ======================================================================


  // ================ systematics for pseudoMELA ==========================
  double slopePS_syst[numPtmpPS] ={-3.32705e-01, -1.90814e-01, -9.77189e-01, -3.81680e-01, 0.0 };
  double yIntrPS_syst[numPtmpPS] ={ 9.05727e-01, 9.95995e-01,  1.40367e+00,  1.12690,      1.0 }; 
  //  ==================================================================

  // ================ systematics for MELA ==========================
  double slope_syst[numPtmp]={1.55836,1.36721,0.846938,-1.72509,-1.821,-1.21998,-1.20542,0.0};
  double yIntr_syst[numPtmp]={0.406822,0.512352,0.733718,1.40255,1.10837,0.941797,0.91357,1.0};
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
      if(newTempValue <= 0.) newTempValue=0.00000001;
      tempUp->SetBinContent(i,j,newTempValue);
      newTempValue = oldTempValue*(-slope[point]*(double)j/30.+2.-yIntr[point]);
      if(newTempValue <= 0.) newTempValue=0.00000001;
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
/*
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
    reweightFunc =new TF1("reweightFunc","gaus",100,1600);
    
    reweightFunc->SetParameter(0,0.4053);
    reweightFunc->SetParameter(1,110.6);
    reweightFunc->SetParameter(2,22.57);
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
    slope=reweightFunc->Eval((double)((i-1)*2+101));

     ==============================================
    // for reweighting pseudo-MELA
    //slope = reweightFunc->Eval((double)((i-1)*2+101));
    ============================================== 

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
*/


void buildChain(TChain* bkgMC, TString channel, int sampleIndex=0) {

  //  TString sample[4]={"H*","ZZTo*","ggZZ*","H*Pse"};
  //  TString sampleName[4]={"signal","qqZZ","ggZZ","signal_PS"};

  TString chPath = (channel=="2e2mu"?"2mu2e":channel); // Adapt to different naming convention...

  //An error is issued on missing files; if a single file is missing in one set it can be safely ignored.
  if(recompute_ && (withPt_||withY_) &&useSqrts==0){
    cout<<"ERROR: use either 7 or 8 TeV to run Pt/Y templates"<<endl;
    if(sqrts==7)useSqrts=1;
    else useSqrts=2;
    cout<<"Setting data sets to "<<sqrts<<"TeV. Please check this is correct!"<<endl;
  }
  if(sampleIndex==0){
    //7TeV
    if(useSqrts<2){
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H115.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H120.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H122.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H124.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H125.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H126.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H128.root");
      //bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H130.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H140.root");
      //bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H150.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H160.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H170.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H180.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H185.root");
      //bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H190.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H200.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H210.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H220.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H250.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H275.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_powheg15H300.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H325.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_powheg15H350.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_powheg15H400.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H425.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_powheg15H450.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H475.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_powheg15H500.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H525.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_powheg15H550.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H575.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_powheg15H600.root");
      if (extendToHighMass) {
	bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_powheg15H650.root");
	bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_powheg15H700.root");
	bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H750.root");
	bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_powheg15H800.root");
	// bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H850.root");
	bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_powheg15H900.root");
	bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H950.root");
	bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_powheg15H1000.root");
	}
    }
    if(useSqrts%2==0){
      //8TeV
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H115.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H120.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H122.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H124.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H125.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H126.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H128.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H130.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H135.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H140.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H145.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H150.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H160.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H175.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H170.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H180.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H185.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H190.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15jhuGenV3H200.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15H225.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15H250.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15H275.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15H300.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H325.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15H350.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H375.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15H400.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H425.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15H450.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H475.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15H500.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H525.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15H550.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H575.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15H600.root");
      if (extendToHighMass) {
	bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15H650.root");
	bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15H700.root");
	bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H750.root");
	bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15H800.root");
	bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H850.root");
	bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15H900.root");
	bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_H950.root");
	bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_powheg15H1000.root");
      }
    }
    
  } else if (sampleIndex==1){
    //7TeV
    if(useSqrts<2){
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ZZTo2e2mu.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ZZTo2e2tau.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ZZTo2mu2tau.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ZZTo4e.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ZZTo4mu.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ZZTo4tau.root");
    }
    //8TeV
    if(useSqrts%2==0){
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
    if(useSqrts%2==0){
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_ggZZ2l2l.root");
      bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_ggZZ4l.root");
    }
    
  } else if(sampleIndex==3){ //this is for alternative signal samples
    abort(); // Standard location of these files still being arranged.
    //       sprintf(temp,"CJLSTtree_Jun25_2012/JHUsignal/HZZ%sTree_%s.root",channel,sample[sampleIndex].c_str());
    //       bkgMC->Add(temp);
  }
    
}


//=======================================================================

TH2F* fillTemplate(TString channel="4mu", int sampleIndex=0,bool isLowMass=true){
  TChain* bkgMC = new TChain("SelectedTree");

  if (isLowMass) {
    buildChain(bkgMC, channel, sampleIndex);
  } else {
    buildChain(bkgMC, "2e2mu", sampleIndex);
    buildChain(bkgMC, "4e", sampleIndex);
    buildChain(bkgMC, "4mu", sampleIndex);
  }

  cout << "Chain for " << channel << " " << sampleIndex << " " << isLowMass << " " << bkgMC->GetEntries() << endl;
  bkgMC->ls();

  float mzz,KD,KD_cut,w;
  //  float interfw=0.;
  float m1=0, m2=0, costheta1=0, costheta2=0, costhetastar=0, phi=0, phi1=0;
  float pt4l=0, Y4l=0;
  float psig=0, pbkg=0;
  
  
  //distinction btw LD and mela needed because we might want 
  //both psMELA (for 2D template) and MELA (for cut)

  TString melaCutName = melaName;  
  if(makePSTemplate) {
    melaName = "ZZpseudoLD";
  }
  
  bkgMC->SetBranchAddress(SigName.Data(),&psig);
  bkgMC->SetBranchAddress(BkgName.Data(),&pbkg);
  if (melaName!=melaCutName) {
    bkgMC->SetBranchAddress(melaCutName.Data(),&KD_cut);
  } else {
    //    KD_cut=KD; //FIXME: This is broken, since KD is not initialized at this point! 
    KD_cut= 1.;      //       Set it to one, just to fix it to something meaningful (it's unused anyhow)
  }
  
  bkgMC->SetBranchAddress("ZZMass",&mzz);
  bkgMC->SetBranchAddress("MC_weight_noxsec",&w);
  /*
  if(sampleIndex == 0 && isLowMass && (channel == "4mu" || channel == "4e"))
    {
      bkgMC->SetBranchAddress("LIWeight",&interfw);
    }
  */
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
  
  TH2F* bkgHist;
  if(!isLowMass)
    bkgHist = new TH2F("bkgHisto","bkgHisto",int((highMzz-180.)/mBinSize+0.5),180,highMzz,30,0,1);
  else
    bkgHist = new TH2F("bkgHisto","bkgHisto",int((180-100)/mBinSize+0.5),100,180,30,0,1);

  bkgHist->Sumw2();

  // fill histogram
	
  for(int i=0; i<bkgMC->GetEntries(); i++){

    bkgMC->GetEntry(i);
    float weight = 0.;

    if(i%100000==0) cout << "event: " << i << "/" << bkgMC->GetEntries() << endl;

    bool cutPassed= (melaCut>0.0) ? (KD_cut>melaCut) : true;
    if(w<.0015 && cutPassed){
      

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
/*
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
*/	
      }

      KD = psig/(psig+pbkg);
      //if(sampleIndex == 0 && isLowMass && (channel == "4mu" || channel == "4e")) weight = w*interfw;
      weight = w;
      bkgHist->Fill(mzz,KD,weight);

    }

  }


  int nXbins=bkgHist->GetNbinsX();
  int nYbins=bkgHist->GetNbinsY();
    
  // normalize slices

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

  
  // average 

  TH2F* notSmooth = new TH2F(*bkgHist);

  if(!isLowMass){
    
    int effectiveArea=1;
    double average=0,binsUsed=0;

    for(int i=1; i<=nXbins; i++){
      for(int j=1; j<=nYbins; j++){
	
	//	binMzz=(i-1)*2+181;
	float binMzz = bkgHist->GetBinCenter(i);

	if( binMzz<300 ) continue;
	if( binMzz>=300 && binMzz<350 ) effectiveArea=1;
	if( binMzz>=350 && binMzz<500 ) effectiveArea=3;
	if( binMzz>=500 && binMzz<600 ) effectiveArea=5;
	if( binMzz>=600 && binMzz<800 ) effectiveArea=7;
	if( binMzz>=800 && binMzz<1000) effectiveArea=11;
	if( binMzz>=1000 && binMzz<1200)effectiveArea=15;
	if( binMzz>=1200 ) effectiveArea=25;
	
	for(int a=-effectiveArea; a<=effectiveArea; a++){
	  if(a+i<1 || a+i>nXbins || j>nYbins || j<1) continue;
	  average+= notSmooth->GetBinContent(a+i,j);
	  binsUsed++;
	}
	
	bkgHist->SetBinContent(i,j,average/binsUsed);

	average=0;
	binsUsed=0;
	
      } // end loop over D
    } // end loop over mZZ
  } // end of horizontal averaging
  
  // smooth

  bkgHist->Smooth();
  if(!isLowMass)
    bkgHist->Smooth();
  
  for(int i=1; i<=nXbins; i++){
    for(int j=1; j<=nYbins; j++){
      if(bkgHist->GetBinContent(i,j)==0)
	bkgHist->SetBinContent(i,j,.00001);
    }// for(int j=1; j<=nYbins; j++){
  }// for(int i=1; i<=nXbins; i++){

  return bkgHist;
  
}

//=======================================================================

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

//=======================================================================

void makeTemplate(TString channel="4mu"){

  //  sprintf(temp,"../datafiles/Dsignal_%s.root",channel.Data());
  TFile* fsig = new TFile(destDir + "Dsignal_" + channel + ".root","RECREATE");
  TFile* fAltsig = 0;
  if (makeAltSignal) {
    //    sprintf(temp,"../datafiles/Dsignal_ALT_%s.root",channel.Data());
    fAltsig = new TFile(destDir + "Dsignal_ALT_" + channel + ".root","RECREATE");
  }
  //  sprintf(temp,"../datafiles/Dbackground_qqZZ_%s.root",channel.Data());
  TFile* fqqZZ = new TFile(destDir + "Dbackground_qqZZ_" + channel + ".root","RECREATE");
  //  sprintf(temp,"../datafiles/Dbackground_ggZZ_%s.root",channel.Data());
  TFile* fggZZ = new TFile(destDir + "Dbackground_ggZZ_" + channel + ".root","RECREATE");
  TFile* fZX = new TFile(destDir + "Dbackground_ZX_" + channel + ".root","RECREATE");
  TH2F* oldTemp;

  pair<TH2F*,TH2F*> histoPair;

  TH2F* low,*high,*h_mzzD;
  
  // ========================================
  // SM Higgs template

  low = fillTemplate(channel,0,true);
  high = fillTemplate(channel,0,false);
  h_mzzD = mergeTemplates(low,high);


  // ---------- apply interference reweighting --------
  
  /*oldTemp = new TH2F(*h_mzzD);
  oldTemp->SetName("oldTemp");

  cout << "correcting for interference and adding syst" << endl;
  if(channel!="2e2mu")
    h_mzzD = reweightForInterference(h_mzzD);  // correct templates for lepton interference

    cout << "h_mzzD: " << h_mzzD << endl;*/

  // --------------------------------------------------

  fsig->cd();

  h_mzzD->Write("h_mzzD");
  //oldTemp->Write("oldTemp");
  h_mzzD->Write("h_mzzD_up");
  h_mzzD->Write("h_mzzD_dn");
  fsig->Close();

  // ========================================
  // alternative signal template

  if (makeAltSignal) {
    low = fillTemplate(channel,3,true);
    high = fillTemplate(channel,3,false);
    h_mzzD = mergeTemplates(low,high);

    // ---------- apply interference reweighting --------
  
    /*oldTemp = new TH2F(*h_mzzD);
    oldTemp->SetName("oldTemp");

    cout << "correcting for interference and adding syst" << endl;
    if(channel!="2e2mu")
      h_mzzD = reweightForInterference(h_mzzD);  // correct templates for lepton interference

      cout << "h_mzzD: " << h_mzzD << endl;*/

  // --------------------------------------------------
    
    fAltsig->cd();
    h_mzzD->Write("h_mzzD");
    //oldTemp->Write("oldTemp");
    histoPair.first->Write("h_mzzD_up");
    histoPair.second->Write("h_mzzD_dn");
    fAltsig->Close();
  }//end if makeAltSignal
  
  // =======================================
  // qqZZ template

  low = fillTemplate(channel,1,true);
  high = fillTemplate(channel,1,false);
  h_mzzD = mergeTemplates(low,high);

  // ---------- apply interference reweighting --------
  
  oldTemp = new TH2F(*h_mzzD);
  oldTemp->SetName("oldTemp");

  cout << "apply systematics for zjets control region" << endl;
  
  histoPair = reweightForCRunc(h_mzzD);

  // --------------------------------------------------

  fqqZZ->cd();
  h_mzzD->Write("h_mzzD");
  oldTemp->Write("oldTemp");
  histoPair.first->Write("h_mzzD_up");
  histoPair.second->Write("h_mzzD_dn");
  fqqZZ->Close();

  // ==========================
  // ggZZ templates
  
  low = fillTemplate(channel,2,true);
  high = fillTemplate(channel,2,false);
  h_mzzD = mergeTemplates(low,high);

  // ---------- apply interference reweighting --------
  
  oldTemp = new TH2F(*h_mzzD);
  oldTemp->SetName("oldTemp");

  cout << "apply systematics for zjets control region" << endl;
  
  histoPair = reweightForCRunc(h_mzzD);

  // --------------------------------------------------

  fggZZ->cd();
  h_mzzD->Write("h_mzzD");
  oldTemp->Write("oldTemp");
  histoPair.first->Write("h_mzzD_up");
  histoPair.second->Write("h_mzzD_dn");
  fggZZ->Close();

  // =======================================
  // ZX template

  //make as if ZZ
  low = fillTemplate(channel,1,true);
  high = fillTemplate(channel,1,false);
  h_mzzD = mergeTemplates(low,high);

  // ---------- apply interference reweighting --------
  
  oldTemp = new TH2F(*h_mzzD);
  oldTemp->SetName("oldTemp");

  //cout << "apply systematics for zjets control region" << endl;
  
  histoPair = reweightForCRunc(h_mzzD);

  // --------------------------------------------------

  fZX->cd();
  //Write ZZ up template as ZX default
  histoPair.first->Write("h_mzzD");
  oldTemp->Write("oldTemp");
  //Write ZX + (ZX-ZZ) as ZX up
  TH2F* ZXup = new TH2F(*h_mzzD);
  for(int ZXbinx = 1; ZXbinx < ZXup->GetNbinsX(); ZXbinx++)
    {
      for(int ZXbiny=1; ZXbiny<=ZXup->GetNbinsY(); ZXbiny++){

      float oldTempValue = histoPair.first->GetBinContent(ZXbinx,ZXbiny);
      float diffTempValue = oldTempValue - (h_mzzD->GetBinContent(ZXbinx,ZXbiny));
      float newTempValue = oldTempValue + diffTempValue;
      if(newTempValue <= 0.) newTempValue=0.00000001;
      ZXup->SetBinContent(ZXbinx,ZXbiny,newTempValue);

    }// end loop over Y bins

      // -------------- normalize mZZ slice ----------------

      double norm_up=(ZXup->ProjectionY("temp",ZXbinx,ZXbinx))->Integral();
      for(int ZXbiny=1; ZXbiny<=ZXup->GetNbinsY(); ZXbiny++){
      
	ZXup->SetBinContent(ZXbinx,ZXbiny,ZXup->GetBinContent(ZXbinx,ZXbiny)/norm_up);

      }

    // ---------------------------------------------------

  }// end loop over X bins
  ZXup->Write("h_mzzD_up");
  //Write ZZ default tempalte as ZX down
  h_mzzD->Write("h_mzzD_dn");
  fZX->Close();


}

//=======================================================================

void storeLDDistribution(){

  makeTemplate("4mu");
  makeTemplate("4e");
  makeTemplate("2e2mu");

}


void generateTemplates() {
  
  //myMELA=0;
  //if (recompute_) myMELA = new Mela(usePowhegTemplate, sqrts); // this is safely leaked
  storeLDDistribution();
}

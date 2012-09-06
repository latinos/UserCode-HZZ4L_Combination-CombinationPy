/* 
 * Create 2D (mass, LD) templates. Script imported from: http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/JHU/MELA/scripts/generateTemplates.C?revision=1.1.2.1&view=markup&pathrev=post_unblinding
 * usage: 
 * -set input paths variables in Config.h
 * -run with:
 * root -q -b generateTemplates.C+
 * 2D templates are written to 
 *
 */

#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include <sstream>
#include <vector>


//----------> SET INPUT VARIABLES in Config.h
#include "Config.h"
//<----------

bool makePSTemplate = false;

bool extendToHighMass = false; // Include signal samples above 600 GeV

float highMzz=(extendToHighMass?1000:800);
float mBinSize=2.;


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

  const int numPoints=8;

  double low[numPoints]   ={100.,        120.,        140.,         160.,     180.  };
  double high[numPoints]  ={120.,        140.,        160.,         180.,     1000. };

  /* ================ systematics for pseudoMELA ==========================
  double slope[numPoints] ={-3.32705e-01, -1.90814e-01, -9.77189e-01, -3.81680e-01, 0.0 };
  double yIntr[numPoints] ={ 9.05727e-01, 9.95995e-01,  1.40367e+00,  1.12690,      1.0 }; 
  ==================================================================*/

  // ================ systematics for MELA ==========================
  double slope[numPoints] ={4.71836e-01, 1.17671e-01, -3.81680e-01, -1.20481, -1.21944, -2.06928, -1.35337, 0.0 };
  double yIntr[numPoints] ={6.83860e-01, 9.38454e-01, 1.12690,      1.24502,  1.72764,  2.11050,  1.52771,  1.0 }; 
  //==================================================================


  for(int i=1; i<=temp->GetNbinsX(); i++){
    point = -1;

    // choose correct scale factor
    for(int p=0; p<numPoints; p++){
      if( (i*2.+101.)>=low[p] && (i*2.+101.)<high[p] ){
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
  TF1* reweightFunc = new TF1("reweightFunc","gaus",100,1000);

  reweightFunc->SetParameter(0,0.354258);
  reweightFunc->SetParameter(1,114.909);
  reweightFunc->SetParameter(2,17.1512);

  /* ===================================================
  // for interference reweighting of pseudo-MELA
  TF1* reweightFunc = new TF1("reweightFunc","([0]+[1]*(x-110) )*0.5*(1 + TMath::Erf([2]*(x -[3]) ))*0.5*(1 + TMath::Erf([4]*([5]-x) ))  ",100,200);

  reweightFunc->SetParameter(0,-5.66409e-01);
  reweightFunc->SetParameter(1, 1.22591e-02);
  reweightFunc->SetParameter(2, 1.64942e+00);
  reweightFunc->SetParameter(3, 1.10080e+02);
  reweightFunc->SetParameter(4, 2.10905e+00);
  reweightFunc->SetParameter(5, 1.78529e+02);
  ==================================================== */

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
    if(i<8){
      slope=.354;
    }else{
      slope=reweightFunc->Eval((double)((i-1)*2+101));
    }

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



void buildChain(TChain* bkgMC, TString channel, int sampleIndex=0) {

  //  TString sample[4]={"H*","ZZTo*","ggZZ*","H*Pse"};
  //  TString sampleName[4]={"signal","qqZZ","ggZZ","signal_PS"};

  TString chPath = (channel=="2e2mu"?"2mu2e":channel); // Adapt to different naming convention...

  //An error is issued on missing files; if a single file is missing in one set it can be safely ignored.

  if(sampleIndex==0){
    //7TeV
    bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H120.root");
    bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H125.root");
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
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H600.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H650.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H700.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H750.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H800.root");
      //     bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H850.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H900.root");
      bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H950.root");
      //     bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_H1000.root");
    }
    

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
    
  } else if (sampleIndex==1){
    //7TeV
    bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ZZTo2e2mu.root");
    bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ZZTo2e2tau.root");
    bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ZZTo2mu2tau.root");
    bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ZZTo4e.root");
    bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ZZTo4mu.root");
    bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ZZTo4tau.root");

    //8TeV
    bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_ZZTo2e2mu.root");
    bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_ZZTo2e2tau.root");
    bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_ZZTo2mu2tau.root");
    bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_ZZTo4e.root");
    bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_ZZTo4mu.root");
    bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_ZZTo4tau.root");

  } else if (sampleIndex==2){
    //7TeV
    bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ggZZ2l2l.root");
    bkgMC->Add(filePath7TeV + "/" + chPath +"/HZZ4lTree_ggZZ4l.root");

    //8TeV
    bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_ggZZ2l2l.root");
    bkgMC->Add(filePath8TeV + "/" + chPath +"/HZZ4lTree_ggZZ4l.root");

    
  } else if(sampleIndex==3){
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

  float mzz,mela,w;
  
  bkgMC->SetBranchAddress("ZZMass",&mzz);
  bkgMC->SetBranchAddress("ZZLD",&mela);
  bkgMC->SetBranchAddress("MC_weight_noxsec",&w);
  
  TH2F* bkgHist;
  if(!isLowMass)
    bkgHist = new TH2F("bkgHisto","bkgHisto",int((highMzz-180.)/mBinSize+0.5),180,highMzz,30,0,1);
  else
    bkgHist = new TH2F("bkgHisto","bkgHisto",int((180-100)/mBinSize+0.5),100,180,30,0,1);

  bkgHist->Sumw2();

  // fill histogram

  for(int i=0; i<bkgMC->GetEntries(); i++){

    bkgMC->GetEntry(i);

    if(w<.0015 /*&& mela>.5 */ ){
      
      bkgHist->Fill(mzz,mela,w);

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

    for(int j=1; j<=nYbins; j++){
      bkgHist->SetBinContent(i,j, bkgHist->GetBinContent(i,j)/norm   );
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
	if( binMzz>=600 ) effectiveArea=7;
	
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
  TString destDir = "../CreateDatacards/templates2D/";

  //  sprintf(temp,"../datafiles/Dsignal_%s.root",channel.Data());
  TFile* fsig = new TFile(destDir + "Dsignal_" + channel + ".root","RECREATE");
  TFile* fpssig = 0;
  if (makePSTemplate) {
    //    sprintf(temp,"../datafiles/Dsignal_ALT_%s.root",channel.Data());
    fpssig = new TFile(destDir + "Dsignal_ALT_" + channel + ".root","RECREATE");
  }
  //  sprintf(temp,"../datafiles/Dbackground_qqZZ_%s.root",channel.Data());
  TFile* fqqZZ = new TFile(destDir + "Dbackground_qqZZ_" + channel + ".root","RECREATE");
  //  sprintf(temp,"../datafiles/Dbackground_ggZZ_%s.root",channel.Data());
  TFile* fggZZ = new TFile(destDir + "Dbackground_ggZZ_" + channel + ".root","RECREATE");
  TH2F* oldTemp;

  pair<TH2F*,TH2F*> histoPair;

  TH2F* low,*high,*h_mzzD;
  
  // ========================================
  // SM Higgs template

  low = fillTemplate(channel,0,true);
  high = fillTemplate(channel,0,false);
  h_mzzD = mergeTemplates(low,high);

  // ---------- apply interference reweighting --------
  
  oldTemp = new TH2F(*h_mzzD);
  oldTemp->SetName("oldTemp");

  cout << "correcting for interference and adding syst" << endl;
  if(channel!="2e2mu")
    h_mzzD = reweightForInterference(h_mzzD);  // correct templates for lepton interference

  cout << "h_mzzD: " << h_mzzD << endl;

  // --------------------------------------------------

  fsig->cd();

  h_mzzD->Write("h_mzzD");
  oldTemp->Write("oldTemp");
  h_mzzD->Write("h_mzzD_up");
  h_mzzD->Write("h_mzzD_dn");
  fsig->Close();

  // ========================================
  // pseudo scalar template

  if (makePSTemplate) {
    low = fillTemplate(channel,3,true);
    high = fillTemplate(channel,3,false);
    h_mzzD = mergeTemplates(low,high);

    // ---------- apply interference reweighting --------
  
    oldTemp = new TH2F(*h_mzzD);
    oldTemp->SetName("oldTemp");

    cout << "correcting for interference and adding syst" << endl;
    if(channel!="2e2mu")
      h_mzzD = reweightForInterference(h_mzzD);  // correct templates for lepton interference

    cout << "h_mzzD: " << h_mzzD << endl;

  // --------------------------------------------------
    
    fpssig->cd();
    h_mzzD->Write("h_mzzD");
    oldTemp->Write("oldTemp");
    histoPair.first->Write("h_mzzD_up");
    histoPair.second->Write("h_mzzD_dn");
    fpssig->Close();
  }
  
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

}

//=======================================================================

void storeLDDistribution(){

  makeTemplate("4mu");
  makeTemplate("4e");
  makeTemplate("2e2mu");

}


void generateTemplates() {
  storeLDDistribution();
}

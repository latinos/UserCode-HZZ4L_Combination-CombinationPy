/* 
 * Create 2D (mass, LD) templates. Script imported from: http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/JHU/MELA/scripts/generateTemplates.C?revision=1.1.2.1&view=markup&pathrev=post_unblinding
 * usage: 
 * -set input paths variables in Config.h
 * -run with:
 * root -q -b generateDummyObsData.C+
 * 2D templates are written to "destDir"
 *
 */

#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include <sstream>
#include <fstream>
#include <vector>


//----------> SET INPUT VARIABLES in Config.h
#include "ConfigSMD.h"
//<----------

bool makePSTemplate = true;

const float melaCut=-0.5;
bool extendToHighMass = false; // Include signal samples above 600 GeV
const bool verbose_=true;
float highMzz=(extendToHighMass?1000:800);
float mBinSize=2.;
const float fracPS=0.50;
string strLumi="XX.YY";

const TString destDir = "./";//it must already exist


void buildChain(TChain* bkgMC, TString channel, int sampleIndex=0) {

  //  TString sample[4]={"H*","ZZTo*","ggZZ*","H*Pse"};
  //  TString sampleName[4]={"signal","qqZZ","ggZZ","signal_PS"};

  TString chPath = channel;//(channel=="2e2mu"?"2mu2e":channel); // Adapt to different naming convention...

  TString filePath(filePath7TeV);
  TString filePathPS(filePath7TeVPS);
    
if (sampleIndex==0){
    //8TeV
  cout<<"Adding Signal-MC tree called "<<(filePath + "/" + chPath +"/HZZ4lTree_H125_withSMD_doubleCBonly.root").Data()<<endl;
    bkgMC->Add(filePath + "/" + chPath +"/HZZ4lTree_H125_withSMD_doubleCBonly.root");


 
  } else if (sampleIndex==1){

    bkgMC->Add(filePath + "/" + chPath +"/HZZ4lTree_ZZTo2e2mu_withSMD_doubleCBonly.root");
    bkgMC->Add(filePath + "/" + chPath +"/HZZ4lTree_ZZTo2e2tau_withSMD_doubleCBonly.root");
    bkgMC->Add(filePath + "/" + chPath +"/HZZ4lTree_ZZTo2mu2tau_withSMD_doubleCBonly.root");
    bkgMC->Add(filePath + "/" + chPath +"/HZZ4lTree_ZZTo4e_withSMD_doubleCBonly.root");
    bkgMC->Add(filePath + "/" + chPath +"/HZZ4lTree_ZZTo4mu_withSMD_doubleCBonly.root");
    bkgMC->Add(filePath + "/" + chPath +"/HZZ4lTree_ZZTo4tau_withSMD_doubleCBonly.root");

    //  } else if (sampleIndex==2){
    //7TeV
    //   bkgMC->Add(filePath + "/" + chPath +"/HZZ4lTree_ggZZ2l2l.root");
    // bkgMC->Add(filePath + "/" + chPath +"/HZZ4lTree_ggZZ4l.root");

    //8TeV
    bkgMC->Add(filePath + "/" + chPath +"/HZZ4lTree_ggZZ2l2l_withSMD_doubleCBonly.root");
    bkgMC->Add(filePath + "/" + chPath +"/HZZ4lTree_ggZZ4l_withSMD_doubleCBonly.root");

    
  }
 else if (sampleIndex==3){
    bkgMC->Add(filePathPS + "/" + chPath +"/HZZ4lTree_jhuPseH125_withSMD_doubleCBonly.root");
 }
 else cout<<"Wrong sample "<<sampleIndex<<endl;
}


//=======================================================================

void fillTemplate(TString channel="4mu", int sampleIndex=0,float lumi=1.0,bool isLowMass=true){


  TChain* treeMCSig = new TChain("SelectedTree");
  TChain* treeMCPS = new TChain("SelectedTree");
  TChain* treeMCBkg = new TChain("SelectedTree");
 TChain* treeMC = new TChain("SelectedTree");
  if(sampleIndex<0){
    buildChain(treeMCSig, channel, 0);
    buildChain(treeMCPS, channel, 3);
    buildChain(treeMCBkg, channel, 1);
    cout << "Chain for " << channel << " " << "ALL 3 SAMPLE TYPES " << " " <<( isLowMass? "LowMass":"HighMass") << "  Total entries in tree-Sig-MC=" << treeMCSig->GetEntries() << endl;
    treeMCSig->ls();
  }
  else {
    cout<<"SMAPLEINDEX=="<<sampleIndex<<endl;
    buildChain(treeMC, channel, sampleIndex);
    cout << "Chain for " << channel << " " << (sampleIndex==0? "SIGNAL" : "BACKGROUND") << " " <<( isLowMass? "LowMass":"HighMass") << "  Total entries in tree=" << treeMC->GetEntries() << endl;
    treeMC->ls();
  }
 


 
  float mzz,w,LD,psLD;
   double mela,SMD;
  char yVarName[32];


  //distinction btw LD and mela needed because we might want 
  //both psMELA (for 2D template) and MELA (for cut)
  //  if(makePSTemplate)sprintf(yVarName,"pseudoLD");
  //  else   sprintf(yVarName,"ZZLD");
  sprintf(yVarName,"pseudoLD");
  treeMCSig->SetBranchAddress("ZZMass",&mzz);
  treeMCSig->SetBranchAddress("ZZLD",&mela);
  //  treeMCSig->SetBranchAddress("pseudoLD",&psLD);
  treeMCSig->SetBranchAddress("superLD",&SMD);
  treeMCSig->SetBranchAddress(yVarName,&LD);
  //  treeMC->SetBranchAddress("MC_weight_noxsec",&w);
  treeMCSig->SetBranchAddress("MC_weight",&w);
 
  
  // fill histogram
  TRandom3 *myR=new TRandom3(6811);
  double dumM,dumLD,dumSuperLD,dumPsLD, dumErr=0.0,dumType;
  ofstream of((destDir + "hzzDummyTMP_" + channel+"_" + strLumi + ".txt").Data(),ios::out);


  int storedEvts=0,storedEvtsSM=0,storedEvtsPS=0,storedEvtsBkg=0;
  double totExpEvt=0.0;
  for(int i=0; i<treeMCSig->GetEntries(); i++){
 
    treeMCSig->GetEntry(i);
    psLD=LD;//only if (yVarName,"pseudoLD")
    float rnum=myR->Rndm();
    if(mzz>105&&mzz<140&&SMD>=0.0)totExpEvt+=w;
    if(verbose_&&i%5000==0)cout<<"SMHiggs: Weight= "<<w<<" RandmNmbr="<< rnum <<"  LD="<<LD<<"  PseudoLD="<< psLD<<"  SuperLD="<<SMD<<endl; 
    bool cutPassed= (melaCut>0.0) ? (mela>melaCut) : true;
    if( rnum < w*lumi && cutPassed && myR->Rndm()<(1.0-fracPS)&&mzz>105&&mzz<140&&SMD>=0.0 ){
      dumM=mzz;
      dumLD=LD;
      dumSuperLD=SMD;
      dumPsLD=psLD;
      dumErr=0.0;
      dumType=0;
      cout<<"Signal-MC SMD="<<SMD<<"  DummySMD="<<dumSuperLD<<endl;
      storedEvts++;
      storedEvtsSM++;
      of<<dumM <<"  "<<dumLD<<" "<<dumErr<<" "<<dumSuperLD<<" "<<dumPsLD <<"  "<<dumType<<endl;
    }

  }//end loop on entries
  cout<<"SMHiggs: total expected events for 1.0 inv fb in [105, 140] (no other cuts): "<<totExpEvt<<endl;



  float mzzPS,wPS,psLDPS,LDPS;
  double melaPS,SMDPS;
  double totExpEvtPS=0.0;
  treeMCPS->SetBranchAddress("ZZMass",&mzzPS);
  treeMCPS->SetBranchAddress("ZZLD",&melaPS);
  //  treeMCPS->SetBranchAddress("pseudoLD",&psLDPS);
  treeMCPS->SetBranchAddress("superLD",&SMDPS);
  treeMCPS->SetBranchAddress(yVarName,&LDPS);
  //  treeMC->SetBranchAddress("MC_weight_noxsec",&w);
  treeMCPS->SetBranchAddress("MC_weight",&wPS);

  for(int i=0; i<treeMCPS->GetEntries(); i++){
 
    treeMCPS->GetEntry(i);
    psLDPS=LDPS;//only if (yVarName,"pseudoLD")

    float rnum=myR->Rndm();
    if(verbose_&&i%5000==0)cout<<"PSHiggs: Weight= "<<wPS<<" RandmNmbr="<< rnum <<"  MELA="<<LDPS<<endl; 
    if(mzzPS>105&&mzzPS<140&&SMDPS>=0.0)totExpEvtPS+=wPS;
 
    bool cutPassed= (melaCut>0.0) ? (melaPS>melaCut) : true;
    if( rnum < wPS*lumi/2.0 && cutPassed && myR->Rndm()<(fracPS) &&mzzPS>105&&mzzPS<140&&SMDPS>=0.0){
     dumM=mzzPS;
      dumLD=LDPS;
      dumSuperLD=SMDPS;
      dumPsLD=psLDPS;
      dumErr=0.0; 
      dumType=-1;
   
      storedEvts++;
      storedEvtsPS++;
      of<<dumM <<"  "<<dumLD<<" "<<dumErr<<" "<<dumSuperLD<<" "<<dumPsLD <<"  "<<dumType<<endl;
    }

  }//end loop on entries
  cout<<"PSHiggs: total expected events for 1.0 inv fb in [105, 140] (no other cuts): "<<totExpEvtPS<<endl;


  treeMCBkg->SetBranchAddress("ZZMass",&mzz);
  treeMCBkg->SetBranchAddress("ZZLD",&mela);
  treeMCBkg->SetBranchAddress(yVarName,&LD);
  // treeMCBkg->SetBranchAddress("pseudoLD",&psLD);
  treeMCBkg->SetBranchAddress("superLD",&SMD);
   treeMCBkg->SetBranchAddress("MC_weight",&w);
  //  treeMC->SetBranchAddress("MC_weight_noxsec",&w);
   //  treeMCBkg->SetBranchAddress("MC_weight",&w);
  for(int i=0; i<treeMCBkg->GetEntries(); i++){
 
    treeMCBkg->GetEntry(i);
   psLD=LD;//only if (yVarName,"pseudoLD")
    float rnum=myR->Rndm();
    if(verbose_&&i%5000==0)cout<<"Bkg qqZZ: Weight= "<<w<<" RandmNmbr="<< rnum <<"  MELA="<<LD<<endl; 
    bool cutPassed= (melaCut>0.0) ? (mela>melaCut) : true;
    if( rnum < w*lumi && cutPassed &&mzz>105&&mzz<140&&SMD>=0.0){
      dumM=mzz;
      dumLD=LD;
      dumErr=0.0;
       dumSuperLD=SMD;
      dumPsLD=psLD;
      dumType=1;
      storedEvts++;
      storedEvtsBkg++;
      of<<dumM <<"  "<<dumLD<<" "<<dumErr<<" "<<dumSuperLD<<" "<<dumPsLD <<"  "<<dumType<<endl;
    }
  }
  cout<<"Stored in outptu tree "<<storedEvts<<" events (SM="<<storedEvtsSM<<"  PS="<<storedEvtsPS<<"  Bkg="<<storedEvtsBkg<<endl;
  //  delete treeMC;
  delete myR;

 
    of.close();
 
  
}


//=======================================================================

void makeDummyObs(TString channel="4mu",float lumi=1.0){


  // ========================================
  // SM Higgs template
  fillTemplate(channel, -1,lumi,true);


}

//=======================================================================

void dumptToTree(TString channel="4mu"){
 ifstream inf((destDir + "hzzDummyTMP_"  + channel+"_" + strLumi + ".txt").Data(),ios::out);
 TFile* ftmp = new TFile(destDir + "hzzDummyTMP_" + channel+"_"+ strLumi  + "_withSMD_doubleCBonly.root","RECREATE");
 double dumM,dumLD,dumSMD,dumPsLD, dumErr=0.0,dumType;
 TTree* dummyT = new TTree("data_obs","data_obs");
  dummyT->Branch("CMS_zz4l_mass",&dumM,"CMS_zz4l_mass/D");
  dummyT->Branch("melaLD",&dumLD,"melaLD/D");
  //  dummyT->Branch("pseudoLD",&dumPsLD,"pseudoLD/D");
  dummyT->Branch("CMS_zz4l_KD",&dumPsLD,"CMS_zz4l_KD/D");
  dummyT->Branch("CMS_zz4l_smd",&dumSMD,"CMS_zz4l_smd/D");
  dummyT->Branch("CMS_zz4l_massErr",&dumErr,"CMS_zz4l_massErr/D");
  dummyT->Branch("dummyEvtType",&dumType,"dummyEvtType/D");

  double oldSMD=-1.0;
 while (!inf.eof() ){
   inf>>dumM>>dumLD>>dumErr>>dumSMD>>dumPsLD>>dumType;
   if(dumSMD!=oldSMD){
     dummyT->Fill();
     oldSMD=dumSMD;
   }
 }
 dummyT->Write();
 delete ftmp;
}

void storeLDDistribution(float lumi=1.0){

  makeDummyObs("4mu",lumi);
  makeDummyObs("4e",lumi);
  makeDummyObs("2e2mu",lumi);

 dumptToTree("4mu");
  dumptToTree("4e");
  dumptToTree("2e2mu");

}


void generateDummyObsData(float lumi=1.0) {

  stringstream ssL;
  ssL<<lumi;
  strLumi=ssL.str();

  storeLDDistribution(lumi);
}

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
#include "TPaveText.h"
#include "TText.h"
#include <sstream>
#include <Riostream.h>
#include <vector>
#include <utility>

/*#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TROOT.h>
#include "ZZMatrixElement/MELA/interface/Mela.h"
#endif
*/
TString canvasHypName;
//----------> SET INPUT VARIABLES in Config.h
#include "ConfigSMD.h"
//<----------

//--

void buildChainSingleMass(TChain* bkgMC, TString channel, int sampleIndex=0, int mh=125) ;
double calcInterfRew(TH1 *h,double KD );
void makePlot1D( TH1 *h ,TString label );
void makePlot2D( TH2 *h ,TString label );
TH2F* fillTemplate(TString channel="4mu", int sampleIndex=0,TString superMelaName="superLD",TString templateName="bkgHisto",  bool smooth=false);
TH1F *fillKDhisto(TString channel="4mu", int sampleIndex=0,float mzzLow=0.0,float mzzHigh=99999.0,  bool smooth=false);
void makeTemplate(TString channel="4mu");
void generateTemplatesSMD_forCR(int altSignal_=3,TString destDirTag="0-") ;

//=======================================================================
double calcInterfRew(TH1 *h,double KD ){
  return h->GetBinContent(h->FindBin(KD));
}

//=======================================================================
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
  c1D->SaveAs((destDir+"can_template1D_SuperMELA_"+label+".png").Data());
  delete c1D;
}

//=======================================================================
void makePlot2D( TH2 *h ,TString label ){

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  char yAxisTitle[64];
  if(altSignal<3)sprintf(yAxisTitle,"KD");
  else   if(altSignal==3)sprintf(yAxisTitle,"D_{0-}");
  else   if(altSignal==4)sprintf(yAxisTitle,"D_{2^{+}_{m}(gg)}");
  else   if(altSignal==5)sprintf(yAxisTitle,"D_{0^{+}_{h}}");
  else   if(altSignal==6)sprintf(yAxisTitle,"D_{1^{+}}");
  else   if(altSignal==7)sprintf(yAxisTitle,"D_{1^{-}}");
  else   if(altSignal==8)sprintf(yAxisTitle,"D_{2^{+}_{m}(qq)}");
  else sprintf(yAxisTitle,"MYDummyKD");

  TCanvas *c2D=new TCanvas("c2d",("CANVAS "+label).Data());
  c2D->cd();
  c2D->SetFillColor(0);
  c2D->SetBorderMode(0);
  c2D->SetBorderSize(2);
  c2D->SetTickx(1);
  c2D->SetTicky(1);
  c2D->SetLeftMargin(0.15);
  c2D->SetRightMargin(0.05);
  c2D->SetTopMargin(0.05);
  c2D->SetBottomMargin(0.15);
  c2D->SetFrameFillStyle(0);
  c2D->SetFrameBorderMode(0);
  c2D->SetFrameFillStyle(0);
  c2D->SetFrameBorderMode(0);

  h->SetXTitle("D_{bkg}");
  h->SetYTitle(yAxisTitle);

  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.007);
  h->GetXaxis()->SetLabelSize(0.045);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleOffset(1.4);
  h->GetXaxis()->SetTitleFont(42); 

  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.007);
  h->GetYaxis()->SetLabelSize(0.045);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(1.1);
  h->GetYaxis()->SetTitleFont(42); 

  //h->GetYaxis()->SetTitleOffset(1.15);
  h->SetTitle("");
  h->Draw("col");

  TPaveText *pt = new TPaveText(0.1577181,0.9562937,0.9580537,0.9947552,"brNDC");
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetFillStyle(0);
  pt->SetTextFont(42);
  pt->SetTextSize(0.05);
  TText *text = pt->AddText(0.01,0.5,"CMS Simulation");
  pt->Draw();   


  char canNameStr[256];
  sprintf(canNameStr,"%s/can_template_SMDvs%s_%s.png",destDir.Data(),canvasHypName.Data(),label.Data());
  c2D->SaveAs(canNameStr);
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

  string suffix=".root";
  //string suffix="_withSMD_doubleCBonly_withProbabilities.root";

 if (sampleIndex==-1){
     if(useSqrts==1){
      cout<<"Readign in 7 TeV for Z+X data control regions"<<endl;
    //7TeV
      bkgMC->Add(filePath7TeV + "/CR/HZZ4lTree_DoubleEle_CREEEEosTree.root");
      bkgMC->Add(filePath7TeV + "/CR/HZZ4lTree_DoubleEle_CREEEEssTree.root");
      bkgMC->Add(filePath7TeV + "/CR/HZZ4lTree_DoubleMu_CRMMMMosTree.root");
      bkgMC->Add(filePath7TeV + "/CR/HZZ4lTree_DoubleMu_CRMMMMssTree.root");
      bkgMC->Add(filePath7TeV + "/CR/HZZ4lTree_DoubleOr_CRMMEEosTree.root");
      bkgMC->Add(filePath7TeV + "/CR/HZZ4lTree_DoubleOr_CRMMEEssTree.root");
      bkgMC->Add(filePath7TeV + "/CR/HZZ4lTree_DoubleOr_CREEMMssTree.root");

    }
    else{
      cout<<"Readign in 8 TeV for Z+X data control regions"<<endl;
      //8TeV
      bkgMC->Add(filePath8TeV + "/CR/HZZ4lTree_DoubleEle_CREEEEosTree.root");
      bkgMC->Add(filePath8TeV + "/CR/HZZ4lTree_DoubleEle_CREEEEssTree.root");
      bkgMC->Add(filePath8TeV + "/CR/HZZ4lTree_DoubleMu_CRMMMMosTree.root");
      bkgMC->Add(filePath8TeV + "/CR/HZZ4lTree_DoubleMu_CRMMMMssTree.root");
      bkgMC->Add(filePath8TeV + "/CR/HZZ4lTree_DoubleOr_CRMMEEosTree.root");
      bkgMC->Add(filePath8TeV + "/CR/HZZ4lTree_DoubleOr_CRMMEEssTree.root");
      bkgMC->Add(filePath8TeV + "/CR/HZZ4lTree_DoubleOr_CREEMMssTree.root");
      
      cout<<"Finished to add files to chain: #entries is "<<bkgMC->GetEntries()<<endl;
    }

  }//end if sampleindex==-1
    
}


//=======================================================================

TH2F* fillTemplate(TString channel, int sampleIndex,TString superMelaName,TString templateName,  bool smooth){
  TChain* bkgMC = new TChain("SelectedTree");


  //  buildChain(bkgMC, channel, sampleIndex);
  // bool use8TeV=false;
  buildChainSingleMass(bkgMC, channel, sampleIndex,mH);
  cout << "Chain for " << channel << " " << sampleIndex << "  "<<(useSqrts==2 ? "8TeV" : "7TeV") << " ===>" << bkgMC->GetEntries() << endl;
  // // bkgMC->ls();

  float mzz,KD,KD_forSel,w=0;
  double sKD;
  float  pbkg=0;
  float p0plus,p0plusN, p0minus,p2minimal;
  float psigM4l, pbkgM4l;
  float psigM4l_ScaleUp, pbkgM4l_ScaleUp,psigM4l_ScaleDown, pbkgM4l_ScaleDown,psigM4l_ResUp, pbkgM4l_ResUp;
  float p0plusVA,p0hplusVA,p0minusVA,p1plusVA,p1minusVA,p2minimalVA,p2minimalVA_qq,pbkgVA;
  bkgMC->SetBranchAddress("ZZMass",&mzz);
  bkgMC->SetBranchAddress("MC_weight_noxsec",&w);
   

  /*
  char psigm4l_name[32], pbkgm4l_name[32];
  if(superMelaName.Contains("syst1Up")){
    sprintf(psigm4l_name,"p0plus_m4l_ScaleUp");
    sprintf(pbkgm4l_name,"bkg_m4l_ScaleUp");
  }
  else  if(superMelaName.Contains("syst1Down")){
    sprintf(psigm4l_name,"p0plus_m4l_ScaleDown");
    sprintf(pbkgm4l_name,"bkg_m4l_ScaleDown");
  }
  else  if(superMelaName.Contains("syst2Up")){
    sprintf(psigm4l_name,"p0plus_m4l_ResUp");
    sprintf(pbkgm4l_name,"bkg_m4l_ResUp");
  }
  else{//use central value
    sprintf(psigm4l_name,"p0plus_m4l");
    sprintf(pbkgm4l_name,"bkg_m4l");
  }
  cout<<"PsigM4L="<<psigm4l_name<<"  PbkgM4L="<<pbkgm4l_name<<endl;
  */

  //  bkgMC->SetBranchAddress(psigm4l_name,&psigM4l);
  bkgMC->SetBranchAddress("p0plus_m4l",&psigM4l);
  bkgMC->SetBranchAddress("bkg_m4l",&pbkgM4l);
  bkgMC->SetBranchAddress("p0plus_m4l_ScaleUp",&psigM4l_ScaleUp);
  bkgMC->SetBranchAddress("bkg_m4l_ScaleUp",&pbkgM4l_ScaleUp);
  bkgMC->SetBranchAddress("p0plus_m4l_ScaleDown",&psigM4l_ScaleDown);
  bkgMC->SetBranchAddress("bkg_m4l_ScaleDown",&pbkgM4l_ScaleDown);
  bkgMC->SetBranchAddress("p0plus_m4l_ResUp",&psigM4l_ResUp);
  bkgMC->SetBranchAddress("bkg_m4l_ResUp",&pbkgM4l_ResUp);

  bkgMC->SetBranchAddress("p0plus_melaNorm",&p0plusN);
  bkgMC->SetBranchAddress("p0plus_mela",&p0plus);
  bkgMC->SetBranchAddress("p0minus_mela",&p0minus);
  bkgMC->SetBranchAddress("p2_mela",&p2minimal);
  bkgMC->SetBranchAddress("bkg_mela",&pbkg);

  bkgMC->SetBranchAddress("p0plus_VAJHU",&p0plusVA);
  bkgMC->SetBranchAddress("p0minus_VAJHU",&p0minusVA);
  bkgMC->SetBranchAddress("p0hplus_VAJHU",&p0hplusVA);
  bkgMC->SetBranchAddress("p1plus_VAJHU",&p1plusVA);
  bkgMC->SetBranchAddress("p1_VAJHU",&p1minusVA);
 bkgMC->SetBranchAddress("p2_VAJHU",&p2minimalVA);
 bkgMC->SetBranchAddress("p2qqb_VAJHU",&p2minimalVA_qq);
  bkgMC->SetBranchAddress("bkg_VAMCFMNorm",&pbkgVA);


  //  bkgMC->SetBranchAddress("",&);


  //  bool cutSameVar=false;
  // bkgMC->SetBranchAddress(melaName.Data(),&KD);
  //  if (melaName!=melaCutName) {
  //  bkgMC->SetBranchAddress(melaCutName.Data(),&KD_forSel);
  // }
  // else {
  ///////    KD_forSel=KD;
  //  cutSameVar=true;
  ////// // bkgMC->SetBranchAddress(melaName.Data(),&KD_forSel);
  // }


  //  bkgMC->SetBranchAddress("superLD",&sKD);
  //bkgMC->SetBranchAddress(superMelaName,&sKD);

  
  //  bool usePseudoLDbin=melaName.Contains("pseudo");
  const int nbinsY=(altSignal==3? nbinsYps : nbinsYgrav);
  float binsY[nbinsY+1];
  for(int ib=0;ib<=nbinsY;ib++){
    if(altSignal==3) binsY[ib]=binsYps[ib];
    else binsY[ib]=binsYgrav[ib];
  }
  TH2F* bkgHist = new TH2F(templateName,templateName,nbinsX,binsX,nbinsY,binsY);
  // TH2F* bkgHist = new TH2F(templateName,templateName,50,0.0,1.0,25,0.0,1.0);

  //  const int nBinsFine=100;
  // float xfine[nBinsFine+1];
  // for(int i=0;i<=nBinsFine;i++)xfine[i]=i*(1.0/nBinsFine);
   //  TH2F* bkgHist = new TH2F(templateName,templateName,nBinsFine,xfine,nbinsY,binsY);

  bkgHist->Sumw2();


  // fill histogram
	
  for(int i=0; i<bkgMC->GetEntries(); i++){

    bkgMC->GetEntry(i);

    //calculate discriminants from individual probabilities (ANALYTICAL APPROACH):
    KD= p0plusN / (p0plusN+pbkg); // traditional MELA analytical
    //  sKD = p0plusN*psigM4l / (p0plusN*psigM4l + pbkg*pbkgM4l)  ;
    // float pseudoKD = p0plus/ (p0plus + p0minus);
    //float graviKD =  p0plus   / (p0plus + 1.2*p2minimal);


    // // using VA

    sKD = double(p0plusVA)*double(psigM4l)    / (double(p0plusVA)*double(psigM4l) + double(pbkgVA)*double(pbkgM4l))  ;
    float sKD_ScaleUp= double(p0plusVA)*psigM4l_ScaleUp    / (double(p0plusVA)*psigM4l_ScaleUp + double(pbkgVA)*pbkgM4l_ScaleUp)  ;
    float sKD_ScaleDown= double(p0plusVA)*psigM4l_ScaleDown    / (double(p0plusVA)*psigM4l_ScaleDown + double(pbkgVA)*pbkgM4l_ScaleDown)  ;
    float sKD_ResUp= double(p0plusVA)*psigM4l_ResUp    / (double(p0plusVA)*psigM4l_ResUp + double(pbkgVA)*pbkgM4l_ResUp)  ;
    if(superMelaName=="superLD_syst1Up")sKD = sKD_ScaleUp;
    if(superMelaName=="superLD_syst1Down")sKD = sKD_ScaleDown;
    if(superMelaName=="superLD_syst2Up")sKD = sKD_ResUp;
    

    float pseudoKD = p0plusVA / (p0plusVA   + p0minusVA);
    float p0hKD = p0plusVA / (p0plusVA   + p0hplusVA);
    float p1plusKD = p0plusVA / (p0plusVA   + p1plusVA);
    float p1minusKD = p0plusVA / (p0plusVA   + p1minusVA);
    float graviKD =  p0plusVA   / (p0plusVA + p2minimalVA);
    float qqgraviKD =  p0plusVA   / (p0plusVA + p2minimalVA_qq);
    if(altSignal==3) KD=pseudoKD;
    else if(altSignal==4)KD=graviKD;
    else if(altSignal==5)KD=p0hKD;
    else if(altSignal==6)KD=p1plusKD;
    else if(altSignal==7)KD=p1minusKD;
    else if(altSignal==8)KD=qqgraviKD;
    

    //    if(cutSameVar)
    KD_forSel=KD;

    bool cutPassed= (kdCut>0.0) ? (KD_forSel>kdCut) : true;
    if( cutPassed &&sKD>=0.0&& mzz>mzzCutLow&&mzz<mzzCutHigh){

      bkgHist->Fill(sKD,KD,w);
      //   bkgHist->Fill(mzz,KD,w);

  //    if(channel=="4mu"&&sampleIndex==0&&useSqrts==2&&superMelaName=="superLD"){
      // cout<<"MZZ="<<mzz<<" pseudoKD="<<KD<<"  SuperMELA="<<sKD <<"  PSigm4l="<<psigM4l<<"  PBkgm4l="<<pbkgM4l<<"   PSigVA="<< p0plusVA<<"  PBkgVA="<<pbkgVA <<endl;
      //}

    }//end if cuts passed

    // if(i%5000==0) cout << "event: " << i << "/" << bkgMC->GetEntries() << endl;
  

  }//end loop on entries


// smooth 

  if(smooth)   bkgHist->Smooth(1,"k5b"); //options:  "k3a", "k5a" , "k5b" 
  
  // normalize TH2
  double totArea=bkgHist->Integral();
  bkgHist->Scale(1.0/totArea);

  //Old Floor
  /*
  // bkgHist->Smooth();
  for(int i=1; i<=bkgHist->GetNbinsX(); i++){
    for(int j=1; j<=bkgHist->GetNbinsY(); j++){
      if(bkgHist->GetBinContent(i,j)<0.00000001)
	bkgHist->SetBinContent(i,j,0.00000001);
    }// for(int j=1; j<=nYbins; j++){
  }// for(int i=1; i<=nXbins; i++){
  */

  //New Floor
  double floor = ((bkgHist->Integral())/(bkgHist->GetNbinsX()*bkgHist->GetNbinsY()))*(0.1/100);
  for(int i = 1; i <= bkgHist->GetNbinsX(); i++)
    {
      for(int j = 1; j <= bkgHist->GetNbinsY(); j++)
	{
	  double orig = bkgHist->GetBinContent(i,j);
	  bkgHist->SetBinContent(i,j,(orig+floor));
	}
    }

  // normalize TH2
  totArea=bkgHist->Integral();
  bkgHist->Scale(1.0/totArea);
  
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


  cout<<"Finishing fillTemplate for sample "<<sampleIndex<<"   "<<channel.Data()<<endl;
  return bkgHist;
  
}


//===================================================
TH1F *fillKDhisto(TString channel, int sampleIndex,float mzzLow,float mzzHigh,bool smooth){


  TChain* bkgMC = new TChain("SelectedTree");
  buildChainSingleMass(bkgMC, channel, sampleIndex,mH);
  float mzz=-1.0,mela=-444.0,w=0;
  double KD=-999., sKD=-999.;
  //char yVarName[32];

  float psig=0, pbkg=0;
  float p0plus, p0minus,p2minimal, psigM4l, pbkgM4l;
 float p0plusVA,pbkgVA;
  bkgMC->SetBranchAddress("ZZMass",&mzz);
  bkgMC->SetBranchAddress("MC_weight_noxsec",&w);
  bkgMC->SetBranchAddress("ZZLD",&mela);
  bkgMC->SetBranchAddress("superLD",&KD);

  bkgMC->SetBranchAddress("p0plus_m4l",&psigM4l);
  bkgMC->SetBranchAddress("bkg_m4l",&pbkgM4l);

  bkgMC->SetBranchAddress("p0plus_mela",&p0plus);
  bkgMC->SetBranchAddress("p0minus_mela",&p0minus);
  // bkgMC->SetBranchAddress("p2_mela",&p2minimal);
  bkgMC->SetBranchAddress("bkg_mela",&pbkg);
  //  bkgMC->SetBranchAddress("",&);
  bkgMC->SetBranchAddress("p0plus_VAJHU",&p0plusVA);
  bkgMC->SetBranchAddress("bkg_VAMCFMNorm",&pbkgVA);
  bkgMC->SetBranchAddress("p2_mela",&p2minimal);
  char hTitle[128];
  sprintf(hTitle,"Distribution of superMELA KD with M_{4l} in [%d, %d]",int(mzzLow),int(mzzHigh));
  TH1F* outHist=new TH1F("finHisto",hTitle,200,0.0,1.0);
  //TH1F* outHist=new TH1F("finHisto",hTitle,nbinsX,binsX);
  outHist->Sumw2();
 // fill histogram
 // cout<<"Looping on tree netries (fillKDHisto) "<<bkgMC->GetEntries()<<endl;
  for(int i=0; i<bkgMC->GetEntries(); i++){
    bkgMC->GetEntry(i);

    //calculate discriminants from individual probabilities
    KD= p0plus / (p0plus+pbkg);
    //    sKD = p0plus*psigM4l / (p0plus*psigM4l + pbkg*pbkgM4l)  ;
    sKD = double(p0plusVA)*double(psigM4l)    / (double(p0plusVA)*double(psigM4l) + double(pbkgVA)*double(pbkgM4l))  ;
    //float pseudoKD = p0plus / (p0plus + 6.0*p0minus);
    //float graviKD =  p0plus / (p0plus + 1.2*p2minimal);
    //if(altSignal==3) KD=pseudoKD;
    //if(altSignal==4)KD=graviKD;

    //   if(i%200==0)cout<<"Entry "<<i<<"  mzz="<<mzz<<"  weight="<<w<<"  LD="<<LD<<endl;
    bool cutPassed= (kdCut>0.0) ? (mela>kdCut) : true;
    if( cutPassed &&mzz>mzzLow&&mzz<mzzHigh){
      outHist->Fill(sKD,w);
    }
  }


// smooth 

    for(int j=1; j<=outHist->GetNbinsX(); j++){
      if(outHist->GetBinContent(j)<0.000001)
	outHist->SetBinContent(j,0.000001);
    }// for(int j=1; j<=nYbins; j++){
    
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
  
  low = fillTemplate(channel,-1,"superLD","zjetsHist",true);
  //  high = fillTemplate(channel,1,false);
  h_mzzD =   (TH2F*)low->Clone("h_mzzD");// mergeTemplates(low,high);
  h_D = fillKDhisto(channel,-1,mzzCutLow,mzzCutHigh,true);//last two are cuts on mZZ
 
  // ---------- apply interference reweighting --------
  
  oldTemp = new TH2F(*h_mzzD);
  oldTemp->SetName("oldTemp");

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

void generateTemplatesSMD_forCR(int altSignal_, TString destDirTag){

  canvasHypName=destDirTag;

  altSignal = altSignal_;

  useSqrts=1;
  destDir=destDirBase+"_"+destDirTag+"_7TeV/";  
    
  makeTemplate("4mu");

  useSqrts=2;
  destDir=destDirBase+"_"+destDirTag+"_8TeV/";  
    
  makeTemplate("4mu");

}



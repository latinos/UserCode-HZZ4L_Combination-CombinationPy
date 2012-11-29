/* 
 * Fit qqZZ background shapes and write parameters in a card fragment.
 * usage: 
 * -set all input variables in Config.h
 * -run with:
 * root -q -b backgroundFits_qqzz_1Dw.C
 * This runs on the 3 final states for 7 and 8 TeV and writes the output in a file (see stdout).
 *
 */

/*
  #ifndef __CINT__
  #include "RooGlobalFunc.h"
  #endif
  #include "RooRealVar.h"
  #include "RooDataSet.h"
  #include "RooGaussian.h"
  #include "RooConstVar.h"
  #include "RooChebychev.h"
  #include "RooAddPdf.h"
  #include "RooWorkspace.h"
  #include "RooPlot.h"
  #include "TCanvas.h"
  #include "TAxis.h"
  #include "TFile.h"
  #include "TH1.h"
*/


#include <iostream>
#include <iomanip>

using namespace RooFit ;


//----------> SET INPUT VARIABLES in Config.h
#include "Config.h"
//<----------

void backgroundFits_qqzz_1Dw(int channel, int sqrts, bool VBFtag);

// Run all final states and sqrts in one go
void backgroundFits_qqzz_1Dw() {

  gSystem->Exec("mkdir -p bkgFigs7TeV");
  gSystem->Exec("mkdir -p bkgFigs8TeV");

  backgroundFits_qqzz_1Dw(1,7,true);
  backgroundFits_qqzz_1Dw(2,7,true);
  backgroundFits_qqzz_1Dw(3,7,true);  
  backgroundFits_qqzz_1Dw(1,8,true);
  backgroundFits_qqzz_1Dw(2,8,true);
  backgroundFits_qqzz_1Dw(3,8,true);

  backgroundFits_qqzz_1Dw(1,7,false);
  backgroundFits_qqzz_1Dw(2,7,false);
  backgroundFits_qqzz_1Dw(3,7,false);  
  backgroundFits_qqzz_1Dw(1,8,false);
  backgroundFits_qqzz_1Dw(2,8,false);
  backgroundFits_qqzz_1Dw(3,8,false);
}

// The actual job
void backgroundFits_qqzz_1Dw(int channel, int sqrts, bool VBFtag)
{
  TString schannel;
  if      (channel == 1) schannel = "4mu";
  else if (channel == 2) schannel = "4e";
  else if (channel == 3) schannel = "2e2mu";
  else cout << "Not a valid channel: " << schannel << endl;

  TString ssqrts = (long) sqrts + TString("TeV");

  cout << "schannel = " << schannel << "  sqrts = " << sqrts << " VBFtag = " << VBFtag << endl;

  TString outfile = "CardFragments/qqzzBackgroundFit_" + ssqrts + "_" + schannel + "_" + Form("%d",int(VBFtag)) + ".txt";
  ofstream of(outfile,ios_base::out);
  of << "### background functions ###" << endl;


  gSystem->AddIncludePath("-I$ROOFITSYS/include");
  gROOT->ProcessLine(".L ../CreateDatacards/include/tdrstyle.cc");
  setTDRStyle(false);
  gStyle->SetPadLeftMargin(0.16);

  TString filepath;
  if (sqrts==7) {
    filepath = filePath7TeV;
  } else if (sqrts==8) {
    filepath = filePath8TeV;
  }

  TChain* tree = new TChain("SelectedTree");
  tree->Add( filepath+ "/" + (schannel=="2e2mu"?"2mu2e":schannel) + "/HZZ4lTree_ZZTo*.root");


  RooRealVar* MC_weight = new RooRealVar("MC_weight","MC_weight",0.,2.) ; 
  RooRealVar* ZZMass = new RooRealVar("ZZMass","ZZMass",100.,1000.);
  RooRealVar* ZZLD = new RooRealVar("ZZLD","ZZLD",0.,1.);
  RooRealVar* NJets = new RooRealVar("NJets","NJets",0.,100.);
  
  char cut[10];
  sprintf(cut,"ZZLD>0.5");
  //RooDataSet* set = new RooDataSet("set","set",tree,RooArgSet(*ZZMass,*MC_weight,*ZZLD),cut,"MC_weight");
  if (VBFtag==true){
    RooDataSet* set = new RooDataSet("set","set",tree,RooArgSet(*ZZMass,*MC_weight,*NJets),"NJets>1","MC_weight");
  }
  else{
    RooDataSet* set = new RooDataSet("set","set",tree,RooArgSet(*ZZMass,*MC_weight,*NJets),"NJets<2","MC_weight");
  }

  double totalweight = 0.;
  double totalweight_z = 0.;
  for (int i=0 ; i<set->numEntries() ; i++) { 
    //set->get(i) ; 
    RooArgSet* row = set->get(i) ;
    //row->Print("v");
    totalweight += set->weight();
    if (row->getRealValue("ZZMass") < 200) totalweight_z += set->weight();
  } 
  cout << "nEntries: " << set->numEntries() << ", totalweight: " << totalweight << ", totalweight_z: " << totalweight_z << endl;
	
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");
	
  //// ---------------------------------------
  //Background
  RooRealVar CMS_qqzzbkg_a0("CMS_qqzzbkg_a0","CMS_qqzzbkg_a0",115.3,0.,200.);
  RooRealVar CMS_qqzzbkg_a1("CMS_qqzzbkg_a1","CMS_qqzzbkg_a1",21.96,0.,200.);
  RooRealVar CMS_qqzzbkg_a2("CMS_qqzzbkg_a2","CMS_qqzzbkg_a2",122.8,0.,200.);
  RooRealVar CMS_qqzzbkg_a3("CMS_qqzzbkg_a3","CMS_qqzzbkg_a3",0.03479,0.,1.);
  RooRealVar CMS_qqzzbkg_a4("CMS_qqzzbkg_a4","CMS_qqzzbkg_a4",185.5,0.,200.);
  RooRealVar CMS_qqzzbkg_a5("CMS_qqzzbkg_a5","CMS_qqzzbkg_a5",12.67,0.,200.);
  RooRealVar CMS_qqzzbkg_a6("CMS_qqzzbkg_a6","CMS_qqzzbkg_a6",34.81,0.,100.);
  RooRealVar CMS_qqzzbkg_a7("CMS_qqzzbkg_a7","CMS_qqzzbkg_a7",0.1393,0.,1.);
  RooRealVar CMS_qqzzbkg_a8("CMS_qqzzbkg_a8","CMS_qqzzbkg_a8",66.,0.,200.);
  RooRealVar CMS_qqzzbkg_a9("CMS_qqzzbkg_a9","CMS_qqzzbkg_a9",0.07191,0.,1.);
  RooRealVar CMS_qqzzbkg_a10("CMS_qqzzbkg_a10","CMS_qqzzbkg_a10",94.11,0.,200.);
  RooRealVar CMS_qqzzbkg_a11("CMS_qqzzbkg_a11","CMS_qqzzbkg_a11",-5.111,-100.,100.);
  RooRealVar CMS_qqzzbkg_a12("CMS_qqzzbkg_a12","CMS_qqzzbkg_a12",4834,0.,10000.);
  RooRealVar CMS_qqzzbkg_a13("CMS_qqzzbkg_a13","CMS_qqzzbkg_a13",0.2543,0.,1.);
	
  if (channel == 1){
    ///* 4mu
    CMS_qqzzbkg_a0.setVal(103.854);
    CMS_qqzzbkg_a1.setVal(10.0718);
    CMS_qqzzbkg_a2.setVal(117.551);
    CMS_qqzzbkg_a3.setVal(0.0450287);
    CMS_qqzzbkg_a4.setVal(185.262);
    CMS_qqzzbkg_a5.setVal(7.99428);
    CMS_qqzzbkg_a6.setVal(39.7813);
    CMS_qqzzbkg_a7.setVal(0.0986891);
    CMS_qqzzbkg_a8.setVal(49.1325);
    CMS_qqzzbkg_a9.setVal(0.0389984);
    CMS_qqzzbkg_a10.setVal(98.6645);
    CMS_qqzzbkg_a11.setVal(-7.02043);
    CMS_qqzzbkg_a12.setVal(5694.66);
    CMS_qqzzbkg_a13.setVal(0.0774525);
    //*/
  }
  else if (channel == 2){
    ///* 4e
    CMS_qqzzbkg_a0.setVal(111.165);
    CMS_qqzzbkg_a1.setVal(19.8178);
    CMS_qqzzbkg_a2.setVal(120.89);
    CMS_qqzzbkg_a3.setVal(0.0546639);
    CMS_qqzzbkg_a4.setVal(184.878);
    CMS_qqzzbkg_a5.setVal(11.7041);
    CMS_qqzzbkg_a6.setVal(33.2659);
    CMS_qqzzbkg_a7.setVal(0.140858);
    CMS_qqzzbkg_a8.setVal(56.1226);
    CMS_qqzzbkg_a9.setVal(0.0957699);
    CMS_qqzzbkg_a10.setVal(98.3662);
    CMS_qqzzbkg_a11.setVal(-6.98701);
    CMS_qqzzbkg_a12.setVal(10.0536);
    CMS_qqzzbkg_a13.setVal(0.110576);
    //*/
  }
  else if (channel == 3){
    ///* 2e2mu
    CMS_qqzzbkg_a0.setVal(110.293);
    CMS_qqzzbkg_a1.setVal(11.8334);
    CMS_qqzzbkg_a2.setVal(116.91);
    CMS_qqzzbkg_a3.setVal(0.0433151);
    CMS_qqzzbkg_a4.setVal(185.817);
    CMS_qqzzbkg_a5.setVal(10.5945);
    CMS_qqzzbkg_a6.setVal(29.6208);
    CMS_qqzzbkg_a7.setVal(0.0826);
    CMS_qqzzbkg_a8.setVal(53.1346);
    CMS_qqzzbkg_a9.setVal(0.0882081);
    CMS_qqzzbkg_a10.setVal(85.3776);
    CMS_qqzzbkg_a11.setVal(-13.3836);
    CMS_qqzzbkg_a12.setVal(7587.95);
    CMS_qqzzbkg_a13.setVal(0.325621);
    //*/
  }
  else {
    cout << "disaster" << endl;
  }
    
  RooqqZZPdf_v2* bkg_qqzz = new RooqqZZPdf_v2("bkg_qqzz","bkg_qqzz",*ZZMass,
					      CMS_qqzzbkg_a0,CMS_qqzzbkg_a1,CMS_qqzzbkg_a2,CMS_qqzzbkg_a3,CMS_qqzzbkg_a4,
					      CMS_qqzzbkg_a5,CMS_qqzzbkg_a6,CMS_qqzzbkg_a7,CMS_qqzzbkg_a8,
					      CMS_qqzzbkg_a9,CMS_qqzzbkg_a10,CMS_qqzzbkg_a11,CMS_qqzzbkg_a12,CMS_qqzzbkg_a13);
  RooArgSet myASet(*ZZMass, CMS_qqzzbkg_a0,CMS_qqzzbkg_a1,CMS_qqzzbkg_a2,CMS_qqzzbkg_a3,CMS_qqzzbkg_a4,
		   CMS_qqzzbkg_a5,CMS_qqzzbkg_a6,CMS_qqzzbkg_a7);
  myASet.add(CMS_qqzzbkg_a8);
  myASet.add(CMS_qqzzbkg_a9);
  myASet.add(CMS_qqzzbkg_a10);
  myASet.add(CMS_qqzzbkg_a11);
  myASet.add(CMS_qqzzbkg_a12);
  myASet.add(CMS_qqzzbkg_a13);
 
  RooFitResult *r1 = bkg_qqzz->fitTo( *set, Save(kTRUE), SumW2Error(kTRUE) );//, Save(kTRUE), SumW2Error(kTRUE)) ;

  cout << endl;
  cout << "------- Parameters for " << schannel << " sqrts=" << sqrts << endl;
  cout << "  a0_bkgd = " << CMS_qqzzbkg_a0.getVal() << endl;
  cout << "  a1_bkgd = " << CMS_qqzzbkg_a1.getVal() << endl;
  cout << "  a2_bkgd = " << CMS_qqzzbkg_a2.getVal() << endl;
  cout << "  a3_bkgd = " << CMS_qqzzbkg_a3.getVal() << endl;
  cout << "  a4_bkgd = " << CMS_qqzzbkg_a4.getVal() << endl;
  cout << "  a5_bkgd = " << CMS_qqzzbkg_a5.getVal() << endl;
  cout << "  a6_bkgd = " << CMS_qqzzbkg_a6.getVal() << endl;
  cout << "  a7_bkgd = " << CMS_qqzzbkg_a7.getVal() << endl;
  cout << "  a8_bkgd = " << CMS_qqzzbkg_a8.getVal() << endl;
  cout << "  a9_bkgd = " << CMS_qqzzbkg_a9.getVal() << endl;
  cout << "  a10_bkgd = " << CMS_qqzzbkg_a10.getVal() << endl;
  cout << "  a11_bkgd = " << CMS_qqzzbkg_a11.getVal() << endl;
  cout << "  a12_bkgd = " << CMS_qqzzbkg_a12.getVal() << endl;
  cout << "  a13_bkgd = " << CMS_qqzzbkg_a13.getVal() << endl;
  cout << "}" << endl;
  cout << "---------------------------" << endl;


  of << "qqZZshape a0_bkgd   " << CMS_qqzzbkg_a0.getVal() << endl;
  of << "qqZZshape a1_bkgd   " << CMS_qqzzbkg_a1.getVal() << endl;
  of << "qqZZshape a2_bkgd   " << CMS_qqzzbkg_a2.getVal() << endl;
  of << "qqZZshape a3_bkgd   " << CMS_qqzzbkg_a3.getVal() << endl;
  of << "qqZZshape a4_bkgd   " << CMS_qqzzbkg_a4.getVal() << endl;
  of << "qqZZshape a5_bkgd   " << CMS_qqzzbkg_a5.getVal() << endl;
  of << "qqZZshape a6_bkgd   " << CMS_qqzzbkg_a6.getVal() << endl;
  of << "qqZZshape a7_bkgd   " << CMS_qqzzbkg_a7.getVal() << endl;
  of << "qqZZshape a8_bkgd   " << CMS_qqzzbkg_a8.getVal() << endl;
  of << "qqZZshape a9_bkgd   " << CMS_qqzzbkg_a9.getVal() << endl;
  of << "qqZZshape a10_bkgd  " << CMS_qqzzbkg_a10.getVal() << endl;
  of << "qqZZshape a11_bkgd  " << CMS_qqzzbkg_a11.getVal() << endl;
  of << "qqZZshape a12_bkgd  " << CMS_qqzzbkg_a12.getVal() << endl;
  of << "qqZZshape a13_bkgd  " << CMS_qqzzbkg_a13.getVal() << endl;
  of << endl << endl;
  of.close();

  cout << endl << "Output written to: " << outfile << endl;
  
    
  double qqzznorm;
  if (channel == 1) qqzznorm = 20.5836;
  else if (channel == 2) qqzznorm = 13.8871;
  else if (channel == 3) qqzznorm = 32.9883;
  else { cout << "disaster!" << endl; }

  ZZMass->setRange("fullrange",100.,1000.);
  ZZMass->setRange("largerange",100.,600.);
  ZZMass->setRange("zoomrange",100.,200.);
    
  double rescale = qqzznorm/totalweight;
  double rescale_z = qqzznorm/totalweight_z;
  cout << "rescale: " << rescale << ", rescale_z: " << rescale_z << endl;


  // Plot m4l and
  RooPlot* frameM4l = ZZMass->frame(Title("M4L"),Range(100,600),Bins(250)) ;
  set->plotOn(frameM4l, MarkerStyle(20), Rescale(rescale)) ;
  
  //set->plotOn(frameM4l) ;
  RooPlot* frameM4lz = ZZMass->frame(Title("M4L"),Range(100,200),Bins(100)) ;
  set->plotOn(frameM4lz, MarkerStyle(20), Rescale(rescale)) ;


  int iLineColor = 1;
  string lab = "blah";
  if (channel == 1) { iLineColor = 2; lab = "4#mu"; }
  if (channel == 3) { iLineColor = 4; lab = "2e2#mu"; }
  if (channel == 2) { iLineColor = 6; lab = "4e"; }

  bkg_qqzz->plotOn(frameM4l,LineColor(iLineColor),NormRange("largerange")) ;
  bkg_qqzz->plotOn(frameM4lz,LineColor(iLineColor),NormRange("zoomrange")) ;
    
//second shape to compare with (if previous comparison code unceommented)
  //bkg_qqzz_bkgd->plotOn(frameM4l,LineColor(1),NormRange("largerange")) ;
  //bkg_qqzz_bkgd->plotOn(frameM4lz,LineColor(1),NormRange("zoomrange")) ;
    
  
  double normalizationBackground_qqzz = bkg_qqzz->createIntegral( RooArgSet(*ZZMass), Range("fullrange") )->getVal();
  cout << "Norm all = " << normalizationBackground_qqzz << endl;
    
  frameM4l->GetXaxis()->SetTitle("m_{4l} [GeV]");
  frameM4l->GetYaxis()->SetTitle("a.u.");
  frameM4lz->GetXaxis()->SetTitle("m_{4l} [GeV]");
  frameM4lz->GetYaxis()->SetTitle("a.u.");

  char lname[192];
  sprintf(lname,"qq #rightarrow ZZ #rightarrow %s", lab.c_str() );
  char lname2[192];
  sprintf(lname2,"Shape Model, %s", lab.c_str() );
  // dummy!
  TF1* dummyF = new TF1("dummyF","1",0.,1.);
  TH1F* dummyH = new TH1F("dummyH","",1, 0.,1.);
  dummyF->SetLineColor( iLineColor );
  dummyF->SetLineWidth( 2 );

  dummyH->SetLineColor( kBlue );
  TLegend * box2 = new TLegend(0.4,0.70,0.80,0.90);
  box2->SetFillColor(0);
  box2->SetBorderSize(0);
  box2->AddEntry(dummyH,"Simulation (POWHEG+Pythia)  ","pe");
  box2->AddEntry(dummyH,lname,"");
  box2->AddEntry(dummyH,"","");
  box2->AddEntry(dummyF,lname2,"l");
    
  TPaveText *pt = new TPaveText(0.15,0.955,0.4,0.99,"NDC");
  pt->SetFillColor(0);
  pt->SetBorderSize(0);
  pt->AddText("CMS Preliminary 2012");
  TPaveText *pt2 = new TPaveText(0.84,0.955,0.99,0.99,"NDC");
  pt2->SetFillColor(0);
  pt2->SetBorderSize(0);
  pt2->AddText("#sqrt{s} = 7 TeV"); 		

  TCanvas *c = new TCanvas("c","c",800,600);
  c->cd();
  frameM4l->Draw();
  frameM4l->GetYaxis()->SetRangeUser(0,0.4);
  if(channel == 3)frameM4l->GetYaxis()->SetRangeUser(0,0.6);
  box2->Draw();
  pt->Draw();
  pt2->Draw();
  TString outputPath = "bkgFigs";
  outputPath = outputPath+ (long) sqrts + "TeV/";
  outputName =  outputPath + "bkgqqzz_" + schannel + "_" + Form("%d",int(VBFtag));
  c->SaveAs(outputName + ".eps");
  c->SaveAs(outputName + ".png");
    
  TCanvas *c2 = new TCanvas("c2","c2",1000,500);
  c2->Divide(2,1);
  c2->cd(1);
  frameM4l->Draw();
  box2->Draw("same");
  c2->cd(2);
  frameM4lz->Draw();
  box2->Draw("same");
  
  outputName = outputPath + "bkgqqzz_" + schannel + "_z" + "_" + Form("%d",int(VBFtag));
  c2->SaveAs(outputName + ".eps");
  c2->SaveAs(outputName + ".png");

  /* TO make the ratio btw 2 shapes, if needed for compairson
  TCanvas *c3 = new TCanvas("c3","c3",1000,500);
   if(sqrts==7)
    sprintf(outputName, "bkgFigs7TeV/bkgqqzz_%s_ratio.eps",schannel.c_str());
  else if(sqrts==8)
    sprintf(outputName, "bkgFigs8TeV/bkgqqzz_%s_ratio.eps",schannel.c_str());

   const int nPoints = 501.;
  double masses[nPoints] ;
  int j=0;
  for (int i=100; i<601; i++){
    masses[j] = i;
    j++;
  }
  cout<<j<<endl;
  double effDiff[nPoints];
  for (int i = 0; i < nPoints; i++){
    ZZMass->setVal(masses[i]);
    double eval = (bkg_qqzz_bkgd->getVal(otherASet)-bkg_qqzz->getVal(myASet))/(bkg_qqzz->getVal(myASet));
    //cout<<bkg_qqzz_bkgd->getVal(otherASet)<<" "<<bkg_qqzz->getVal(myASet)<<" "<<eval<<endl;
    effDiff[i]=eval;
  }
  TGraph* grEffDiff = new TGraph( nPoints, masses, effDiff );
  grEffDiff->SetMarkerStyle(20);
  grEffDiff->Draw("AL");

  //c3->SaveAs(outputName);
  */

  outputName = outputPath + "bkgqqzz_" + schannel + "_z" + "_" + Form("%d",int(VBFtag)) + ".root";
  TFile* outF = new TFile(outputName,"RECREATE");
  outF->cd();
  c2->Write();
  frameM4l->Write();
  frameM4lz->Write();	
  outF->Close();


  delete c;
  delete c2;
}


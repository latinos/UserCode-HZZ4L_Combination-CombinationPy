#include <Riostream.h>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#include "TArrow.h"
#include "TGaxis.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TArrow.h"
#include "TColor.h"
#include "TROOT.h"
#include "Math/DistFunc.h"

void setTDRStyle();

int extractSignificanceStats(bool unblind=false, TString legALT="0^{-}", TString nameALT="0m"){

  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  setTDRStyle();


  const float lumi7TeV=5.051;
  const float lumi8TeV=12.21;

  //char fileName[128];
  //sprintf(fileName,"qmu.root");
  //TFile *fq=new TFile(fileName,"READ");
  //TTree *t=(TTree*)fq->Get("q");
  TChain* t = new TChain("q");
  t->Add("fixedMu/qmu_*.root");
  //t->Add("prefitMu/qmu_*.root");

  float q,m,w;
  int type;
  t->SetBranchAddress("q",&q);
  t->SetBranchAddress("mh",&m);
  t->SetBranchAddress("weight",&w);
  t->SetBranchAddress("type",&type);


  TH1F *hSM=new TH1F("hSM;S = -2 #times ln(L_{PS}/L_{SM});Number of Toys","",8000,-40,40);
  TH1F *hPS=new TH1F("hPS;S = -2 #times ln(L_{PS}/L_{SM});Number of Toys","",8000,-40,40);
  TH1F *hObs=new TH1F("hObserved","",8000,-40,40);

  std::vector<float> v_SM, v_PS,v_Obs;

  for(int i=0;i<t->GetEntries();i++){
    t->GetEntry(i);

    if(i==0)cout<<"MASS in the TREE = "<<m<<endl<<endl;

    q*=2.0;
    if(type<0){ //SM hypothesis 
      hSM->Fill(-q);
      v_SM.push_back(-q);
    }
    else if(type>0){//ALT hypothesis
      hPS->Fill(-q);
      v_PS.push_back(-q);
    }
    else{
      hObs->Fill(-q);
      cout<<"DATA -q ="<<(-q)<<endl;
      v_Obs.push_back(-q);
    }
  }//end loop on tree entries
  cout<<"Finished to loop, sorting vectors "<<v_SM.size()<<" "<<v_PS.size()<<" "<<v_Obs.size()<<endl;
  sort(v_SM.begin(),v_SM.end());//sort in ascending order
  sort(v_PS.begin(),v_PS.end()); 
  sort(v_Obs.begin(),v_Obs.end());
  int ntoysSM= hSM->GetEntries();
  int ntoysPS= hPS->GetEntries();

  //we assume that SM is on the right and PS on the left of zero
  if(v_PS.at(0)>v_SM.at(ntoysSM-1)){
    cout<<"Swapped distributions !!! The alternative model shouldstay on the negative side of the significance."<<endl;
    cout<<"Please edit the code and change the sign of q when filling histos and vectors in the loop on tree entries"<<endl;
    return 1;
  }

  if((int(v_SM.size())!= ntoysSM)||(int(v_PS.size())!= ntoysPS)){
    cout<<"Mismatch in size of vectors and #entries of historgams ! v_SM.size()="<< v_SM.size() <<"  ntoysSM="<<ntoysSM<<endl;
    return 1;
  }

  float medianSM=v_SM.at(int(ntoysSM/2));
  float medianPS=v_PS.at(int(ntoysPS/2));
  cout<<"Toys generated "<<ntoysSM<<"\t"<<ntoysPS<<endl;
  cout<<"Mean of SM/PS hypothesis: "<<hSM->GetMean()<<"\t"<<hPS->GetMean()<<endl;
  cout<<"RMS  of SM/PS hypothesis: "<<hSM->GetRMS()<<"\t"<<hPS->GetRMS()<<endl;
  cout<<"Median of SM/PS hypothesis: "<<medianSM<<"\t"<<medianPS<<endl;

  const float step=0.05;
  float coverage=0.0;
  float diff=10.0;
  float cut=v_PS.at(0)-step;
  float crosspoint=-99.0;
  int startSM=ntoysSM-1, startPS=0;
  cout<<"Starting to loop with cut at "<<cut<<endl;

  //apply a cut on the vectors with the results of toys,
  //for each cut check the area in the tail for PS and SM
  //and calculate the difference.
  //Find the value of cut that minimizes the difference.
 

  while(cut<=v_SM.at(ntoysSM-1)+step){
    //    if(int(cut*100)%100==0)
    // cout<<"Cutting at "<<cut<<endl;
    float cutSM=-1.0,cutPS=-1.0;

    for(int iSM=startSM;iSM>=0;iSM--){      
      //entries in v_SM and v_PS are sorted
      if(v_SM.at(iSM)<cut){//gotcha
	cutSM=ntoysSM-iSM;
	break;
      }
      //      //      else cout<<"SM "<<v_SM.at(iSM)<<" > "<<cut<<endl;
    }

    for(int iPS=startPS;iPS<ntoysPS;iPS++){
      if(v_PS.at(iPS)>cut){//gotcha
	cutPS=iPS;
	break;
      }
      ////      else cout<<v_PS.at(iPS)<<" < "<<cut<<endl;
  
    }

    if(cutSM>=0&&cutPS>=0){
      float fracSM=(ntoysSM-cutSM)/ntoysSM;
      float fracPS=(ntoysPS-cutPS)/ntoysPS;
      // //   cout<<"FracSM: "<<fracSM<<"   FracPS: "<<fracPS<<endl;
      if(fabs(fracSM-fracPS)<diff){
	diff=fabs(fracSM-fracPS);
	coverage=fabs(fracSM+fracPS)/2.0;
	crosspoint=cut;
	//cout<<"New coverage="<<coverage<<" at xpoint="<<crosspoint<<" with diff "<<diff<<"  FracSM="<<fracSM<<"  FracPS="<<fracPS<<endl;
      }
      ////      else cout<<"Diff is too high: "<<fabs(fracSM-fracPS)<<"   fracSM="<<fracSM<<"  fracPS="<<fracPS<<endl;
    }//end if both cuutSM and cutPS are >=0
    //// else cout<<"For cut="<<cut <<" Negative cutSM/cutPS: "<<cutSM<<"  "<<cutPS<<endl;

    cut+=step;
  }//end while loop
 
  cout<<"Finished loop on vector elements, min diff is "<<diff<<", looped until cut_fin="<<cut<<endl;
  cout<<"q value where SM and ALT distributions have same area on opposite sides: "<<crosspoint<<"  Coverage="<<coverage<<endl;
  float separation=2*ROOT::Math::normal_quantile_c(1.0 - coverage, 1.0);
  cout<<"Separation from tail prob: "<<separation<<endl<<endl<<endl;
  

  float integralSM=hSM->Integral();
  float integralPS=hPS->Integral();
 
  float tailSM=hSM->Integral(1,hSM->FindBin(medianPS))/integralSM;
  float tailPS=hPS->Integral(hPS->FindBin(medianSM),hPS->GetNbinsX())/integralPS;
  //cout<<"Median point prob SM: "<<tailSM<<"  ("<<ROOT::Math::normal_quantile_c(tailSM,1.0) <<" sigma)"<<endl;
  //cout<<"Median point prob PS: "<<tailPS<<"  ("<<ROOT::Math::normal_quantile_c(tailPS,1.0) <<" sigma)"<<endl;

  diff=10.0;
  coverage=0.0;
  for(int i=1;i<hSM->GetNbinsX();i++){
    
    float fracSM=hSM->Integral(1,i) / integralSM;
    float fracPS=hPS->Integral(i,hPS->GetNbinsX()) / integralPS;
    if(fabs(fracSM-fracPS)<diff){
      diff=fabs(fracSM-fracPS);
      coverage=(fracSM+fracPS)/2.0;
    }

  }

  float sepH= 2*ROOT::Math::normal_quantile_c(1.0 - coverage, 1.0);
  cout<<"Separation from histograms = "<<sepH<<" with coverage "<<coverage<<endl;

  if(unblind){
    if(v_Obs.size()<1){
      cout<<"Ooops ! The size of the vector with the observed separation is not 1 but "<<v_Obs.size()<<" ! I am not going to plot the observed results."<<endl;
      unblind=false;
    }
    else{
      float obsTailSM=hSM->Integral(1,hSM->FindBin(v_Obs.at(0)))/integralSM;
      float obsTailPS=hPS->Integral(hPS->FindBin(v_Obs.at(0)),hPS->GetNbinsX())/integralPS;
      cout<<"P(SM < Obs): "<<obsTailSM<<"  ("<<ROOT::Math::normal_quantile_c(obsTailSM,1.0) <<" sigma)"<<endl;
      cout<<"P(PS > Obs): "<<obsTailPS<<"  ("<<ROOT::Math::normal_quantile_c(obsTailPS,1.0) <<" sigma)"<<endl;

      float obsCLsRatio = obsTailPS / (1.0 - obsTailSM);
      cout<<"CLs criterion P(PS > Obs) / P(SM > Obs) : "<<obsCLsRatio<<endl;//"  ("<<ROOT::Math::normal_quantile_c(obsCLsRatio,1.0) <<" sigma)"<<endl;
    }
  

    /*
    cout << "\n\nOBSERVED SIGNIFICANCE" << endl;
    cout << "observation: " << v_Obs[0] << endl;
    cout << "bin: " << hObs->GetMaximumBin() << endl;
    cout << " --------------- " << endl;
    double obsPval_SM = 1-hSM->Integral(0,hObs->GetMaximumBin())/integralSM;
    cout << "pvalue SM: " << obsPval_SM << endl;
    cout << "signif SM: " << ROOT::Math::normal_quantile_c(obsPval_SM,1.0) << endl;
    double obsPval_PS =  hPS->Integral(0,hObs->GetMaximumBin())/integralPS;
    cout << "pvalue PS: " << obsPval_PS << endl;
    cout << "signif PS: " << ROOT::Math::normal_quantile_c(obsPval_PS,1.0) << endl<<endl<<endl;
    */

  }//end if unblinding

  //Plotting
  gStyle->SetOptStat(0);
  TCanvas *c1=new TCanvas("c","c",500,500);
  c1->cd();
  hSM->Rebin(50);
  hPS->Rebin(50);
  hSM->SetXTitle(" -2 #times ln(L_{ALT}/L_{0+})");
  hSM->SetYTitle("Pseudoexperiments");
  hPS->SetXTitle(" -2 #times ln(L_{ALT}/L_{0+})");
  hPS->SetYTitle("Pseudoexperiments");
  hSM->SetLineColor(kRed+2);
  hSM->SetLineStyle(2);
  hSM->SetFillColor(798);
  hSM->SetLineWidth(2);
  //hSM->SetFillStyle(3605);
  //TColor *col = gROOT->GetColor(927);
  //col->SetAlpha(0.01);
  //TColor *col=new TColor(1000,72,164,255,"bluish",1);
  //hPS->SetLineColor(927);
  //hPS->SetFillColor(1000);
  hPS->SetFillColor(kAzure+7);
  hPS->SetLineColor(kBlue);
  hPS->SetLineWidth(1);
  hPS->SetFillStyle(3001);
  //hPS->SetFillStyle(3695);

  hObs->SetLineColor(kRed);
  hObs->SetLineWidth(2);

  TGraph *grObs=new TGraph();//dummy, just for the legend
  grObs->SetLineColor(kRed);
  grObs->SetLineWidth(1);
  
  hSM->GetXaxis()->SetRangeUser(-30.0,30.0);
  hSM->GetXaxis()->SetLabelFont(42);
  hSM->GetXaxis()->SetLabelOffset(0.007);
  hSM->GetXaxis()->SetLabelSize(0.045);
  hSM->GetXaxis()->SetTitleSize(0.05);
  hSM->GetXaxis()->SetTitleOffset(1.15);
  hSM->GetXaxis()->SetTitleFont(42);
  hSM->GetYaxis()->SetLabelFont(42);
  hSM->GetYaxis()->SetLabelOffset(0.007);
  hSM->GetYaxis()->SetLabelSize(0.045);
  hSM->GetYaxis()->SetTitleSize(0.05);
  hSM->GetYaxis()->SetTitleOffset(1.8);
  hSM->GetYaxis()->SetTitleFont(42); 
  //TGaxis::SetMaxDigits(2); 
  hSM->Scale(1./hSM->Integral("width"));
  hPS->Scale(1./hPS->Integral("width"));
  float maxhSM=hSM->GetBinContent(hSM->GetMaximumBin());
  float maxhPS=hPS->GetBinContent(hPS->GetMaximumBin());
  if(maxhPS>maxhSM){
    hSM->SetMaximum(maxhPS*1.3);
    hPS->SetMaximum(maxhPS*1.3);
  }
  else{
    hSM->SetMaximum(maxhSM*1.3);
    hPS->SetMaximum(maxhSM*1.3);
  }


  hSM->Draw();
  hPS->Draw("sames");

  TArrow *obsArrow=0;
  if(unblind)obsArrow=new TArrow(v_Obs.at(0),hSM->GetMaximum()/5.0,v_Obs.at(0),0.0,0.05,"|->");
  else obsArrow=new TArrow(0.0,hSM->GetMaximum()/5.0,0.0,0.0,0.05,"|->");
  obsArrow->SetLineColor(kRed);
  obsArrow->SetLineWidth(2.0);
  if(unblind)  obsArrow->Draw();

  //TLegend *leg = new TLegend(0.63,0.73,0.92,0.93);
  TLegend *leg = new TLegend(0.63,0.73,0.88,0.93);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);

  leg->AddEntry(hSM, "0^{+} CJLST","f");
  leg->AddEntry(hPS, legALT+" CJLST","f");
  if(unblind) leg->AddEntry(hObs,"Observed CJLST","L");
  //  if(unblind) leg->AddEntry(hObs,"  Simulated data","L");

  c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.18);
   c1->SetRightMargin(0.05);
   c1->SetTopMargin(0.05);
   c1->SetBottomMargin(0.15);
   c1->SetFrameFillStyle(0);
   c1->SetFrameBorderMode(0);
   c1->SetFrameFillStyle(0);
   c1->SetFrameBorderMode(0);

 TPaveText *pt = new TPaveText(0.1577181,0.9562937,0.9580537,0.9947552,"brNDC");
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetFillStyle(0);
  pt->SetTextFont(42);
  pt->SetTextSize(0.03);
  TText *text = pt->AddText(0.01,0.5,"CMS preliminary");
  text = pt->AddText(0.3,0.6,"#sqrt{s} = 7 TeV, L = 5.1 fb^{-1}  #sqrt{s} = 8 TeV, L = 19.6 fb^{-1}");
  //text = pt->AddText(0.3,0.6,"#sqrt{s} = 8 TeV, L = 19.6 fb^{-1}");
  pt->Draw();   


  TFile *UFLresults= new TFile("/afs/cern.ch/user/s/sbologne/workspace/superME_alessioCode/UFLresults/mu1/sigsep_sepExample_ggSpin0M_plot.root","READ");
  TCanvas *UFLcanv= (TCanvas*) UFLresults->Get("c1");
  TH1F* hUFLPS = (TH1F*)UFLcanv->GetPrimitive("hPS");
  TH1F* hUFLSM = (TH1F*)UFLcanv->GetPrimitive("hSM");
  //TArrow* UFLArrow=new TArrow(10.974,hSM->GetMaximum()/5.0,10.974,0.0,0.05,"|->");
  TArrow* UFLArrow=new TArrow(10.955,hSM->GetMaximum()/5.0,10.955,0.0,0.05,"|->");

  hUFLPS->SetLineColor(kBlack);
  hUFLPS->SetFillStyle(0);
  hUFLSM->SetLineColor(kBlack);
  hUFLSM->SetLineStyle(2);
  hUFLSM->SetFillStyle(0);
  UFLArrow->SetLineColor(kGreen);
  UFLArrow->SetLineWidth(3);
  c1->cd();
  hUFLPS->Scale(1./hUFLPS->Integral("width"));
  hUFLSM->Scale(1./hUFLSM->Integral("width"));
  hUFLPS->Draw("same");
  hUFLSM->Draw("same");
  UFLArrow->Draw("same");

  TH1F* hObs2=new TH1F();
  hObs2->SetLineColor(kGreen);
  hObs2->SetLineWidth(3);
  leg->AddEntry(hUFLSM, "0^{+} UFL","l");
  leg->AddEntry(hUFLPS, legALT+" UFL","l");
  if(unblind) leg->AddEntry(hObs2,"Observed UFL","L");
  leg->Draw();

  c1->SaveAs("sigsep_"+nameALT+".eps");
  c1->SaveAs("sigsep_"+nameALT+".png");
  c1->SaveAs("sigsep_"+nameALT+".root");
  c1->SaveAs("sigsep_"+nameALT+".C");

  return 0;
}//end main

void setTDRStyle() {
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(750); //Height of canvas
  tdrStyle->SetCanvasDefW(1050); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  //  tdrStyle->SetErrorMarker(20);
  tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);

//For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.010);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.25);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.14);
  tdrStyle->SetPadRightMargin(0.04);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.005);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.25);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  //tdrStyle->SetNdivisions(505, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

// Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  Int_t n=35;
  Int_t *colors = new Int_t[n];
  for (Int_t i =0; i < n; i++) {
    //colors[i] = i+61;
    //colors[i] = 70-i;
    colors[i] = i+63;
    }
  /*colors[0]=61
  colors[1]=62;
  colors[2]=63;
  colors[3]=;
  colors[4]=;
  colors[5]=;
  colors[6]=;
  colors[7]=;
  colors[8]=;
  colors[9]=;
  colors[10]=;*/
  gStyle->SetPalette(n, colors);
  
  //gStyle->SetPalette(1, 0);
  tdrStyle->cd();
}
/*
void damnedPlot(){
TCanvas *c1 =  new TCanvas("paddy","paddy",700,700)
  c1->Divide(1,2,0,0);
 c1->cd(1).SetBottomMargin(0.001);
 c1->cd(1).SetTopMargin(0.01);
 c1->cd(1).SetRightMargin(0.01);
 c1.cd(2);
 c1.cd(2).SetTopMargin(0.001);
 c1.cd(2).SetRightMargin(0.01);

}
*/

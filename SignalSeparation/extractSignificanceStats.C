#include <Riostream.h>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#include "TArrow.h"
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
#include "Math/DistFunc.h"

int extractSignificanceStats(){

  const float lumi7TeV=5.051;
  const float lumi8TeV=5.261;

  char fileName[128];
  sprintf(fileName,"qmu.root");
  TFile *fq=new TFile(fileName,"READ");
  TTree *t=(TTree*)fq->Get("q");

  float q,m,w;
  int type;
  t->SetBranchAddress("q",&q);
  t->SetBranchAddress("mh",&m);
  t->SetBranchAddress("weight",&w);
  t->SetBranchAddress("type",&type);

  TH1F *hSM=new TH1F("hSM;S = -2 #times ln(L_{1}/L_{2});Number of Toys","",8000,-15,15);
  TH1F *hPS=new TH1F("hPS;S = -2 #times ln(L_{1}/L_{2});Number of Toys","",8000,-15,15);
  TH1F *hObs=new TH1F("hObserved","",8000,-15,15);
  cout<<"Start to lopp on tree in file "<<fileName<<endl;

  std::vector<float> v_SM, v_PS,v_Obs;

  for(int i=0;i<t->GetEntries();i++){
    t->GetEntry(i);
    if(i==0)cout<<"MASS in the TREE = "<<m<<endl<<endl;

    if(type<0){ //SM hypothesis 
      hSM->Fill(-q);
      v_SM.push_back(-q);
    }
    else if(type>0){//ALT hypothesis
      hPS->Fill(-q);
      v_PS.push_back(-q);
    }
    else{
      hObs->Fill(q);
      v_Obs.push_back(q);
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

  /*
  while(cut<=v_SM.at(ntoysSM-1)+step){
    //    if(int(cut*100)%100==0)
cout<<"Cutting at "<<cut<<endl;
    for(int iSM=startSM;iSM>=0;iSM--){
      
      if(v_SM.at(iSM)<cut){
	startSM=ntoysSM-iSM;
	//break;
      }
      else cout<<"SM "<<v_SM.at(iSM)<<" > "<<cut<<endl;
    }

    for(int iPS=startPS;iPS<ntoysPS;iPS++){
      if(v_PS.at(iPS)>cut){
	startPS=iPS;
	//break;
      }
      else cout<<v_PS.at(iPS)<<" < "<<cut<<endl;
  
    }
    float fracSM=(ntoysSM-startSM)/ntoysSM;
    float fracPS=startPS/ntoysPS;
    cout<<"Frac "<<fracSM<<" "<<fracPS<<endl;
    if(fabs(fracSM-fracPS)<diff){
      diff=fabs(fracSM-fracPS);
      coverage=fabs(fracSM-fracPS)/2.0;
      crosspoint=cut;
      cout<<"New coverage="<<coverage<<" at xpoint="<<crosspoint<<"  "<<startSM<<endl;
    }
    cut+=step;
  }//end while loop
  */
  cout<<"Finished loop on vector elements, min is "<<diff<<" cut is at "<<cut<<endl;
  cout<<"q value where SM and ALT distributions have same area on opposite sides: "<<crosspoint<<endl;
  cout<<"Coverage "<<coverage<<endl;
  float separation=2*ROOT::Math::normal_quantile_c(1.0 - coverage, 1.0);
  cout<<"Separation: "<<separation<<endl<<endl<<endl;

  float integralSM=hSM->Integral();
  float integralPS=hPS->Integral();
 
  float tailSM=hSM->Integral(1,hSM->FindBin(medianPS))/integralSM;
  float tailPS=hPS->Integral(hPS->FindBin(medianSM),hPS->GetNbinsX())/integralPS;
  cout<<"Tail prob SM: "<<tailSM<<"  ("<<ROOT::Math::normal_quantile_c(tailSM,1.0) <<" sigma)"<<endl;
  cout<<"Tail prob PS: "<<tailPS<<"  ("<<ROOT::Math::normal_quantile_c(tailPS,1.0) <<" sigma)"<<endl;

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


  cout << "OBSERVED SIGNIFICANCE" << endl;

  cout << "observation: " << v_Obs[0] << endl;
  cout << "bin: " << hObs->GetMaximumBin() << endl;
  
  cout << " --------------- " << endl;
  double obsPval_SM = 1-hSM->Integral(0,hObs->GetMaximumBin())/hSM->Integral();
  cout << "pvalue SM: " << obsPval_SM << endl;
  cout << "signif SM: " << ROOT::Math::normal_quantile_c(obsPval_SM,1.0) << endl;;
  double obsPval_PS =  hPS->Integral(0,hObs->GetMaximumBin())/hPS->Integral();
  cout << "pvalue PS: " << obsPval_PS << endl;
  cout << "signif PS: " << ROOT::Math::normal_quantile_c(obsPval_PS,1.0) << endl;;


  gStyle->SetOptStat(0);
  TCanvas *c1=new TCanvas("c1","c1",500,500);
  c1->cd();
  hSM->Rebin(50);
  hPS->Rebin(50);
  float maxhSM=hSM->GetBinContent(hSM->GetMaximumBin());
  float maxhPS=hPS->GetBinContent(hPS->GetMaximumBin());
  if(maxhPS>maxhSM){
    hSM->SetMaximum(maxhPS*1.15);
    hPS->SetMaximum(maxhPS*1.15);
  }
  else{
    hSM->SetMaximum(maxhSM*1.15);
    hPS->SetMaximum(maxhSM*1.15);
  }
  hSM->SetXTitle("S = -2 #times ln(L_{1}/L_{2})");
  hSM->SetYTitle("Generated experiments");
  hPS->SetXTitle("S = -2 #times ln(L_{1}/L_{2})");
  hPS->SetYTitle("Generated experiments");
  hSM->SetLineColor(kMagenta-3);
  hSM->SetFillColor(kMagenta-3);
  hSM->SetLineWidth(2);
  hSM->SetFillStyle(3605);
  hPS->SetLineColor(kBlue+1);
  hPS->SetFillColor(kBlue+1);
  hPS->SetLineWidth(2);
  hPS->SetFillStyle(3695);
  hObs->SetLineColor(kGreen+3);
  hObs->SetLineWidth(2);
  hSM->Draw();
  hPS->Draw("sames");
  
  TArrow *obsArrow  = new TArrow(v_Obs[0],hSM->GetMaximum()/2.,v_Obs[0],0.0,.05,"|>");
  obsArrow->SetLineWidth(3);
  obsArrow->Draw("same");
  //hObs->Draw("sames");

  TLegend *leg = new TLegend(0.7,0.6,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hSM,"  PS, 0-","f");
  leg->AddEntry(hPS,"  SM, 0+","f");
  leg->Draw();


  TPaveText pt(0.16,0.95,0.45,0.99,"NDC");
  pt.SetFillColor(0);
  pt.AddText("CMS Expected");
  pt.SetBorderSize(0);
  TPaveText pt2(0.55,0.95,0.99,0.99,"NDC");
  pt2.SetFillColor(0);
  pt2.AddText(Form(" #sqrt{s} = 7 TeV, L = %.3f fb^{-1}; #sqrt{s} = 8 TeV, L = %.3f fb^{-1}",lumi7TeV,lumi8TeV));
  pt2.SetBorderSize(0);
  pt.Draw();
  pt2.Draw();
  c1->SaveAs("sigsep_combine.eps");
  c1->SaveAs("sigsep_combine.root");

  return 0;
}//end main

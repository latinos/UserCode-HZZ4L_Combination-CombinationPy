#include <algorithm>

void computeObservation(){

  //l1 = [Y( Alt hyp. )*pdf( x=event1 | x~Alt hyp. ) + Y( ZZ )*pdf( x=event1 | x~ZZ )] / [Y( H )*pdf( x=event1 | x~H ) + Y( ZZ )*pdf( x=event1 | x~ZZ )] 

  //0m
  double corr2e2mu = 1.12188;
  double corr4mu = 0.912885;
  double corr4e = 0.858691;
  /*
  //gg->2mp
  double corr2e2mu = 1.12826;
  double corr4mu = 0.900717;
  double corr4e = 0.865112;
  //0h+
  double corr2e2mu = 1.06477;
  double corr4mu = .941477;
  double corr4e = 0.947087;
  //1+
  double corr2e2mu = 1.07114;
  double corr4mu = 0.968514;
  double corr4e = 0.882397;
  //1-
  //double corr2e2mu = 1.1036;
  //double corr4mu = 0.939779;
  //double corr4e = 0.854791;
  //qq->2mp
  double corr2e2mu = 1.1417;
  double corr4mu = 0.884544;
  double corr4e = 0.861437;
  */
  

  TString dataLoc="./dataTemplates";

  TString templateLoc_8TeV = "templates2D/templates_8TeV";
  TString templateLoc_7TeV = "templates2D";
	
  TFile *f_data_2e2mu_7TeV= new TFile( dataLoc+"/Hist_Scatter_Data_2e2m_7TeV.dat.root","READ");
  TH2F* dataHist_2e2mu_7TeV= f_data_2e2mu_7TeV->Get("hist");
  dataHist_2e2mu_7TeV->SetName("data_1m_2e2mu_7TeV");

  TFile *f_templ_SM_2e2mu_7TeV= new TFile(templateLoc_7TeV+"/Dsignal_2e2mu.root","READ");
  TH2F* templ_SM_2e2mu_7TeV = f_templ_SM_2e2mu_7TeV->Get("h_superDpsD");
  templ_SM_2e2mu_7TeV->SetName("h_superDpsD_1m_SM_2e2mu_7TeV");

  TFile *f_templ_ALT_2e2mu_7TeV= new TFile(templateLoc_7TeV+"/Dsignal_ALT_2e2mu.root","READ");
  TH2F* templ_ALT_2e2mu_7TeV = f_templ_ALT_2e2mu_7TeV->Get("h_superDpsD");
  templ_ALT_2e2mu_7TeV->SetName("h_superDpsD_1m_ALT_2e2mu_7TeV");

  TFile *f_templ_bkg_2e2mu_7TeV= new TFile(templateLoc_7TeV+"/Dbackground_qqZZ_2e2mu.root","READ");
  TH2F* templ_bkg_2e2mu_7TeV = f_templ_bkg_2e2mu_7TeV->Get("h_superDpsD");
  templ_bkg_2e2mu_7TeV->SetName("h_superDpsD_1m_bkg_2e2mu_7TeV");
 
  TFile *f_data_4e_7TeV= new TFile( dataLoc+"/Hist_Scatter_Data_4e_7TeV.dat.root","READ");
  TH2F* dataHist_4e_7TeV= f_data_4e_7TeV->Get("hist");
  dataHist_4e_7TeV->SetName("data_1m_4e_7TeV");

  TFile *f_templ_SM_4e_7TeV= new TFile(templateLoc_7TeV+"/Dsignal_4e.root","READ");
  TH2F* templ_SM_4e_7TeV = f_templ_SM_4e_7TeV->Get("h_superDpsD");
  templ_SM_4e_7TeV->SetName("h_superDpsD_1m_SM_4e_7TeV");

  TFile *f_templ_ALT_4e_7TeV= new TFile(templateLoc_7TeV+"/Dsignal_ALT_4e.root","READ");
  TH2F* templ_ALT_4e_7TeV = f_templ_ALT_4e_7TeV->Get("h_superDpsD");
  templ_ALT_4e_7TeV->SetName("h_superDpsD_1m_ALT_4e_7TeV");

  TFile *f_templ_bkg_4e_7TeV= new TFile(templateLoc_7TeV+"/Dbackground_qqZZ_4e.root","READ");
  TH2F* templ_bkg_4e_7TeV = f_templ_bkg_4e_7TeV->Get("h_superDpsD");
  templ_bkg_4e_7TeV->SetName("h_superDpsD_1m_bkg_4e_7TeV");
 
  TFile *f_data_4mu_7TeV= new TFile( dataLoc+"/Hist_Scatter_Data_4m_7TeV.dat.root","READ");
  TH2F* dataHist_4mu_7TeV= f_data_4mu_7TeV->Get("hist");
  dataHist_4mu_7TeV->SetName("data_1m_4mu_7TeV");

  TFile *f_templ_SM_4mu_7TeV= new TFile(templateLoc_7TeV+"/Dsignal_4mu.root","READ");
  TH2F* templ_SM_4mu_7TeV = f_templ_SM_4mu_7TeV->Get("h_superDpsD");
  templ_SM_4mu_7TeV->SetName("h_superDpsD_1m_SM_4mu_7TeV");

  TFile *f_templ_ALT_4mu_7TeV= new TFile(templateLoc_7TeV+"/Dsignal_ALT_4mu.root","READ");
  TH2F* templ_ALT_4mu_7TeV = f_templ_ALT_4mu_7TeV->Get("h_superDpsD");
  templ_ALT_4mu_7TeV->SetName("h_superDpsD_1m_ALT_4mu_7TeV");

  TFile *f_templ_bkg_4mu_7TeV= new TFile(templateLoc_7TeV+"/Dbackground_qqZZ_4mu.root","READ");
  TH2F* templ_bkg_4mu_7TeV = f_templ_bkg_4mu_7TeV->Get("h_superDpsD");
  templ_bkg_4mu_7TeV->SetName("h_superDpsD_1m_bkg_4mu_7TeV");
 

  cout<<"2e2mu 7TeV -----------"<<endl;
  double l_2e2mu_7TeV=0;
  for (int i=1; i<=dataHist_2e2mu_7TeV->GetXaxis()->GetNbins(); i++){
    for (int j=1; j<=dataHist_2e2mu_7TeV->GetYaxis()->GetNbins(); j++){
      int evNumber=dataHist_2e2mu_7TeV->GetBinContent(i,j);
      if (evNumber<1)
	continue;

      double Y_SM= 1.78015;
      double Y_ALT=1.78015*corr2e2mu;
      double Y_bkg=2.3083+0.7711;
    
      double pdf_SM = templ_SM_2e2mu_7TeV->GetBinContent(i,j)/templ_SM_2e2mu_7TeV->Integral();
      double pdf_bkg = templ_bkg_2e2mu_7TeV->GetBinContent(i,j)/templ_bkg_2e2mu_7TeV->Integral();
      double pdf_ALT = templ_ALT_2e2mu_7TeV->GetBinContent(i,j)/templ_ALT_2e2mu_7TeV->Integral();
    
      double l = (Y_ALT*pdf_ALT + Y_bkg*pdf_bkg)/(Y_SM*pdf_SM + Y_bkg*pdf_bkg);
      cout<<"l "<<l<<endl;
      double l_bin = pow(l,evNumber);
      cout<<"l_bin "<<l_bin<<endl;

      //l_2e2mu_7TeV= l_2e2mu_7TeV*l;
      l_2e2mu_7TeV= l_2e2mu_7TeV + log(l);
      cout<<"l_2e2mu_7TeV "<<l_2e2mu_7TeV<<endl;
    }
  }
  cout<<"---------------"<<endl<<endl;

  cout<<"4e 7TeV -----------"<<endl;
  double l_4e_7TeV=0;
  for (int i=1; i<=dataHist_4e_7TeV->GetXaxis()->GetNbins(); i++){
    for (int j=1; j<=dataHist_4e_7TeV->GetYaxis()->GetNbins(); j++){
      int evNumber=dataHist_4e_7TeV->GetBinContent(i,j);
      if (evNumber<1)
	continue;

      double Y_SM= 0.723492;
      double Y_ALT=0.723492*corr4e;
      double Y_bkg=0.8552+0.3653;
    
      double pdf_SM = templ_SM_4e_7TeV->GetBinContent(i,j)/templ_SM_4e_7TeV->Integral();
      double pdf_bkg = templ_bkg_4e_7TeV->GetBinContent(i,j)/templ_bkg_4e_7TeV->Integral();
      double pdf_ALT = templ_ALT_4e_7TeV->GetBinContent(i,j)/templ_ALT_4e_7TeV->Integral();
    
      double l = (Y_ALT*pdf_ALT + Y_bkg*pdf_bkg)/(Y_SM*pdf_SM + Y_bkg*pdf_bkg);
      cout<<"l "<<l<<endl;
      double l_bin = pow(l,evNumber);
      cout<<"l_bin "<<l_bin<<endl;

      //l_4e_7TeV= l_4e_7TeV*l;
      l_4e_7TeV= l_4e_7TeV+log(l);
      cout<<"l_4e_7TeV "<<l_4e_7TeV<<endl;
    }
  }
  cout<<"---------------"<<endl<<endl;

  cout<<"4mu 7TeV -----------"<<endl;
  double l_4mu_7TeV=0;
  for (int i=1; i<=dataHist_4mu_7TeV->GetXaxis()->GetNbins(); i++){
    for (int j=1; j<=dataHist_4mu_7TeV->GetYaxis()->GetNbins(); j++){
      int evNumber=dataHist_4mu_7TeV->GetBinContent(i,j);
      if (evNumber<1)
	continue;

      double Y_SM= 1.30233;
      double Y_ALT=1.30233*corr4mu;
      double Y_bkg=1.8182+0.5872;
    
      double pdf_SM = templ_SM_4mu_7TeV->GetBinContent(i,j)/templ_SM_4mu_7TeV->Integral();
      double pdf_bkg = templ_bkg_4mu_7TeV->GetBinContent(i,j)/templ_bkg_4mu_7TeV->Integral();
      double pdf_ALT = templ_ALT_4mu_7TeV->GetBinContent(i,j)/templ_ALT_4mu_7TeV->Integral();
    
      double l = (Y_ALT*pdf_ALT + Y_bkg*pdf_bkg)/(Y_SM*pdf_SM + Y_bkg*pdf_bkg);
      cout<<"l "<<l<<endl;
      double l_bin = pow(l,evNumber);
      cout<<"l_bin "<<l_bin<<endl;

      //l_4mu_7TeV= l_4mu_7TeV*l;
      l_4mu_7TeV= l_4mu_7TeV+log(l);
      cout<<"l_4mu_7TeV "<<l_4mu_7TeV<<endl;
    }
  }
  cout<<"---------------"<<endl<<endl;


  TFile *f_data_2e2mu_8TeV= new TFile( dataLoc+"/Hist_Scatter_Data_2e2m_8TeV.dat.root","READ");
  TH2F* dataHist_2e2mu_8TeV= f_data_2e2mu_8TeV->Get("hist");
  dataHist_2e2mu_8TeV->SetName("data_1m_2e2mu_8TeV");

  TFile *f_templ_SM_2e2mu_8TeV= new TFile(templateLoc_8TeV+"/Dsignal_2e2mu.root","READ");
  TH2F* templ_SM_2e2mu_8TeV = f_templ_SM_2e2mu_8TeV->Get("h_superDpsD");
  templ_SM_2e2mu_8TeV->SetName("h_superDpsD_1m_SM_2e2mu_8TeV");

  TFile *f_templ_ALT_2e2mu_8TeV= new TFile(templateLoc_8TeV+"/Dsignal_ALT_2e2mu.root","READ");
  TH2F* templ_ALT_2e2mu_8TeV = f_templ_ALT_2e2mu_8TeV->Get("h_superDpsD");
  templ_ALT_2e2mu_8TeV->SetName("h_superDpsD_1m_ALT_2e2mu_8TeV");

  TFile *f_templ_bkg_2e2mu_8TeV= new TFile(templateLoc_8TeV+"/Dbackground_qqZZ_2e2mu.root","READ");
  TH2F* templ_bkg_2e2mu_8TeV = f_templ_bkg_2e2mu_8TeV->Get("h_superDpsD");
  templ_bkg_2e2mu_8TeV->SetName("h_superDpsD_1m_bkg_2e2mu_8TeV");
 
  TFile *f_data_4e_8TeV= new TFile( dataLoc+"/Hist_Scatter_Data_4e_8TeV.dat.root","READ");
  TH2F* dataHist_4e_8TeV= f_data_4e_8TeV->Get("hist");
  dataHist_4e_8TeV->SetName("data_1m_4e_8TeV");

  TFile *f_templ_SM_4e_8TeV= new TFile(templateLoc_8TeV+"/Dsignal_4e.root","READ");
  TH2F* templ_SM_4e_8TeV = f_templ_SM_4e_8TeV->Get("h_superDpsD");
  templ_SM_4e_8TeV->SetName("h_superDpsD_1m_SM_4e_8TeV");

  TFile *f_templ_ALT_4e_8TeV= new TFile(templateLoc_8TeV+"/Dsignal_ALT_4e.root","READ");
  TH2F* templ_ALT_4e_8TeV = f_templ_ALT_4e_8TeV->Get("h_superDpsD");
  templ_ALT_4e_8TeV->SetName("h_superDpsD_1m_ALT_4e_8TeV");

  TFile *f_templ_bkg_4e_8TeV= new TFile(templateLoc_8TeV+"/Dbackground_qqZZ_4e.root","READ");
  TH2F* templ_bkg_4e_8TeV = f_templ_bkg_4e_8TeV->Get("h_superDpsD");
  templ_bkg_4e_8TeV->SetName("h_superDpsD_1m_bkg_4e_8TeV");
 
  TFile *f_data_4mu_8TeV= new TFile( dataLoc+"/Hist_Scatter_Data_4m_8TeV.dat.root","READ");
  TH2F* dataHist_4mu_8TeV= f_data_4mu_8TeV->Get("hist");
  dataHist_4mu_8TeV->SetName("data_1m_4mu_8TeV");

  TFile *f_templ_SM_4mu_8TeV= new TFile(templateLoc_8TeV+"/Dsignal_4mu.root","READ");
  TH2F* templ_SM_4mu_8TeV = f_templ_SM_4mu_8TeV->Get("h_superDpsD");
  templ_SM_4mu_8TeV->SetName("h_superDpsD_1m_SM_4mu_8TeV");

  TFile *f_templ_ALT_4mu_8TeV= new TFile(templateLoc_8TeV+"/Dsignal_ALT_4mu.root","READ");
  TH2F* templ_ALT_4mu_8TeV = f_templ_ALT_4mu_8TeV->Get("h_superDpsD");
  templ_ALT_4mu_8TeV->SetName("h_superDpsD_1m_ALT_4mu_8TeV");

  TFile *f_templ_bkg_4mu_8TeV= new TFile(templateLoc_8TeV+"/Dbackground_qqZZ_4mu.root","READ");
  TH2F* templ_bkg_4mu_8TeV = f_templ_bkg_4mu_8TeV->Get("h_superDpsD");
  templ_bkg_4mu_8TeV->SetName("h_superDpsD_1m_bkg_4mu_8TeV");
 

  cout<<"2e2mu 8TeV -----------"<<endl;
  double l_2e2mu_8TeV=0;
  for (int i=1; i<=dataHist_2e2mu_8TeV->GetXaxis()->GetNbins(); i++){
    for (int j=1; j<=dataHist_2e2mu_8TeV->GetYaxis()->GetNbins(); j++){
      int evNumber=dataHist_2e2mu_8TeV->GetBinContent(i,j);
      if (evNumber<1)
	continue;

      double Y_SM= 7.99256;
      double Y_ALT=7.99256*corr2e2mu;
      double Y_bkg=9.4170+1.7876;
    
      double pdf_SM = templ_SM_2e2mu_8TeV->GetBinContent(i,j)/templ_SM_2e2mu_8TeV->Integral();
      double pdf_bkg = templ_bkg_2e2mu_8TeV->GetBinContent(i,j)/templ_bkg_2e2mu_8TeV->Integral();
      double pdf_ALT = templ_ALT_2e2mu_8TeV->GetBinContent(i,j)/templ_ALT_2e2mu_8TeV->Integral();
    
      double l = (Y_ALT*pdf_ALT + Y_bkg*pdf_bkg)/(Y_SM*pdf_SM + Y_bkg*pdf_bkg);
      cout<<"l "<<l<<endl;
      double l_bin = pow(l,evNumber);
      cout<<"l_bin "<<l_bin<<endl;

      //l_2e2mu_8TeV= l_2e2mu_8TeV*l;
      l_2e2mu_8TeV= l_2e2mu_8TeV+log(l);
      cout<<"l_2e2mu_8TeV "<<l_2e2mu_8TeV<<endl;
    }
  }
  cout<<"---------------"<<endl<<endl;



  cout<<"4e 8TeV -----------"<<endl;
  double l_4e_8TeV=0;
  for (int i=1; i<=dataHist_4e_8TeV->GetXaxis()->GetNbins(); i++){
    for (int j=1; j<=dataHist_4e_8TeV->GetYaxis()->GetNbins(); j++){
      int evNumber=dataHist_4e_8TeV->GetBinContent(i,j);
      if (evNumber<1)
	continue;

      double Y_SM= 3.16754;
      double Y_ALT= 3.16754*corr4e;
      double Y_bkg= 3.2180+1.1435;
    
      double pdf_SM = templ_SM_4e_8TeV->GetBinContent(i,j)/templ_SM_4e_8TeV->Integral();
      double pdf_bkg = templ_bkg_4e_8TeV->GetBinContent(i,j)/templ_bkg_4e_8TeV->Integral();
      double pdf_ALT = templ_ALT_4e_8TeV->GetBinContent(i,j)/templ_ALT_4e_8TeV->Integral();
    
      double l = (Y_ALT*pdf_ALT + Y_bkg*pdf_bkg)/(Y_SM*pdf_SM + Y_bkg*pdf_bkg);
      cout<<"l "<<l<<endl;
      double l_bin = pow(l,evNumber);
      cout<<"l_bin "<<l_bin<<endl;

      //l_4e_8TeV= l_4e_8TeV*l;
      l_4e_8TeV= l_4e_8TeV+log(l);
      cout<<"l_4e_8TeV "<<l_4e_8TeV<<endl;
    }
  }
  cout<<"---------------"<<endl<<endl;


  cout<<"4mu 8TeV -----------"<<endl;
  double l_4mu_8TeV=0;
  for (int i=1; i<=dataHist_4mu_8TeV->GetXaxis()->GetNbins(); i++){
    for (int j=1; j<=dataHist_4mu_8TeV->GetYaxis()->GetNbins(); j++){
      int evNumber=dataHist_4mu_8TeV->GetBinContent(i,j);
      if (evNumber<1)
	continue;

      double Y_SM= 6.07817;
      double Y_ALT= 6.07817*corr4mu;
      double Y_bkg= 7.6762 + 0.9318;
    
      double pdf_SM = templ_SM_4mu_8TeV->GetBinContent(i,j)/templ_SM_4mu_8TeV->Integral();
      double pdf_bkg = templ_bkg_4mu_8TeV->GetBinContent(i,j)/templ_bkg_4mu_8TeV->Integral();
      double pdf_ALT = templ_ALT_4mu_8TeV->GetBinContent(i,j)/templ_ALT_4mu_8TeV->Integral();
    
      double l = (Y_ALT*pdf_ALT + Y_bkg*pdf_bkg)/(Y_SM*pdf_SM + Y_bkg*pdf_bkg);
      cout<<"l "<<l<<endl;
      double l_bin = pow(l,evNumber);
      cout<<"l_bin "<<l_bin<<endl;

      //l_4mu_8TeV= l_4mu_8TeV*l;
      l_4mu_8TeV= l_4mu_8TeV+log(l);
      cout<<"l_4mu_8TeV "<<l_4mu_8TeV<<endl;
    }
  }
  cout<<"---------------"<<endl<<endl;


  //double l_7TeV = l_2e2mu_7TeV*l_4e_7TeV*l_4mu_7TeV  ;
  //double l_8TeV = l_2e2mu_8TeV*l_4e_8TeV*l_4mu_8TeV  ;
  //double l_tot = l_7TeV*l_8TeV;
  double l_7TeV = l_2e2mu_7TeV+l_4e_7TeV+l_4mu_7TeV  ;
  double l_8TeV = l_2e2mu_8TeV+l_4e_8TeV+l_4mu_8TeV  ;
  double l_tot = l_7TeV+l_8TeV;
  
  cout<<"l_7TeV "<<l_7TeV<<endl;
  cout<<"l_8TeV "<<l_8TeV<<endl;
  cout<<"l_tot "<<l_tot<<endl;

  cout << "Final result -2xln(q) = "<<(-2*(l_tot)) << endl;
}

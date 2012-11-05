void plotScan1D(){


  TFile *f=new TFile("higgsCombine1D.MultiDimFit.mH125.root");
  TTree *t=(TTree*)f->Get("limit");
  t->Draw("2*deltaNLL:x", "deltaNLL > 0","PL");

  TGraph *gr0 = (TGraph*) gROOT->FindObject("Graph")->Clone();
  gr0->SetName("scan1D");

  TCanvas *c1=new TCanvas("can1","CANVAS-SCAN1D",800,800);
  c1->cd();
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  gr0->GetXaxis()->SetTitle("f_{a3}");
  gr0->GetYaxis()->SetTitle("2 #times #Delta(NLL)");
  gr0->GetXaxis()->SetLabelSize(0.04);
  gr0->GetYaxis()->SetLabelSize(0.04);
  gr0->Draw("APL");

 float lumi7TeV=5.051;
  float lumi8TeV=12.21;
  TPaveText pt(0.16,0.95,0.40,0.99,"NDC");
  pt.SetFillColor(0);
  pt.SetTextAlign(12);
  pt.SetTextSize(0.027);
  pt.AddText("CMS Preliminary");
  pt.SetBorderSize(0);
  TPaveText pt2(0.53,0.95,0.98,0.99,"NDC");
  pt2.SetFillColor(0);
  pt2.SetTextAlign(32);
  pt2.SetTextSize(0.027);
  pt2.AddText(Form(" #sqrt{s} = 7 TeV, L = %.3f fb^{-1}; #sqrt{s} = 8 TeV, L = %.2f fb^{-1}",lumi7TeV,lumi8TeV));
  pt2.SetBorderSize(0);
  pt.Draw();
  pt2.Draw();
  c1->SaveAs("can_scan1D_x_HCPlumi.root");
  c1->SaveAs("can_scan1D_x_HCPlumi.eps");
}

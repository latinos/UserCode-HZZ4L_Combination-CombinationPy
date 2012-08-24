bool cms_excluded = false;
bool lep_excluded = false;
bool tev_excluded = false, tev_excluded_alsolow = false;
TH1 *fakeLEP = 0, *fakeTEV = 0, *fakeCMS = 0;
void drawExclusions(double xmin0, double xmax0, double ymin, double ymax) {
    if (cms_excluded) {
        double xmin = 127.5, xmax = TMath::Min(xmax0,600.);
        TBox box(xmin, ymin, xmax, ymax);
        box.SetLineStyle(0);
        box.SetFillStyle(3345);
        box.SetFillColor(95);
        box.DrawClone();
        TLine line(xmin, ymin, xmin, ymax);
        line.SetLineWidth(2);
        line.SetLineColor(95);
        line.DrawClone();
        if (xmax0 > 600) line.DrawLine(xmax, ymin, xmax, ymax);
        if (fakeCMS == 0) {
            fakeCMS = new TH1F("fake_CMS","fake_CMS",1,0.,1.);
            fakeCMS->SetLineStyle(0);
            fakeCMS->SetFillStyle(3345);
            fakeCMS->SetFillColor(95);
            fakeCMS->SetLineWidth(2);
            fakeCMS->SetLineColor(95);
        }
    } 
    if (tev_excluded && xmax0 > 147) {
        double xmin = 147, xmax = TMath::Min(179., xmax0);
        TBox box(xmin, ymin, xmax, ymax);
        box.SetLineStyle(0);
        box.SetFillStyle(3354);
        box.SetFillColor(213);
        box.DrawClone();
        TLine line(xmin, ymin, xmin, ymax);
        line.SetLineWidth(2);
        line.SetLineColor(213);
        line.DrawClone();
        if (xmax0 > 179) line.DrawLine(xmax, ymin, xmax, ymax);
        xmin = 100, xmax = 109;
        TBox box2(xmin, ymin, xmax, ymax);
        box2.SetLineStyle(0);
        box2.SetFillStyle(3354);
        box2.SetFillColor(213);
        if (xmin0 < 109 && tev_excluded_alsolow) box2.DrawClone();
        TLine line2(xmin, ymin, xmin, ymax);
        line2.SetLineWidth(2);
        line2.SetLineColor(213);
        if (xmin0 < 109 && tev_excluded_alsolow) line2.DrawClone();
        if (xmin0 < 109 && tev_excluded_alsolow) line2.DrawLine(xmax, ymin, xmax, ymax);
        if (fakeTEV == 0) {
            fakeTEV = new TH1F("fake_tev","fake_tev",1,0.,1.);
            fakeTEV->SetLineStyle(0);
            fakeTEV->SetFillStyle(3354);
            fakeTEV->SetFillColor(213);
            fakeTEV->SetLineWidth(2);
            fakeTEV->SetLineColor(213);
        }
    } 
    if (lep_excluded) {
        TBox box(xmin0, ymin, 114.4, ymax);
        box.SetLineStyle(0);
        box.SetFillStyle(3345);
        box.SetFillColor(209);
        box.DrawClone();
        TLine line(114.4, ymin, 114.4, ymax);
        line.SetLineWidth(2);
        line.SetLineColor(209);
        line.DrawClone();
        if (fakeLEP == 0) {
            fakeLEP = new TH1F("fake_lep","fake_lep",1,0.,1.);
            fakeLEP->SetLineStyle(0);
            fakeLEP->SetFillStyle(3345);
            fakeLEP->SetFillColor(209);
            fakeLEP->SetLineWidth(2);
            fakeLEP->SetLineColor(209);
        }
    } 

}

void setCanvas(TH1 *first, TString title, double min=0.1, double max=30, TString ytitle="Limit (#sigma_{95%}/#sigma_{SM})") { 
    if (first == 0) std::cerr << "Error for " << title << std::endl;
    first->SetTitle(title);
    first->GetYaxis()->SetTitle(ytitle);
    c1->SetLogy(min > 0 || title.Contains("pval"));    
    first->GetYaxis()->SetRangeUser(min,max);
    if (min <= 0) first->GetYaxis()->SetDecimals(true);
    first->GetYaxis()->SetTitleOffset(1.00+0.2*isSquareCanvas);
    first->GetXaxis()->SetTitle("Higgs boson mass ("+massUnits+")");
    if (SM == "SM4") {
        first->GetXaxis()->SetTitle("SM4 Higgs boson mass ("+massUnits+")");
    } else if (SM == "FP") {
        first->GetXaxis()->SetTitle("FP Higgs boson mass ("+massUnits+")");
    }
    c1->SetTickx(1);    
    c1->SetTicky(1);    
    double xmin = first->GetXaxis()->GetXmin();
    double xmax = first->GetXaxis()->GetXmax(); 
    drawExclusions(xmin, xmax, min, max);
    if (ytitle.Contains("p-value")) {
        //TBox twoSig(xmin, 1, xmax, ROOT::Math::normal_cdf_c(2)); twoSig.SetFillColor(220); twoSig.DrawClone();
        //TBox oneSig(xmin, 1, xmax, ROOT::Math::normal_cdf_c(1)); oneSig.SetFillColor(211); oneSig.DrawClone();
        TLine threeSig(xmin, ROOT::Math::normal_cdf_c(3), xmax, ROOT::Math::normal_cdf_c(3)); 
        threeSig.SetLineColor(lineAt1Color); threeSig.SetLineWidth(2); threeSig.SetLineStyle(c1->GetGridy() ? 7 : lineAt1Style); 
        TLatex latex; latex.SetTextFont(42); latex.SetTextSize(0.045); latex.SetTextColor(2);
        bool protrude = c1->GetGridy(); 
        for (double i = 1; i <= 8; ++i) {
            if (min >  ROOT::Math::normal_cdf_c(i)) break;
            threeSig.DrawLine(xmin, ROOT::Math::normal_cdf_c(i), xmax+(xmax-xmin)*0.015*protrude, ROOT::Math::normal_cdf_c(i));
            TString sspam = TString::Format("%d#sigma", int(i));
            latex.DrawLatex(xmax+(xmax-xmin)*0.01, ROOT::Math::normal_cdf_c(i)*1.1, sspam.Data());
        }
    }
    if (ytitle.Contains("CL_{S} of")) {
        double xmin = first->GetXaxis()->GetXmin();
        double xmax = first->GetXaxis()->GetXmax(); 
        //if (max > 1) {
        //    TLine oneLine(xmin,1,xmax,1); oneLine.SetLineWidth(2);
        //    oneLine.DrawClone();
        //}
        double lines[3] = { 1-0.90, 1-0.95, 1-0.99 };
        TLine hline(xmin, ROOT::Math::normal_cdf_c(3), xmax, ROOT::Math::normal_cdf_c(3)); 
        hline.SetLineColor(lineAt1Color); 
        TLatex latex; latex.SetTextFont(42); latex.SetTextSize(0.045 - 0.008*isSquareCanvas); latex.SetTextColor(2);
        for (int i = 0; i < 3; ++i) {
            hline.SetLineWidth(4); hline.SetLineStyle(lineAt1Style);
            //if (i == 1) { hline.SetLineWidth(4); hline.SetLineStyle(1); }
            //else        { hline.SetLineWidth(2); hline.SetLineStyle(2); }
            bool protrude = c1->GetGridy(); 
            hline.DrawLine(xmin, lines[i], xmax+(xmax-xmin)*0.02*protrude, lines[i]);
            TString sspam = TString::Format("%.0f%%", 100*(1-lines[i]));
            latex.DrawLatex(xmax+(xmax-xmin)*0.01, lines[i]*1.1, sspam.Data());
        }
    }

}

TLegend *newLegend(double x1, double y1, double x2, double y2) {
    TLegend *leg = new TLegend(x1,y1,x2,y2); 
    leg->SetFillColor(0);
    leg->SetShadowColor(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    return leg;
}

TPaveText *cmsprel;
void spam(const char *text=0, double x1=0.17, double y1=0.89, double x2=0.58, double y2=0.94, int textAlign=-1, bool fill=true, float fontSize=0) {
   //std::cout << "SPAM " << text << " x1 = " << x1  << ", y1 = " << y1 << ", x2 = " << x2 << ", y2 = " << y2 << std::endl;
   if (TString(text).Contains("#splitline")) { 
       if (y1 > 0.5) y1 -= 0.065+0.04*isTiny; else y2 += 0.065+0.04*isTiny; 
   } else if (TString(text).Contains("\n")) {
      char buff[20480]; sprintf(buff,text);
      double stride = (y1 > 0.5 ? -0.0475 : +0.0475);  if (fontSize) stride *= fontSize/(isSquareCanvas ? 0.04 : 0.045);
      int lines = 0; sprintf(buff,text);
      for (char *line = strtok(buff,"\n"); line != 0; line = strtok((char*)0,"\n")) {
        lines++;
      }
      if (stride > 0) { lines--; y1 += lines*stride; y2 += lines*stride; stride *= -1; }
      sprintf(buff,text);
      for (char *line = strtok(buff,"\n"); line != 0; line = strtok((char*)0,"\n")) {
         spam(line,x1,y1,x2,y2,textAlign==-1 ? 22 : textAlign,fill,fontSize);
         y1 += stride; y2 += stride;
      }
      return;
   }
   if (textAlign == -1) textAlign=12;
   //cmsprel = new TPaveText(x1,y1,x2,y2,"brtlNDC");
   cmsprel = new TPaveText(x1,y1,x2,y2,"brtlNDC");
   if (fontSize == 0) fontSize = (isSquareCanvas ? 0.04 : 0.045) + (isTiny ? 0.02 : 0);
   cmsprel->SetTextSize(fontSize);
   cmsprel->SetFillColor(0);
   cmsprel->SetFillStyle((fill || text == 0) ? 1001 : 0);
   cmsprel->SetLineStyle(0);
   cmsprel->SetLineColor(0);
   cmsprel->SetLineWidth(0);
   cmsprel->SetTextAlign(textAlign);
   cmsprel->SetTextFont(42);
   //cmsprel->SetTextColor(text.Contains("PRIVATE") ? 205 : 1);
   cmsprel->SetTextColor( TString(text ? text : SPAM.Data()).Contains("PRIVATE") ? 205 : 1 );
   cmsprel->AddText(text ? text : SPAM.Data());
   cmsprel->SetBorderSize(0);
   cmsprel->Draw("same");
}

void finalizeNoSave(TString name, double xmin, double xmax, double ymin, double ymax, const char *tspam=0, bool spamLow=false) {
    TLine line(xmin,1,xmax,1); 
    line.SetLineColor(lineAt1Color); 
    line.SetLineStyle(lineAt1Style); 
    line.SetLineWidth(4);
    if (name.Contains("smcls") || name.Contains("smacls")) { 
        line.SetY1(0.05); line.SetY2(0.05); 
        line.DrawClone();
        //line.SetLineStyle(7); line.SetLineWidth(2);
        //line.DrawLine(xmin,1-0.68,xmax,1-0.68);
        line.DrawLine(xmin,1-0.90,xmax,1-0.90);
        line.DrawLine(xmin,0.01,  xmax,0.01);
        //line.DrawLine(xmin,0.003, xmax,0.003);
    } else if (name.Contains("pval")) {
        if (name.Contains("_band")) {
            TLine threeSig(xmin, ROOT::Math::normal_cdf_c(3), xmax, ROOT::Math::normal_cdf_c(3)); 
            threeSig.SetLineColor(lineAt1Color); threeSig.SetLineWidth(2); threeSig.SetLineStyle(c1->GetGridy() ? 7 : lineAt1Style); 
            TLatex latex; latex.SetTextFont(42); latex.SetTextSize(0.045); latex.SetTextColor(2);
            bool protrude = c1->GetGridy(); 
            for (double z = 1; z <= 8; ++z) {
                if (ymin >  ROOT::Math::normal_cdf_c(z)) break;
                threeSig.DrawLine(xmin, ROOT::Math::normal_cdf_c(z), xmax+(xmax-xmin)*0.015*protrude, ROOT::Math::normal_cdf_c(z));
                TString sspam = TString::Format("%d#sigma", int(z));
                latex.DrawLatex(xmax+(xmax-xmin)*0.01, ROOT::Math::normal_cdf_c(z)*1.1, sspam.Data());
            }
        }
    } else {
        if (lineAt1Style) line.DrawClone();
        if (gPad->GetGridy() == 0 && ymin < 0) {
            line.SetLineColor(kBlack); line.SetLineWidth(2); line.SetLineStyle(lineAt1Style);
            line.DrawLine(xmin,0,xmax,0);
        }
    }
    if (gPad->GetLogx() && xmin <= 100 && xmax >= 200) {    
        TLine tick; tick.SetLineWidth(1); tick.SetLineColor(1);
        double dyh = ymax * 0.08;
        double dyl = ymin * 0.08; //fabs(c1->PixeltoY(c1->VtoPixel(0.95)) - c1->PixeltoY(c1->VtoPixel(0.94)));
        if (gPad->GetLogy() && log(ymax/ymin) > log(1e6)) { dyh *= 2; dyl *= 2; }
        if (gPad->GetLogy() == 0) { dyh = dyl = 0.01*(ymax-ymin); }
        if (isTiny) { dyh *= 2; dyl *= 2; }
        for (int i = 100; i < xmax; i += 10)  {
            if (i > 400 && i % 20 == 10) continue;
            tick.DrawLine(i, ymin, i, ymin+(i % 100 == 0 ? 2*dyl : dyl)); 
            tick.DrawLine(i, ymax, i, ymax-(i % 100 == 0 ? 2*dyh : dyh)); 
        }
    }
    if (leg) leg->Draw();
    double xstart = 0.17 - 0.05*isTiny;
    double xend   = (isSquareCanvas ? .60 : .64);
    if (name.Contains("smcls") || name.Contains("smacls")) { xend -= 0.08; }
    if (TString(tspam).Contains("\n")) xend -= 0.005;
    if (tspam) spam(tspam, xstart,
                           0.87-0.71*spamLow + 0.04*isTiny, 
                           xend, 
                           0.93-0.71*spamLow + 0.04*isTiny);
}

void justSave(TString name, double xmin=0., double xmax=0., double ymin=0., double ymax=0., const char *tspam=0, bool spamLow=false) {
    c1->Print(globalPrefix+name+".eps");
    //c1->Print(globalPrefix+name+".png");
    //gSystem->Exec("convert "+globalPrefix+name+".eps "+globalPrefix+name+".png");
    TString convOpt = "-q  -dBATCH -dSAFER  -dNOPAUSE  -dAlignToPixels=0 -dEPSCrop  -dPrinted -dTextAlphaBits=4 -dGraphicsAlphaBits=4 -sDEVICE=png16m";
    TString convCmd = Form("gs %s -sOutputFile=%s.png -q \"%s.eps\" -c showpage -c quit", convOpt.Data(), (globalPrefix+name).Data(), (globalPrefix+name).Data());
    gSystem->Exec(convCmd);
    if (name.Contains("pval_ml_")) {
        gSystem->Exec("epstopdf "+globalPrefix+name+".eps --outfile="+globalPrefix+name+".pdf");
    } else {
        c1->Print(globalPrefix+name+".pdf");
    }
}

void finalize(TString name, double xmin, double xmax, double ymin, double ymax, const char *tspam=0, bool spamLow=false) {
    finalizeNoSave(name,xmin, xmax, ymin, ymax, tspam, spamLow);
    justSave(name,xmin, xmax, ymin, ymax, tspam, spamLow);
}

void minMaxY(TGraphAsymmErrors *a, double &ymin, double &ymax, double xmax=999, double hardymin=-999) {
    if (hardymin == -999) hardymin = (SM == "SM" ? 0.08 : 0.02); 
    if (a == 0) return;
    ymin = 0.6; ymax = (hardymin < 0 ? 1 : 12); double yavg = 1.; double npoints = 0;
    double xhalf = 0.5*TMath::Min(a->GetX()[a->GetN()-1],xmax);
    for (int i = 0, n = a->GetN(); i < n; ++i) {
        double yhi = a->GetY()[i] + a->GetErrorYhigh(i);
        double ylo = a->GetY()[i] - a->GetErrorYlow(i);
        if (a->GetX()[i] > xmax) continue;
        npoints++;
        if (ylo * yhi > 0) yavg *= (ylo*yhi);
        if (yhi*3   > ymax) ymax = yhi * 3;
        if (a->GetX()[i] > xhalf && yhi*6 > ymax) ymax = yhi * 6; 
        if (ylo*0.6 < ymin) ymin = ylo * 0.6;
    }
    if (ymin < hardymin) ymin = hardymin;
    yavg = pow(yavg, 0.5/double(npoints));
    if (ymax < 10*yavg) ymax = 10*yavg;
    if (forceYmin != 0) ymin = forceYmin;
}


void squareCanvas(bool gridx=1,bool gridy=1) {
    gStyle->SetCanvasDefW(600); //Width of canvas
    gStyle->SetPaperSize(20.,20.);
    c1->Close();
    c1 = new TCanvas("c1","c1");
    c1->cd();
    c1->SetWindowSize(600 + (600 - c1->GetWw()), 600 + (600 - c1->GetWh()));
    c1->SetRightMargin(0.05);
    c1->SetGridy(gridy); c1->SetGridx(gridx);
    isSquareCanvas = true;
}

void rectangleCanvas(bool gridx=1,bool gridy=1) {
    gStyle->SetCanvasDefW(850); //Width of canvas
    gStyle->SetPaperSize(950./500.*20.,20.);
    c1->Close();
    c1 = new TCanvas("c1","c1");
    c1->cd();
    c1->SetWindowSize(950 + (950 - c1->GetWw()), 500 + (500 - c1->GetWh()));
    c1->SetRightMargin(0.04);
    c1->SetGridy(gridy); c1->SetGridx(gridx);
    isSquareCanvas = false;
}


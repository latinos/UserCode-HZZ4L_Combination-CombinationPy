TCanvas *c1 = new TCanvas("c1","c1");
TLegend *leg = 0;
TString globalPrefix = "";

TString SPAM = "CMS Preliminary,  #sqrt{s} = 7 TeV";
TString SPAM2L = "CMS Preliminary\n#sqrt{s} = 7 TeV";
//TString SPAM = "CMS Priv,  #sqrt{s} = 7 TeV";
//TString SPAM = "PRIVATE,  #sqrt{s} = 7 TeV";
TString SM = "SM";
TString oneSigmaText="(68%)"; //"#pm 1#sigma";
TString twoSigmaText="(95%)"; //"#pm 2#sigma";
TString oneSigmaFitText="68% CL band";
TString oneSigmaFitCCC=" (68%)";
TString massUnits="GeV";
TString lumiSymbol="L"; // L_{int}
bool justLumiForCombined = true;
bool doSquares = true;
bool noLineStyles = false;
bool channelSpamOnRightHandSide = false;
bool doLEESPAM = false;
TString LEESPAM_1L = "NOT USED";
TString LEESPAM_2L = "NOT USED";
TString LEESPAM_3L = "Global significance\n0.8#sigma for 110-600 GeV range\n2.1#sigma for 110-145 GeV range";

bool isSquareCanvas = false, isTiny = false, lessSpacesInLegends = false;
bool track_missing = false;

double x_zoom = 145.;
bool doZoom2 = false;
double x_zoom2_min = 160, x_zoom2_max = 300;
double forceYmin = 0; //0.08;
double forceYmax = 0; //12.8;

const int nchann_all = 1+5+3+3+7;
const int nchann = 14, nchann2 = 6, nchann3 = 3, nchann4 = 13;
const char *chann[nchann_all] = { "comb",
                                  "vhbb", "htt", "hgg", "hww", "hzz", 
                                  "combp", "combs", "combl",
                                  "httm", "vhtt", "vhww3l", 
                                  "htt0", "hww2l", "hzz4l", "hzz2l2t", "hzz2l2q_low", "hzz2l2q", "hzz2l2nu"
                                 };
const char *chann2[nchann2] = { "comb", "vhbb", "htt", "hgg", "hww", "hzz" };
const char *chann3[nchann3] = { "comb", "combp_full", "combl_full" };
const char *chann4[nchann4] = { "comb", "vhbb", "htt0", "httm", "vhtt", "hgg", "hww2l", "vhww3l", "hzz4l", "hzz2l2t", "hzz2l2q_low", "hzz2l2q", "hzz2l2nu" };
const int color68 = 80, color95 = 90, color50 = TColor::GetColor(20,20,20);
const int colorFit68 = 80;
int lineAt1Style = 1;
int lineAt1Color = 2;


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
    if (tev_excluded && xmax0 > 156) {
        double xmin = 156, xmax = TMath::Min(177., xmax0);
        TBox box(xmin, ymin, xmax, ymax);
        box.SetLineStyle(0);
        box.SetFillStyle(3354);
        box.SetFillColor(213);
        box.DrawClone();
        TLine line(xmin, ymin, xmin, ymax);
        line.SetLineWidth(2);
        line.SetLineColor(213);
        line.DrawClone();
        if (xmax0 > 177) line.DrawLine(xmax, ymin, xmax, ymax);
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
    if (name.Contains("smcls")) { 
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
        line.DrawClone();
    }
    if (gPad->GetLogx() && xmin <= 100 && xmax >= 600) {    
        TLine tick; tick.SetLineWidth(1); tick.SetLineColor(1);
        double dyh = ymax * 0.08;
        double dyl = ymin * 0.08; //fabs(c1->PixeltoY(c1->VtoPixel(0.95)) - c1->PixeltoY(c1->VtoPixel(0.94)));
        if (gPad->GetLogy() && log(ymax/ymin) > log(1e6)) { dyh *= 2; dyl *= 2; }
        if (gPad->GetLogy() == 0) { dyh = dyl = 0.01*(ymax-ymin); }
        if (isTiny) { dyh *= 2; dyl *= 2; }
        for (int i = 100; i < 600; i += 10)  {
            if (i > 400 && i % 20 == 10) continue;
            tick.DrawLine(i, ymin, i, ymin+(i % 100 == 0 ? 2*dyl : dyl)); 
            tick.DrawLine(i, ymax, i, ymax-(i % 100 == 0 ? 2*dyh : dyh)); 
        }
    }
    if (leg) leg->Draw();
    if (name.Contains("pval")) {   
        if (name.Contains("_all")) {
            //spam("#splitline{Look-elsewhere effect}{not included}", 0.17 - 0.05*isTiny, 0.21 + 0.10*isTiny, 0.47 - 0.05*isTiny, 0.26 + 0.10*isTiny);
            bool isZoom = name.Contains("zoom");
            if (isSquareCanvas) {
                if (doLEESPAM) {
                //if (SM == "SM") spam(LEESPAM_2L, 0.19, 0.35, 0.59, 0.40, 22, true);
                if (SM == "SM") spam(LEESPAM_3L, 0.19, 0.35, 0.59, 0.40, 22, true, 0.034);
                if (SM == "FP") spam(LEESPAM_2L, 0.19, 0.35, 0.59, 0.40, 22, true);
                if (SM == "SM4") spam(LEESPAM_2L, 0.19, 0.35, 0.59, 0.40, 22, true);
                }
                spam(tspam, .17, .15-0.00*(!doLEESPAM), .49, .20-0.00*(!doLEESPAM), 22);
            } else {
                //if (doLEESPAM) spam(LEESPAM_1L, 0.175 - 0.05*isTiny, 0.15 + 0.10*isTiny, 0.59 - 0.05*isTiny, 0.20 + 0.10*isTiny);
                //if (doLEESPAM) spam(LEESPAM_1L, 0.175, 0.35 , 0.59, 0.40);
                if (SM == "SM") spam(LEESPAM_3L, 0.19, 0.35, 0.51, 0.40, 22, true, 0.04);
                if (SM == "FP") spam(LEESPAM_2L, 0.19, 0.35, 0.51, 0.40, 22, true);
                if (SM == "SM4") spam(LEESPAM_2L, 0.19, 0.35, 0.51, 0.40, 22, true);
                spam(tspam, .17, .15-0.000*(!doLEESPAM), .59, .20-0.000*(!doLEESPAM), 22);
            }
            tspam = 0;
        /*} else if (name.Contains("_comb")) {
            spam(tspam, .575, .415, .945, .465, 22 );
            spam("Interpretation requires look-elsewhere effect correction", 0.17 - 0.05*isTiny, 0.06 + 0.00*isTiny, 0.67 - 0.05*isTiny, 0.21 + 0.00*isTiny);
            tspam = 0;*/
        } else if (tspam) {
            if (isTiny && isSquareCanvas) {
                spam(tspam, 0.17, 0.18-0.05*(!doLEESPAM), 0.57, 0.21-0.05*(!doLEESPAM));
                if (doLEESPAM) spam("Interpretation requires look-elsewhere effect correction", 0.17 + 0.05*isTiny, 0.18 - 0.12*isTiny, 0.94 - 0.12*isTiny, 0.23 - 0.10*isTiny, 22);
            } else if (TString(tspam).Contains("#splitline")) {
                spam(tspam, 0.17, 0.22-0.05*(!doLEESPAM), 0.57, 0.25-0.05*(!doLEESPAM));
                if (doLEESPAM) spam("#splitline{Interpretation requires look-}{elsewhere effect correction}", 0.55, 0.34, 0.93, 0.38);
            } else {
                spam(tspam, 0.17 + 0.00*isTiny, 0.23 - 0.05*isTiny-0.05*(!doLEESPAM), 0.94 - 0.12*isTiny, 0.29 - 0.05*isTiny, 22);
                if (doLEESPAM) spam("Interpretation requires look-elsewhere effect correction", 0.17 + 0.00*isTiny, 0.18 - 0.12*isTiny, 0.94 - 0.12*isTiny, 0.23 - 0.10*isTiny, 22);
                if (name.Contains("_combe")) {
                    spam("#splitline{Using event-by-event mass}{resolution in H #rightarrow ZZ #rightarrow 4l}", 0.55, 0.34, 0.93, 0.38);
                }
            }
            tspam = 0;
        }
    }
    /*
    if (SM != "SM") {
        if (name.Contains("cls")) {
            spam(SM == "SM4" ? "#splitline{Standard model with 4}{generations of fermions}" : 
                               "Fermiophobic model",
                 0.17, 0.75, 0.46, 0.82);
        } else if (name.Contains("pval")) {
            spam(SM == "SM4" ? "SM with 4 generations of fermions" :
                               "Fermiophobic model", 
                 0.17, 0.15, 0.63, 0.20);
        }
    }
    */ 
    bool isFullComb = name.Contains("_all") || (channelFromName(name,false,true) == "Combined");
    if (tspam) spam(tspam, 0.17 - 0.05*isTiny,
                           0.89-0.74*spamLow + 0.04*isTiny-0.005*TString(tspam).Contains("\n"), 
                           (isSquareCanvas ? .62 : .64)+(name.Contains("pvala") || name.Contains("smcls") ? -0.07: 0)- 0.15*isFullComb-0.005*TString(tspam).Contains("\n"), 
                           0.94-0.74*spamLow + 0.04*isTiny);
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

TGraphAsymmErrors *slidingWindowAverage(TGraphAsymmErrors *input, int slidingWindow) {
    TGraphAsymmErrors *out  = (TGraphAsymmErrors *) input->Clone();
    bool isLogPlot = TString(input->GetName()).Contains("smcls");
    for (int i = 0, n = input->GetN(); i < n; ++i) {
        /*if (i >= 1 && i < n-1 && (input->GetX()[i] < input->GetX()[i-1] || input->GetX()[i] > input->GetX()[i+1])) {
            out->GetX()[i] = 0.5*(input->GetX()[i-1] + input->GetX()[i+1]);
            out->GetY()[i] = 0.5*(input->GetY()[i-1] + input->GetY()[i+1]);
            out->GetEYlow()[i] = 0.5*(input->GetEYlow()[i-1] + input->GetEYlow()[i+1]);
            out->GetEYhigh()[i] = 0.5*(input->GetEYhigh()[i-1] + input->GetEYhigh()[i+1]);
            continue;
        }*/
        double y0 = input->GetY()[i];
        double sum = 0, sumhi = 0, sumlo = 0, sumw = 0;
        for (int j = i-slidingWindow; j <= i+slidingWindow; ++j) {
            if (j < 0 || j >= n) continue;
            double y = input->GetY()[j], w = 1.0; // /(1.0 + abs(i-j));
            if (isLogPlot) {
                if (fabs(log(y0/y)) > log(2)) continue;    
            } else {
                if (fabs(y0-y) > 0.1*y0) continue;    
            }
            if (isLogPlot) {
                sum   += w*log(y);
                sumlo += w*log(input->GetEYlow()[j]);
                sumhi += w*log(input->GetEYhigh()[j]);
            } else {
                sum   += w*y;
                sumlo += w*input->GetEYlow()[j];
                sumhi += w*input->GetEYhigh()[j];
            }
            sumw  += w;
        }
        if (sumw == 0) continue;
        if (isLogPlot) {
            out->GetY()[i] = exp(sum/sumw);
            out->GetEYlow()[i] = exp(sumlo/sumw);
            out->GetEYhigh()[i] = exp(sumhi/sumw);
        } else {
            out->GetY()[i] = sum/sumw;
            out->GetEYlow()[i] = sumlo/sumw;
            out->GetEYhigh()[i] = sumhi/sumw;
        }
    }
    return out;
}

TGraphAsymmErrors *smoothSMCLs(TGraphAsymmErrors *input, int slidingWindow, int order) {
    double xi[99], yi[99];
    TString who = TString(input->GetName());
    bool nonzero   = true; //TString(input->GetName()).Contains("cls");
    bool isLogPlot = (who.Contains("smcls") || who.Contains("pval"));
    bool absolute  = (who.Contains("smcls") || who.Contains("pval"));
    double maxdel  = (who.Contains("smcls") || who.Contains("pval")) ? 0 : log(1.25); 
    TGraphAsymmErrors *inp = (TGraphAsymmErrors *) input->Clone();
    TGraphAsymmErrors *out = (TGraphAsymmErrors *) input->Clone();
    if (absolute) {
        for (int i = 0, n = input->GetN(); i < n; ++i) {
            inp->GetEYlow() [i] = inp->GetY()[i] - inp->GetEYlow()[i];
            inp->GetEYhigh()[i] = inp->GetY()[i] + inp->GetEYhigh()[i];
            out->GetEYlow() [i] = out->GetY()[i] - out->GetEYlow()[i];
            out->GetEYhigh()[i] = out->GetY()[i] + out->GetEYhigh()[i];
        }
    } 
    for (int which = -1; which <= +1; ++which) {
        double *inY = (which == 0 ? inp->GetY() : (which > 0 ? inp->GetEYhigh() : inp->GetEYlow()));
        double *ouY = (which == 0 ? out->GetY() : (which > 0 ? out->GetEYhigh() : out->GetEYlow()));
        for (int i = 0, n = input->GetN(); i < n; ++i) {
            double y0 = inY[i], x0 = input->GetX()[i]; int points = 0;
            if (nonzero && y0 == 0) continue;
            for (int j = i-slidingWindow; j <= i+slidingWindow; ++j) {
                if (j < 0 || j >= n) continue;
                if (nonzero && inY[j] == 0) continue;
                if (maxdel > 0) {
                    double del = inp->GetY()[j]/inp->GetY()[i]; 
                    if (fabs(log(del)) > maxdel) continue;
                }
                xi[points] = input->GetX()[j];
                yi[points] = isLogPlot ? log(inY[j]/y0) : inY[j];
                points++;
            }
            double ynew = isLogPlot ? 1 : y0;
            if (points > order+2) {
                ynew = smoothWithPolyFit(x0, order+1, points, xi, yi);
            } else if (points > 2) {
                ynew = smoothWithPolyFit(x0, 1, points, xi, yi);
            } else if (maxdel == 0 && nonzero && points == 1) {
                ouY[i] = 0; continue; // kill the blip
            } else continue; // nothing to do
            ouY[i] = isLogPlot ? y0 * exp(ynew) : ynew;
        }
    }
    if (absolute) {
        for (int i = 0, n = input->GetN(); i < n; ++i) {
            out->GetEYlow() [i] = +out->GetY()[i] - out->GetEYlow()[i];
            out->GetEYhigh()[i] = -out->GetY()[i] + out->GetEYhigh()[i];
        }
    }
    delete inp;
    return out;
}

TGraphAsymmErrors * removeGlitches(TGraphAsymmErrors *out) {
    do {
        int bad = -1; bool hasgood = false;
        for (int i = 0; i < out->GetN(); ++i) {
            if (out->GetEYlow()[i] == 0 && out->GetEYhigh()[i] == 0) {
                bad = i;
            } else { 
                hasgood = true; 
            }
        }
        if (!hasgood) return out;
        if (bad == -1) return out;
        out->RemovePoint(bad); 
    } while (1);
}

int nmasses = 0;
double masses[200];
void loadMasses(const char *file="masses.txt") {
    FILE *f = fopen(file,"r");
    float mass;
    while (fscanf(f,"%f", &mass) == 1) {
        masses[nmasses++] = mass;
    }
}

TGraphAsymmErrors *missingPoints(TGraphAsymmErrors *points) {
    if (!track_missing) return 0;
    if (points == 0) return 0;
    int n = points->GetN(); Double_t *xj = points->GetX(); Double_t *yj = points->GetY();
    double xmin = xj[0], xmax = xj[n-1];
    bool logint = true;
    TGraphAsymmErrors *ret = new TGraphAsymmErrors(); ret->SetName(TString("missing_")+points->GetName());
    // check if we have half-integer points
    bool halfint = false;
    for (int i = 0; i < n; ++i) { 
        if (points->GetX()[i] - floor(points->GetX()[i]) > 0.4) { halfint = true; break; }
    } 
    int nmiss = 0;
    for (int i = 0; i < nmasses; ++i) {
        double x = masses[i];
        if (masses[i] < xmin || masses[i] > xmax) continue;
        if ((x - floor(x)) > 0.4 && !halfint) continue;
        bool found = false; double xlo = 0, ylo = 0, xhi = 0, yhi = 0;
        for (int j = 0; j < n; ++j) {
            if (xj[j] < x) { xlo = xj[j]; ylo = yj[j]; }
            else if (xj[j] == x) { found = true; break; }
            else if (xhi == 0) { xhi = xj[j]; yhi = yj[j]; break; }
        }
        if (!found) {
            double y = ( yhi * (x-xlo) + ylo * (xhi - x) ) / (xhi - xlo);
            if (yhi > 0 && ylo > 0 && logint) {
                y = exp( ( log(yhi) * (x-xlo) + log(ylo) * (xhi - x) ) / (xhi - xlo) );
            }
            //printf("Interpolated missing point %3d from (%4.0f, %7.3f) + (%4.0f, %7.3f) --> (%4.0f, %7.3f)\n", masses[i], xlo, ylo, xhi, yhi, x, y);
            nmiss++;
            ret->Set(nmiss);
            ret->SetPoint(nmiss-1, x, y);
        }
    }
    ret->SetMarkerStyle(20);
    ret->SetMarkerSize(0.7);
    ret->SetMarkerColor(100);
    if (nmiss == 0) { delete ret; return 0; }
    return ret;
}


TGraph *draw2(TString who, int fillColor68, int fillColor95, int lineColor, bool same=true, bool mean=false, int slidingWindow=0, int smoothorder=2) {
    TGraphAsymmErrors *mean68 = (TGraphAsymmErrors*) gROOT->FindObject(who+(mean?"_mean":"_median"));
    TGraphAsymmErrors *mean95 = (TGraphAsymmErrors*) gROOT->FindObject(who+(mean?"_mean":"_median")+"_95");
    if (mean68 == 0) { std::cerr << "MISSING " << who << (mean ? "_mean"    : "_median")    << std::endl; return 0; }
    if (mean95 == 0) { std::cerr << "MISSING " << who << (mean ? "_mean_95" : "_median_95") << std::endl; return 0; }
    mean68 = removeGlitches(mean68);
    mean95 = removeGlitches(mean95);
    if (mean68->GetN() == 1) {
        mean68->SetPointError(0, 10, 10, mean68->GetErrorYlow(0), mean68->GetErrorYhigh(0));
        mean95->SetPointError(0, 10, 10, mean95->GetErrorYlow(0), mean95->GetErrorYhigh(0));
    }
    TGraphAsymmErrors *meanL = (TGraphAsymmErrors*) mean68->Clone();
    for (int i= 0; i < meanL->GetN(); ++i) { meanL->SetPointError(i, 0,0,0,0); }

    if (slidingWindow != 0 && mean68->GetN() > 5) {
        if (slidingWindow > 0) {
            mean68 = slidingWindowAverage(mean68, slidingWindow);
            mean95 = slidingWindowAverage(mean95, slidingWindow);
            meanL  = slidingWindowAverage(meanL,  slidingWindow);
        } else {
            mean68 = smoothSMCLs(mean68, -slidingWindow, smoothorder);
            mean95 = smoothSMCLs(mean95, -slidingWindow, smoothorder);
            meanL  = smoothSMCLs(meanL,  -slidingWindow, smoothorder);
        }
    }

    meanL->SetLineColor(lineColor); mean68->SetLineColor(lineColor); mean95->SetLineColor(lineColor);
    meanL->SetMarkerColor(lineColor); mean68->SetMarkerColor(lineColor); mean95->SetMarkerColor(lineColor);
    meanL->SetLineWidth(3); mean68->SetLineWidth(3); mean95->SetLineWidth(3);
    meanL->SetMarkerSize(1.6); mean68->SetMarkerSize(1.6); mean95->SetMarkerSize(1.6);
    meanL->SetLineStyle(7);     mean68->SetLineStyle(7);     mean95->SetLineStyle(7); 
    mean68->SetLineColor(fillColor68); mean95->SetLineColor(fillColor95);
    mean68->SetLineWidth(1);           mean95->SetLineWidth(1);
    mean68->SetFillColor(fillColor68);  
    mean95->SetFillColor(fillColor95);
    mean95->Draw(same ? "E3 SAME" : "AE3");
    mean68->Draw("E3 SAME");
    meanL->Draw("LX SAME");
    return mean95;
}

TGraph *draw1(TString who, int lineColor, TString option="L") {
    TGraph *g = (TGraph*) gROOT->FindObject(who);
    if (g == 0) { std::cerr << "Graph " << who << " empty." << std::endl; return 0; }
    if (g->GetN() == 0) { std::cerr << "Graph " << who << " empty." << std::endl; return 0; }
    g->SetLineColor(lineColor);
    g->SetMarkerColor(lineColor);
    g->SetLineWidth(3); 
    g->Draw(option+" SAME");
    return g;
}

void minMaxY(TGraphAsymmErrors *a, double &ymin, double &ymax, double xmax=999, double hardymin=-999) {
    if (hardymin == -999) hardymin = (SM == "SM" ? 0.08 : 0.02); 
    if (a == 0) return;
    ymin = 0.6; ymax = (hardymin < 0 ? 1 : 12); double yavg = 1.; double npoints = 0;
    for (int i = 0, n = a->GetN(); i < n; ++i) {
        double yhi = a->GetY()[i] + a->GetErrorYhigh(i);
        double ylo = a->GetY()[i] - a->GetErrorYlow(i);
        if (a->GetX()[i] > xmax) continue;
        npoints++;
        if (ylo * yhi > 0) yavg *= (ylo*yhi);
        if (yhi*3   > ymax) ymax = yhi * 3;
        if (ylo*0.6 < ymin) ymin = ylo * 0.6;
    }
    if (ymin < hardymin) ymin = hardymin;
    yavg = pow(yavg, 0.5/double(npoints));
    if (ymax < 10*yavg) ymax = 10*yavg;
    if (forceYmin != 0) ymin = forceYmin;
}

TString channelFromName(TString who, bool withSpaces=false, bool noLumi=false) {
    TString name = who, space, lumi;
    if (who.Contains("comb"))  { name = "Combined"; space = "", lumi= "4.6-4.8 fb^{-1}"; }
    if (who.Contains("hzz"))  { 
        name = withSpaces ? "H #rightarrow ZZ" : "H #rightarrow ZZ Comb."; 
        space = "              "; lumi= "4.7 fb^{-1}"; 
    }
    if (who.Contains("combp"))  { 
        name = isSquareCanvas ? "ZZ + #gamma#gamma" : "H #rightarrow ZZ + #gamma#gamma"; 
        space = "     "; lumi= "4.8 fb^{-1}"; 
        if (withSpaces) { space = ""; lumi = ""; }
    }
    if (who.Contains("combl"))  { 
        name = isSquareCanvas ? "bb + #tau#tau + WW" : "H #rightarrow bb + #tau#tau + WW" ; 
        space = "   "; lumi= "4.6 fb^{-1}"; 
        if (withSpaces) { space = ""; lumi = ""; }
    }
    if (who.Contains("combs"))  { 
        name = isSquareCanvas ? "VBF+VH excl." : "VBF+VH exclusive" ; 
        space = "   "; lumi= "4.6-4.8 fb^{-1}"; 
        if (withSpaces) { space = ""; lumi = ""; }
    }
    if (who.Contains("hww"))   { name = "H #rightarrow WW"; space = "           "; lumi="4.6 fb^{-1}"; }
    if (who.Contains("hww2l"))   { name = "H #rightarrow WW #rightarrow 2l 2#nu"; space = ""; lumi="4.6 fb^{-1}"; }
    if (who.Contains("vhww3l"))   { name = "WH #rightarrow 3l 3#nu"; space = "     "; lumi="4.6 fb^{-1}"; }
    if (who.Contains("vhbb"))  { name = "H #rightarrow bb"; space = "             "; lumi="4.7 fb^{-1}"; }
    if (who.Contains("hgg_vbfonly"))   { name = "VBF H #rightarrow #gamma#gamma"; space = "        "; lumi="4.8 fb^{-1}"; }
    if (who.Contains("hgg"))   { name = "H #rightarrow #gamma#gamma"; space = "              "; lumi="4.8 fb^{-1}"; }
    if (who.Contains("htt"))   { name = "H #rightarrow #tau#tau"; space = "              "; lumi="4.6 fb^{-1}"; }
    if (who.Contains("vhtt"))   { name = "VH #rightarrow #tau_{h} 2l"; space = "         "; lumi="4.6 fb^{-1}"; }
    if (who.Contains("httm"))   { name = "H #rightarrow #tau#tau #rightarrow #mu#mu"; space = "    "; lumi="4.6 fb^{-1}"; }
    if (who.Contains("hzz4l")) { name = "H #rightarrow ZZ #rightarrow 4l"; space = "     "; lumi="4.7 fb^{-1}"; }
    if (who.Contains("hzz2l2q"))  { name = "H #rightarrow ZZ #rightarrow 2l 2q"; space=""; lumi="4.6 fb^{-1}"; }
    if (who.Contains("hzz2l2nu")) { name = "H #rightarrow ZZ #rightarrow 2l 2#nu"; space=""; lumi="4.6 fb^{-1}"; }
    if (who.Contains("hzz2l2t")) { name = "H #rightarrow ZZ #rightarrow 2l 2#tau"; space=""; lumi="4.6 fb^{-1}"; }
    if (lessSpacesInLegends) {
        int nsp = space.Length()-12; space = "";
        for (int i = 0; i < nsp; ++i) space += " ";
    }
    if (noLumi) return name;
    if (withSpaces) return (name == "Combined" || lumi == "") ? name : (name + space + "  ("+lumi+")");
    if (name == "Combined" && justLumiForCombined || channelSpamOnRightHandSide) return lumiSymbol+" = "+lumi;
    return name + ", "+lumiSymbol+" = "+lumi;
}

int colorFromName(TString who, bool dark=false) {
    if (who.Contains("combp")) return (dark ? 205 : 2);
    if (who.Contains("combl")) return (dark ? 213 : 4);
    if (who.Contains("comb"))  return (dark ? 19 : 1);
    if (who.Contains("hww2l"))   return (dark ? 215 : 215);
    if (who.Contains("vhww"))   return (dark ? 213 : 213);
    if (who.Contains("hww"))   return (dark ? 4 : 4);
    if (who.Contains("hgg"))   return (dark ? 209 : 209);
    if (who.Contains("vhbb"))  return (dark ? 67 : 67 );
    if (who.Contains("vhtt"))   return (dark ? 51 : 51);
    if (who.Contains("httm"))   return (dark ? 223 : 223);
    if (who.Contains("htt"))   return (dark ? 221 : 221);
    if (who.Contains("hzz4l"))  return (dark ? 2 : 2);
    if (who.Contains("hzz2l2q"))  return (dark ? 93 : 93);
    if (who.Contains("hzz2l2nu")) return (dark ? 28 : 28);
    if (who.Contains("hzz2l2t"))  return (dark ? 223 : 223);
    if (who.Contains("hzz")) return (dark ? 205 : 2);
    return 39;
}
int lineStyleFromName(TString who) {
    if (noLineStyles) return 1;
    if (who.Contains("combzz"))   return 1;
    if (who.Contains("comb"))     return 1;
    if (who.Contains("vhbb"))     return 3;
    if (who.Contains("hww"))      return 7;
    if (who.Contains("hgg"))      return 5;
    if (who.Contains("htt"))      return 2;
    if (who.Contains("hzz4l"))    return 1;
    if (who.Contains("hzz2l2q"))  return 5;
    if (who.Contains("hzz2l2nu")) return 2;
    if (who.Contains("hzz2l2t"))  return 3;
    return 39;
}

void drawOnePlot(TString who, TString what="auto") {
    TGraphAsymmErrors *obs = (TGraphAsymmErrors *) gFile->Get(who+"_obs");
    if (obs == 0) { std::cout << "Missing "+who << std::endl; return; }
    TGraphAsymmErrors *miss = missingPoints(obs);
    double xmin = obs->GetX()[0]             - obs->GetErrorXlow(0), xmin0 = xmin;
    double xmax = obs->GetX()[obs->GetN()-1] + obs->GetErrorXhigh(obs->GetN()-1);
    if (xmin <= 120 && xmax > 200) { xmin = 99.98; xmax = (SM == "FP" ? 300.1 : 600.1); }
    int col = 1; //colorFromName(who);
    obs->SetLineWidth(2);
    obs->SetLineColor(col); obs->SetMarkerColor(col); 
    obs->SetMarkerStyle(21); obs->SetMarkerSize(0.8);
    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    double ymin, ymax; minMaxY(obs, ymin, ymax);
    bool isML = who.Contains("ml_") || who.Contains("mlz_");
    if (who.Contains("ml_"))  ymin = 0.02;
    if (who.Contains("mlz_")) { 
        ymin = -2.5; ymax = 5; 
        if (who.Contains("comb") && !who.Contains("combzz")) { ymin = -1; ymax = 2.5; }
    }
    if (isML) { obs->SetFillColor(colorFit68); obs->Draw("E3"); frame0.Draw("AXIGSAME");}
    obs->Draw("LPX");    
    if (miss) miss->Draw("P");
    if (what == "auto") {
        if (who.Contains("pla_")) what = "P.L. Approx limit #sigma_{95%}/#sigma_{"+SM+"}";
        if (who.Contains("bayes_")) what = "Bayesian limit #sigma_{95%}/#sigma_{"+SM+"}";
        if (isML) what = "Best fit #sigma/#sigma_{"+SM+"}";
    }
    leg = 0;
    TH1F dummyBox("dummyBox","dummyBox",1,0.,1.); dummyBox.SetFillColor(colorFit68); dummyBox.SetLineColor(colorFit68);
    if (isML) {
        bool notJustComb = (who != "mlz_comb");
        leg = newLegend(.65-0.09*isSquareCanvas+0.07*isSquareCanvas*notJustComb,.85,.90+0.01*isSquareCanvas,.92); 
        leg->SetTextSize(0.045-0.000*isSquareCanvas-0.007*isSquareCanvas*notJustComb);
        leg->AddEntry(&dummyBox, oneSigmaFitText, "F");
        leg->Draw();
    }
    setCanvas(&frame0, "", ymin, ymax, what);
    TString myspam = "#splitline{"+SPAM+"}{"+channelFromName(who)+"}";
    if (isSquareCanvas && SPAM.Contains("Preliminary")) myspam = SPAM2L + "\n" + channelFromName(who);
    finalize(who,xmin,xmax,ymin,ymax, myspam);
    if (xmax > x_zoom && xmin < x_zoom) {
        c1->SetLogx(1);
        frame0.GetXaxis()->SetMoreLogLabels(); frame0.GetXaxis()->SetNoExponent();
        finalize(who+"_logx",xmin,xmax,ymin,ymax, myspam);
        c1->SetLogx(0);
        xmin = xmin0;
        TH1D frame2("frame2","frame2", 1, xmin, x_zoom); frame2.Draw(); 
        obs->SetFillColor(colorFit68);
        if (isML) { obs->SetFillColor(colorFit68); obs->Draw("E3"); frame2.Draw("AXIGSAME"); }
        obs->Draw("LXP");  
        if (miss) miss->Draw("P");
        if (leg) leg->Draw();
        if (!who.Contains("mlz_")) minMaxY(obs, ymin, ymax, x_zoom); 
        if (who.Contains("ml_")) ymin = 0.02; 
        setCanvas(&frame2, "", ymin, ymax, what);
        finalize(who+"_zoom",xmin,x_zoom,ymin,ymax, myspam);
        if (doZoom2) {
            TH1D frame3("frame3","frame3", 1, x_zoom2_min,x_zoom2_max); frame3.Draw(); 
            obs->SetFillColor(colorFit68);
            if (isML) { obs->SetFillColor(colorFit68); obs->Draw("E3"); frame3.Draw("AXIGSAME"); }
            obs->Draw("LXP");  
            if (miss) miss->Draw("P");
            if (leg) leg->Draw();
            setCanvas(&frame3, "", ymin, ymax, what);
            finalize(who+"_zoom2", x_zoom2_min,x_zoom2_max, ymin,ymax, myspam);
        }
    }
    leg = 0;
}

bool allow_halfint_pvals = true;
double minPValue(TGraph *g, double x_min, double x_max, int minPoints=0, double spikeKiller=0) {
    int points = 0, imin = 0; double pmin = 1.0;
    for (int i = 0, n = g->GetN(); i < n; ++i) {
        if (g->GetX()[i] < x_min || g->GetX()[i] > x_max) continue;
        if (!allow_halfint_pvals && float(int(g->GetX()[i])) != g->GetX()[i]) continue;
        if (spikeKiller > 0) {
            double Z  = ROOT::Math::normal_quantile_c(g->GetY()[i],1);
            double ZP = (i > 0)   ? (ROOT::Math::normal_quantile_c(g->GetY()[i-1],1.)) : 0.0;
            double ZN = (i < n-1) ? (ROOT::Math::normal_quantile_c(g->GetY()[i+1],1.)) : 0.0;
            if (Z - TMath::Min(ZP,ZN) > spikeKiller) continue;
        }
        points++;
        if (g->GetY()[i] < pmin) {
            pmin = g->GetY()[i]; imin = i;
        }
    }
    if (points < minPoints) { 
        std::cout << "Plot " << g->GetName() << " has " << points << " points, less than " << minPoints << std::endl; 
        return -1; 
    }
    return pmin;
}
void drawOnePVal(TString who, TString what="Local p-value", bool drawToys=false) {
    TGraphAsymmErrors *obs = (TGraphAsymmErrors *) gFile->Get(who+"_obs");
    if (obs == 0) return;
    TGraphAsymmErrors *exp = (TGraphAsymmErrors *) gFile->Get(who+"_median");
    if (exp == 0) exp = (TGraphAsymmErrors *) gFile->Get(who+"_asimov");

    TGraphAsymmErrors *toys = 0, *toysExp = 0;
    if (drawToys && who.Index("pvala") == 0) {
        TString twho = who; twho.ReplaceAll("pvala","pval");
        toys = (TGraphAsymmErrors *) gFile->Get(twho+"_obs");
        if (toys == 0) return;
        toys->SetLineWidth(4);
        toys->SetMarkerStyle(20);
        toys->SetMarkerSize(1.4);
        toys->SetLineColor(62);
        toys->SetMarkerColor(62);
        toysExp = (TGraphAsymmErrors *) gFile->Get(twho+"_median");
        if (toysExp != 0) {
            toysExp->SetLineWidth(4);
            toysExp->SetMarkerStyle(27);
            toysExp->SetMarkerSize(1.2);
            toysExp->SetLineColor(51);
            toysExp->SetMarkerColor(51);
        }
    }

    TGraphAsymmErrors *miss = toys ? 0 : missingPoints(obs);
    double xmin = obs->GetX()[0]             - obs->GetErrorXlow(0), xmin0 = xmin;
    double xmax = obs->GetX()[obs->GetN()-1] + obs->GetErrorXhigh(obs->GetN()-1);
    if (xmin <= 120 && xmax > 200) { xmin = 99.98; xmax = (SM == "FP" ? 300.1 : 600.1); }
    int col = 1; //colorFromName(who);
    obs->SetLineWidth(toys ? 3 : 2);
    obs->SetLineColor(col); obs->SetMarkerColor(col); 
    obs->SetMarkerStyle(21); obs->SetMarkerSize(0.9);
    if (exp) {
        TGraphAsymmErrors *smooth = smoothSMCLs(exp, 7, 2);
        if (smooth == 0) std::cout << "Smoothing of " << expi[i]->GetName() << " returned ZERO" << std::endl;
        else if (SM != "SM4") exp = smooth;
        exp->SetLineWidth(2); exp->SetLineStyle(1); exp->SetLineColor(64);
    }
    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    double ymin, ymax; minMaxY(obs, ymin, ymax, 999, 0); 
    ymax = 1.0; ymin = minPValue(obs,xmin,xmax)/10; if (ymin <= 0 || ymin > 1e-5) ymin = 1e-5;
    
    setCanvas(&frame0, "", ymin, ymax, what);
    frame0.GetYaxis()->SetTitleOffset(1.10+0.25*isSquareCanvas);
    if (toys) toys->Draw("P");
    if (toysExp) toysExp->Draw("P");
    obs->Draw(toys ? "LX" : "LXP");    
    if (exp) exp->Draw("LX");
    if (miss) miss->Draw("P");
    TString myspam = (toys ? "#splitline{"+SPAM+"}{"+channelFromName(who)+"}" : SPAM+", "+channelFromName(who)) ;
    if (toys) {
        leg = newLegend(.66-0.065*isSquareCanvas,.18,.93,.32+((exp!=0)+(toysExp!=0))*0.04); leg->SetTextSize(0.04-0.003*isSquareCanvas);
        leg->AddEntry(obs,  "Asymptotic Obs.", "L");
        leg->AddEntry(toys, "Ensemble Obs.",   "P");
        if (toysExp) leg->AddEntry(toysExp, "Ensemble Exp.",   "P");
        if (exp) leg->AddEntry(exp, "Asymptotic Exp.", "L");
    } else {
        leg = 0;
    }
    
    finalize(who+(drawToys?"_comp":""),xmin,xmax,ymin,ymax,myspam,true);
    if (xmax > x_zoom && xmin < x_zoom) {
        c1->SetLogx(1);
        frame0.GetXaxis()->SetMoreLogLabels(); frame0.GetXaxis()->SetNoExponent();
        finalize(who+(drawToys?"_comp_logx":"_logx"),xmin,xmax,ymin,ymax,myspam,true);
        c1->SetLogx(0);
        xmin = xmin0;
        TH1D frame2("frame2","frame2", 1, xmin, x_zoom); frame2.Draw(); 
        setCanvas(&frame2, "", ymin, ymax, what);
        if (toys) toys->Draw("P");
        if (toysExp) toysExp->Draw("P");
        obs->Draw(toys ? "LX" : "LXP");    
        if (exp) exp->Draw("LX");
        if (miss) miss->Draw("P");
        finalize(who+(drawToys?"_comp_zoom":"_zoom"),xmin,x_zoom,ymin,ymax,myspam,true);
        if (doZoom2) {
            TH1D frame3("frame3","frame3", 1, x_zoom2_min,x_zoom2_max); frame3.Draw(); 
            setCanvas(&frame3, "", ymin, ymax, "Local p-value");
            obs->Draw("LXP");
            if (exp) exp->Draw("L");
            if (miss) miss->Draw("P");
            finalize(who+"_zoom2",x_zoom2_min,x_zoom2_max,ymin,ymax,myspam,true);
        }
    }
    leg = 0;
}

void drawOnePValBand(TString who, TString who2="") {
    TGraphAsymmErrors *obs = (TGraphAsymmErrors *) gFile->Get(who+"_obs");
    TGraphAsymmErrors *m68 = (TGraphAsymmErrors *) gFile->Get(who+"_median");
    TGraphAsymmErrors *m95 = (TGraphAsymmErrors *) gFile->Get(who+"_median_95"); 
    TGraphAsymmErrors *obs2 = (TGraphAsymmErrors *) gFile->Get(who2+"_obs");
    if (obs == 0 || m68 == 0 || m95 == 0) return;

    TString label1="", label2="", postfix="_band", style2="";
    if (obs2 && who.Index("pvala") == 0 && who2.Index("pval_") == 0) {
        label1="Asympt. "; label2 = "Ensemble "; postfix="_band_as"; style2="P";
        obs2->SetLineWidth(4);
        obs2->SetMarkerStyle(20);
        obs2->SetMarkerSize(1.4);
        obs2->SetLineColor(62);
        obs2->SetMarkerColor(62);
    } else if (obs2 && who.Index("pval_") == 0 && who2.Index("pvala_") == 0) {
        label1="Ensemble "; label2 = "Asympt. "; postfix="_band_t"; style2="L";
        obs2->SetLineColor(2);
        obs2->SetLineStyle(1);
        obs2->SetLineWidth(3);
    }

    double xmin = obs->GetX()[0]             - obs->GetErrorXlow(0), xmin0 = xmin;
    double xmax = obs->GetX()[obs->GetN()-1] + obs->GetErrorXhigh(obs->GetN()-1);
    if (xmin <= 120 && xmax > 200) { xmin = 99.98; xmax = (SM == "FP" ? 300.1 : 600.1); }
    int col = 1; //colorFromName(who);
    obs->SetLineWidth(4);
    obs->SetLineColor(1); obs->SetMarkerColor(1); 
    obs->SetMarkerStyle(21); obs->SetMarkerSize(0.8);
    m68->SetFillColor(color68); m68->SetLineColor(color50); m68->SetLineStyle(2); m68->SetLineWidth(2);
    m95->SetFillColor(color95); m95->SetLineColor(color50); m95->SetLineStyle(2); m95->SetLineWidth(2);

    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    double ymin, ymax; minMaxY(obs, ymin, ymax, 999, 0); 
    ymax = 1.0; ymin = minPValue(obs,xmin,xmax)/10; if (ymin <= 0 || ymin > 1e-5) ymin = 1e-5;
    
    setCanvas(&frame0, "", ymin, ymax, "Local p-value");
    frame0.GetYaxis()->SetTitleOffset(1.10+0.25*isSquareCanvas);
    draw2(who, color68, color95, color50, true, false, -5);
    frame0.Draw("AXIGSAME");
    if (style2 == "L") {
        obs->Draw("LXP");    
        if (obs2) obs2->Draw("LX");
    } else {
        if (obs2) obs2->Draw("PX");
        obs->Draw("LX");    
    }
    leg = newLegend(.66-0.065*isSquareCanvas,.18,.93,.38); leg->SetTextSize(0.04-0.003*isSquareCanvas);
    leg->AddEntry(obs, obs2 ? (label1+"Obs.") : "Observed", (style2 == "L" ? "LP" : "L"));
    leg->AddEntry(m68, "Expected "+oneSigmaText,   "LF");
    leg->AddEntry(m95, "Expected "+twoSigmaText,   "LF");
    if (obs2) leg->AddEntry(obs2, label2+"Obs.",  style2);
    
    leg->Draw();

    TString myspam =  "#splitline{"+SPAM+"}{"+channelFromName(who)+"}";
    if (isSquareCanvas && SPAM.Contains("Preliminary")) myspam = SPAM2L + "\n" + channelFromName(who);
    finalize(who+postfix,xmin,xmax,ymin,ymax,myspam,true);
    if (xmax > x_zoom && xmin < x_zoom) {
        c1->SetLogx(1);
        frame0.GetXaxis()->SetMoreLogLabels(); frame0.GetXaxis()->SetNoExponent();
        finalize(who+postfix+"_logx",xmin,xmax,ymin,ymax,myspam,true);
        c1->SetLogx(0);
        xmin = xmin0;
        TH1D frame2("frame2","frame2", 1, xmin, x_zoom); frame2.Draw(); 
        setCanvas(&frame2, "", ymin, ymax, "Local p-value");
        draw2(who, color68, color95, color50, true, false, -5);
        frame0.Draw("AXIGSAME");
        if (style2 == "L") {
            obs->Draw("LXP");    
            if (obs2) obs2->Draw("LX");
        } else {
            if (obs2) obs2->Draw("PX");
            obs->Draw("LX");    
        }
        leg->Draw();
        finalize(who+postfix+"_zoom",xmin,xmax,ymin,ymax,myspam,true);
        }
    leg = 0;
}

void drawOnePValFitCollage(TString who, TString whichFit, double xmin, double xmax, bool logx, TString postfix) {
    TGraphAsymmErrors *obs = (TGraphAsymmErrors *) gFile->Get("pvala_"+who+"_obs");
    TGraphAsymmErrors *exp = (TGraphAsymmErrors *) gFile->Get("pvala_"+who+"_median");
    TGraphAsymmErrors *mlf = (TGraphAsymmErrors *) gFile->Get(whichFit+"_"+who+"_obs");
    if (obs == 0 || mlf == 0) return;

    int col = 1; 
    mlf->SetLineWidth(2);
    mlf->SetLineColor(col); mlf->SetMarkerColor(col); 
    mlf->SetMarkerStyle(21); mlf->SetMarkerSize(0.8);

    bool zero = whichFit.Contains("mlz");
    isTiny = true;
    TPad *c1_1 = new TPad("pad1", "The pad 80% of the height",0.0,0.45,1.0,1.00,-1); c1_1->Draw();
    TPad *c1_2 = new TPad("pad2", "The pad 20% of the height",0.0,0.00,1.0,0.45,-1); c1_2->Draw();
    c1_1->SetLogx(logx); c1_2->SetLogx(logx); 
    c1_1->SetLogy(1); c1_2->SetLogy(!zero);
    c1_1->SetBottomMargin(0);
    c1_2->SetTopMargin(0.0);
    c1_2->SetBottomMargin(0.25);
    c1_1->SetLeftMargin(0.13);
    c1_2->SetLeftMargin(0.13);
    c1_1->SetRightMargin(0.04);
    c1_2->SetRightMargin(0.04);
    c1_1->cd();

    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    double ymin, ymax; minMaxY(obs, ymin, ymax, 999, 0); ymax = 1.0; ymin = 2e-5;
    setCanvas(&frame0, "", ymin, ymax, "Local p-value");
    frame0.GetYaxis()->SetLabelSize(isSquareCanvas ? 0.065 : 0.07);
    frame0.GetXaxis()->SetLabelSize(isSquareCanvas ? 0.065 : 0.07);
    frame0.GetYaxis()->SetTitleSize(isSquareCanvas ? 0.075 : 0.08);
    frame0.GetXaxis()->SetTitleSize(isSquareCanvas ? 0.075 : 0.08);
    frame0.GetYaxis()->SetTitleOffset(isSquareCanvas ? 0.68 : 0.6);
    frame0.GetXaxis()->SetTitle("");
    //if (exp) exp->Draw("L");
    obs->Draw("LXP");
    TString myspam = SPAM+", "+channelFromName(who);
    if (isSquareCanvas) myspam = "#splitline{"+SPAM+"}{"+channelFromName(who)+"}";
    finalizeNoSave("pvala_"+who,xmin,xmax,ymin,ymax,myspam,true);

    c1_2->cd();
    TH1D frameF("frame","frame", 1, xmin, xmax); frameF.Draw(); gStyle->SetOptStat(0);
    double yminF = 0.02, ymaxF = 20; 
    if (zero) { 
        minMaxY(mlf, yminF, ymaxF, 999, -1);
        ymaxF /= 2.5; yminF = -1;
    }
    setCanvas(&frameF, "", yminF, ymaxF, "Best fit #sigma/#sigma_{"+SM+"} ");
    frameF.GetYaxis()->SetTitleSize(isSquareCanvas ? 0.09 : 0.10);
    frameF.GetYaxis()->SetLabelSize(isSquareCanvas ?  0.07 : 0.08);
    frameF.GetYaxis()->SetTitleOffset(isSquareCanvas ? 0.44 : 0.4);
    frameF.GetYaxis()->SetNoExponent(1);
    frameF.GetYaxis()->SetNdivisions(505);
    frameF.GetXaxis()->SetTitleSize(0.10);
    frameF.GetXaxis()->SetLabelSize(0.08);
    frameF.GetXaxis()->SetMoreLogLabels(); frameF.GetXaxis()->SetNoExponent();
    mlf->SetFillColor(colorFit68);
    mlf->Draw("E3");    
    frameF.Draw("AXIGSAME");
    mlf->Draw("LXP");    
    leg = newLegend(.74-0.05*isSquareCanvas,.27+.50*zero,.92+0.00*isSquareCanvas,.40+.50*zero); leg->SetTextSize(0.08);
    leg->AddEntry(mlf, oneSigmaFitText+" from fit", "F");
    leg->Draw();
    finalizeNoSave(whichFit+"_"+who,xmin,xmax, yminF, ymaxF);

    c1->cd();
    justSave("pval_"+whichFit+"_"+who+postfix,xmin,xmax,ymin,ymax);
    leg = 0;

    if (!zero) {
        c1_2->cd();
        c1_2->SetLogy(0);
        c1_2->Clear(); 
        frameF.Draw();
        minMaxY(mlf, yminF, ymaxF, 999, -1); ymaxF /= 3;
        setCanvas(&frameF, "", 0, ymaxF, "Best fit #sigma/#sigma_{"+SM+"} ");
        frameF.GetYaxis()->SetTitle("Best fit #sigma/#sigma_{"+SM+"} ");
        frameF.GetYaxis()->SetTitleOffset(0.3+0.05*isSquareCanvas);
        mlf->Draw("E3");    
        frameF.Draw("AXIGSAME");
        mlf->Draw("LXP");    
        leg = newLegend(.74-0.05*isSquareCanvas,.77,.92+0.00*isSquareCanvas,.90); leg->SetTextSize(0.08);
        leg->AddEntry(mlf, oneSigmaFitText+" from fit", "F");
        leg->Draw();
        finalizeNoSave(whichFit+"_"+who,xmin,xmax, 0.0, ymaxF);

        c1->cd();
        justSave("pval_"+whichFit+"_"+who+postfix+"_liny",xmin,xmax,ymin,ymax);
    }

    c1_2->Delete();
    c1_1->Delete();
    c1->Clear();
    isTiny = false;
    leg = 0;
}

void drawOnePValFitCollage(TString who, TString whichFit="ml") {
    TGraphAsymmErrors *obs = (TGraphAsymmErrors *) gFile->Get("pvala_"+who+"_obs");
    TGraphAsymmErrors *mlf = (TGraphAsymmErrors *) gFile->Get(whichFit+"_"+who+"_obs");
    if (obs == 0) return;
    double xmin = obs->GetX()[0]             - obs->GetErrorXlow(0), xmin0 = xmin;
    double xmax = obs->GetX()[obs->GetN()-1] + obs->GetErrorXhigh(obs->GetN()-1);
    if (xmin <= 120 && xmax > 200) { 
        drawOnePValFitCollage(who, whichFit, 99.98, 600.1, false, "");
        drawOnePValFitCollage(who, whichFit, 99.98, 600.1, true,  "_logx");
        drawOnePValFitCollage(who, whichFit, xmin,  x_zoom, false, "_zoom");
    } else {
        drawOnePValFitCollage(who, whichFit, xmin, xmax, false, "");
        if (xmax > x_zoom) drawOnePValFitCollage(who, whichFit, xmin, x_zoom, false, "_zoom");
    }
}


void drawOnePValFitCollageOld(TString who, TString what="auto") {
    TGraphAsymmErrors *obs = (TGraphAsymmErrors *) gFile->Get("pvala_"+who+"_obs");
    TGraphAsymmErrors *mlf = (TGraphAsymmErrors *) gFile->Get("ml_"+who+"_obs");
    if (obs == 0 || mlf == 0) return;

    double xmin = obs->GetX()[0]             - obs->GetErrorXlow(0), xmin0 = xmin;
    double xmax = obs->GetX()[obs->GetN()-1] + obs->GetErrorXhigh(obs->GetN()-1);
    if (xmin <= 120 && xmax > 200) { xmin = 99.98; xmax = 600.1; }

    int col = 1; 
    obs->SetLineWidth(2);
    obs->SetLineColor(col); obs->SetMarkerColor(col); 
    obs->SetMarkerStyle(21); obs->SetMarkerSize(0.8);
    mlf->SetLineWidth(2);
    mlf->SetLineColor(col); mlf->SetMarkerColor(col); 
    mlf->SetMarkerStyle(21); mlf->SetMarkerSize(0.8);

    isTiny = true;
    TPad *c1_1 = new TPad("pad1", "The pad 80% of the height",0.0,0.65,1.0,1.00,-1); c1_1->Draw();
    TPad *c1_2 = new TPad("pad2", "The pad 20% of the height",0.0,0.00,1.0,0.65,-1); c1_2->Draw();
    if (xmax > 500 && xmin < 200) { c1_1->SetLogx(1); c1_2->SetLogx(1); }
    c1_1->SetLogy(1); c1_2->SetLogy(1);
    c1_1->SetBottomMargin(0);
    c1_2->SetTopMargin(0.02);
    c1_1->SetLeftMargin(0.10);
    c1_2->SetLeftMargin(0.10);
    c1_1->SetRightMargin(0.04);
    c1_2->SetRightMargin(0.04);

    c1_1->cd();
    TH1D frameF("frame","frame", 1, xmin, xmax); frameF.Draw(); gStyle->SetOptStat(0);
    setCanvas(&frameF, "", 0.02, 20, "Best fit #sigma/#sigma_{"+SM+"}");
    frameF.GetYaxis()->SetTitleSize(0.12);
    frameF.GetYaxis()->SetLabelSize(0.10);
    frameF.GetYaxis()->SetTitleOffset(isSquareCanvas ? 0.42 : 0.4);
    frameF.GetXaxis()->SetTitle("");
    frameF.GetYaxis()->SetNoExponent(1);
    mlf->SetFillColor(colorFit68);
    mlf->Draw("E3");    
    frameF.Draw("AXIGSAME");
    mlf->Draw("LXP");    
    leg = newLegend(.74-0.05*isSquareCanvas,.75,.90+0.04*isSquareCanvas,.90); leg->SetTextSize(0.10);
    leg->AddEntry(mlf, oneSigmaFitText+" from fit", "F");
    leg->Draw();
    leg = 0;
    finalizeNoSave("ml_"+who,xmin,xmax,0.02, 20);

    c1_2->cd();
    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    double ymin, ymax; minMaxY(obs, ymin, ymax, 999, 0); ymax = 1.0; if (ymin > 2e-6) ymin = 2e-6;
    setCanvas(&frame0, "", ymin, ymax, "Local p-value");
    frame0.GetYaxis()->SetLabelSize(isSquareCanvas ? 0.055 : 0.06);
    frame0.GetXaxis()->SetLabelSize(isSquareCanvas ? 0.055 : 0.06);
    frame0.GetYaxis()->SetTitleSize(isSquareCanvas ? 0.065 : 0.07);
    frame0.GetXaxis()->SetTitleSize(isSquareCanvas ? 0.065 : 0.07);
    frame0.GetYaxis()->SetTitleOffset(isSquareCanvas ? 0.78 : 0.7);
    frame0.GetXaxis()->SetMoreLogLabels(); frame0.GetXaxis()->SetNoExponent();
    obs->Draw("LP");    
    TString myspam = "#splitline{"+SPAM+"}{"+channelFromName(who)+"}";
    finalizeNoSave("pvala_"+who,xmin,xmax,ymin,ymax,myspam,true);
    /*if (xmax >= 200 && xmin < 200) {
        c1->SetLogx(1);
        frame0.GetXaxis()->SetMoreLogLabels(); frame0.GetXaxis()->SetNoExponent();
        finalize(who+"_logx",xmin,xmax,ymin,ymax,myspam,true);
        c1->SetLogx(0);
        xmin = xmin0;
        TH1D frame2("frame2","frame2", 1, xmin, x_zoom); frame2.Draw(); 
        setCanvas(&frame2, "", ymin, ymax, "Local p-value");
        obs->Draw("LXP");
        finalize(who+"_zoom",xmin,x_zoom,ymin,ymax,myspam,true);
    }*/

    c1->cd();
    justSave("pval_fit_"+who,xmin,xmax,ymin,ymax);
    c1_2->Delete();
    c1_1->Delete();
    c1->Clear();
    isTiny = false;
    leg = 0;
}

void drawCompLimit(TString who1, TString who2, TString label="One", TString oldlabel="Two", TString what="auto") {
    TGraphAsymmErrors *obs = (TGraphAsymmErrors *) gFile->Get(who1);
    TGraphAsymmErrors *old = (TGraphAsymmErrors *) gFile->Get(who2);
    if (obs == 0 || old == 0) return;
    double xmin = obs->GetX()[0]             - obs->GetErrorXlow(0), xmin0 = xmin;
    double xmax = obs->GetX()[obs->GetN()-1] + obs->GetErrorXhigh(obs->GetN()-1);
    if (xmin <= 120 && xmax > 200) { xmin = 99.98; xmax = (SM == "FP" ? 300.1 : 600.1); }
    int col = 2, dcol = 4;
    obs->SetLineWidth(2);
    obs->SetLineColor(col); obs->SetMarkerColor(col); 
    obs->SetMarkerStyle(21); obs->SetMarkerSize(0.8);
    old->SetLineWidth(2);
    old->SetLineColor(dcol); old->SetMarkerColor(dcol); 
    old->SetMarkerStyle(25); old->SetMarkerSize(0.8);
    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    double ymin, ymax; minMaxY(obs, ymin, ymax);
    double ymin2, ymax2; minMaxY(old, ymin2, ymax2);
    if (ymin2 < ymin) ymin = ymin2; if (ymax2 > ymax) ymax = ymax2;
    if (who1.Contains("mlz_")) { ymin = -2.5; ymax = 5; }
    if ((who1.Contains("median") && who2.Contains("median")) ||
        (who1.Contains("mlz") && who2.Contains("mlz"))) {
        old->SetFillStyle(1001); old->SetFillColor(216);
        obs->SetFillStyle(3444); obs->SetFillColor(2); obs->SetLineColor(216);
        old->Draw("E3"); obs->Draw("E3");    
        old->Draw("LPX"); obs->Draw("LPX");    
    } else {
        old->Draw("LP"); obs->Draw("LP");    
    }
    leg = newLegend(.69,.75,.84,.92); leg->SetTextSize(0.05);
    leg->AddEntry(obs,    label, "LP");
    leg->AddEntry(old, oldlabel, "LP");
    leg->Draw();
    setCanvas(&frame0, "", ymin, ymax, "Limit");
    TString myspam("#splitline{PRIVATE FOR COMB.}{"+channelFromName(who1)+"}");
    finalize(who1+"_vs_"+who2,xmin,xmax,ymin,ymax,myspam);
    if (xmax >= 200 && xmin < x_zoom) {
        xmin = xmin0;
        TH1D frame2("frame2","frame2", 1, xmin, x_zoom); frame2.Draw(); 
        if ((who1.Contains("median") && who2.Contains("median")) ||
            (who1.Contains("mlz") && who2.Contains("mlz"))) {
            old->Draw("E3"); obs->Draw("E3");    
            old->Draw("LPX"); obs->Draw("LPX");    
        } else {
            old->Draw("LP"); obs->Draw("LP");    
        }
        leg->Draw();
        minMaxY(obs, ymin, ymax, x_zoom); minMaxY(old, ymin2, ymax2, x_zoom);
        if (ymin2 < ymin) ymin = ymin2; if (ymax2 > ymax) ymax = ymax2;
        if (who1.Contains("mlz_")) { ymin = -2.5; ymax = 5; }
        setCanvas(&frame2, "", ymin, ymax, what);
        finalize(who1+"_vs_"+who2+"_zoom",xmin,x_zoom,ymin,ymax,myspam);
        if (doZoom2) {
            TH1D frame3("frame3","frame3", 1,x_zoom2_min, x_zoom2_max); frame3.Draw(); 
            if ((who1.Contains("median") && who2.Contains("median")) ||
                (who1.Contains("mlz") && who2.Contains("mlz"))) {
                old->Draw("E3"); obs->Draw("E3");    
                old->Draw("LPX"); obs->Draw("LPX");    
            } else {
                old->Draw("LP"); obs->Draw("LP");    
            }
            leg->Draw();
            setCanvas(&frame3, "", ymin, ymax);
            finalize(who1+"_vs_"+who2+"_zoom2",x_zoom2_min,x_zoom2_max,ymin,ymax,myspam);
        }
    }
    leg = 0;
}
void drawCompPVal(TString who1, TString who2, TString label="One", TString oldlabel="Two", TString what="auto", double ymin=1e-4) {
    TGraphAsymmErrors *obs = (TGraphAsymmErrors *) gFile->Get(who1); if (obs == 0) std::cout << "Cannot find " << who1 << std::endl;
    TGraphAsymmErrors *old = (TGraphAsymmErrors *) gFile->Get(who2); if (old == 0) std::cout << "Cannot find " << who2 << std::endl;
    if (obs == 0 || old == 0) return;
    double xmin = obs->GetX()[0]             - obs->GetErrorXlow(0), xmin0 = xmin;
    double xmax = obs->GetX()[obs->GetN()-1] + obs->GetErrorXhigh(obs->GetN()-1);
    if (xmin <= 120 && xmax > 200) { xmin = 99.98; xmax = (SM == "FP" ? 300.1 : 600.1); }
    int col = 2, dcol = 4;
    obs->SetLineWidth(2);
    obs->SetLineColor(col); obs->SetMarkerColor(col); 
    obs->SetMarkerStyle(21); obs->SetMarkerSize(0.8);
    old->SetLineWidth(2);
    old->SetLineColor(dcol); old->SetMarkerColor(dcol); 
    old->SetMarkerStyle(25); old->SetMarkerSize(0.8);
    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    double ymax = 1.0;
    old->Draw("LP");    
    obs->Draw("LP");    
    leg = newLegend(.75,.14,.92,.28); leg->SetTextSize(0.05);
    leg->AddEntry(obs,    label, "LP");
    leg->AddEntry(old, oldlabel, "LP");
    leg->Draw();
    setCanvas(&frame0, "", ymin, ymax, "Local p-value");
    TString myspam("#splitline{PRIVATE FOR COMB.}{"+channelFromName(who1)+"}");
    finalize(who1+"_vs_"+who2,xmin,xmax,ymin,ymax,myspam,true);
    if (xmax >= 200 && xmin < x_zoom) {
        c1->SetLogx(1);
        frame0.GetXaxis()->SetMoreLogLabels(); frame0.GetXaxis()->SetNoExponent();
        finalize(who1+"_vs_"+who2+"_logx",xmin,xmax,ymin,ymax,myspam,true);
        c1->SetLogx(0);
        xmin = xmin0;
        TH1D frame2("frame2","frame2", 1, xmin, x_zoom); frame2.Draw(); 
        old->Draw("LP");    
        obs->Draw("LP");  
        leg->Draw();
        setCanvas(&frame2, "", ymin, ymax, "Local p-value");
        finalize(who1+"_vs_"+who2+"_zoom",xmin,200,ymin,ymax,myspam,true);
        if (doZoom2) {
            TH1D frame3("frame3","frame3", 1,x_zoom2_min, x_zoom2_max); frame3.Draw(); 
            old->Draw("LP");    
            obs->Draw("LP");  
            leg->Draw();
            setCanvas(&frame3, "", ymin, ymax, "Local p-value");
            finalize(who1+"_vs_"+who2+"_zoom2",x_zoom2_min,x_zoom2_max,ymin,ymax,myspam,true);
        }
    }
    leg = 0;
}

TGraphAsymmErrors *makeGrid(TString who, TString who2, double lowfactor=0.2, double highfactor=4) {
    TGraphAsymmErrors *g1 = (TGraphAsymmErrors *) gFile->Get(who),  *m1 = g1 ? missingPoints(g1) : 0; 
    TGraphAsymmErrors *g2 = (TGraphAsymmErrors *) gFile->Get(who2), *m2 = g2 ? missingPoints(g2) : 0;
    if (g1 == 0 && g2 == 0) { std::cerr << "Missing both '" << who << "' and '" << who2 << "'" << std::endl; return 0; }
    TGraphAsymmErrors *ret = new TGraphAsymmErrors(nmasses); 
    for (int i = 0; i < nmasses; ++i) {
        double x = masses[i];
        int i1 = findBin(g1, x), j1 = findBin(m1, x);
        int i2 = findBin(g2, x), j2 = findBin(m2, x);
        double ymin = 999, ymax = -999;
        if (i1 != -1 && g1->GetY()[i1] != 0) { ymin = TMath::Min(ymin, g1->GetY()[i1]); ymax = TMath::Max(ymax, g1->GetY()[i1]); }
        if (i2 != -1 && g2->GetY()[i2] != 0) { ymin = TMath::Min(ymin, g2->GetY()[i2]); ymax = TMath::Max(ymax, g2->GetY()[i2]); }
        if (j1 != -1 && m1->GetY()[j1] != 0) { ymin = TMath::Min(ymin, m1->GetY()[j1]); ymax = TMath::Max(ymax, m1->GetY()[j1]); }
        if (j2 != -1 && m2->GetY()[j2] != 0) { ymin = TMath::Min(ymin, m2->GetY()[j2]); ymax = TMath::Max(ymax, m2->GetY()[j2]); }
        if (SM != "SM4") {
            if (ymin < 0.2) ymin = 0.2; if (ymin > 3) ymin = 3;
            if (who.Contains("comb")) {
                if (ymax < 0.2) ymax = 0.2; 
                if (ymax > 8)   ymax = 8;
            } else {
                if (ymax < 0.5) ymax = 0.5; 
            }
        } else {
            if (ymin < 0.02) ymin = 0.02; if (ymin > 3) ymin = 3;
            if (who.Contains("comb")) {
                if (ymax < 0.02) ymax = 0.02; 
                if (ymax > 8)   ymax = 8;
            } else {
                if (ymax < 0.05) ymax = 0.05; 
            }
        }
        double y = sqrt(ymin*ymax);
        ret->SetPoint(i, x, y);
        ret->SetPointError(i, 0, 0, y - lowfactor*ymin, highfactor*ymax - y);
    }
    //ret->SetLineStyle(1);
    //ret->SetMarkerStyle(0);
    ret->SetFillColor(16);
    ret->SetFillStyle(3244);
    ret->SetLineColor(223);
    ret->SetLineWidth(4);
    ret->SetLineStyle(3);
    return ret;
}


bool CLs_debug_apriori_grid=true;
TList CLs_grids;
TGraphAsymmErrors *makeAPrioriGrid(TString who) {
    std::cout << "Asked for grid " << who << std::endl; //CLs_grids.ls();
    if (CLs_grids.FindObject("grid_"+who)) return (TGraphAsymmErrors *) CLs_grids.FindObject("grid_"+who);
    return makeGrid(who+"_obs", who+"_median", 0.4, 2.5);
}


bool Draw_TEV=false;
void drawOneCLs(TString who) {
    TGraphAsymmErrors *obs   = (TGraphAsymmErrors *) gFile->Get(who+"_obs");
    if (obs == 0) return;
    TGraphAsymmErrors *apriori = CLs_debug_apriori_grid ? makeAPrioriGrid(who) : 0;
    TGraphAsymmErrors *miss = missingPoints(obs);
    TGraphAsymmErrors *m68 = (TGraphAsymmErrors *) gFile->Get(who+"_median");
    TGraphAsymmErrors *m95 = (TGraphAsymmErrors *) gFile->Get(who+"_median_95"); 
    TGraphAsymmErrors *miss68 = m68 ? missingPoints(m68) : 0;
    if (miss68) { miss68->SetMarkerStyle(24); miss68->SetMarkerSize(0.9); miss68->SetMarkerColor(4); }

    TGraphAsymmErrors *obsTEV = Draw_TEV ? (TGraphAsymmErrors *) gFile->Get("tevatron_obs")    : 0;
    TGraphAsymmErrors *expTEV = Draw_TEV ? (TGraphAsymmErrors *) gFile->Get("tevatron_median") : 0;
    if (Draw_TEV) {
        obsTEV->SetLineColor(kBlue); obsTEV->SetLineWidth(5); obsTEV->SetLineStyle(1);
        expTEV->SetLineColor(kBlue); expTEV->SetLineWidth(5); expTEV->SetLineStyle(2);
    }

    double xmin = obs->GetX()[0]             - obs->GetErrorXlow(0), xmin0 = xmin;
    double xmax = obs->GetX()[obs->GetN()-1] + obs->GetErrorXhigh(obs->GetN()-1);
    if (xmin <= 120 && xmax > 200) { xmin = 99.98; xmax = (SM == "FP" ? 300.1 : 600.1); }
    int col = 1, dcol = 1, smooth = 5, smoothorder = 0;
    obs->SetLineWidth(2); 
    obs->SetLineColor(col); obs->SetMarkerColor(col); 
    obs->SetMarkerStyle(21); obs->SetMarkerSize(0.8);
    double ymin, ymax;   minMaxY(obs, ymin, ymax);
    double ymin2, ymax2; if (m68) minMaxY(m68, ymin2, ymax2); 
    if (ymax2 > ymax) ymax = ymax2; if (ymin2 < ymin) ymin = ymin2;
    if (Draw_TEV) ymax *=2;
    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    frame0.GetYaxis()->SetTitleOffset(1.10+0.25*isSquareCanvas);
    TString what = "95% CL limit on #sigma/#sigma_{"+SM+"}";
    if (who.Contains("cls99_")) TString what = "99% CL limit on #sigma/#sigma_{"+SM+"}";
    if (who.Contains("cls90_")) TString what = "90% CL limit on #sigma/#sigma_{"+SM+"}";
    //if (who.Contains("acls_")) what = "Asymptotic "+what;
    setCanvas(&frame0, "", ymin, ymax, what);
    if (apriori) apriori->Draw("E3");
    draw2(who, color68, color95, color50, true, false, smooth, smoothorder);
    if (m68) { m68->SetFillColor(color68); m68->SetLineColor(color50); m68->SetLineStyle(2); m68->SetLineWidth(2); }
    if (m95) { m95->SetFillColor(color95); m95->SetLineColor(color50); m95->SetLineStyle(2); m95->SetLineWidth(2); }
    frame0.Draw("AXIGSAME");
    double leg_y_hi = 0.94;
    double leg_y_lo = (isSquareCanvas ? 0.78 : .75)- 0.05*(lep_excluded+tev_excluded+1.6*Draw_TEV);
    double leg_x_lo = isSquareCanvas ? .645-0.06*(Draw_TEV||tev_excluded) : .65;
    leg = newLegend(leg_x_lo, leg_y_lo,.93+0.005*isSquareCanvas, leg_y_hi); leg->SetTextSize(isSquareCanvas || lep_excluded?  0.034 : 0.037);
    obs->Draw("LP");    
    if (miss68) miss68->Draw("P");
    if (miss) miss->Draw("P");
    if (apriori) apriori->Draw("LX");
    leg->AddEntry(obs, "Observed", "LP");
    leg->AddEntry(m68, "Expected "+oneSigmaText, "LF");
    leg->AddEntry(m95, "Expected "+twoSigmaText, "LF");
    if (Draw_TEV) {
        leg->AddEntry(obsTEV, "Tevatron Observed", "L");
        leg->AddEntry(expTEV, "Tevatron Expected", "L");
    }
    if (lep_excluded) leg->AddEntry(fakeLEP, "LEP excluded", "F");
    if (tev_excluded) leg->AddEntry(fakeTEV, "Tevatron excluded", "F");
    if (cms_excluded) leg->AddEntry(fakeCMS, "CMS excluded", "F");
    if (Draw_TEV) { expTEV->Draw("L SAME"); obsTEV->Draw("L SAME"); }
    leg->Draw();
    TString myspam = "#splitline{"+SPAM+"}{"+channelFromName(who)+"}";
    if (isSquareCanvas && SPAM.Contains("Preliminary")) myspam = SPAM2L + "\n" + channelFromName(who);
    if (!c1->GetGridx() && !c1->GetGridy()) frame0.Draw("AXIGSAME");
    if (channelSpamOnRightHandSide) {
        spam(channelFromName(who,false,true)+" only", .65, leg_y_lo-0.01, 0.93, leg_y_lo-0.06, 22);
    }
    finalize(who,xmin,xmax,ymin,ymax,myspam);
    if (xmax > x_zoom && xmin < x_zoom) {
        c1->SetLogx(1);
        frame0.GetXaxis()->SetMoreLogLabels(); frame0.GetXaxis()->SetNoExponent();
        finalize(who+"_logx",xmin,xmax,ymin,ymax,myspam);
        c1->SetLogx(0);
        xmin = xmin0;
        TH1D frame2("frame2","frame2", 1, xmin, x_zoom); frame2.Draw(); 
        setCanvas(&frame2, "", ymin, ymax, what);
        if (apriori) apriori->Draw("E3");
        draw2(who, color68, color95, color50, true, false, smooth, smoothorder);
        frame2.Draw("AXIGSAME");
        obs->Draw("LP");    
        bool noMoreTev = tev_excluded && (x_zoom <= 156);
        leg = newLegend(leg_x_lo+0.06*noMoreTev*isSquareCanvas, leg_y_lo+0.05*noMoreTev,.93+0.005*isSquareCanvas, leg_y_hi); 
        leg->SetTextSize(isSquareCanvas || lep_excluded?  0.034 : 0.037);
        leg->AddEntry(obs, "Observed", "LP");
        leg->AddEntry(m68, "Expected "+oneSigmaText, "LF");
        leg->AddEntry(m95, "Expected "+twoSigmaText, "LF");
        if (Draw_TEV) {
            leg->AddEntry(obsTEV, "Tevatron Observed", "L");
            leg->AddEntry(expTEV, "Tevatron Expected", "L");
        }
        if (lep_excluded) leg->AddEntry(fakeLEP, "LEP excluded", "F");
        if (tev_excluded && !noMoreTev) leg->AddEntry(fakeTEV, "Tevatron excluded", "F");
        if (cms_excluded) leg->AddEntry(fakeCMS, "CMS excluded", "F");

        if (miss68) miss68->Draw("P");
        if (miss) miss->Draw("P");
        if (Draw_TEV) { expTEV->Draw("L SAME"); obsTEV->Draw("L SAME"); }
        if (apriori) apriori->Draw("LX");
        leg->Draw();
        minMaxY(obs, ymin, ymax, x_zoom);
        if (m68) minMaxY(m68, ymin2, ymax2, x_zoom); 
        if (ymax2 > ymax) ymax = ymax2; if (ymin2 < ymin) ymin = ymin2;
        if (Draw_TEV) ymax *=2;
        if (!c1->GetGridx() && !c1->GetGridy()) frame2.Draw("AXIGSAME");
        if (channelSpamOnRightHandSide) {
            spam(channelFromName(who,false,true)+" only", .65, leg_y_lo-0.01, 0.93, leg_y_lo-0.06, 22);
        }
        finalize(who+"_zoom",xmin,x_zoom,ymin,ymax,myspam);
        if (doZoom2) {
            TH1D frame3("frame3","frame3", 1, x_zoom2_min,x_zoom2_max); frame3.Draw(); 
            setCanvas(&frame3, "", ymin, ymax, what);
            if (apriori) apriori->Draw("E3");
            draw2(who, color68, color95, color50, true, false, 5);
            frame3.Draw("AXIGSAME");
            obs->Draw("LP");    
            if (miss68) miss68->Draw("P");
            if (miss) miss->Draw("P");
            if (apriori) apriori->Draw("LX");
            leg->Draw();
            finalize(who+"_zoom2",x_zoom2_min,x_zoom2_max,ymin,ymax,myspam);
        }
    }
    leg = 0;
}
void drawTwoCLs(TString who, TString who2, TString name2="Other", TString postfix="", int color2=kBlue, bool doObs2=true) {
    TGraphAsymmErrors *obs   = (TGraphAsymmErrors *) gFile->Get(who+"_obs");
    if (obs == 0) return;
    TGraphAsymmErrors *m68 = (TGraphAsymmErrors *) gFile->Get(who+"_median");
    TGraphAsymmErrors *m95 = (TGraphAsymmErrors *) gFile->Get(who+"_median_95"); 

    TGraphAsymmErrors *obs2 = (TGraphAsymmErrors *) gFile->Get(who2+"_obs");
    TGraphAsymmErrors *exp2 = (TGraphAsymmErrors *) gFile->Get(who2+"_median");
    if (obs2 == 0 || exp2 == 0) { std::cerr << "Missing " << who2 << " to compare with " << who << std::endl; return; }
    obs2->SetLineColor(color2); obs2->SetLineWidth(5); obs2->SetLineStyle(1);
    exp2->SetLineColor(color2); exp2->SetLineWidth(5); exp2->SetLineStyle(2);

    double xmin = obs->GetX()[0]             - obs->GetErrorXlow(0), xmin0 = xmin;
    double xmax = obs->GetX()[obs->GetN()-1] + obs->GetErrorXhigh(obs->GetN()-1);
    if (xmin <= 120 && xmax > 200) { xmin = 99.98; xmax = (SM == "FP" ? 300.1 : 600.1); }
    int col = 1, dcol = 1;
    obs->SetLineWidth(2); 
    obs->SetLineColor(col); obs->SetMarkerColor(col); 
    obs->SetMarkerStyle(21); obs->SetMarkerSize(0.8);
    double ymin, ymax;   minMaxY(obs, ymin, ymax); 
    double ymin2, ymax2; minMaxY(obs2, ymin2, ymax2); 
    if (ymax2 > ymax) ymax = ymax2; if (ymin2 < ymin) ymin = ymin2;
    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    frame0.GetYaxis()->SetTitleOffset(1.10+0.25*isSquareCanvas);
    TString what = "95% CL limit on #sigma/#sigma_{"+SM+"}";
    //if (who.Contains("acls_")) what = "Asymptotic "+what;
    setCanvas(&frame0, "", ymin, ymax, what);
    draw2(who, color68, color95, color50, true, false, 5);
    if (m68) { m68->SetFillColor(color68); m68->SetLineColor(color50); m68->SetLineStyle(2); m68->SetLineWidth(2); }
    if (m95) { m95->SetFillColor(color95); m95->SetLineColor(color50); m95->SetLineStyle(2); m95->SetLineWidth(2); }
    frame0.Draw("AXIGSAME");
    double leg_y_hi = 0.94;
    double leg_y_lo = (isSquareCanvas ? 0.78 : .75)- 0.05*(1+doObs2);
    double leg_x_lo = isSquareCanvas ? .645-0.06 : .65;
    leg = newLegend(leg_x_lo, leg_y_lo,.93+0.005*isSquareCanvas, leg_y_hi); leg->SetTextSize(isSquareCanvas?  0.034 : 0.037);
    obs->Draw("LP");    
    exp2->Draw("LX SAME"); if (doObs2) obs2->Draw("LX SAME");
    leg->AddEntry(obs, "Observed", "LP");
    leg->AddEntry(m68, "Expected "+oneSigmaText, "LF");
    leg->AddEntry(m95, "Expected "+twoSigmaText, "LF");
    if (doObs2) leg->AddEntry(obs2, name2+" Obs.", "L");
    leg->AddEntry(exp2, name2+" Exp.", "L");
    leg->Draw();
    TString myspam = "#splitline{"+SPAM+"}{"+channelFromName(who)+"}";
    if (isSquareCanvas && SPAM.Contains("Preliminary")) myspam = SPAM2L + "\n" + channelFromName(who);
    finalize(who+postfix,xmin,xmax,ymin,ymax,myspam);
    if (xmax > x_zoom && xmin < x_zoom) {
        c1->SetLogx(1);
        frame0.GetXaxis()->SetMoreLogLabels(); frame0.GetXaxis()->SetNoExponent();
        finalize(who+postfix+"_logx",xmin,xmax,ymin,ymax,myspam);
        c1->SetLogx(0);
        xmin = xmin0;
        TH1D frame2("frame2","frame2", 1, xmin, x_zoom); frame2.Draw(); 
        setCanvas(&frame2, "", ymin, ymax, what);
        draw2(who, color68, color95, color50, true, false, 5);
        frame2.Draw("AXIGSAME");
        obs->Draw("LP");    
        exp2->Draw("LX SAME"); if (doObs2) obs2->Draw("LX SAME");
        leg->Draw();
        minMaxY(obs, ymin, ymax, x_zoom); minMaxY(exp2, ymin, ymax, x_zoom);
        if (ymax2 > ymax) ymax = ymax2; if (ymin2 < ymin) ymin = ymin2;
        finalize(who+postfix+"_zoom",xmin,x_zoom,ymin,ymax,myspam);
    }
    leg = 0;
}

void drawMethods(TString who, bool include_pla=true, bool bands=true) {
    TGraphAsymmErrors *obs = (TGraphAsymmErrors *) gFile->Get(who+"_obs");
    if (obs == 0) return;
    TGraphAsymmErrors *m68 = (TGraphAsymmErrors *) gFile->Get(who+"_median");    m68->SetFillColor(211); m68->SetLineColor(39);
    TGraphAsymmErrors *m95 = (TGraphAsymmErrors *) gFile->Get(who+"_median_95"); m95->SetFillColor(220); m95->SetLineColor(39);
    bool isToy = !who.Contains("acls");

    TString what = "95% CL limit on #sigma/#sigma_{"+SM+"}";

    TString pla = who; if (isToy) { pla.ReplaceAll("cls_","acls_"); } else { pla.ReplaceAll("acls_","pla_"); }
    std::cout << "PLA is " << pla << std::endl;
    TString bayes = who; bayes.ReplaceAll(isToy ? "cls_" : "acls_","bayes_"); 
    TGraphAsymmErrors *pla_obs  = (TGraphAsymmErrors *) gFile->Get(pla+"_obs");
    TGraphAsymmErrors *bayes_obs = (TGraphAsymmErrors *) gFile->Get(bayes+"_obs");
    double xmin = obs->GetX()[0]             - obs->GetErrorXlow(0), xmin0 = xmin;
    double xmax = obs->GetX()[obs->GetN()-1] + obs->GetErrorXhigh(obs->GetN()-1);
    if (xmin <= 120 && xmax > 200) { xmin = 99.98; xmax = (SM == "FP" ? 300.1 : 600.1); }
    obs->SetLineWidth(2); 
    obs->SetLineColor(1); obs->SetMarkerColor(1); 
    obs->SetMarkerStyle(21); obs->SetMarkerSize(0.8);
    double ymin, ymax; minMaxY(obs, ymin, ymax);
    double ymin2, ymax2; minMaxY(m68, ymin2, ymax2); 
    if (ymax2 > ymax) ymax = ymax2; if (ymin2 < ymin) ymin = ymin2;
    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    setCanvas(&frame0, "", ymin, ymax, what);
    leg = newLegend(.65-0.13*isSquareCanvas,.75 - 0.04*bands+0.01*isSquareCanvas,.95-0.02*isSquareCanvas,.94); leg->SetTextSize(0.037-0.003*isSquareCanvas);
    if (bands) {
        m68->SetFillColor(color68); m68->SetLineColor(color50); m68->SetLineStyle(2); m68->SetLineWidth(2);
        m95->SetFillColor(color95); m95->SetLineColor(color50); m95->SetLineStyle(2); m95->SetLineWidth(2);
        draw2(who, color68, color95, color50, true, false, 5);
        frame0.Draw("AXIGSAME");
    } else {
        m68->SetLineColor(222); m68->SetLineStyle(9); m68->SetLineWidth(4); 
        m68->Draw("LX");
    } 
    obs->Draw("LP"); 
    leg->AddEntry(obs, isToy ? "CLs Observed" : "Asym. CLs Obs.", "LP");
    if (bands) {
        leg->AddEntry(m68, isToy ? "CLs Expected "+oneSigmaText : "CLs Exp. "+oneSigmaText, "LF");
        leg->AddEntry(m95, isToy ? "CLs Expected "+twoSigmaText : "CLs Exp. "+twoSigmaText, "LF");
    } else {
        leg->AddEntry(m68, isToy ? "CLs Expected" : "Asym. CLs Exp.", "L");
    }
    if (bayes_obs) {
        bayes_obs->SetLineColor(215);
        bayes_obs->SetMarkerColor(215);
        bayes_obs->SetLineStyle(2);
        bayes_obs->SetMarkerStyle(24);
        bayes_obs->SetMarkerSize(0.7);
        bayes_obs->SetLineWidth(3);
        bayes_obs->Draw("LP");
        leg->AddEntry(bayes_obs, "Bayesian Observed", "LP");
    }
    if (include_pla && pla_obs != 0) { 
        pla_obs->SetLineColor(2);
        pla_obs->SetLineStyle(1);
        pla_obs->SetLineWidth(2);
        pla_obs->Draw("LX");
        leg->AddEntry(pla_obs, isToy ? "Asymptotic CLs Obs." : "PL Approx. Observed", "L");
    }
    leg->Draw();
    TString myspam = "#splitline{"+SPAM+"}{"+channelFromName(who)+"}";
    if (isSquareCanvas && SPAM.Contains("Preliminary")) myspam = SPAM2L + "\n" + channelFromName(who);
    TString xname = (include_pla ? "_comp3" : "_comp");
    finalize(who+xname,xmin,xmax,ymin,ymax,myspam);
    if (xmax > x_zoom && xmin < x_zoom) {
        double xcut = (SM == "SM" ? x_zoom : 140);
        c1->SetLogx(1);
        frame0.GetXaxis()->SetMoreLogLabels(); frame0.GetXaxis()->SetNoExponent();
        finalize(who+xname+"_logx",xmin,xmax,ymin,ymax,myspam);
        c1->SetLogx(0);
        xmin = xmin0;
        TH1D frame2("frame2","frame2", 1, xmin, xcut); frame2.Draw(); 
        setCanvas(&frame2, "", ymin, ymax, what);
        draw2(who, color68, color95, color50, true, false, 5);
        frame2.Draw("AXIGSAME");
        obs->Draw("LP");
        if (bayes_obs) bayes_obs->Draw("LP");
        if (include_pla && pla_obs != 0) pla_obs->Draw("LX");
        minMaxY(obs, ymin, ymax, xcut);
        minMaxY(m68, ymin2, ymax2, xcut); 
        if (ymax2 > ymax) ymax = ymax2; if (ymin2 < ymin) ymin = ymin2;
        leg->Draw();
        finalize(who+xname+"_zoom",xmin,xcut,ymin,ymax,myspam);
    }
    leg = 0;
}


void drawSMCLs(TString who, bool showAs=false) {
    TGraphAsymmErrors *obs   = (TGraphAsymmErrors *) gFile->Get(who+"_obs");
    TGraphAsymmErrors *m68 = (TGraphAsymmErrors *) gFile->Get(who+"_median");    
    TGraphAsymmErrors *m95 = (TGraphAsymmErrors *) gFile->Get(who+"_median_95"); 
    if (obs == 0 || m68 == 0 || m95 == 0) return;
    m68->SetFillColor(211); m68->SetLineColor(39);
    m95->SetFillColor(220); m95->SetLineColor(39);
    //bool savGridY = c1->GetGridy(); c1->SetGridy(0); 
    double savRightMargin = c1->GetRightMargin(); c1->SetRightMargin(0.08);

    TString aswho = who; aswho.ReplaceAll("smcls_", "smacls_");
    TGraphAsymmErrors *as = showAs ? (TGraphAsymmErrors *) gFile->Get(aswho+"_obs") : 0;
    if (showAs && as == 0) std::cerr << "ERROR: cannot find " << aswho+"_obs" << std::endl;

    double xmin = obs->GetX()[0]             - obs->GetErrorXlow(0), xmin0 = xmin;
    double xmax = obs->GetX()[obs->GetN()-1] + obs->GetErrorXhigh(obs->GetN()-1);
    if (xmin <= 120 && xmax > 200) { xmin = 99.98; xmax = (SM == "FP" ? 300.1 : 600.1); }
    int col = colorFromName(who), dcol = colorFromName(who, /*dark=*/true);
    obs->SetLineWidth(2); 
    obs->SetLineColor(col); obs->SetMarkerColor(col); 
    obs->SetMarkerStyle(21); obs->SetMarkerSize(0.8);
    if (as) {
        as->SetLineColor(215);
        as->SetMarkerColor(215);
        as->SetLineStyle(2);
        as->SetMarkerStyle(24);
        as->SetMarkerSize(0.7);
        as->SetLineWidth(3);
    }
    int smooth = -5;
    double ymin = 5e-4, ymax = 9;
    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    TString what = "CL_{S} of "+SM+" Higgs hypothesis";
    setCanvas(&frame0, "", ymin, ymax, what);
    draw2(who, color68, color95, color50, true, false, smooth);
    frame0.Draw("AXIGSAME");
    obs->Draw("LP");
    if (as) as->Draw("LP");
    //spam("#splitline{"+SPAM+"}{"+channelFromName(who)+"}", 0.17,.15,.54,.20,32);
    TString myspam = "#splitline{"+SPAM+"}{"+channelFromName(who)+"}";
    if (isSquareCanvas && SPAM.Contains("Preliminary")) myspam = SPAM2L + "\n" + channelFromName(who);
    spam(myspam, 0.19,.88,.58-0.1*isSquareCanvas,.93);
    //leg = newLegend(.64,.15,.88,.31+showAs*0.04); leg->SetTextSize(0.04);
    leg = newLegend(.67-0.11*isSquareCanvas,.93-.16-showAs*0.04,.91,.93); leg->SetTextSize(0.04);
    m68->SetFillColor(color68); m68->SetLineColor(color50); m68->SetLineStyle(2); m68->SetLineWidth(2);
    m95->SetFillColor(color95); m95->SetLineColor(color50); m95->SetLineStyle(2); m95->SetLineWidth(2);
    leg->AddEntry(obs, "Observed", "LP");
    leg->AddEntry(m68, "Expected "+oneSigmaText, "LF");
    leg->AddEntry(m95, "Expected "+twoSigmaText, "LF");
    if (showAs) leg->AddEntry(as, "Asymptotic Obs.", "LP");
    frame0.Draw("AXIS SAME"); // re-draw ticks
    finalize(who+(showAs?"_comp":""),xmin,xmax,ymin,ymax);
    //finalize(who,xmin,xmax,ymin,ymax,SPAM,true);
    if (xmax > x_zoom && xmin < x_zoom) {
        c1->SetLogx(1);
        frame0.GetXaxis()->SetMoreLogLabels(); frame0.GetXaxis()->SetNoExponent();
        //finalize(who+"_logx",xmin,xmax,ymin,ymax,SPAM,true);
        finalize(who+(showAs?"_comp_logx":"_logx"),xmin,xmax,ymin,ymax);
        c1->SetLogx(0);
        xmin = xmin0;
        TH1D frame2("frame2","frame2", 1, xmin, x_zoom); frame2.Draw(); 
        setCanvas(&frame2, "", ymin, ymax, what);
        draw2(who, color68, color95, color50, true, false, smooth);
        frame2.Draw("AXIGSAME"); // grid
        obs->Draw("LP");
        if (as) as->Draw("LP");
        //spam("#splitline{"+SPAM+"}{"+channelFromName(who)+"}", 0.17,.15,.54,.20,32);
        spam(myspam, 0.19,.88,.58-0.1*isSquareCanvas,.93);
        frame2.Draw("AXIS SAME"); // ticks
        finalize(who+(showAs ? "_comp_zoom" : "_zoom"),xmin,x_zoom,ymin,ymax);
    }
    c1->SetRightMargin(savRightMargin); //c1->SetGridy(savGridY); 
}

void drawTwoSMCLs(TString who, TString who2, TString name2="Other", TString postfix="", int color2=kBlue) {
    TGraphAsymmErrors *obs   = (TGraphAsymmErrors *) gFile->Get(who+"_obs");
    TGraphAsymmErrors *m68 = (TGraphAsymmErrors *) gFile->Get(who+"_median");    
    TGraphAsymmErrors *m95 = (TGraphAsymmErrors *) gFile->Get(who+"_median_95"); 
    if (obs == 0 || m68 == 0 || m95 == 0) return;
    m68->SetFillColor(211); m68->SetLineColor(39);
    m95->SetFillColor(220); m95->SetLineColor(39);
    c1->SetRightMargin(0.08);

    TGraphAsymmErrors *obs2 = (TGraphAsymmErrors *) gFile->Get(who2+"_obs");
    TGraphAsymmErrors *exp2 = (TGraphAsymmErrors *) gFile->Get(who2+"_median");
    if (obs2 == 0 || exp2 == 0) { std::cerr << "Missing " << who2 << " to compare with " << who << std::endl; return; }
    obs2->SetLineColor(color2); obs2->SetLineWidth(5); obs2->SetLineStyle(1);
    exp2->SetLineColor(color2); exp2->SetLineWidth(5); exp2->SetLineStyle(2);

    double xmin = obs->GetX()[0]             - obs->GetErrorXlow(0), xmin0 = xmin;
    double xmax = obs->GetX()[obs->GetN()-1] + obs->GetErrorXhigh(obs->GetN()-1);
    if (xmin <= 120 && xmax > 200) { xmin = 99.98; xmax = (SM == "FP" ? 300.1 : 600.1); }
    int col = colorFromName(who), dcol = colorFromName(who, /*dark=*/true);
    obs->SetLineWidth(2); 
    obs->SetLineColor(col); obs->SetMarkerColor(col); 
    obs->SetMarkerStyle(21); obs->SetMarkerSize(0.8);
    double ymin = 8e-5, ymax = 1;
    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    TString what = "CL_{S} of "+SM+" Higgs hypothesis";
    setCanvas(&frame0, "", ymin, ymax, what);
    draw2(who, color68, color95, color50, true, false, 5);
    frame0.Draw("AXIGSAME");
    obs->Draw("LP");
    exp2->Draw("LX SAME"); obs2->Draw("LX SAME");
    spam("#splitline{"+SPAM+"}{"+channelFromName(who)+"}", 0.17,.15,.54,.20,32);
    leg = newLegend(.64,.15,.88,.37); leg->SetTextSize(0.04);
    m68->SetFillColor(color68); m68->SetLineColor(color50); m68->SetLineStyle(2); m68->SetLineWidth(2);
    m95->SetFillColor(color95); m95->SetLineColor(color50); m95->SetLineStyle(2); m95->SetLineWidth(2);
    leg->AddEntry(obs, "Observed", "LP");
    leg->AddEntry(m68, "Expected "+oneSigmaText, "LF");
    leg->AddEntry(m95, "Expected "+twoSigmaText, "LF");
    leg->AddEntry(obs2, name2+" Observed", "L");
    leg->AddEntry(exp2, name2+" Expected", "L");
    finalize(who+postfix,xmin,xmax,ymin,ymax);
    //finalize(who,xmin,xmax,ymin,ymax,SPAM,true);
    if (xmax >= 200 && xmin < x_zoom) {
        c1->SetLogx(1);
        frame0.GetXaxis()->SetMoreLogLabels(); frame0.GetXaxis()->SetNoExponent();
        //finalize(who+"_logx",xmin,xmax,ymin,ymax,SPAM,true);
        finalize(who+postfix+"_logx",xmin,xmax,ymin,ymax);
        c1->SetLogx(0);
        xmin = xmin0;
        TH1D frame2("frame2","frame2", 1, xmin, x_zoom); frame2.Draw(); 
        setCanvas(&frame2, "", ymin, ymax, what);
        draw2(who, color68, color95, color50, true, false, 5);
        frame2.Draw("AXIGSAME");
        obs->Draw("LP");
        exp2->Draw("LX SAME"); obs2->Draw("LX SAME");
        spam("#splitline{"+SPAM+"}{"+channelFromName(who)+"}", 0.17,.15,.54,.20,32);
        finalize(who+postfix+"_zoom",xmin,x_zoom,ymin,ymax);
    }
    c1->SetRightMargin(0.04);
}


void findCrossings(TString who, TString xname, double threshold, double xmin, double xmax, TString what="95% CL limit on #sigma/#sigma_{SM}") {
    TGraphAsymmErrors *obs = (TGraphAsymmErrors *) gFile->Get(who);
    if (who == 0) return;
    int ihigh = -1, ilow = -1;
    for (int i = 0, n = obs->GetN(); i < n; ++i) {
        if (xmin <= obs->GetX()[i] && obs->GetX()[i] <= xmax) {
            if (xname.Contains("low")) {
                if (obs->GetY()[i] < threshold) {
                    ilow = i;
                    break;
                } else { 
                    ihigh = i;
                }
            } else  {
                if (obs->GetY()[i] < threshold) {
                    ilow = i;
                } else { 
                    ihigh = i;
                    break;
                }
            }
        }
    }
    if (ilow == -1 || ihigh == -1) { std::cout << "didn't find points." << std::endl; return; }
    double x1 = obs->GetX()[ilow], x2 = obs->GetX()[ihigh];
    double y1 = obs->GetY()[ilow], y2 = obs->GetY()[ihigh];
    TF1 *linear  = new TF1("linear", "[0] * (x - [1]) + [2]", xmin, xmax); 
    linear->SetParameters((y2-y1)/(x2-x1), x1, y1);
    TF1 *linlog = new TF1("linlog", "[0] * pow([1], (x - [2])/[3])", xmin, xmax);
    linlog->SetParameters(y1, y2/y1, x1, x2-x1);
    TF1 *loglog  = new TF1("linear", "[0] * pow([1], log(x/[2])/[3])", xmin, xmax);
    loglog->SetParameters(y1, y2/y1, x1, log(x2/x1));
    loglog->SetLineWidth(3); loglog->SetLineColor( 63); 
    linlog->SetLineWidth(5); linlog->SetLineColor(210); linlog->SetLineStyle(2);
    linear->SetLineWidth(3); linear->SetLineColor(  1);  linear->SetLineStyle(9);
    obs->SetMarkerStyle(21); obs->SetMarkerSize(1.3);
    double x_linear = linear->GetX(threshold, xmin, xmax);
    double x_linlog = linlog->GetX(threshold, xmin, xmax);
    double x_loglog = loglog->GetX(threshold, xmin, xmax);
    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    linear->Draw("SAME");  linlog->Draw("SAME");     loglog->Draw("SAME");
    obs->Draw("PX");    
    double ymin = threshold/2, ymax = threshold * 3;
    double xleg = 0.67, yleg = 0.18;
    if (who.Contains("smcls")) { 
        ymax = 3*threshold; ymin = 0;  
        xleg = 0.27; yleg = 0.55;
    }
    if (xname.Contains("low")) {
        xleg = 0.67 ;
        yleg = 0.66;
    }
    leg = newLegend(xleg,yleg,xleg+.28,yleg+.25); leg->SetTextSize(0.04);
    leg->AddEntry(linear, Form("Linear,  m_{H} = %.1f", x_linear), "L");
    leg->AddEntry(linlog, Form("Lin-log, m_{H} = %.1f", x_linlog), "L");
    leg->AddEntry(loglog, Form("Log-log, m_{H} = %.1f", x_loglog), "L");
    setCanvas(&frame0, "", ymin, ymax, what);
    TString myspam = "#splitline{"+SPAM+"}{"+channelFromName(who)+"}";
    finalize(who+"_crossing_"+xname,xmin,xmax,ymin,ymax, myspam);
}
 

void printInfernalTable(TString fileName, double xmin, double xmax) {
    TGraphAsymmErrors *pvala_obs = (TGraphAsymmErrors *) gFile->Get("pvala_comb_obs");
    TGraphAsymmErrors *smcls_obs = (TGraphAsymmErrors *) gFile->Get("smcls_comb_obs");
    TGraphAsymmErrors *smcls_exp = (TGraphAsymmErrors *) gFile->Get("smcls_comb_median");
    TGraphAsymmErrors *cls_obs   = (TGraphAsymmErrors *) gFile->Get("cls_comb_obs");
    TGraphAsymmErrors *cls_exp   = (TGraphAsymmErrors *) gFile->Get("cls_comb_median");
    TGraphAsymmErrors *bayes_obs = (TGraphAsymmErrors *) gFile->Get("bayes_comb_obs");
    if (cls_obs == 0) return;
    
    FILE *fout = fopen(fileName.Data(), "w");
    fprintf(fout, 
        "\tHiggs boson &   observed        &   observed (expected)         &   \\multicolumn {2} {c|} {observed (expected) 95\\%% C.L. limits on $\\mu$}      \\\\\n");
    fprintf(fout, 
        "\tmass (GeV)  & $\\tilde p$-value  & $\\mathrm{CL_s}$ for $\\mu=1$   &   $\\mathrm{CL_s}$-method  &  Bayesian                  \\\\ \\hline \n");
    //printf( 
    //    "\tHiggs boson &   observed        &   observed (expected)         &   \\multicolumn {2} {c|} {observed (expected) 95\\%% C.L. limits on $\\mu$}      \\\\\n");
    //printf( 
    //    "\tmass (GeV)  & $\\tilde p$-value  & $\\mathrm{CL_s}$ for $\\mu=1$   &   $\\mathrm{CL_s}$-method  &  Bayesian                  \\\\ \\hline \n");
    for (int i = 0, n = cls_obs->GetN(); i < n; ++i) {
        int mass = cls_obs->GetX()[i]; float y_cls_obs = cls_obs->GetY()[i];
        if (mass < xmin || mass > xmax) continue;
        int ipvala_obs = findBin(pvala_obs, mass); float y_pvala_obs = (ipvala_obs != -1 ? pvala_obs->GetY()[ipvala_obs] : -1.);
        int ismcls_obs = findBin(smcls_obs, mass); float y_smcls_obs = (ismcls_obs != -1 ? smcls_obs->GetY()[ismcls_obs] : -1.);
        int ismcls_exp = findBin(smcls_exp, mass); float y_smcls_exp = (ismcls_exp != -1 ? smcls_exp->GetY()[ismcls_exp] : -1.);
        int icls_exp   = findBin(cls_exp, mass);   float y_cls_exp = (icls_exp != -1 ? cls_exp->GetY()[icls_exp] : -1.);
        int ibayes_obs = findBin(bayes_obs, mass); float y_bayes_obs = (ibayes_obs != -1 ? bayes_obs->GetY()[ibayes_obs] : -1.);
    fprintf(fout,
        "\t%3d         &    {\\bf{%5.3f}}   &     {\\bf{%5.3f}} (%5.3f)        &     {\\bf{%5.2f}} (%4.2f)    &    {\\bf{%5.2f}}       \\\\ \n",
        mass, y_pvala_obs, y_smcls_obs, y_smcls_exp, y_cls_obs, y_cls_exp, y_bayes_obs);
    //printf(
    //    "\t%3d         &    {\\bf{%5.3f}}   &     {\\bf{%5.3f}} (%5.3f)        &     {\\bf{%5.2f}} (%4.2f)    &    {\\bf{%5.2f}}       \\\\ \n",
    //    mass, y_pvala_obs, y_smcls_obs, y_smcls_exp, y_cls_obs, y_cls_exp, y_bayes_obs);
    }
    fclose(fout);
}

bool PLC_debug = false;
void drawCombObs(TString who, TString what="P.L. Approx limit #sigma_{95%}/#sigma_{"+SM+"}", const char **chann, int nchann, bool combined, TString postfix="") {
    TGraphAsymmErrors *obs = (TGraphAsymmErrors *) gFile->Get(who+"_comb_obs");
    if (obs == 0) return;
    double xmin = obs->GetX()[0]             - obs->GetErrorXlow(0), xmin0 = xmin;
    double xmax = obs->GetX()[obs->GetN()-1] + obs->GetErrorXhigh(obs->GetN()-1);
    if (xmin <= 120 && xmax > 200) { xmin = 99.98; xmax = (SM == "FP" ? 300.1 : 600.1); }
    bool isC3 = (postfix == "_c3");
    double ymin, ymax; minMaxY(obs, ymin, ymax); ymax *= (isC3 ? 2 : 8);
    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    TGraphAsymmErrors *obsi[nchann_all];
    float legxstart = isSquareCanvas ? (lessSpacesInLegends ? 0.60 : 0.50) 
                                     : (lessSpacesInLegends ? 0.72 : 0.62);
    leg = newLegend(legxstart,.65+0.10*isC3,.95-0.02*isSquareCanvas,.94); leg->SetTextSize(0.032);
    int ntot = 0, nlow = 0;
    for (int i = 0; i < nchann; ++i) {
        if (i == 0 && !combined) continue;
        obsi[i] = (TGraphAsymmErrors *) gFile->Get(who+"_"+chann[i]+"_obs");
        if (obsi[i] == 0) continue;
        int col = colorFromName(chann[i]);
        //obsi[i]->SetLineWidth(2);
        obsi[i]->SetLineWidth(i == 0 ? 4 : 3);
        obsi[i]->SetLineColor(col); obsi[i]->SetMarkerColor(col); 
        obsi[i]->SetMarkerStyle(21); obsi[i]->SetMarkerSize(0.8);
        obsi[i]->SetLineStyle(lineStyleFromName(chann[i]));
        if (!TString(chann[i]).Contains("_low")) {
            leg->AddEntry(obsi[i], channelFromName(chann[i],/*withspaces=*/true), who.Contains("cls") ? "L" : "LP");
            ntot++;
        }
        if (obsi[i]->GetX()[0] < x_zoom) nlow++;
    }
    for (int i = nchann-1; i >= 0; --i) {
        if (i == 0 && !combined) continue;
        if (obsi[i] == 0) continue;
        obsi[i]->Draw(who.Contains("cls") ? "LX" : "LPX");
    }
    setCanvas(&frame0, "", ymin, ymax, what);
    TString myspam = "#splitline{"+SPAM+"}{"+channelFromName(chann[0])+"}";
    if (isSquareCanvas && SPAM.Contains("Preliminary")) myspam = SPAM2L + "\n" + channelFromName(chann[0]);
    spam(myspam, 0.17,.89,.58,.94);
    leg->Draw();
    finalize(who+"_all"+postfix,xmin,xmax,ymin,ymax);
    if (xmax >= 200 && xmin < x_zoom) {
        c1->SetLogx(1);
        frame0.GetXaxis()->SetMoreLogLabels(); frame0.GetXaxis()->SetNoExponent();
        finalize(who+"_all"+postfix+"_logx",xmin,xmax,ymin,ymax);
        c1->SetLogx(0);
        xmin = xmin0;
        leg = newLegend(legxstart,.65+0.10*isC3+(ntot-nlow)*0.032,.95-0.02*isSquareCanvas,.94); leg->SetTextSize(0.032);
        TH1D frame2("frame2","frame2", 1, xmin, x_zoom); frame2.Draw(); 
        minMaxY(obs, ymin, ymax, 200.); ymax *= 15;
        for (int i = nchann-1; i >= 0; --i) {
            if (i == 0 && !combined) continue;
            if (obsi[i] == 0) continue;
            if (obsi[i]->GetX()[1] > x_zoom) continue;
            obsi[i]->Draw(who.Contains("cls") ? "LX" : "LPX");
        }
        for (int i = 0; i < nchann; ++i) {
            if (i == 0 && !combined) continue;
            if (obsi[i] == 0) continue;
            if (obsi[i]->GetX()[1] > x_zoom) continue;
            leg->AddEntry(obsi[i], channelFromName(chann[i],/*withspaces=*/true), who.Contains("cls") ? "L" : "LP");
        }
        leg->Draw();
        //spam(SPAM, 0.16,.85,.56,.91);
        spam(myspam, 0.17,.89,.58,.94);
        setCanvas(&frame2, "", ymin, ymax, what);
        finalize(who+"_all"+postfix+"_zoom",xmin,x_zoom,ymin,ymax);
    }
    leg = 0;
}

void drawCombBoth(TString who, TString what="P.L. Approx limit #sigma_{95%}/#sigma_{"+SM+"}", const char **chann, int nchann, bool observed=true, bool combined=true, TString postfix="") {
    TGraphAsymmErrors *obs = (TGraphAsymmErrors *) gFile->Get(who+"_comb_obs");
    if (obs == 0) return;
    double xmin = obs->GetX()[0]             - obs->GetErrorXlow(0), xmin0 = xmin;
    double xmax = obs->GetX()[obs->GetN()-1] + obs->GetErrorXhigh(obs->GetN()-1);
    bool isC3 = (postfix == "_c3");
    if (xmin <= 120 && xmax > 200) { xmin = 99.98; xmax = (SM == "FP" ? 300.1 : 600.1); }
    double ymin, ymax; minMaxY(obs, ymin, ymax); ymax *= (isC3?2:8);
    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    TGraphAsymmErrors *obsi[nchann_all], *expi[nchann_all];
    float legxstart = isSquareCanvas ? (lessSpacesInLegends ? 0.60 : 0.50) 
                                     : (lessSpacesInLegends ? 0.72 : 0.62);
    float legystart = .65-0.05*combined+0.14*isC3+0.03*lessSpacesInLegends+0.06*(SM=="FP");
    leg = newLegend(legxstart, legystart, .95-0.02*isSquareCanvas, .94); leg->SetTextSize(0.032);
    if (!observed) leg->SetHeader("Expected limits");
    int ntot = 0, nlow = 0;
    for (int i = 0; i < nchann; ++i) {
        if (i == 0 && !combined) continue;
        obsi[i] = (TGraphAsymmErrors *) gFile->Get(who+"_"+chann[i]+"_obs");
        expi[i] = (TGraphAsymmErrors *) gFile->Get(who+"_"+chann[i]+"_median");
        if (expi[i] == 0 || obsi[i] == 0) continue;
        if (TString(chann[i]).Contains("hzz2l2q") || 
            (noLineStyles == true && TString(chann[i]).Contains("comb")) ||
            (noLineStyles == true && TString(chann[i]) == "hzz")) {
            TGraphAsymmErrors *smooth = smoothSMCLs(expi[i], 7, 2);
            if (smooth == 0) std::cout << "Smoothing of " << expi[i]->GetName() << " returned ZERO" << std::endl;
            else expi[i] = smooth;
        }
        int col = colorFromName(chann[i]); //if (col == 96) col = 67;
        obsi[i]->SetLineWidth(i == 0 ? 4 : 3);
        obsi[i]->SetLineColor(col);  obsi[i]->SetMarkerColor(col); 
        obsi[i]->SetMarkerStyle(21); obsi[i]->SetMarkerSize(0.8);
        expi[i]->SetLineWidth(i == 0 ? 4 : (observed ? 3 : 4));
        expi[i]->SetLineColor(col);  expi[i]->SetMarkerColor(col); 
        expi[i]->SetMarkerStyle(21); expi[i]->SetMarkerSize(0.8);
        obsi[i]->SetLineStyle(observed ||  noLineStyles ? 1 : 2);
        expi[i]->SetLineStyle(observed || !noLineStyles ? 2 : 1);
        if (i == 0 && expi[i] != 0) {
            if (observed) {
                leg->AddEntry(obsi[i], channelFromName(chann[i],/*withspaces=*/true)+(lessSpacesInLegends?" obs.":" observed"), who.Contains("cls") ? "L" : "LP");
                leg->AddEntry(expi[i], channelFromName(chann[i],/*withspaces=*/true)+(lessSpacesInLegends?" exp.":" expected"), who.Contains("cls") ? "L" : "LP");
            } else {
                leg->AddEntry(expi[i], channelFromName(chann[i],/*withspaces=*/true), who.Contains("cls") ? "L" : "LP");
            }
            ntot++;
        } else {
            if (!TString(chann[i]).Contains("_low")) {
                leg->AddEntry(obsi[i], channelFromName(chann[i],/*withspaces=*/true), who.Contains("cls") ? "L" : "LP");
                ntot++;
            }
        }
        if (obsi[i]->GetX()[0] < x_zoom) nlow++;
    }
    for (int i = nchann-1; i >= 0; --i) {
        if (expi[i] == 0 || obsi[i] == 0) continue;
        expi[i]->Draw(who.Contains("cls") ? "LX" : "LPX");
        if (observed) obsi[i]->Draw(who.Contains("cls") ? "LX" : "LPX");
    }
    setCanvas(&frame0, "", ymin, ymax, what);
    TString myspam = "#splitline{"+SPAM+"}{"+channelFromName(chann[0])+"}";
    if (isSquareCanvas && SPAM.Contains("Preliminary")) myspam = SPAM2L + "\n" + channelFromName(chann[0]);
    spam(myspam, 0.17,.89,.58,.94);
    leg->Draw();
    finalize(who+(observed?"_all2":"_allexp")+postfix,xmin,xmax,ymin,ymax);
    if (xmax >= 200 && xmin < x_zoom) {
        c1->SetLogx(1);
        frame0.GetXaxis()->SetMoreLogLabels(); frame0.GetXaxis()->SetNoExponent();
        finalize(who+(observed?"_all2"+postfix+"_logx":"_allexp"+postfix+"_logx"),xmin,xmax,ymin,ymax);
        c1->SetLogx(0);
        xmin = xmin0;
        leg = newLegend(legxstart,legystart+(ntot-nlow)*0.032,.95-0.02*isSquareCanvas,.94); leg->SetTextSize(0.032);
        if (!observed) leg->SetHeader("Expected limits");
        TH1D frame2("frame2","frame2", 1, xmin, x_zoom); frame2.Draw(); 
        minMaxY(obs, ymin, ymax, 200.); ymax *= (isC3?2:15);
        for (int i = nchann-1; i >= 0; --i) {
            if (i == 0 && !combined) continue;
            if (expi[i] == 0 || obsi[i] == 0) continue;
            if (obsi[i]->GetX()[1] > x_zoom) continue;
            expi[i]->Draw(who.Contains("cls") ? "LX" : "LPX");
            if (observed) obsi[i]->Draw(who.Contains("cls") ? "LX" : "LPX");
        }
        for (int i = 0; i < nchann; ++i) {
            if (i == 0 && !combined) continue;
            if (expi[i] == 0 || obsi[i] == 0) continue;
            if (obsi[i]->GetX()[1] > x_zoom) continue;
            if (i == 0 && expi[i] != 0) {
                if (observed) {
                    leg->AddEntry(obsi[i], channelFromName(chann[i],/*withspaces=*/true)+(lessSpacesInLegends?" obs.":" observed"), who.Contains("cls") ? "L" : "LP");
                    leg->AddEntry(expi[i], channelFromName(chann[i],/*withspaces=*/true)+(lessSpacesInLegends?" exp.":" expected"), who.Contains("cls") ? "L" : "LP");
                } else {
                    leg->AddEntry(expi[i], channelFromName(chann[i],/*withspaces=*/true), who.Contains("cls") ? "L" : "LP");
                }
            } else {
                leg->AddEntry(obsi[i], channelFromName(chann[i],/*withspaces=*/true), who.Contains("cls") ? "L" : "LP");
            }
        }
        leg->Draw();
        spam(myspam, 0.17,.89,.58,.94);
        setCanvas(&frame2, "", ymin, ymax, what);
        finalize(who+(observed?"_all2"+postfix+"_zoom":"_allexp"+postfix+"_zoom"),xmin,x_zoom,ymin,ymax);
    }
    leg = 0;
}

void drawCombPVal(TString who, const char **chann, int nchann, TString postfix="", TString toywho="") {
    TGraphAsymmErrors *obs = (TGraphAsymmErrors *) gFile->Get(who+"_"+chann[0]+"_obs");
    if (obs == 0) return;
    TGraphAsymmErrors *exp = (TGraphAsymmErrors *) gFile->Get(who+"_"+chann[0]+"_median");
    if (exp == 0) exp = (TGraphAsymmErrors *) gFile->Get(who+"_"+chann[0]+"_asimov");

    TGraphAsymmErrors *toys = toywho != "" ? (TGraphAsymmErrors *) gFile->Get(toywho+"_obs") : 0;
    if (toys) {
        toys->SetLineWidth(4);
        toys->SetMarkerStyle(20);
        toys->SetMarkerSize(toys->GetN() > 3 && SM == "SM"? 1.3 : 1.8);
        toys->SetLineColor(92);
        toys->SetMarkerColor(92);
    }


    double xmin = obs->GetX()[0]             - obs->GetErrorXlow(0), xmin0 = xmin;
    double xmax = obs->GetX()[obs->GetN()-1] + obs->GetErrorXhigh(obs->GetN()-1);
    bool isC3 = (postfix == "_c3");
    if (xmin <= 120 && xmax > 200) { xmin = 99.98; xmax = (SM == "FP" ? 300.1 : 600.1); }
    if (TString(chann[0]).Contains("_toy_")) { xmin = x_zoom2_min; xmax = x_zoom2_max; }
    double ymin, ymax; minMaxY(obs, ymin, ymax, 999, 0); ymax = 1.0; ymin = 8e-7;
    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    setCanvas(&frame0, "", ymin, ymax, "Local p-value");
    TGraphAsymmErrors *obsi[nchann_all];
    float legxstart = isSquareCanvas ? (lessSpacesInLegends ? 0.60 : 0.50)
                                     : (lessSpacesInLegends ? 0.72 : 0.62);
    float legyend = .43-0.10*isC3; 
    if (toys) legyend += 0.04;
    TLegend *leg = newLegend(legxstart,.15,.95-0.02*isSquareCanvas,legyend); leg->SetTextSize(0.032);
    for (int i = 0; i < nchann; ++i) {
        obsi[i] = (TGraphAsymmErrors *) gFile->Get(who+"_"+chann[i]+"_obs");
        if (obsi[i] == 0) continue;
        int col = colorFromName(chann[i]);
        //obsi[i]->SetLineWidth(2);
        obsi[i]->SetLineWidth(i == 0 ? 4 : 4);
        obsi[i]->SetLineColor(col); obsi[i]->SetMarkerColor(col); 
        obsi[i]->SetMarkerStyle(21); obsi[i]->SetMarkerSize(0.8);
        obsi[i]->SetLineStyle(lineStyleFromName(chann[i]));
        if (!TString(chann[i]).Contains("_low")) {
            if (i == 0 && exp != 0) {
                leg->AddEntry(obsi[i], channelFromName(chann[i],/*withspaces=*/true)+(lessSpacesInLegends?" obs.":" observed"), "L");
                leg->AddEntry(exp,     (lessSpacesInLegends?"Exp. for SM Higgs":"Expected for SM Higgs"), "L");
                if (toys) leg->AddEntry(toys, "Comb. ensemble", "P");
            } else {
                leg->AddEntry(obsi[i], channelFromName(chann[i],/*withspaces=*/true), "L");
            }
        }
    }
    if (exp) {
        exp->SetLineWidth(4); exp->SetLineStyle(7); exp->SetLineColor(1);
        TGraphAsymmErrors *smooth = smoothSMCLs(exp, 7, 2);
        if (smooth == 0) std::cout << "Smoothing of " << expi[i]->GetName() << " returned ZERO" << std::endl;
        else exp = smooth;
        exp->Draw("LX");
    }
    if (toys) toys->Draw("P");
    for (int i = nchann-1; i >= 0; --i) {
        if (obsi[i] == 0) continue;
        obsi[i]->Draw("LX");
    }
    leg->Draw();
    TString myspam = "#splitline{"+SPAM+"}{"+channelFromName(chann[0])+"}";
    if (isSquareCanvas && SPAM.Contains("Preliminary")) myspam = SPAM2L + "\n" + channelFromName(chann[0]);
    finalize(who+"_all"+postfix,xmin,xmax,ymin,ymax,myspam,true);
    if (xmax >= x_zoom && xmin <= x_zoom) {
        c1->SetLogx(1);
        frame0.GetXaxis()->SetMoreLogLabels(); frame0.GetXaxis()->SetNoExponent();
        finalize(who+"_all"+postfix+"_logx",xmin,xmax,ymin,ymax,myspam,true);
        c1->SetLogx(0);
        leg = newLegend(legxstart,.153,.95-0.02*isSquareCanvas,legyend); leg->SetTextSize(0.031);
        xmin = xmin0;
        TH1D frame2("frame2","frame2", 1, xmin, x_zoom); frame2.Draw(); 
        setCanvas(&frame2, "", ymin, ymax, "Local p-value");
        if (exp) exp->Draw("LX");
        if (toys) toys->Draw("P");
        for (int i = nchann-1; i >= 0; --i) {
            if (obsi[i] == 0) continue;
            if (obsi[i]->GetX()[1] > x_zoom) continue;
            obsi[i]->Draw("LX");
        }
        for (int i = 0; i < nchann; ++i) {
            if (obsi[i] == 0) continue;
            if (obsi[i]->GetX()[1] > x_zoom) continue;
            if (i == 0 && exp != 0) {
                leg->AddEntry(obsi[i], channelFromName(chann[i],/*withspaces=*/true)+(lessSpacesInLegends?" obs.":" observed"), "L");
                leg->AddEntry(exp,     (lessSpacesInLegends?"Exp. for SM Higgs":"Expected for SM Higgs"), "L");
                if (toys) leg->AddEntry(toys, "Comb. ensemble", "P");
            } else {
                leg->AddEntry(obsi[i], channelFromName(chann[i],/*withspaces=*/true), "L");
            }
        }
        leg->Draw();
        finalize(who+"_all"+postfix+"_zoom",xmin,x_zoom,ymin,ymax,myspam,true);
    }
    leg = 0;
}

TGraphAsymmErrors *getNSigmaLine(TGraphAsymmErrors *g, int sigma/* = +/-1 */) {
    if (g == 0) return 0;
    int n = g->GetN();
    TGraphAsymmErrors *ret = new TGraphAsymmErrors(n);
    for (int i = 0; i < n; ++i) {
        ret->SetPoint(i, g->GetX()[i], g->GetY()[i] + (sigma > 0 ? g->GetErrorYhigh(i) : - g->GetErrorYlow(i)));
        ret->SetPointError(i, 0, 0, 0, 0);
    }
    ret->SetLineWidth(g->GetLineWidth() > 1 ? g->GetLineWidth() - 1 : 1);
    ret->SetLineColor(g->GetLineColor());
    ret->SetLineStyle(2);
    return ret;
}

void drawCombMuHat(TString who, double xmin, double xmax, double ymin, double ymax, TString label="", bool fill=true) {
    TGraphAsymmErrors *obs = (TGraphAsymmErrors *) gFile->Get(who+"_comb_obs");
    if (obs == 0) return;
    if (xmax == 600) xmax = 600.1;
    bool isZ = who.Contains("mlz");
    const int nch = 6;
    const char *muhats[nch] = { "vhbb", "htt",  "hzz4l", "hww", "hgg","comb" };
    const int    linec[nch] = {  209,     14,     215,    206,   217,  1     };
    const int    fillc[nch] = {   82,     17,      64,    208,    90,  1     };
    const int    fills[nch] = { 1001,   1001,    1001,   1001,  1001,  3244  };
    TGraphAsymmErrors *obsi[nchann], *obsi_h[nch], *obsi_l[nch];
    for (int i = 0; i < nch; ++i) {
        obsi[i] = (TGraphAsymmErrors *) gFile->Get(who+"_"+muhats[i]+"_obs"); 
        if (obsi[i] == 0) continue;
        //if (TString(muhats[i]) == "htt") continue;
        obsi[i]->SetLineColor(linec[i]); obsi[i]->SetLineWidth((i == nch-1 ? 5 : 3)+(fill?0:2));
        if (fill) { obsi[i]->SetFillColor(fillc[i]); obsi[i]->SetFillStyle(fills[i]); }
        obsi_h[i] = getNSigmaLine(obsi[i], +1);
        obsi_l[i] = getNSigmaLine(obsi[i], -1);
        if (i == nch-1) { obsi_h[i]->SetLineStyle(1); obsi_l[i]->SetLineStyle(1); }
    }
    leg = newLegend(.62,.73,.95,.94); leg->SetTextSize(0.032);
    leg->AddEntry(obsi[nch-1], "Combined", (fill ? "LF" : "L"));
    for (int i = 0; i < nch-1; ++i) {
        if (obsi[i]) leg->AddEntry(obsi[i], channelFromName(muhats[i],/*withspaces=*/true), (fill ? "LF" : "L"));
    }
    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    for (int i = 0; i < nch; ++i) { if (obsi[i]) { if (fill) obsi[i]->Draw("E3"); } }
    for (int i = 0; i < nch; ++i) { if (obsi[i]) { obsi_h[i]->Draw("L"); obsi_l[i]->Draw("L");}}
    for (int i = 0; i < nch; ++i) { if (obsi[i]) { obsi[i]->Draw("LX"); } }
    setCanvas(&frame0, "", ymin, ymax, "Best fit #sigma/#sigma_{SM}");
    frame0.Draw("AXIG SAME");
    spam(SPAM, 0.17,.89,.58,.94);
    leg->Draw();
    finalize(who+"_all"+label,xmin,xmax,ymin,ymax);
    leg = 0;
}

void cccPlot(TString who, float where, const char **mychann, int nchann, TString postfix="", double rMin=-2.0, double rMax=5.0) {
    TGraphAsymmErrors *obs = (TGraphAsymmErrors *) gFile->Get(who+"_"+mychann[0]+"_obs");
    if (obs == 0) return;
    int j = findBin(obs, where); if (j == -1) return;
    double r0 = obs->GetY()[j], r0min = r0 - obs->GetErrorYlow(j), r0max = r0 + obs->GetErrorYhigh(j);
    bool halfint = (where - floor(where) > 0.3);
    c1->SetLogx(0); c1->SetLogy(0); c1->SetTicky(0);  c1->SetTickx(0);
    TGraphAsymmErrors *points = new TGraphAsymmErrors(); int npoints = 0;
    TString names[nchann_all];
    for (int i = nchann-1; i > 0; --i) {
       TGraphAsymmErrors *obsi = (TGraphAsymmErrors *) gFile->Get(who+"_"+mychann[i]+"_obs"); if (obsi == 0) continue;
       int ji = findBin(obsi, where); 
       if (ji != -1) {
           points->Set(++npoints);
           points->SetPoint(npoints-1, obsi->GetY()[ji], npoints - 0.5);
           points->SetPointError(npoints-1, obsi->GetErrorYlow(ji), obsi->GetErrorYhigh(ji), 0., 0.);
           names[npoints-1] = channelFromName(mychann[i],true,true);
       } else if (halfint && ((ji = findBin(obsi, floor(where))) != -1) && findBin(obsi, ceil(where)) != -1) {
           points->Set(++npoints);
           points->SetPoint(npoints-1, 0.5*(obsi->GetY()[ji]+obsi->GetY()[ji+1]), npoints - 0.5);
           points->SetPointError(npoints-1, 0.5*(obsi->GetErrorYlow(ji) + obsi->GetErrorYlow(ji+1)), 
                                            0.5*(obsi->GetErrorYhigh(ji) + obsi->GetErrorYhigh(ji+1)), 
                                            0., 0.);
           names[npoints-1] = channelFromName(mychann[i],true,true);
       }
    }
    TH2F frame("frame",";Best fit #sigma/#sigma_{SM};",1, rMin, rMax, npoints+2, 0., npoints+2);
    for (int k = 1; k <= npoints; ++k) {
       frame.GetYaxis()->SetBinLabel(k, names[k-1]);
    }
    frame.GetYaxis()->SetTickLength(0);

    double leftMargin = c1->GetRightMargin(); c1->SetLeftMargin(0.24+0.05*(where > 130));
    points->SetLineColor(kRed);
    points->SetLineWidth(5);
    points->SetMarkerStyle(21);
    points->SetMarkerSize(1.7);
    frame.GetXaxis()->SetTitleSize(0.05);
    frame.GetXaxis()->SetLabelSize(0.04);
    frame.GetYaxis()->SetLabelSize(0.06);
    frame.Draw(); gStyle->SetOptStat(0);
    TBox globalFitBand(r0min, 0, r0max, npoints); TH1F fake("fake","fake",1,0,1); 
    //globalFitBand.SetFillStyle(3344); fake.SetFillStyle(3344);
    //globalFitBand.SetFillColor(65);   fake.SetFillColor(65);
    globalFitBand.SetFillStyle(1001); fake.SetFillStyle(1001);
    globalFitBand.SetFillColor(colorFit68);   fake.SetFillColor(colorFit68);
    globalFitBand.SetLineStyle(0);    
    globalFitBand.DrawClone();
    TLine globalFitLine(r0, 0, r0, npoints);
    globalFitLine.SetLineWidth(4);     fake.SetLineWidth(0);
    //globalFitLine.SetLineColor(214); fake.SetLineColor(214); 
    globalFitLine.SetLineColor(1); fake.SetLineColor(colorFit68); 
    globalFitLine.DrawClone();
    globalFitLine.SetLineStyle(2);  
    points->Draw("P SAME");
    frame.Draw("AXIS SAME");
    double xoff = 0; //-0.1*SPAM.Contains("Preliminary");
    leg = newLegend(.26+xoff,.74,.63+xoff,.935); leg->SetTextSize(0.037); // was TextSize 0.4, and goind up to 0.65+xoffs
    if (halfint) {
        leg->SetHeader(Form("   m_{H} = %.1f %s", where, massUnits.Data()));
    } else {
        leg->SetHeader(Form("   m_{H} = %.0f %s", where, massUnits.Data()));
    }
    leg->AddEntry(&fake,  "Combined "+oneSigmaFitCCC, "F"); //LEF
    //leg->AddEntry(points, "Single channel "+oneSigmaFitCCC, "LP");
    leg->AddEntry(points, "Single channel", "LP");
    if (SPAM.Contains("Preliminary")) {
        spam("CMS Preliminary\n#sqrt{s} = 7 TeV\n"+channelFromName(mychann[0]), .65+xoff, .88, .94, .93, 22);
    } else {
        spam("#splitline{"+SPAM+"}{"+channelFromName(mychann[0])+"}", .65+xoff, .88, .94, .93, 22+10*SPAM.Contains("P") );
    }
    leg->Draw();
    justSave(Form("%s_ccc_mH%.1f%s", who.Data(), where, postfix.Data()));
    c1->SetLeftMargin(leftMargin); c1->SetTicky(1); c1->SetTickx(1);
}

void writeGrid(TString who, TString inputPrefix) {
    TGraphAsymmErrors *g = makeAPrioriGrid(TString::Format("acls_%s", who.Data()));
    if (g == 0) continue;
    g = slidingWindowAverage(g, 3);
    FILE *grid = fopen(TString::Format("grids/%s%s.txt", inputPrefix.Data(), who.Data()).Data(), "w");
    if (grid == 0) { std::cout << "Cannot write to grid " << TString::Format("grids/%s%s.txt", inputPrefix.Data(), who.Data()).Data() << std::endl; continue; }
    for (int j = 0, n = g->GetN(); j < n; ++j) {
        double x  = g->GetX()[j];
        double y  = g->GetY()[j];
        double yl = y-g->GetEYlow()[j];
        double yh = y+g->GetEYhigh()[j];
        fprintf(grid, "%3.0f %5.2f  %5.2f  %5.2f\n", x, yl, y, yh);
    }
    fclose(grid);
    CLs_grids.Add(g);
}
void writeGrids(TString inputPrefix) {
    for (int i = 0; i < nchann_all; ++i) {
        writeGrid(chann[i], inputPrefix);
        break;
    }
}
void readGrids(TString inputPrefix) {
    for (int i = 0; i < nchann; ++i) {
        FILE *fgrid = fopen(TString::Format("grids/%s%s.txt", inputPrefix.Data(), chann[i]).Data(), "r");
        if (fgrid == 0) continues;
        float x,y,yl,yh;
        TGraphAsymmErrors *grid = new TGraphAsymmErrors(); int points = 0;
        while(fscanf(fgrid,"%f %f %f %f", &x, &yl, &y, &yh) == 4) {
            grid->Set(points+1);
            grid->SetPoint(points,x,y);
            grid->SetPointError(points,0,0,y-yl,yh-y);
            points++;
        }
        std::cout << "Input grid for " << chann[i] << " contains " << points << " points" << std::endl;
        fclose(fgrid);
        grid->SetName(TString::Format("grid_cls_%s", chann[i]));
        grid->SetFillColor(16);
        grid->SetFillStyle(3244);
        CLs_grids.Add(grid);
    }
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

void plots(int which) {
  std::cout << 1 << std::endl;
    gROOT->LoadMacro("$CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/plotting/bandUtils.cxx+");
    std::cout << 2 << std::endl;
    gROOT->ProcessLine(".x tdrstyle.cc");
    std::cout << 3 << std::endl;
   TString inputPrefix = "";
    switch (which)  {
        case 1: inputPrefix = "SM4_"; SM="SM4"; 
                LEESPAM_2L = "#splitline{Global p-value ~0.5}{for 110-600 GeV range}";
                break;
        case 2: inputPrefix = "FF_";  SM="FP";  
                LEESPAM_2L = "#splitline{Global significance 1.1#sigma}{for 110-300 GeV range}";
                break;
    }
    std::cout << 4 << std::endl;

    loadMasses("masses.txt");
    std::cout << 5 << std::endl;
    TString globalPrefix0 = inputPrefix+"plots/";
    globalPrefix =  globalPrefix0+"/";
    gStyle->SetOptTitle(0);
    gSystem->Exec("cp results/"+inputPrefix+"bands.root "+inputPrefix+"bands.work.root");
    std::cout << 6 << std::endl;
  
  TFile::Open(inputPrefix+"bands.work.root", "UPDATE");
  std::cout << 7 << std::endl;

  std::cout << inputPrefix << std::endl;
    // WRITE OUT THE GRIDS TO RUN CLS AND BAYESIAN
  //writeGrids(inputPrefix); return; 
  writeGrid("hzz4l",""); return; 
  switch (which) {
        case 0:
            selectedPointsBands(gFile, "pval_comb", "pval_comb125", 125.0);
            selectedPointsBands(gFile, "pval_comb", "pval_comb12x", 125.0, 124.5, 124.0, 125.5, 126.0);
            break;
        case 1:
            selectedPointsBands(gFile, "pval_comb", "pval_comb125", 119.0);
            selectedPointsBands(gFile, "pval_comb", "pval_comb12x", 119.0, 120.0, 118.0, 121.0, 117.0);
            break;
        case 2:
            selectedPointsBands(gFile, "pval_comb", "pval_comb125", 126.0);
            selectedPointsBands(gFile, "pval_comb", "pval_comb12x", 126.0, 124.0, 125.0, 127.0, 128.0);
            break;
    }
    std::cout << 8 << std::endl;

    c1 = new TCanvas("c1","c1");
    rectangleCanvas();
    std::cout << 9 << std::endl;

#if 0
    // PAPER PLOTS, keep it simple
    squareCanvas(/*grid x,y=*/0,0); x_zoom = 145;
    globalPrefix = globalPrefix0+"/pre-paper/sqr_";
    noLineStyles = true; // always use solid lines

    // fig 1
    //drawSMCLs(TString::Format("smcls_%s", "comb"));

    // fig 2
    //drawOneCLs(TString::Format("cls_%s",  "comb"));

    // fig 3
    drawCombObs("acls", "95% CL limit on #sigma/#sigma_{"+SM+"}", chann, nchann);
    lineAt1Style = 2;
    drawCombBoth("acls", "95% CL limit on #sigma/#sigma_{"+SM+"}", chann, nchann, /*obs=*/false, /*combined=*/false);
    lineAt1Style = 1;

    // fig 4
    forceYmin = 0.08; forceYmax = 12;
    channelSpamOnRightHandSide = true;
    drawOneCLs(TString::Format("acls_%s",  "hzz4l"));
    //drawOneCLs(TString::Format("acls_%s",  "combl_full"));
    channelSpamOnRightHandSide = false;
    forceYmin = 0; forceYmax = 0;
    //fig 5a
    drawCombPVal("pvala", chann, nchann);

    //fig 5b
    drawOnePlot(TString::Format("mlz_%s", "hzz4l")); 

    // fig 6
    globalPrefix = globalPrefix0+"/pre-paper/";
    cccPlot("mlz", 119.5,  chann, nchann, "", -1.0, 4.0);
    cccPlot("mlz", 124.0,  chann, nchann, "", -1.0, 3.0);
    return;
#endif

#if 0
    // PAS PLOTS, keep it simple
    squareCanvas(/*grid x,y=*/0,0); 
    globalPrefix = globalPrefix0+"/sqr_";
    noLineStyles = true; // always use solid lines
    lessSpacesInLegends = true; // compact legends (we only use decay modes combined)

    if (which == 0) {
        drawSMCLs(TString::Format("smcls_%s", "comb"));

        drawOneCLs(TString::Format("cls_%s",  "comb"));

        drawCombBoth("acls", "95% CL limit on #sigma/#sigma_{"+SM+"}", chann2, nchann2, true,  true, "_zzc");
        
        drawCombPVal("pvala", chann2, nchann2, "_zzc");
        drawOnePlot(TString::Format("mlz_%s", "comb")); 
        
        globalPrefix = globalPrefix0+"/";
        cccPlot("mlz", 119.5,  chann, nchann, "", -2.0, 4.0);
        cccPlot("mlz", 125.0,  chann, nchann, "", -1.0, 4.0);
    } else if (which == 2) {
        drawOneCLs(TString::Format("acls_%s",  "comb"));

        drawCombBoth("acls", "95% CL limit on #sigma/#sigma_{"+SM+"}", chann2, nchann2, true,  true, "_zzc");
        
        drawCombPVal("pvala", chann2, nchann2, "_zzc");
        drawOnePlot(TString::Format("mlz_%s", "comb")); 
    }
    return;
#endif


#if 0
    squareCanvas(); globalPrefix =  globalPrefix0+"/"; c1->SetGridy(0); 
    switch (which) {
        case 0:
        for (int gridx = 0; gridx < 2; ++gridx) {
            globalPrefix =  globalPrefix0+(gridx ? "/grid_" : "/");
            c1->SetGridx(gridx);
            cccPlot("mlz", 119.5,  chann2, nchann2, "", -3.0, 4.0);
            cccPlot("mlz", 124.0,  chann2, nchann2, "", -1.0, 4.0);
            cccPlot("mlz", 125.0,  chann2, nchann2, "", -1.0, 4.0);
            cccPlot("mlz", 136.,  chann2, nchann2, "", -1.0, 7.0);
            cccPlot("mlz", 144.,  chann2, nchann2, "", -1.5, 5.0);
        }
        break;
        case 2:
        for (int gridx = 0; gridx < 2; ++gridx) {
            globalPrefix =  globalPrefix0+(gridx ? "/grid_" : "/");
            c1->SetGridx(gridx);
            cccPlot("mlz", 119.5,  chann2, nchann2, "", -1.0, 5.0);
            cccPlot("mlz", 124.0,  chann2, nchann2, "", -1.0, 4.0);
            cccPlot("mlz", 125.0,  chann2, nchann2, "", -1.0, 4.0);
            cccPlot("mlz", 126.0,  chann2, nchann2, "", -1.0, 4.0);
            cccPlot("mlz", 130.5,  chann2, nchann2, "", -1.0, 7.0);
        }
        break;
    }
    return;
#endif

    const char *prefixes[4] = { "", "sqr_", "sqr_grid_", "grid_" };
    const int   grids   [4] = {  0, 0,      1,           1,      };
    const int   squared [4] = {  0, 1,      1,           0,      };
    for (int version = 0; version < 4; ++version) {
        //if (version != 1) continue;
        if (squared[version]) squareCanvas(grids[version],grids[version]);
        else rectangleCanvas(grids[version],grids[version]);

        globalPrefix =  globalPrefix0+"/"+prefixes[version];
#if 1   
        noLineStyles = true;
        //if (gFile->Get("acls_comb_obs")) drawCombObs("acls",  "95% CL limit on #sigma/#sigma_{"+SM+"}", chann4, nchann4, true);
        //if (gFile->Get("acls_comb_obs")) drawCombBoth("acls", "95% CL limit on #sigma/#sigma_{"+SM+"}", chann4, nchann4);
        //if (gFile->Get("acls_comb_obs")) drawCombBoth("acls", "95% CL limit on #sigma/#sigma_{"+SM+"}", chann4, nchann4, false);
        noLineStyles = false;
	std::cout << 10 << std::endl;

        // ----- ZZ combined ------
        lessSpacesInLegends = true; 
        noLineStyles = true;
        if (gFile->Get("acls_comb_obs")) drawCombObs("acls",  "95% CL limit on #sigma/#sigma_{"+SM+"}", chann2, nchann2, true,        "_zzc"); 
        if (gFile->Get("acls_comb_obs")) drawCombBoth("acls", "95% CL limit on #sigma/#sigma_{"+SM+"}", chann2, nchann2, true,  true, "_zzc");
        if (gFile->Get("acls_comb_obs")) drawCombBoth("acls", "95% CL limit on #sigma/#sigma_{"+SM+"}", chann2, nchann2, false, true, "_zzc");
	std::cout << 11 << std::endl;
	noLineStyles = true;
        //if (gFile->Get("acls_comb_obs")) drawCombObs("acls",  "95% CL limit on #sigma/#sigma_{"+SM+"}", chann2, nchann2, true,        "_zzc_solid"); 
        //if (gFile->Get("acls_comb_obs")) drawCombBoth("acls", "95% CL limit on #sigma/#sigma_{"+SM+"}", chann2, nchann2, false, true, "_zzc_solid");
        noLineStyles = false;
        lessSpacesInLegends = false;
	std::cout << 12 << std::endl;

        // ----- Two SubCombs -----
        //if (gFile->Get("acls_comb_obs")) drawCombBoth("acls", "95% CL limit on #sigma/#sigma_{"+SM+"}", chann3, nchann3, true, true, "_c3");
        // ====== PVALS ======
#endif
        noLineStyles = true;
        if (gFile->Get("pvala_comb_obs")) drawCombPVal("pvala", chann, nchann);
        if (gFile->Get("pvala_comb_obs")) { lessSpacesInLegends = true; drawCombPVal("pvala", chann2, nchann2, "_zzc"); lessSpacesInLegends = false; }
        if (gFile->Get("pvala_comb_obs")) { lessSpacesInLegends = true; drawCombPVal("pvala", chann2, nchann2, "_zzc_t",  "pval_comb125"); lessSpacesInLegends = false; }
        if (gFile->Get("pvala_comb_obs")) { lessSpacesInLegends = true; drawCombPVal("pvala", chann2, nchann2, "_zzc_t2", "pval_comb12x"); lessSpacesInLegends = false; }
        doLEESPAM = true;
        if (gFile->Get("pvala_comb_obs")) { lessSpacesInLegends = true; drawCombPVal("pvala", chann2, nchann2, "_zzc_gz"); lessSpacesInLegends = false; }
        if (gFile->Get("pvala_comb_obs")) { lessSpacesInLegends = true; drawCombPVal("pvala", chann2, nchann2, "_zzc_t_gz", "pval_comb125"); lessSpacesInLegends = false; }
        if (gFile->Get("pvala_comb_obs")) { lessSpacesInLegends = true; drawCombPVal("pvala", chann2, nchann2, "_zzc_t2_gz", "pval_comb12x"); lessSpacesInLegends = false; }
        doLEESPAM = false;
        //if (gFile->Get("pvala_comb_obs")) drawCombPVal("pvala", chann3, nchann3, "_c3");
        noLineStyles = false;
        //drawCombMuHat("mlz", 110, 150, -2.5, 5, "");  
        //drawCombMuHat("mlz", 110, 150, -2.5, 5, "_nofill", 0);  
        //drawCombMuHat("mlz", 114, 132, -2, 4, "_zoom");  
        //drawCombMuHat("mlz", 114, 132, -2, 4, "_zoom_nofill", 0);  


        for (int i = 0; i < nchann; ++i) {
            if (!TString(chann[i]).Contains("comb")) { 
                globalPrefix =  globalPrefix0+"/subchannels/"+prefixes[version];
            } else { 
                globalPrefix =  globalPrefix0+"/"+prefixes[version];
            }

            //if (!TString(chann[i]).Contains("comb")) continue;
            //globalPrefix =  globalPrefix0+"/private/";

            //if (which <= 1 && (TString(chann[i]) != "comb" && TString(chann[i]) != "hgg")) continue;
            //if (i < 1+5+3+3) continue;
            if (TString(chann[i]) != "comb") continue;
            //if (!TString(chann[i]).Contains("comb")) continue;
            //if (TString(chann[i]) != "vhtt") continue;
#if 0       //---- Asymptotics ----
            //CLs_debug_apriori_grid = TString(chann[i]).Contains("comb");
	    std::cout << 13 << std::endl;
	    std::cout << chann[i] << std::endl;
            drawOneCLs(TString::Format("acls90_%s", "hzz4l"));
            drawOneCLs(TString::Format("acls99_%s", "hzz4l"));
	    
            //drawOnePVal(TString::Format("pvala_%s", chann[i]), "Local p-value", /*include_toys=*/false);
	    std::cout << 14 << std::endl;
	    //drawOnePlot(TString::Format("mlz_%s", chann[i])); 
            ////drawOnePValFitCollage(TString::Format("%s", chann[i]), "mlz");
#endif
            drawOneCLs(TString::Format("acls_%s", "hzz4l"));

        
#if 0       //---- TOYS -----
            drawOneCLs(TString::Format("cls_%s",  chann[i]));
            drawSMCLs(TString::Format("smcls_%s", chann[i]));
            drawMethods(TString::Format("cls_%s", chann[i]), /*include_pla=*/true);
            //drawSMCLs(TString::Format("smcls_%s", chann[i]), true);
            ////drawOnePVal(TString::Format("pvala_%s", chann[i]), "Local p-value", /*toys=*/true);
            //if (gFile->Get(Form("bayes_%s_obs", chann[i]))) drawMethods(TString::Format("acls_%s", chann[i]), /*include_pla=*/false);
#endif

#if 0       // P-values with bands
            if (i == 0) {
                drawOnePValBand("pvala_comb");
                drawOnePValBand("pval_comb");
                drawOnePValBand("pvala_comb","pval_comb");
                drawOnePValBand("pval_comb","pvala_comb");
            }
#endif

#if 0       // P.R. PLOTS
            if (i == 0) {
                globalPrefix =  globalPrefix0+"/lte_"+prefixes[version]; lep_excluded = true; tev_excluded = true; 
                drawOneCLs(TString::Format("cls_%s", chann[i]));
                globalPrefix =  globalPrefix0+"/ltc_"+prefixes[version]; lep_excluded = true; tev_excluded = true; cms_excluded = true;
                drawOneCLs(TString::Format("cls_%s", chann[i]));
                globalPrefix =  globalPrefix0+"/ltl_"+prefixes[version]; lep_excluded = true; tev_excluded = false; cms_excluded = false; Draw_TEV = true; 
                drawOneCLs(TString::Format("cls_%s", chann[i]));
                globalPrefix =  globalPrefix0+"/"+prefixes[version]; lep_excluded = false; tev_excluded = false; Draw_TEV = false; cms_excluded = false;
            }
#endif
#if 0       // NOW VS PAPER, COOL
            if (which == 0 && i == 0) {
                drawTwoCLs(TString::Format("cls_%s", chann[i]),  TString::Format("cls_%s_paper", chann[i]), "Previous", "_hist", /*color=*/59, /*doObs2=*/false);
            }
#endif

#if 0       // NOW VS PAPER
            if (which == 0) {
                TString herePref = globalPrefix;
                globalPrefix =  globalPrefix0+"/history/"+prefixes[version];
                TString channLP = chann[i];
                if (gFile->Get(TString::Format("cls_%s_obs", chann[i])) && gFile->Get(TString::Format("cls_%s_paper_obs", channLP.Data()))) {
                    std::cout << "Doing 1-1 comparison of CLs plots for " << chann[i] << std::endl;
                    drawTwoCLs(TString::Format("cls_%s", chann[i]),  TString::Format("cls_%s_paper", channLP.Data()), "HIG-11", "", /*color=*/59);
                } else if (gFile->Get(TString::Format("acls_%s_obs", chann[i])) && gFile->Get(TString::Format("acls_%s_paper_obs", channLP.Data()))) {
                    std::cout << "Doing 1-1 comparison of Asym CLs plots for " << chann[i] << std::endl;
                    drawTwoCLs(TString::Format("acls_%s", chann[i]), TString::Format("acls_%s_paper", channLP.Data()), "HIG-11", "", /*color=*/59);
                }
                if (gFile->Get(TString::Format("smcls_%s_obs", chann[i])) && gFile->Get(TString::Format("smcls_%s_paper_obs", channLP.Data()))) {
                    std::cout << "Doing 1-1 comparison of SM CLs plots for " << chann[i] << std::endl;
                    drawTwoSMCLs(TString::Format("smcls_%s", chann[i]),  TString::Format("smcls_%s_paper", channLP.Data()), "HIG-11", "", /*color=*/59);
                }
                globalPrefix = herePref;
        }
#endif

        } // loop over plots
    } // loop over versions

#if 0   // TIMING 
        globalPrefix = globalPrefix0+"/timing/";
        drawOnePlot(TString::Format("timecls_%s", chann[i]), "Time per CLs toy [s]");
#endif

#if 0 // Interpolation stuff
    globalPrefix =  globalPrefix0+"/interpolation/";
    //findCrossings("smcls_comb_obs", "low", 0.001, 120, 135, "CL_{S} of SM Higgs boson hypothesis");
    findCrossings("smcls_comb_obs", "low99", 0.01,  120, 135, "CL_{S} of SM Higgs boson hypothesis");
    findCrossings("smcls_comb_obs", "low95", 0.05,  120, 135, "CL_{S} of SM Higgs boson hypothesis");
    //findCrossings("smcls_comb_obs", "high", 0.001, 120, 135, "CL_{S} of SM Higgs boson hypothesis");
    findCrossings("smcls_comb_obs", "high99", 0.01,  450, 600, "CL_{S} of SM Higgs boson hypothesis");
    //findCrossings("smcls_comb_obs", "high95", 0.05,  500, 600, "CL_{S} of SM Higgs boson hypothesis");
    //findCrossings("smcls_comb_median", "low", 0.001, 120, 135, "CL_{S} of SM Higgs boson hypothesis");
    findCrossings("smcls_comb_median", "low99", 0.01, 120, 135, "CL_{S} of SM Higgs boson hypothesis");
    findCrossings("smcls_comb_median", "low95", 0.05, 110, 125, "CL_{S} of SM Higgs boson hypothesis");
    //findCrossings("smcls_comb_median", "high", 0.001, 120, 135, "CL_{S} of SM Higgs boson hypothesis");
    findCrossings("smcls_comb_median", "high99", 0.01, 400, 550, "CL_{S} of SM Higgs boson hypothesis");
    findCrossings("smcls_comb_median", "high95", 0.05, 500, 600, "CL_{S} of SM Higgs boson hypothesis");
    findCrossings("cls_comb_obs",      "low95", 1.0,  120, 135, "95% CL limit on #sigma/#sigma_{SM}");
    findCrossings("cls_comb_median",   "low95", 1.0,  110, 125, "95% CL limit on #sigma/#sigma_{SM}");
    findCrossings("cls_comb_obs",      "low95", 1.0,  120, 135, "95% CL limit on #sigma/#sigma_{SM}");
    findCrossings("cls_comb_median",   "low95", 1.0,  110, 125, "95% CL limit on #sigma/#sigma_{SM}");
    return;
#endif


#if 0
    /// ASYMPTOTICS VS TOYS
    rectangleCanvas(); 
    globalPrefix =  globalPrefix0+"/asymptotic-valid/";
    for (int i = 0; i < nchann; ++i) {
        drawCompLimit( TString::Format("acls_%s_obs",       chann[i]), TString::Format("cls_%s_obs",       chann[i]), "Asymp.", "Toys");
        drawCompLimit( TString::Format("acls_%s_median",    chann[i]), TString::Format("cls_%s_median",    chann[i]), "Asymp.", "Toys");
        drawCompLimit( TString::Format("acls_%s_median_95", chann[i]), TString::Format("cls_%s_median_95", chann[i]), "Asymp.", "Toys");
    }
#endif
#if 0
    // ROOSTATS VS LANDS
    globalPrefix =  globalPrefix0+"/validation/";
    for (int i = 0; i < 1; ++i) {
        drawCompLimit( TString::Format("acls_%s_obs",  chann[i]), TString::Format("acls_%s_lands_obs",  chann[i]), "RooStats", "LandS", "Asym CL_{S} obs. limit");
        drawCompLimit( TString::Format("acls_%s_median",  chann[i]), TString::Format("acls_%s_lands_median",  chann[i]), "RooStats", "LandS", "Asym CL_{S} exp. limit");
        drawCompPVal(  TString::Format("pvala_%s_obs",    chann[i]), TString::Format("pvala_%s_lands_obs", chann[i]), "RooStats", "LandS", "Local p-value obs.");
    }
#endif

#if 0
    globalPrefix =  globalPrefix0+"/private/";
    //drawCompPVal("pvala_comb_asimov", "pval_comb_median",  "Asimov", "Med Toys", "Local p-value exp", 1e-8); 
    //drawCompPVal("pvala_comb_asimov", "pvalh_comb_median", "Asimov", "Med Asymp", "Local p-value exp", 1e-8); 
    //drawCompPVal("pvalh_comb_median", "pval_comb_median",  "Med Asymp", "Med Toys", "Local p-value exp", 1e-8); 
    //drawCompPVal("pvala_comb_obs", "pvala_comb_lands_obs", "RooStats", "LandS", "Local p-value obs.");
    drawOnePValBand("pvala_comb");
    drawOnePValBand("pval_comb");
    drawOnePValBand("pvala_comb","pval_comb");
    drawOnePValBand("pval_comb","pvala_comb");
#endif


#if 0
    globalPrefix =  globalPrefix0+"/pvtests/";
    const int npvflav = 6;
    const char *pvflav[npvflav] = { "M2M", "S", "SZ3", "BZ1", "BZ", "TQ0" };
    for (int i = 0; i < npvflav; ++i) {
        std::cout << "Compare with " << Form("pvala_%s_comb_obs", pvflav[i]) << std::endl;
        drawCompPVal(Form("pvala_%s_comb_obs", pvflav[i]), "pvala_comb_obs",   pvflav[i], "Median", "Local p-value"); 
        //drawCompPVal(Form("pvala_%s_comb_obs", pvflav[i]), "pval_comb_obs",    pvflav[i], "Toys", "Local p-value"); 
        //drawCompPVal(Form("pvala_%s_comb_obs", pvflav[i]), "pvala_BZ_comb_obs", pvflav[i], "BZ", "Local p-value"); 
    }
#endif

#if 0
    // EVENT-BY-EVENT
    globalPrefix =  globalPrefix0+"/other/";
    drawCompPVal( "pvala_comb_obs", "pvala_combe_obs", "Default", "Ev-by-ev.", "Local p-value obs."),
#endif

}

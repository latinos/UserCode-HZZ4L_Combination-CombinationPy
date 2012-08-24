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

void cccPlot(TString who, float where, const char **mychann, int nchann, TString srcPostfix="", TString postfix="", double rMin=-2.0, double rMax=5.0, bool individuals=false, bool comb=true) {
    TGraphAsymmErrors *obs = (TGraphAsymmErrors *) gFile->Get(who+"_"+mychann[0]+srcPostfix+"_obs");
    if (obs == 0) return;
    int j = findBin(obs, where); if (j == -1) return;
    double r0 = obs->GetY()[j], r0min = r0 - obs->GetErrorYlow(j), r0max = r0 + obs->GetErrorYhigh(j);
    bool halfint = (where - floor(where) > 0.3);
    c1->SetLogx(0); c1->SetLogy(0); c1->SetTicky(0);  c1->SetTickx(0);
    TGraphAsymmErrors *points = new TGraphAsymmErrors(); int npoints = 0;
    TString names[99];
    for (int i = nchann-1; i > 0; --i) {
       TGraphAsymmErrors *obsi = (TGraphAsymmErrors *) gFile->Get(who+"_"+mychann[i]+srcPostfix+"_obs"); if (obsi == 0) continue;
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
    int morePoints = 2;
    if (npoints >  7) morePoints += 2;
    TH2F frame("frame",";Best fit #sigma/#sigma_{SM};",1, rMin, rMax, npoints+morePoints, 0., npoints+morePoints);
    for (int k = 1; k <= npoints; ++k) {
       frame.GetYaxis()->SetBinLabel(k, names[k-1]);
    }
    frame.GetYaxis()->SetTickLength(0);

    double leftMargin = c1->GetLeftMargin(); c1->SetLeftMargin(individuals ? 0.30 : 0.24);
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
    if (comb) globalFitBand.DrawClone();
    TLine globalFitLine(r0, 0, r0, npoints);
    globalFitLine.SetLineWidth(4);     fake.SetLineWidth(0);
    //globalFitLine.SetLineColor(214); fake.SetLineColor(214); 
    globalFitLine.SetLineColor(1); fake.SetLineColor(colorFit68); 
    if (comb) globalFitLine.DrawClone();
    globalFitLine.SetLineStyle(2);  
    points->Draw("P SAME");
    frame.Draw("AXIS SAME");
    leg = newLegend(.31,.88,.55,.93); leg->SetTextSize(0.037);
    leg->SetHeader(Form(" m_{H} = %.0f %s", where, massUnits.Data()));
    spam(SPAM, .57, .88, .93, .93);
    leg->Draw();
    justSave(Form("%s_ccc_mH%.1f%s", who.Data(), where, postfix.Data()));
    c1->SetLeftMargin(leftMargin); c1->SetTicky(1); c1->SetTickx(1);
    leg = 0;
}

void cccPlot2(TString who1, TString who2, float where, const char **mychann, int nchann, TString srcPostfix1, TString srcPostfix2, TString label1, TString label2, TString postfix="", double rMin=-2.0, double rMax=5.0, bool individuals=false, bool comb=true) {
    c1->SetLogx(0); c1->SetLogy(0); c1->SetTicky(0);  c1->SetTickx(0);
    TGraphAsymmErrors *points1 = new TGraphAsymmErrors(); int npoints1 = 0;
    TGraphAsymmErrors *points2 = new TGraphAsymmErrors(); int npoints2 = 0;
    TString names[99]; int npoints = 0;
    for (int i = nchann-1, j = -1; i >= 0; --i) {
       TGraphAsymmErrors *obsi1 = (TGraphAsymmErrors *) gFile->Get(who1+"_"+mychann[i]+srcPostfix1+"_obs"); 
       TGraphAsymmErrors *obsi2 = (TGraphAsymmErrors *) gFile->Get(who2+"_"+mychann[i]+srcPostfix2+"_obs"); 
       if (obsi1 == 0 && obsi2 == 0) { std::cerr << "Missing " << mychann[i] << std::endl; continue; }
       names[npoints++] = (i == 0 ? TString("Combined") : channelFromName(mychann[i],true,true));
       if (obsi1 && (j = findBin(obsi1, where)) != -1) { 
           points1->Set(++npoints1);
           points1->SetPoint(npoints1-1, obsi1->GetY()[j], npoints - 0.35);
           points1->SetPointError(npoints1-1, obsi1->GetErrorYlow(j), obsi1->GetErrorYhigh(j), 0., 0.);
       }
       if (obsi2 && (j = findBin(obsi2, where)) != -1) { 
           points2->Set(++npoints2);
           points2->SetPoint(npoints2-1, obsi2->GetY()[j], npoints - 0.65);
           points2->SetPointError(npoints2-1, obsi2->GetErrorYlow(j), obsi2->GetErrorYhigh(j), 0., 0.);
       }
    }
    int morePoints = 2;
    if (npoints <= 4) morePoints -= 1;
    if (npoints >  7) morePoints += 2;
    TH2F frame("frame",";Best fit #sigma/#sigma_{SM};",1, rMin, rMax, npoints+morePoints, 0., npoints+morePoints);
    for (int k = 1; k <= npoints; ++k) {
       frame.GetYaxis()->SetBinLabel(k, names[k-1]);
    }
    frame.GetYaxis()->SetTickLength(0);

    double leftMargin = c1->GetLeftMargin(); c1->SetLeftMargin(individuals ? 0.30 : 0.24);
    points1->SetLineColor(kRed);
    points1->SetMarkerColor(205);
    points1->SetLineWidth(4);
    points1->SetMarkerStyle(21);
    points1->SetMarkerSize(1.4);
    points2->SetLineColor(kBlue);
    points2->SetMarkerColor(213);
    points2->SetLineWidth(4);
    points2->SetMarkerStyle(20);
    points2->SetMarkerSize(1.4);
    frame.GetXaxis()->SetTitleSize(0.05);
    frame.GetXaxis()->SetLabelSize(0.04);
    frame.GetYaxis()->SetLabelSize(0.06);
    frame.Draw(); gStyle->SetOptStat(0);
    TLine line(1, 0, 1, npoints); line.SetLineStyle(7); 
    line.DrawLine(0,0,0,npoints);
    line.DrawLine(1,0,1,npoints);
    points1->Draw("P SAME");
    points2->Draw("P SAME");
    frame.Draw("AXIS SAME");
    leg = newLegend(.31,.78,.55,.93); leg->SetTextSize(0.037);
    leg->SetHeader(Form(" m_{H} = %.0f %s", where, massUnits.Data()));
    leg->AddEntry(points1, label1, "LP");
    leg->AddEntry(points2, label2, "LP");
    leg->Draw();
    //leg = newLegend(.72,.15,.94,.25); leg->SetTextSize(0.037);
    //leg->Draw();
    spam(SPAM, .57, .88, .93, .93);
    justSave(Form("%s_ccc_mH%.1f%s", who1.Data(), where, postfix.Data()));
    c1->SetLeftMargin(leftMargin); c1->SetTicky(1); c1->SetTickx(1);
    leg = 0;
}

double muHatFromTGraphScan(TGraph *graph, bool interpolate=false) {
    int iMin = 0, n = graph->GetN(); double *yi = graph->GetY();
    for (int i = 1; i < n; ++i) {
        if (yi[i] < yi[iMin]) iMin = i;
    }
    if (!interpolate || iMin == 0 || iMin == n-1) return graph->GetX()[iMin];
    TSpline3 spline("dummy",graph);
    double x0 = graph->GetX()[iMin-1], x2 = graph->GetX()[iMin+1];
    double xmin = x0, ymin = graph->GetY()[iMin-1];
    for (double x = x0, dx = (x2-x0)*0.02; x < x2; x += dx) {
        double y = spline.Eval(x);
        if (y < ymin) { ymin = y; xmin = x; }
    }
    return x;
}
void makeMLZFromScan(TDirectory *bands, TString outName, float outMass, TString scanName, TString channel, bool interpolate=false) {
    TGraph *scan = gFile->Get(scanName+"_"+channel+"_obs"); if (scan == 0) return;
    TGraphAsymmErrors *ret  = new TGraphAsymmErrors(1);
    double muHat = muHatFromTGraphScan(scan, interpolate);
    double muLow = findCrossingOfScan1D(scan, 1.0, true), muHigh = findCrossingOfScan1D(scan, 1.0, false);
    ret->SetPoint(0, outMass, muHat);
    ret->SetPointError(0, 0., 0., muHat-muLow, muHigh-muHat);
    ret->SetName(outName+"_"+channel+"_obs");
    bands->WriteTObject(ret);
}
void makeMLZFromScans(TDirectory *bands, TString outName, float outMass, TString scanName, const char **mychann, int nchann, TString postfix="", bool interpolate=false) {
    for (int i = 0; i < nchann; ++i) {
        makeMLZFromScan(bands,outName,outMass,scanName,chann[i]+postfix, interpolate); 
    } 
}


void doCompatibilityChi2(TString who, float where, const char **mychann, int nchann, TString postfix="", bool nminus1) {
    // get reference mu or best fit
    double muhat = 1;
    if (nminus1) {
        TGraph *comb = (TGraph*) gFile->Get(who+"_"+mychann[0]+postfix+"_obs");
        if (comb == 0) return;
        muhat = muHatFromTGraphScan(comb);
        double muLow = findCrossingOfScan1D(comb, 1.0, true), muHigh = findCrossingOfScan1D(comb, 1.0, false);
        printf("Will use  %+.3f   %+.3f/%+.3f as best fit value\n", muhat, muLow-muhat, muHigh-muhat);
    } else {
        printf("Will compute chi-square with respect to the mu = 1 hypothesis\n");
    }
    // compute NLL
    int n = nminus1 ? -1 : 0;
    double chi2 = 0;
    for (int i = 1; i < nchann; ++i) {
        TGraph *ch = (TGraph*) gFile->Get(who+"_"+mychann[i]+postfix+"_obs");
        if (ch == 0 || ch->GetN() == 0) { std::cout << "   missing " << who+"_"+mychann[i]+"_obs" << std::endl; continue; }
        double myhat = muHatFromTGraphScan(ch);
        double muLow = findCrossingOfScan1D(ch, 1.0, true), muHigh = findCrossingOfScan1D(ch, 1.0, false);
        double mychi2 = ch->Eval(muhat);
        printf("    channel %-10s with muhat  %+.3f   %+.3f/%+.3f  contributes by %.3f to the global chi2.\n",
                    mychann[i], myhat, muLow-myhat, muHigh-myhat, mychi2);
        chi2 += mychi2; n++;
    }
    std::cout << "Global chi2 " << chi2 << "   ndf " << n << "  probability: " << TMath::Prob(chi2,n) << std::endl;
}

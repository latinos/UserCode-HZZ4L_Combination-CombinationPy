double findCrossingOfScan1D(TGraph *graph, double threshold, bool leftSide) {
    double *x = graph->GetX();
    double *y = graph->GetY();
    int imin = 0, n = graph->GetN();
    for (int i = 1; i < n; ++i) {
        if (y[i] < y[imin]) imin = i;
    }
    int imatch = -1;
    if (leftSide) {
        for (int i = 0; i < imin; ++i) {
            if (y[i] > threshold && y[i+1] < threshold) {
                imatch = i; break;
            }
        }
        if (imatch == -1) return x[0];
    } else {
        for (int i = imin; i < n; ++i) {
            if (y[i-1] < threshold && y[i] > threshold) {
                imatch = i-1; break;
            }
        }
        if (imatch == -1) return x[n-1];
    }
    double d1 = fabs(y[imatch] - threshold), d2 = fabs(y[imatch+1] - threshold);
    return (x[imatch]*d2 + x[imatch+1]*d1)/(d1+d2);
}

void drawOneScan1D(TString who, TString what="auto", double xmin0=-1, double xmax0=-1, bool drawFast=false) { 
    TGraph *obs = (TGraph *) gFile->Get(who+"_obs"), *fast = 0;
    if (obs == 0) return;
    if (drawFast) {
        TString fastwho = who; 
        if (fastwho.Contains("scan_1d")) fastwho("scan_1d") = "scan_1d_fast";
        else fastwho("scan") = "scan_fast";
        fast = (TGraph*) gFile->Get(fastwho+"_obs");
    }
    double xmin = obs->GetX()[0];
    double xmax = obs->GetX()[obs->GetN()-1];
    if (xmin0 < xmax0) { xmin = xmin0; xmax = xmax0; } 
    int col = 1;
    obs->SetLineWidth(4);
    obs->SetLineColor(col); obs->SetMarkerColor(col); 
    if (fast) {
        fast->SetLineWidth(4);   fast->SetLineWidth(2); fast->SetLineStyle(7);
        fast->SetLineColor(col); fast->SetMarkerColor(col); 
        leg = newLegend(.66,.78,.93,.93); leg->SetTextSize(0.04);
        leg->AddEntry(obs,  "with syst.", "L");
        leg->AddEntry(fast, "no syst.",   "L");
    }
    TH1D frame0("frame","frame", 100, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    double ymin = 0, ymax = 10; 
    if (fast) fast->Draw("CX");
    obs->Draw("CX");    
    if (what == "auto") {
        if (who.Contains("mass_scan_1d")) what = "Higgs Boson Mass ("+massUnits+")";
    }
    setCanvas(&frame0, "", ymin, ymax, "- 2 #Delta ln L");
    TLine chi2(xmin, 3.84, xmax, 3.84);
    chi2.SetLineColor(lineAt1Color); chi2.SetLineStyle(lineAt1Style); chi2.SetLineWidth(2);
    chi2.DrawClone();
    chi2.SetLineWidth(2);
    double hi68 = findCrossingOfScan1D(obs, 1.00, false);
    double lo68 = findCrossingOfScan1D(obs, 1.00, true);
    double hi95 = findCrossingOfScan1D(obs, 3.84, false);
    double lo95 = findCrossingOfScan1D(obs, 3.84, true);
    chi2.DrawLine(lo95, 0, lo95, 3.84);    chi2.DrawLine(lo68, 0, lo68, 1.00);
    chi2.DrawLine(hi95, 0, hi95, 3.84);    chi2.DrawLine(hi68, 0, hi68, 1.00);
    frame0.GetXaxis()->SetTitle(what);
    frame0.GetXaxis()->SetNdivisions(1005);
    TString myspam = SPAM; if (!who.Contains("comb")) myspam("\n") = "\n" + channelFromName(who) + "\n";
    finalize(who+(fast?"_comp":""),xmin,xmax,ymin,ymax, myspam);
    leg = 0;
}

void drawCombScan1D(TString who, const char **chann, int nchann, TString srcPostfix="", TString postfix="", TString what="auto", double xmin0=-1, double xmax0=-1) { 
    TGraph *obs = (TGraph*) gFile->Get(who+"_"+chann[0]+srcPostfix+"_obs");
    if (obs == 0) return;
    double xmin = obs->GetX()[0]             - obs->GetErrorXlow(0);
    double xmax = obs->GetX()[obs->GetN()-1] + obs->GetErrorXhigh(obs->GetN()-1);
    if (xmin0 < xmax0) { xmin = xmin0; xmax = xmax0; } 
    TGraph *obsi[10];
    leg = newLegend(.66,.75,.95-0.02*isSquareCanvas,.94); leg->SetTextSize(0.035);
    for (int i = 0; i < nchann; ++i) {
        obsi[i] = (TGraph*) gFile->Get(who+"_"+chann[i]+srcPostfix+"_obs");
        if (obsi[i] == 0) continue;
        int col = (i == 0 ? 1 : colorFromName(chann[i]));
        obsi[i]->SetLineWidth(i == 0 ? 4 : 3);
        obsi[i]->SetLineColor(col); 
        TString label =  "Combined"; if (i) label = channelFromName(chann[i],/*withspaces=*/false,/*nolumi=*/true);
        leg->AddEntry(obsi[i], label, "L");
    }
    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    double ymin = 0, ymax = 10; 
    for (int i = nchann-1; i >= 0; --i) { if (obsi[i]) obsi[i]->Draw("LX"); } 
    if (what == "auto") {
        if (who.Contains("mass_scan_1d")) what = "Higgs Boson Mass ("+massUnits+")";
    }
    leg->Draw();
    setCanvas(&frame0, "", ymin, ymax, "- 2 #Delta ln L");
    TLine chi2(xmin, 3.84, xmax, 3.84);
    chi2.SetLineColor(lineAt1Color); chi2.SetLineStyle(lineAt1Style); chi2.SetLineWidth(2);
    chi2.DrawClone();
    chi2.SetLineWidth(2);
    frame0.GetXaxis()->SetTitle(what);
    frame0.GetXaxis()->SetNdivisions(1005);
    finalize(who+"_all"+postfix,xmin,xmax,ymin,ymax, SPAM);
    leg = 0;
}

void drawCompScan1D(TString who, TString who2, TString label, TString oldlabel, TString what="auto", double xmin0=-1, double xmax0=-1) { 
    TGraph *obs = (TGraph *) gFile->Get(who);
    TGraph *old = (TGraph *) gFile->Get(who2);
    if (obs == 0 || old == 0) return;
    double xmin = obs->GetX()[0]             - obs->GetErrorXlow(0);
    double xmax = obs->GetX()[obs->GetN()-1] + obs->GetErrorXhigh(obs->GetN()-1);
    if (xmin0 < xmax0) { xmin = xmin0; xmax = xmax0; } 
    int col = 2, dcol = 4;
    obs->SetLineWidth(2);
    obs->SetLineColor(col); obs->SetMarkerColor(col); 
    obs->SetMarkerStyle(21); obs->SetMarkerSize(0.8);
    old->SetLineWidth(4);
    old->SetLineColor(dcol); old->SetMarkerColor(dcol); 
    old->SetMarkerStyle(25); old->SetMarkerSize(1.2);
    leg = newLegend(.65,.80,.92,.92); leg->SetTextSize(0.04);
    leg->AddEntry(obs,    label, "LP");
    leg->AddEntry(old, oldlabel, "LP");
    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    leg->Draw();
    double ymin = 0, ymax = 10; 
    if (what == "auto") {
        if (who.Contains("mass_scan_1d")) what = "Higgs Boson Mass ("+massUnits+")";
    }
    setCanvas(&frame0, "", ymin, ymax, "- 2 #Delta ln L");
    TLine chi2(xmin, 3.84, xmax, 3.84);
    chi2.SetLineColor(lineAt1Color); chi2.SetLineStyle(lineAt1Style); chi2.SetLineWidth(2);
    chi2.DrawClone();
    chi2.SetLineWidth(2);
    double hi68 = findCrossingOfScan1D(obs, 1.00, false);
    double lo68 = findCrossingOfScan1D(obs, 1.00, true);
    double hi95 = findCrossingOfScan1D(obs, 3.84, false);
    double lo95 = findCrossingOfScan1D(obs, 3.84, true);
    chi2.DrawLine(lo95, 0, lo95, 3.84);    chi2.DrawLine(lo68, 0, lo68, 1.00);
    chi2.DrawLine(hi95, 0, hi95, 3.84);    chi2.DrawLine(hi68, 0, hi68, 1.00);
    chi2.SetLineColor(dcol);
    hi68 = findCrossingOfScan1D(old, 1.00, false);
    lo68 = findCrossingOfScan1D(old, 1.00, true);
    hi95 = findCrossingOfScan1D(old, 3.84, false);
    lo95 = findCrossingOfScan1D(old, 3.84, true);
    chi2.DrawLine(lo95, 0, lo95, 3.84);    chi2.DrawLine(lo68, 0, lo68, 1.00);
    chi2.DrawLine(hi95, 0, hi95, 3.84);    chi2.DrawLine(hi68, 0, hi68, 1.00);
    frame0.GetXaxis()->SetTitle(what);
    frame0.GetXaxis()->SetNdivisions(1005);
    old->Draw("LP"); obs->Draw("LP");    
    TString myspam = SPAM; if (!who.Contains("comb")) myspam("\n") = "\n" + channelFromName(who) + "\n";
    finalize(who+"_vs_"+who2,xmin,xmax,ymin,ymax,myspam);
    leg = 0;
}

void styleMultiGraph(TList *tmg, int lineColor, int lineWidth, int lineStyle) {
    for (int i = 0; i < tmg->GetSize(); ++i) {
        TGraph *g = (TGraph*) tmg->At(i);
        g->SetLineColor(lineColor); g->SetLineWidth(lineWidth); g->SetLineStyle(lineStyle);
    }
}

void drawOneScan2D(TString who, TString xwhat="auto", TString ywhat="auto", bool hist2d = true, bool contours = true, bool bestfit = true, double xmin0=-1, double xmax0=-1, double ymin0=-1, double ymax0=-1, bool isChi2=true, TString postfix="") { 
    TH2    *obs_2d   = (TH2    *) gFile->Get(who+"_th2");
    TList *obs_c68  = (TList *) gFile->Get(who+"_c68");
    TList *obs_c95  = (TList *) gFile->Get(who+"_c95");
    TGraph *obs_best = (TGraph *) gFile->Get(who+"_best");
    if (hist2d   && obs_2d   == 0) return;
    double rm = c1->GetRightMargin();
    double lc1 = lineAt1Color; lineAt1Color = 1;
    if (isChi2) { c1->SetRightMargin(0.24); }
    //if (contours && obs_c68  == 0) return;
    //if (contours && obs_c95  == 0) return;
    //if (bestfit  && obs_best == 0) return;
    double xmin = xmin0, xmax = xmax0, ymin = ymin0, ymax = ymax0;
    if (hist2d && obs_2d != 0) {
        obs_2d->GetXaxis()->SetRangeUser(xmin, xmax);
        obs_2d->GetYaxis()->SetRangeUser(ymin, ymax);
        /*xmin = obs_2d->GetXaxis()->GetXmin();
        xmax = obs_2d->GetXaxis()->GetXmax();
        ymin = obs_2d->GetYaxis()->GetXmin();
        ymax = obs_2d->GetYaxis()->GetXmax();*/
    } else {
        obs_2d = new TH2F("frame","frame",1,xmin,xmax,1,ymin,ymax);
    }
    if (xwhat == "auto") {
        if (who.Contains("mass_scan_2d") || who.Contains("mass_bayes")) { xwhat = "Higgs Boson Mass ("+massUnits+")"; ywhat = "#sigma/#sigma_{SM}"; }
    }
    c1->SetLogy(0);
    if (isChi2) {
        obs_2d->GetZaxis()->SetRangeUser(0, 20.);
        obs_2d->GetZaxis()->SetTitle("-2 #Delta ln L");
    }
    obs_2d->GetXaxis()->SetTitle(xwhat);
    obs_2d->GetYaxis()->SetTitle(ywhat);
    obs_2d->GetXaxis()->SetDecimals(1);
    obs_2d->GetYaxis()->SetDecimals(1);

    obs_2d->Draw(hist2d ? "COLZ" : ""); gStyle->SetOptStat(0);
    obs_2d->SetContour(100);
    obs_2d->GetXaxis()->SetNdivisions(1005);
    if (contours && obs_c68) {
        styleMultiGraph(obs_c68, /*color=*/1, /*width=*/3, /*style=*/1);
        obs_c68->Draw("L SAME");
    }
    if (contours && obs_c95) {
        styleMultiGraph(obs_c95, /*color=*/1, /*width=*/3, /*style=*/9);
        obs_c95->Draw("L SAME");
    }
    if (bestfit && obs_best) {
        obs_best->SetMarkerStyle(34); obs_best->SetMarkerSize(2.0);
        obs_best->Draw("P SAME");
    }
    TString myspam = SPAM; if (!TString(who).Contains("comb")) myspam("\n") = "\n" + channelFromName(who) + "\n";
    finalize(who+postfix,xmin,xmax,ymin,ymax, myspam);
    if (isChi2) { c1->SetRightMargin(rm); }
    lineAt1Color = lc1;
}


void drawCombScan2D(TString who, const char **chann, int nchann, TString srcPostfix="", TString postfix="", TString xwhat, TString ywhat, bool hist2d = true, bool c68 = true, bool c95 = true, bool bestfit = true, double xmin0=-1, double xmax0=-1, double ymin0=-1, double ymax0=-1, bool isChi2=true) {
    TH2    *obs_2d   = (TH2    *) gFile->Get(who+"_"+chann[0]+srcPostfix+"_th2");
    if (hist2d && obs_2d   == 0) { std::cerr << "Could not find " << who+chann[0]+srcPostfix+"_th2" << std::endl; return; }
    if (c68 && gFile->Get(who+"_"+chann[0]+srcPostfix+"_c68") == 0) return;
    if (bestfit && gFile->Get(who+"_"+chann[0]+srcPostfix+"_best") == 0) return;
    double rm = c1->GetRightMargin();
    double lc1 = lineAt1Color; lineAt1Color = 1;
    if (false && isChi2 && hist2d) { c1->SetRightMargin(0.24); }
    double xmin = xmin0, xmax = xmax0, ymin = ymin0, ymax = ymax0;
    if (hist2d && obs_2d != 0) {
        obs_2d->GetXaxis()->SetRangeUser(xmin, xmax);
        obs_2d->GetYaxis()->SetRangeUser(ymin, ymax);
    } else {
        obs_2d = new TH2F("frame","frame",1,xmin,xmax,1,ymin,ymax);
    }
    if (xwhat == "auto") {
        if (who.Contains("mass_scan_2d") || who.Contains("mass_bayes")) { xwhat = "Higgs Boson Mass ("+massUnits+")"; ywhat = "#sigma/#sigma_{SM}"; }
    }
    c1->SetLogy(0);
    if (isChi2) {
        obs_2d->GetZaxis()->SetRangeUser(0, 20.);
        obs_2d->GetZaxis()->SetTitle("-2 #Delta ln L");
    }
    obs_2d->GetXaxis()->SetTitle(xwhat);
    obs_2d->GetYaxis()->SetTitle(ywhat);
    obs_2d->GetXaxis()->SetDecimals(1);
    obs_2d->GetYaxis()->SetDecimals(1);

    obs_2d->Draw(hist2d ? "COLZ" : ""); gStyle->SetOptStat(0);
    obs_2d->SetContour(100);
    obs_2d->GetXaxis()->SetNdivisions( 510);
    leg = newLegend(.66,.75,.95-0.02*isSquareCanvas,.94); leg->SetTextSize(0.035);
    TList *c68i[99], *c95i[99]; TGraph *besti[99];
    for (int i = 0; i < nchann; ++i) {
        c68i[i] = (TList*) gFile->Get(who+"_"+chann[i]+srcPostfix+"_c68");
        c95i[i] = (TList*) gFile->Get(who+"_"+chann[i]+srcPostfix+"_c95");
        besti[i] = (TGraph*) gFile->Get(who+"_"+chann[i]+srcPostfix+"_best");
        if (c68 && c68i[i] == 0) continue;
        if (c95 && c95i[i] == 0) continue;
        if (bestfit && besti[i] == 0) continue;
        int col = (i == 0 ? 1 : colorFromName(chann[i]));
        int lw  = (i == 0 ? 4 : 2+hist2d);
        if (c68i[i]) styleMultiGraph(c68i[i], col, lw,   1);
        if (c95i[i]) styleMultiGraph(c95i[i], col, lw-1, 7);
        if (besti[i]) { besti[i]->SetLineWidth(lw); besti[i]->SetLineColor(col); besti[i]->SetMarkerColor(col); besti[i]->SetMarkerStyle(34); besti[i]->SetMarkerSize(2.0); }
        TString label =  "Combined"; if (i) label = channelFromName(chann[i],/*withspaces=*/false,/*nolumi=*/true);
        if (bestfit) leg->AddEntry(besti[i], label, "P");
        else if (c68) leg->AddEntry(c68i[i]->At(0), label, "L");
    }
    for (int i = nchann-1; i >= 0; --i) {
        if (c68i[i] && c68) c68i[i]->Draw("C SAME");
        if (c95i[i] && c95) c95i[i]->Draw("C SAME");
        if (bestfit && besti[i]) besti[i]->Draw("P SAME");
    }
    TString myspam = SPAM; if (!TString(chann[0]).Contains("comb")) myspam("\n") = "\n" + channelFromName(chann[0]) + "\n";
    finalize(who+"_all"+postfix,xmin,xmax,ymin,ymax, myspam);
    if (isChi2) { c1->SetRightMargin(rm); }
    if (obs_2d->GetName() == TString("frame")) delete obs_2d;
    lineAt1Color = lc1;
    leg = 0;
}



void drawCompContour2D(TString who, TString who2, TString label, TString oldlabel, double xmin, double xmax, double ymin, double ymax, TString what="auto") { 
    TList *obs = (TList *) gFile->Get(who);
    TList *old = (TList *) gFile->Get(who2);
    if (obs == 0 || old == 0) return;
    int col = 2, dcol = 4;
    styleMultiGraph(obs, /*color=*/ col, /*width=*/2, /*style=*/1);
    styleMultiGraph(old, /*color=*/dcol, /*width=*/4, /*style=*/1);
    leg = newLegend(.65,.80,.92,.92); leg->SetTextSize(0.04);
    leg->AddEntry(obs->At(0),    label, "L");
    leg->AddEntry(old->At(0), oldlabel, "L");
    obs_2d = new TH2F("frame","frame",1,xmin,xmax,1,ymin,ymax);
    TString xwhat, ywhat;
    if (what == "auto") {
        if (who.Contains("mass_scan_2d") || who.Contains("mass_bayes")) { xwhat = "Higgs Boson Mass ("+massUnits+")"; ywhat = "#sigma/#sigma_{SM}"; }
    }
    c1->SetLogy(0);
    obs_2d->GetXaxis()->SetTitle(xwhat);
    obs_2d->GetYaxis()->SetTitle(ywhat);
    obs_2d->Draw(""); gStyle->SetOptStat(0);
    obs_2d->GetXaxis()->SetNdivisions(1005);
    leg->Draw();
    old->Draw("C SAME");
    obs->Draw("C SAME");
    TString myspam = SPAM; if (!TString(who).Contains("comb")) myspam("\n") = "\n" + channelFromName(who) + "\n";
    finalize(who+"_vs_"+who2,xmin,xmax,ymin,ymax,myspam);
    delete obs_2d;
    leg = 0;
}



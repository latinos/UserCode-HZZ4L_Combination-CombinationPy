bool Draw_TEV=false;
bool Draw_CLs_Band = true, Draw_CLs_ExpLine = false, Draw_CLs_Obs = true;

void drawOneCLs(TString who) {
    TGraphAsymmErrors *obs   = (TGraphAsymmErrors *) gFile->Get(who+"_obs");
    TGraphAsymmErrors *apriori = CLs_debug_apriori_grid ? makeAPrioriGrid(who) : 0;
    TGraphAsymmErrors *miss = obs ? missingPoints(obs) : 0;
    TGraphAsymmErrors *m68 = (TGraphAsymmErrors *) gFile->Get(who+"_median");
    TGraphAsymmErrors *m95 = (TGraphAsymmErrors *) gFile->Get(who+"_median_95"); 
    TGraphAsymmErrors *miss68 = m68 ? missingPoints(m68) : 0;
    if (miss68) { miss68->SetMarkerStyle(24); miss68->SetMarkerSize(0.9); miss68->SetMarkerColor(4); }
    if (obs == 0 && m68 == 0) return;

    TGraphAsymmErrors *obsTEV = Draw_TEV ? (TGraphAsymmErrors *) gFile->Get("tevatron_obs")    : 0;
    TGraphAsymmErrors *expTEV = Draw_TEV ? (TGraphAsymmErrors *) gFile->Get("tevatron_median") : 0;
    if (Draw_TEV) {
        obsTEV->SetLineColor(kBlue); obsTEV->SetLineWidth(5); obsTEV->SetLineStyle(1);
        expTEV->SetLineColor(kBlue); expTEV->SetLineWidth(5); expTEV->SetLineStyle(2);
    }

    TGraphAsymmErrors *ref = (obs ? obs : m68);
    double xmin = ref->GetX()[0]             - ref->GetErrorXlow(0), xmin0 = xmin;
    double xmax = ref->GetX()[ref->GetN()-1] + ref->GetErrorXhigh(ref->GetN()-1);
    if (xmin <= 120 && xmax > 200) { xmin = 99.98; xmax = (SM == "FP" ? 300.1 : 600.1); }
    //int col = 1, dcol = 1, smooth = 5, smoothorder = 0;
    int col = 1, dcol = 1, smooth = !who.Contains("acls"), smoothorder = 0;
    if (obs) {
        obs->SetLineWidth(2); 
        obs->SetLineColor(col); obs->SetMarkerColor(col); 
        obs->SetMarkerStyle(21); obs->SetMarkerSize(0.8);
    }
    double ymin, ymax;   minMaxY(ref, ymin, ymax);
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
    if (Draw_CLs_Band) draw2(who, color68, color95, color50, true, false, smooth, smoothorder);
    if (Draw_CLs_ExpLine) { 
        TGraphAsymmErrors *smoothExp = smoothSMCLs(m68, smooth, smoothorder);
        if (smoothExp != 0) m68 = smoothExp;
        m68->Draw("LX"); 
    }
    if (m68) { m68->SetFillColor(color68); m68->SetLineColor(color50); m68->SetLineStyle(2); m68->SetLineWidth(2); }
    if (m95) { m95->SetFillColor(color95); m95->SetLineColor(color50); m95->SetLineStyle(2); m95->SetLineWidth(2); }
    frame0.Draw("AXIGSAME");
    double leg_y_hi = 0.94;
    double leg_y_lo = (isSquareCanvas ? 0.78 : .75)- 0.05*(lep_excluded+tev_excluded+1.6*Draw_TEV);
    double leg_x_lo = isSquareCanvas ? .645-0.06*(Draw_TEV||tev_excluded) : .65;
    leg = newLegend(leg_x_lo, leg_y_lo,.93+0.005*isSquareCanvas, leg_y_hi); leg->SetTextSize(isSquareCanvas || lep_excluded?  0.034 : 0.037);
    if (obs && Draw_CLs_Obs) obs->Draw("LP");    
    if (miss68) miss68->Draw("P");
    if (miss) miss->Draw("P");
    if (apriori) apriori->Draw("LX");
    if (obs) leg->AddEntry(obs, "Observed", "LP");
    if (m68) { leg->AddEntry(m68, "Expected "+oneSigmaText, "LF");
               leg->AddEntry(m95, "Expected "+twoSigmaText, "LF"); }
    if (Draw_TEV) {
        leg->AddEntry(obsTEV, "Tevatron Observed", "L");
        leg->AddEntry(expTEV, "Tevatron Expected", "L");
    }
    if (lep_excluded) leg->AddEntry(fakeLEP, "LEP excluded", "F");
    if (tev_excluded) leg->AddEntry(fakeTEV, "Tevatron excluded", "F");
    if (cms_excluded) leg->AddEntry(fakeCMS, "CMS excluded", "F");
    if (Draw_TEV) { expTEV->Draw("L SAME"); obsTEV->Draw("L SAME"); }
    leg->Draw();
    TString myspam = SPAM; if (!who.Contains("comb")) myspam("\n") = "\n" + channelFromName(who) + "\n";
    if (!c1->GetGridx() && !c1->GetGridy()) frame0.Draw("AXIGSAME");
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
        if (Draw_CLs_Band) draw2(who, color68, color95, color50, true, false, smooth, smoothorder);
        if (Draw_CLs_ExpLine) m68->Draw("LX");
        frame2.Draw("AXIGSAME");
        if (obs && Draw_CLs_Obs) obs->Draw("LP");    
        bool noMoreTev = tev_excluded && (x_zoom <= 156);
        leg = newLegend(leg_x_lo+0.06*noMoreTev*isSquareCanvas, leg_y_lo+0.05*noMoreTev,.93+0.005*isSquareCanvas, leg_y_hi); 
        leg->SetTextSize(isSquareCanvas || lep_excluded?  0.034 : 0.037);
        if (obs) leg->AddEntry(obs, "Observed", "LP");
        if (m68) { leg->AddEntry(m68, "Expected "+oneSigmaText, "LF");
                   leg->AddEntry(m95, "Expected "+twoSigmaText, "LF"); }
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
        minMaxY(ref, ymin, ymax, x_zoom);
        if (m68) minMaxY(m68, ymin2, ymax2, x_zoom); 
        if (ymax2 > ymax) ymax = ymax2; if (ymin2 < ymin) ymin = ymin2;
        if (Draw_TEV) ymax *=2;
        if (!c1->GetGridx() && !c1->GetGridy()) frame2.Draw("AXIGSAME");
        finalize(who+"_zoom",xmin,x_zoom,ymin,ymax,myspam);
        if (doZoom2) {
            TH1D frame3("frame3","frame3", 1, x_zoom2_min,x_zoom2_max); frame3.Draw(); 
            setCanvas(&frame3, "", ymin, ymax, what);
            if (apriori) apriori->Draw("E3");
            draw2(who, color68, color95, color50, true, false, 5);
            frame3.Draw("AXIGSAME");
            if (obs) obs->Draw("LP");    
            if (miss68) miss68->Draw("P");
            if (miss) miss->Draw("P");
            if (apriori) apriori->Draw("LX");
            leg->Draw();
            finalize(who+"_zoom2",x_zoom2_min,x_zoom2_max,ymin,ymax,myspam);
        }
    }
    leg = 0;
}
void drawTwoCLs(TString who, TString who2, TString name2="Other", TString postfix="", int color2=kBlue, bool doObs2=true, bool noExpObsLabels=false) {
    TGraphAsymmErrors *obs   = (TGraphAsymmErrors *) gFile->Get(who+"_obs");
    TGraphAsymmErrors *m68 = (TGraphAsymmErrors *) gFile->Get(who+"_median");
    TGraphAsymmErrors *m95 = (TGraphAsymmErrors *) gFile->Get(who+"_median_95"); 
    if (obs == 0 && m68 == 0) return;

    TGraphAsymmErrors *obs2 = (TGraphAsymmErrors *) gFile->Get(who2+"_obs");
    TGraphAsymmErrors *exp2 = (TGraphAsymmErrors *) gFile->Get(who2+"_median");
    if (obs2 == 0 && exp2 == 0) { std::cerr << "Missing " << who2 << " to compare with " << who << std::endl; return; }
    if (obs2) { obs2->SetLineColor(color2); obs2->SetLineWidth(5); obs2->SetLineStyle(1); }
    if (exp2) { exp2->SetLineColor(color2); exp2->SetLineWidth(5); exp2->SetLineStyle(2); }

    TGraphAsymmErrors *ref = (obs ? obs : m68);
    double xmin = ref->GetX()[0]             - ref->GetErrorXlow(0), xmin0 = xmin;
    double xmax = ref->GetX()[ref->GetN()-1] + ref->GetErrorXhigh(ref->GetN()-1);
    if (xmin <= 120 && xmax > 200) { xmin = 99.98; xmax = (SM == "FP" ? 300.1 : 600.1); }
    int col = 1, dcol = 1;
    if (obs) {
        obs->SetLineWidth(2); 
        obs->SetLineColor(col); obs->SetMarkerColor(col); 
        obs->SetMarkerStyle(21); obs->SetMarkerSize(0.8);
    }
    double ymin, ymax;   minMaxY(ref, ymin, ymax); 
    double ymin2, ymax2; minMaxY(exp2 ? exp2 : obs2, ymin2, ymax2); 
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
    if (obs) obs->Draw("LP");    
    if (exp2) exp2->Draw("LX SAME"); 
    if (obs2 && doObs2) obs2->Draw("LX SAME");
    if (obs) leg->AddEntry(obs, "Observed", "LP");
    leg->AddEntry(m68, "Expected "+oneSigmaText, "LF");
    leg->AddEntry(m95, "Expected "+twoSigmaText, "LF");
    if (obs2 && doObs2) leg->AddEntry(obs2, name2+(noExpObsLabels ? "" : " Obs."), "L");
    if (exp2) leg->AddEntry(exp2, name2+(noExpObsLabels ? "" : " Exp."), "L");
    leg->Draw();
    TString myspam = SPAM; if (!who.Contains("comb")) myspam("\n") = "\n" + channelFromName(who) + "\n";
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
        if (obs) obs->Draw("LP");    
        if (exp2) exp2->Draw("LX SAME"); 
        if (obs2 && doObs2) obs2->Draw("LX SAME");
        leg->Draw();
        minMaxY(ref, ymin, ymax, x_zoom); minMaxY(exp2 ? exp2 : obs2, ymin, ymax, x_zoom);
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
    TString myspam = SPAM; if (!who.Contains("comb")) myspam("\n") = "\n" + channelFromName(who) + "\n";
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
    float legxstart = isSquareCanvas ? 0.60 : 0.72;
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
    TString myspam = SPAM; if (!TString(chann[0]).Contains("comb")) myspam("\n") = "\n" + channelFromName(chann[0]) + "\n";
    setCanvas(&frame0, "", ymin, ymax, what);
    leg->Draw();
    finalize(who+"_all"+postfix,xmin,xmax,ymin,ymax, myspam);
    if (xmax >= 200 && xmin < x_zoom) {
        c1->SetLogx(1);
        frame0.GetXaxis()->SetMoreLogLabels(); frame0.GetXaxis()->SetNoExponent();
        finalize(who+"_all"+postfix+"_logx",xmin,xmax,ymin,ymax);
        c1->SetLogx(0);
        xmin = xmin0;
        leg = newLegend(legxstart,.65+0.10*isC3+(ntot-nlow)*0.032,.95-0.02*isSquareCanvas,.94); leg->SetTextSize(0.032);
        TH1D frame2("frame2","frame2", 1, xmin, x_zoom); frame2.Draw(); 
        //minMaxY(obs, ymin, ymax, 200.); ymax *= 15;
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
        setCanvas(&frame2, "", ymin, ymax, what);
        finalize(who+"_all"+postfix+"_zoom",xmin,x_zoom,ymin,ymax, myspam);
    }
    leg = 0;
}

void drawCombBoth(TString who, TString what="P.L. Approx limit #sigma_{95%}/#sigma_{"+SM+"}", const char **chann, int nchann, TString srcPostfix="", bool observed=true, bool combined=true, TString postfix="") {
    TGraphAsymmErrors *obs = (TGraphAsymmErrors *) gFile->Get(Form("%s_%s%s_obs",who.Data(),chann[0],srcPostfix.Data()));
    if (obs == 0) return;
    double xmin = obs->GetX()[0]             - obs->GetErrorXlow(0), xmin0 = xmin;
    double xmax = obs->GetX()[obs->GetN()-1] + obs->GetErrorXhigh(obs->GetN()-1);
    bool isC3 = (postfix == "_c3" || nchann == 3);
    if (xmin <= 120 && xmax > 200) { xmin = 99.98; xmax = (SM == "FP" ? 300.1 : 600.1); }
    double ymin, ymax; minMaxY(obs, ymin, ymax); ymax *= (isC3?2:8);
    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    TGraphAsymmErrors *obsi[nchann_all], *expi[nchann_all];
    float legxstart = isSquareCanvas ? 0.65 : 0.72;
    float legystart = .65-0.05*combined+0.14*isC3+0.03*lessSpacesInLegends+0.06*(SM=="FP");
    leg = newLegend(legxstart, legystart, .95-0.02*isSquareCanvas, .94); leg->SetTextSize(0.032);
    if (!observed) leg->SetHeader("Expected limits");
    int ntot = 0, nlow = 0;
    for (int i = 0; i < nchann; ++i) {
        if (i == 0 && !combined) continue;
        obsi[i] = (TGraphAsymmErrors *) gFile->Get(who+"_"+chann[i]+srcPostfix+"_obs");
        expi[i] = (TGraphAsymmErrors *) gFile->Get(who+"_"+chann[i]+srcPostfix+"_median");
        if (expi[i] == 0 || obsi[i] == 0) continue;
        /*if (TString(chann[i]).Contains("hzz2l2q") || 
            (noLineStyles == true && TString(chann[i]).Contains("comb")) ||
            (noLineStyles == true && TString(chann[i]) == "hzz")) {
            TGraphAsymmErrors *smooth = smoothSMCLs(expi[i], 7, 2);
            if (smooth == 0) std::cout << "Smoothing of " << expi[i]->GetName() << " returned ZERO" << std::endl;
            else expi[i] = smooth;
        }*/
        int col = (i == 0 ? 1 : colorFromName(chann[i])); //if (col == 96) col = 67;
        obsi[i]->SetLineWidth(i == 0 ? 4 : 3);
        obsi[i]->SetLineColor(col);  obsi[i]->SetMarkerColor(col); 
        obsi[i]->SetMarkerStyle(21); obsi[i]->SetMarkerSize(0.8);
        expi[i]->SetLineWidth(i == 0 ? 4 : (observed ? 3 : 4));
        expi[i]->SetLineColor(col);  expi[i]->SetMarkerColor(col); 
        expi[i]->SetMarkerStyle(21); expi[i]->SetMarkerSize(0.8);
        obsi[i]->SetLineStyle(observed ||  noLineStyles ? 1 : 2);
        expi[i]->SetLineStyle(observed || !noLineStyles ? 2 : 1);
        if (i == 0 && expi[i] != 0) {
            TString whoAmI = "Combined"; //channelFromName(chann[i],/*withspaces=*/true);
            if (observed) {
                leg->AddEntry(obsi[i], whoAmI+(lessSpacesInLegends?" obs.":" observed"), who.Contains("cls") ? "L" : "LP");
                leg->AddEntry(expi[i], whoAmI+(lessSpacesInLegends?" exp.":" expected"), who.Contains("cls") ? "L" : "LP");
            } else {
                leg->AddEntry(expi[i], whoAmI, who.Contains("cls") ? "L" : "LP");
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
    TString myspam = SPAM; if (!TString(chann[0]).Contains("comb")) myspam("\n") = "\n" + channelFromName(chann[0]) + "\n";
    setCanvas(&frame0, "", ymin, ymax, what);
    leg->Draw();
    finalize(who+(observed?"_all2":"_allexp")+postfix,xmin,xmax,ymin,ymax, myspam);
    if (xmax >= 200 && xmin < x_zoom) {
        c1->SetLogx(1);
        frame0.GetXaxis()->SetMoreLogLabels(); frame0.GetXaxis()->SetNoExponent();
        finalize(who+(observed?"_all2"+postfix+"_logx":"_allexp"+postfix+"_logx"),xmin,xmax,ymin,ymax);
        c1->SetLogx(0);
        xmin = xmin0;
        leg = newLegend(legxstart,legystart+(ntot-nlow)*0.032,.95-0.02*isSquareCanvas,.94); leg->SetTextSize(0.032);
        if (!observed) leg->SetHeader("Expected limits");
        TH1D frame2("frame2","frame2", 1, xmin, x_zoom); frame2.Draw(); 
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
        setCanvas(&frame2, "", ymin, ymax, what);
        finalize(who+(observed?"_all2"+postfix+"_zoom":"_allexp"+postfix+"_zoom"),xmin,x_zoom,ymin,ymax, myspam);
    }
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
    old->SetLineWidth(4);
    old->SetLineColor(dcol); old->SetMarkerColor(dcol); 
    old->SetMarkerStyle(25); old->SetMarkerSize(1.2);
    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    double ymin, ymax; minMaxY(obs, ymin, ymax);
    double ymin2, ymax2; minMaxY(old, ymin2, ymax2);
    if (ymin2 < ymin) ymin = ymin2; if (ymax2 > ymax) ymax = ymax2;
    if (who1.Contains("mlz_")) { ymin = -2.5; ymax = 5; }
    if ((who1.Contains("median") && who2.Contains("median")) ||
        (who1.Contains("mlz") && who2.Contains("mlz"))) {
        old->SetFillStyle(1001); old->SetFillColor(216);
        obs->SetFillStyle(3444); obs->SetFillColor(2); obs->SetLineColor(206);
        old->Draw("E3"); obs->Draw("E3");    
        old->Draw("LPX"); obs->Draw("LPX");    
    } else {
        old->Draw("LP"); obs->Draw("LP");    
    }
    leg = newLegend(.65,.80,.92,.92); leg->SetTextSize(0.04);
    leg->AddEntry(obs,    label, "LP");
    leg->AddEntry(old, oldlabel, "LP");
    leg->Draw();
    setCanvas(&frame0, "", ymin, ymax, "Limit");
    TString myspam = SPAM; if (!who1.Contains("comb")) myspam("\n") = "\n" + channelFromName(who1) + "\n";
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



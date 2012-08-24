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
    double ymin = 5e-4, ymax = 60;
    if (who.Contains("smacls")) ymin = 7e-8;
    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    TString what = "CL_{S} of "+SM+" Higgs hypothesis";
    setCanvas(&frame0, "", ymin, ymax, what);
    draw2(who, color68, color95, color50, true, false, smooth);
    frame0.Draw("AXIGSAME");
    obs->Draw("LP");
    if (as) as->Draw("LP");
    //spam("#splitline{"+SPAM+"}{"+channelFromName(who)+"}", 0.17,.15,.54,.20,32);
    TString myspam = SPAM; if (!who.Contains("comb")) myspam("\n") = "\n" + channelFromName(who) + "\n";
    //leg = newLegend(.64,.15,.88,.31+showAs*0.04); leg->SetTextSize(0.04);
    leg = newLegend(.67-0.11*isSquareCanvas,.93-.16-showAs*0.04,.91,.93); leg->SetTextSize(0.04);
    m68->SetFillColor(color68); m68->SetLineColor(color50); m68->SetLineStyle(2); m68->SetLineWidth(2);
    m95->SetFillColor(color95); m95->SetLineColor(color50); m95->SetLineStyle(2); m95->SetLineWidth(2);
    leg->AddEntry(obs, "Observed", "LP");
    leg->AddEntry(m68, "Expected "+oneSigmaText, "LF");
    leg->AddEntry(m95, "Expected "+twoSigmaText, "LF");
    if (showAs) leg->AddEntry(as, "Asymptotic Obs.", "LP");
    frame0.Draw("AXIS SAME"); // re-draw ticks
    finalize(who+(showAs?"_comp":""),xmin,xmax,ymin,ymax,myspam);
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
        frame2.Draw("AXIS SAME"); // ticks
        finalize(who+(showAs ? "_comp_zoom" : "_zoom"),xmin,x_zoom,ymin,ymax,myspam);
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
    TString myspam = SPAM; if (!who.Contains("comb")) myspam("\n") = "\n" + channelFromName(who) + "\n";
    spam(myspam, 0.17,.15,.54,.20,32);
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


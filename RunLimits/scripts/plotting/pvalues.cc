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
void drawOnePVal(TString who, TString what="Local p-value", TString drawToys="", TString toyPostfix="_comp") {
    TGraphAsymmErrors *obs = (TGraphAsymmErrors *) gFile->Get(who+"_obs");
    if (obs == 0) return;
    TGraphAsymmErrors *exp = (TGraphAsymmErrors *) gFile->Get(who+"_median");
    if (exp == 0) exp = (TGraphAsymmErrors *) gFile->Get(who+"_asimov");

    TGraphAsymmErrors *toys = 0, *toysExp = 0;
    if (drawToys != "") {
        if (drawToys == "auto" && who.Index("pvala") == 0) {
            drawToys = who; drawToys.ReplaceAll("pvala","pval");
        }
        toys = (TGraphAsymmErrors *) gFile->Get(drawToys+"_obs");
        if (toys == 0) return;
        std::cout << "Toys are " << toys->GetName() << std::endl;
        toys->SetLineWidth(4);
        toys->SetMarkerStyle(20);
        toys->SetMarkerSize(1.4);
        toys->SetLineColor(92);
        toys->SetMarkerColor(92);
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
        exp->SetLineWidth(2); exp->SetLineStyle(7); exp->SetLineColor(214); 
    }
    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    double ymin, ymax; minMaxY(obs, ymin, ymax, 999, 0); 
    ymax = 1.0; ymin = minPValue(obs,xmin,xmax)/1000; if (ymin <= 0 || ymin > 1e-5) ymin = 1e-5;
    
    setCanvas(&frame0, "", ymin, ymax, what);
    frame0.GetYaxis()->SetTitleOffset(1.10+0.25*isSquareCanvas);
    if (toys) toys->Draw("P");
    if (toysExp) toysExp->Draw("P");
    obs->Draw(toys ? "LX" : "LXP");    
    if (exp) exp->Draw("LX");
    if (miss) miss->Draw("P");
    TString myspam = SPAM; if (!who.Contains("comb")) myspam("\n") = "\n" + channelFromName(who) + "\n";
    if (toys) {
        leg = newLegend(.66-0.065*isSquareCanvas,.16,.93,.30+((exp!=0)+(toysExp!=0))*0.04); leg->SetTextSize(0.04-0.003*isSquareCanvas);
        leg->AddEntry(obs,  "Asymptotic Obs.", "L");
        leg->AddEntry(toys, "Ensemble Obs.",   "P");
        if (toysExp) leg->AddEntry(toysExp, "Ensemble Exp.",   "P");
        if (exp) leg->AddEntry(exp, "Asymptotic Exp.", "L");
    } else {
        leg = 0;
    }
    
    finalize(who+(toys?toyPostfix:""),xmin,xmax,ymin,ymax,myspam,true);
    if (xmax > x_zoom && xmin < x_zoom) {
        c1->SetLogx(1);
        frame0.GetXaxis()->SetMoreLogLabels(); frame0.GetXaxis()->SetNoExponent();
        finalize(who+(toys?toyPostfix+"_logx":"_logx"),xmin,xmax,ymin,ymax,myspam,true);
        c1->SetLogx(0);
        xmin = xmin0;
        TH1D frame2("frame2","frame2", 1, xmin, x_zoom); frame2.Draw(); 
        setCanvas(&frame2, "", ymin, ymax, what);
        if (toys) toys->Draw("P");
        if (toysExp) toysExp->Draw("P");
        obs->Draw(toys ? "LX" : "LXP");    
        if (exp) exp->Draw("LX");
        if (miss) miss->Draw("P");
        finalize(who+(toys?toyPostfix+"_zoom":"_zoom"),xmin,x_zoom,ymin,ymax,myspam,true);
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

    TString myspam = SPAM; if (!TString(chann[0]).Contains("comb")) myspam("\n") = "\n" + channelFromName(chann[0]) + "\n";
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

void drawCombPVal(TString who, const char **chann, int nchann, TString srcPostfix="", TString postfix="", TString toywho="", double xmin0=+1, double xmax0=-1) {
    TGraphAsymmErrors *obs = (TGraphAsymmErrors *) gFile->Get(who+"_"+chann[0]+srcPostfix+"_obs");
    if (obs == 0) return;
    TGraphAsymmErrors *exp = (TGraphAsymmErrors *) gFile->Get(who+"_"+chann[0]+srcPostfix+"_median");
    if (exp == 0) exp = (TGraphAsymmErrors *) gFile->Get(who+"_"+chann[0]+srcPostfix+"_asimov");

    TGraphAsymmErrors *toys = toywho != "" ? (TGraphAsymmErrors *) gFile->Get(toywho+srcPostfix+"_obs") : 0;
    if (toys) {
        toys->SetLineWidth(4);
        toys->SetMarkerStyle(20);
        toys->SetMarkerSize(toys->GetN() > 3 && SM != "FP"? 1.3 : 1.8);
        toys->SetLineColor(92);
        toys->SetMarkerColor(92);
    }


    double xmin = obs->GetX()[0]             - obs->GetErrorXlow(0), xmin2 = xmin;
    double xmax = obs->GetX()[obs->GetN()-1] + obs->GetErrorXhigh(obs->GetN()-1);
    bool isC3 = (postfix.Contains("_c3") || nchann == 3);
    if (xmin <= 120 && xmax > 200) { xmin = 99.98; xmax = (SM == "FP" ? 300.1 : 600.1); }
    if (TString(chann[0]).Contains("_toy_")) { xmin = x_zoom2_min; xmax = x_zoom2_max; }
    if (xmin0 < xmax0) { xmin = xmin0; xmax = xmax0; }
    double ymin, ymax; ymax = 1.0; ymin = 8e-12;
    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    setCanvas(&frame0, "", ymin, ymax, "Local p-value");
    TGraphAsymmErrors *obsi[nchann_all];
    float legxstart = isSquareCanvas ? 0.61 : 0.72;
    float legyend = .43-0.10*isC3; 
    if (toys) legyend += 0.04;
    TLegend *leg = newLegend(legxstart,.15,.95-0.02*isSquareCanvas,legyend); leg->SetTextSize(0.032);
    for (int i = 0; i < nchann; ++i) {
        obsi[i] = (TGraphAsymmErrors *) gFile->Get(who+"_"+chann[i]+srcPostfix+"_obs");
        if (obsi[i] == 0) continue;
        int col = colorFromName(chann[i]);
        //obsi[i]->SetLineWidth(2);
        obsi[i]->SetLineWidth(i == 0 ? 4 : 4);
        obsi[i]->SetLineColor(col); obsi[i]->SetMarkerColor(col); 
        obsi[i]->SetMarkerStyle(21); obsi[i]->SetMarkerSize(0.8);
        obsi[i]->SetLineStyle(lineStyleFromName(chann[i]));
        if (!TString(chann[i]).Contains("_low")) {
            if (i == 0 && exp != 0) {
                leg->AddEntry(obsi[i], TString("Combined ")+(lessSpacesInLegends?" obs.":" observed"), "L");
                leg->AddEntry(exp,     (lessSpacesInLegends?"Exp. for "+SM+" Higgs":"Expected for "+SM+" Higgs"), "L");
                if (toys) leg->AddEntry(toys, "Comb. ensemble", "P");
            } else {
                leg->AddEntry(obsi[i], channelFromName(chann[i],/*withspaces=*/true), "L");
            }
        }
    }
    if (exp) {
        TGraphAsymmErrors *smooth = (SM=="SM4" ? smoothSMCLs(exp, 5, 3) : smoothSMCLs(exp, 7,2));
        if (smooth == 0) std::cout << "Smoothing of " << expi[i]->GetName() << " returned ZERO" << std::endl;
        else exp = smooth;
        exp->SetLineWidth(4); exp->SetLineStyle(7); exp->SetLineColor(1);
        exp->Draw("LX");
    }
    if (toys) toys->Draw("P");
    for (int i = nchann-1; i >= 0; --i) {
        if (obsi[i] == 0) continue;
        obsi[i]->Draw("LX");
    }
    leg->Draw();
    TString myspam = SPAM; if (!TString(chann[0]).Contains("comb")) myspam("\n") = "\n" + channelFromName(chann[0]) + "\n";
    finalize(who+"_all"+postfix,xmin,xmax,ymin,ymax,myspam,true);
    if (xmax >= x_zoom && xmin <= x_zoom) {
        c1->SetLogx(1);
        frame0.GetXaxis()->SetMoreLogLabels(); frame0.GetXaxis()->SetNoExponent();
        finalize(who+"_all"+postfix+"_logx",xmin,xmax,ymin,ymax,myspam,true);
        c1->SetLogx(0);
        leg = newLegend(legxstart,.153,.95-0.02*isSquareCanvas,legyend); leg->SetTextSize(0.031);
        xmin = xmin2;
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
                leg->AddEntry(obsi[i], TString("Combined ")+(lessSpacesInLegends?" obs.":" observed"), "L");
                leg->AddEntry(exp,     (lessSpacesInLegends?"Exp. for "+SM+" Higgs":"Expected for "+SM+" Higgs"), "L");
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


void drawCompPVal(TString who1, TString who2, TString label="One", TString oldlabel="Two", TString what="auto", double ymin=1e-10) {
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
    old->SetLineWidth(4);
    old->SetLineColor(dcol); old->SetMarkerColor(dcol); 
    old->SetMarkerStyle(25); old->SetMarkerSize(1.2);
    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    double ymax = 1.0;
    old->Draw("LP");    
    obs->Draw("LP");    
    leg = newLegend(.65,.14,.92,.28); leg->SetTextSize(0.04);
    leg->AddEntry(obs,    label, "LP");
    leg->AddEntry(old, oldlabel, "LP");
    leg->Draw();
    setCanvas(&frame0, "", ymin, ymax, "Local p-value");
    TString myspam = SPAM; if (!who1.Contains("comb")) myspam("\n") = "\n" + channelFromName(who1) + "\n";
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



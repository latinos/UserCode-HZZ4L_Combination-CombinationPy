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
    if (isML) { obs->SetFillColor(colorFit68); obs->Draw("E3"); frame0.Draw("AXIGSAME"); frame0.Draw("AXIS SAME");}
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
        if (isML) { obs->SetFillColor(colorFit68); obs->Draw("E3"); frame2.Draw("AXIS SAME"); frame2.Draw("AXIGSAME"); }
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



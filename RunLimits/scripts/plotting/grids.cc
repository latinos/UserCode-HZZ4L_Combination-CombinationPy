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


bool CLs_debug_apriori_grid=false;
TList CLs_grids;
TGraphAsymmErrors *makeAPrioriGrid(TString who) {
    std::cout << "Asked for grid " << who << std::endl; //CLs_grids.ls();
    if (CLs_grids.FindObject("grid_"+who)) return (TGraphAsymmErrors *) CLs_grids.FindObject("grid_"+who);
    return makeGrid(who+"_obs", who+"_median", 0.4, 2.5);
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




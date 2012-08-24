TString prefix = "";

void cleanCLs(TGraphAsymmErrors *graph) {
    bool foundBad = false;
    do {
        foundBad = false;
        for (int i = 0, n = graph->GetN(); i < n; ++i) {
            if (graph->GetX()[i] == 0) { 
                graph->RemovePoint(i); 
                foundBad = true; 
                break; 
            } 
        }
    } while (foundBad);
    
}
void doIt(TFile *bands, TString name, TString rootfile) {
    if (/*error=*/gSystem->AccessPathName(rootfile)) return;
    zero_is_valid = (name.Contains("pval") || name.Contains("sig") || name.Contains("smcls") || name.Contains("ml_") || name.Contains("mlz_"));
    do_bands_asimov = name.Contains("pval");
    use_precomputed_quantiles = name.Contains("cls") || name.Contains("ml_") || name.Contains("mlz_");
    do_bands_95 = !(name.Contains("ml_") || name.Contains("mlz_"));
    //if (name.Contains("pval") && name.Contains("comb")) do_bands_cputime = true;
    makeBands(bands, name, rootfile, 0);
    std::cout << "Tryring to do " << name << std::endl;
    if (name.Index("ml") == 0) {
        if (bands->Get(name+"_median") == 0) {  std::cerr << "Missing median??" << std::endl; return; } 
        TGraphAsymmErrors *obs = (TGraphAsymmErrors *) bands->Get(name+"_median")->Clone();
        if (obs == 0) return;
        bands->Delete(name+"_obs;*");
        bands->Delete(name+"_median;*");
        obs->SetName(name+"_obs");
        bands->WriteTObject(obs, name+"_obs");
        printLineAErr(bands, name+"_obs", "results/"+prefix+name+".txt");
    } else {
        printBand(bands, name, "results/"+prefix+name+".txt", /*mean=*/false);
        if (rootfile.Contains("PLPE")) {
            printLine(bands, name+"_asimov", "results/"+prefix+name+"_expected.txt", /*mean=*/false);
        }
    }
    zero_is_valid = false;
    use_precomputed_quantiles = false;
    do_bands_95 = true;
    do_bands_asimov = false;
    do_bands_cputime = false; do_bands_realtime = false;
}

const int    nfcbands = 4;
const char * fcbandnam[nfcbands] = { "68", "95", "99", "9973" };
const double fcbandwid[nfcbands] = { 0.68, 0.95, 0.99, 0.9973 };
void doFc(TFile *bands, TString name, TString rootfile) {
    if (/*error=*/gSystem->AccessPathName(rootfile)) return;
    TFile *rfile = TFile::Open(rootfile);
    for (int i = 0; i < nfcbands; ++i) {
        std::cout << "Make " << name << std::endl;
        TGraphAsymmErrors *fcb = theFcBelt(rfile, 1, 0, Observed, fcbandwid[i]);
        if (fcb == 0) continue;
        fcb->SetName(name+"_"+fcbandnam[i]);
        bands->WriteTObject(fcb);
    }
    printFcBand(bands, name, name+".txt", nfcbands, fcbandnam);
}
void importCloned(TFile *bands, TString oldfile, TString oldname, int nwhat, const char **what, int nwho, const char **who) {
    TFile *oldBands = TFile::Open(oldfile);
    for (int i = -2; i < nwho; ++i) {
        TString whoI, who0;
        switch (i) {
            case -1: whoI = "comb"; who0 = "comb"; break;
            case -2: whoI = "hzz2l2q"; who0 = "hzz2l2q"; break;
            default: whoI = who[i]; who0 = who[i]; break;
        }
        if (oldname == "lp11" && whoI == "hww") whoI = "hwwc";
        for (int j = 0; j < nwhat; ++j) {
            TString name  = TString::Format("%s_%s", what[j], who0.Data());
            TString nameI = TString::Format("%s_%s", what[j], whoI.Data());
            TGraphAsymmErrors *oldg_obs    = (TGraphAsymmErrors *) oldBands->Get(nameI+"_obs");
            TGraphAsymmErrors *oldg_median = (TGraphAsymmErrors *) oldBands->Get(nameI+"_median");
            TGraphAsymmErrors *oldg_medi95 = (TGraphAsymmErrors *) oldBands->Get(nameI+"_median_95");
            if (oldg_obs) {
                cleanCLs(oldg_obs);
                bands->WriteTObject(oldg_obs->Clone(name+"_"+oldname+"_obs"), name+"_"+oldname+"_obs");
                std::cout  << "Imported " << oldname << " result " << name << "_obs" << std::endl;
            }
            if (oldg_median) {
                cleanCLs(oldg_median);
                bands->WriteTObject(oldg_median->Clone(name+"_"+oldname+"_median"), name+"_"+oldname+"_median");
                std::cout  << "Imported " << oldname << " result " << name << "_median" << std::endl;
            }
            if (oldg_medi95) {
                cleanCLs(oldg_medi95);
                bands->WriteTObject(oldg_medi95->Clone(name+"_"+oldname+"_median_95"), name+"_"+oldname+"_median_95");
                std::cout  << "Imported " << oldname << " result " << name << "_median_95" << std::endl;
            }
        }
    }
    oldBands->Close();
}

void importLandS(TDirectory *bands) {
    importLandS(bands, "mlz_comb_lands",    "lands/comb_maxllfit.root", 1, 0);
    importLandS(bands, "mlz_combe_lands",   "lands/comb_maxllfit_ebereso.root", 1, 0);
    importLandS(bands, "pvala_comb_lands",  "lands/comb_pvalue.root",   1, 0);
    importLandS(bands, "pvala_comb_lands",  "lands/comb_pvalueSB.root", 0, 1);
    importLandS(bands, "pvala_combe_lands", "lands/comb_pvalue_ebereso.root", 1, 0);
}
void makeBands(int which=0) {
    gROOT->LoadMacro("$CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/plotting/bandUtils.cxx+");

    switch (which)  {
        case 1: prefix = "SM4_"; break;
        case 2: prefix = "FF_";  break;
    }

    TFile *bands = new TFile(Form("results/%sbands.root",prefix.Data()),"RECREATE");

    std::cout << "Will save bands to " << bands->GetName() << std::endl;
    
    obs_avg_mode = MedianObs;
    do_bands_nosyst = false;
    do_bands_mean   = false;
    do_bands_ntoys  = false;
    do_bands_asimov = false;
    halfint_masses = true;

    const int nwhat = 13;
    const char *what[nwhat] = { "pla", "pvala", "pvala", "pval", "acls" , "ml", "mlz", "cls",  "bayes", "smcls", "smacls", "acls90", "acls99"};
    const char *WHAT[nwhat] = { "PLC", "PLP",   "PLPE",  "PVAL", "ASCLS", "ML", "MLZ", "FREQ", "BAYES", "SMCLS", "SMASCLS", "ASCLS90", "ASCLS99"};

    
    const int nwho = 20;
    const char *who[nwho] = { "comb", "htt", "vhtt", "httm", "htt0", "hzz", "hzz2l2t", "hzz2l2nu", "hzz2l2q_all", "hzz4l", "hww", "vhww3l", "hww2l", "hgg_novbf", "hgg_vbf", "hgg", "vhbb", "combs", "combp_low", "combl_low" };
    const char *WHO[nwho] = { "COMB", "HTT", "VHTT", "HTTM", "HTT0", "HZZ", "HZZ2L2T", "HZZ2L2NU", "HZZ2L2Q",     "HZZ4L", "HWW", "VHWW3L", "HWW2L", "HGG_NOVBF", "HGG_VBF", "HGG", "VHBB", "COMBS", "COMBP", "COMBL" };

    for (int i = 0; i < nwho; ++i) {
        //if (!TString(who[i]).Contains("comb")) continue;
        //doFc(bands, TString::Format("fc_%s",who[i]), TString::Format("higgsCombine%s_FC.root", WHO[i]));
        //continue;
        for (int j = 0; j < nwhat; ++j) {
            //if (TString(what[j]) != "acls") continue;
            doIt(bands, TString::Format("%s_%s", what[j], who[i]),
                        TString::Format("results/higgsCombine%s%s_%s.root", prefix.Data(), WHO[i], WHAT[j]));
        }
        //break;
        #if 0
        if (!/*error=*/gSystem->AccessPathName(TString::Format("timecls_%s.txt", who[i]))) {
            TGraphErrors *time = new TGraphErrors(TString::Format("timecls_%s.txt", who[i]), "Mass: %lg Time: %lg +/- %lg s");
            time->SetName(TString::Format("timecls_%s_obs", who[i]));
            bands->WriteTObject(time);
        }
        #endif
    }

    const int npvflav = 9;
    const char *pvflav[npvflav] = { "M2M", "S", "SZ3", "BZ1", "BZ", "TQ0", "TQ03", "TQ05N", "TQ05" };
    for (int i = 0; i < npvflav; ++i) {
        doIt(bands, TString::Format("pvala_%s_comb", pvflav[i]),
                    TString::Format("results/higgsCombine%sCOMB_PLP.%s.root", prefix.Data(), pvflav[i]));
    }
#if 0
    if (which == 0) {
        importLine(bands, "tevatron_obs",    "tevatron.obs.txt");    
        importLine(bands, "tevatron_median", "tevatron.exp.txt");    
        std::cout  << "Imported Tevatron results" << std::endl;
        importCloned(bands, "results/bands.paper.root", "paper", nwhat, what, nwho, who);
    }
#endif

    for (int j = 0; j < nwhat; ++j) {
        TString w = what[j];
        cutBands(bands, w+"_hzz2l2q_all", w+"_hzz2l2q_low", 130, 164);
        cutBands(bands, w+"_hzz2l2q_all", w+"_hzz2l2q",     200, 600);
        cutBands(bands, w+"_hzz", w+"_combp_high", 150.5, 600);
        cutBands(bands, w+"_hww", w+"_combl_high", 145.5, 600);
        pasteBands(bands, w+"_combl_low", w+"_combl_high", w+"_combl");
        pasteBands(bands, w+"_combp_low", w+"_combp_high", w+"_combp");
    }
    if (which == 0) {
        cutBands(bands, "smcls_comb",       "smcls_comb_new", 110, 200);
        cutBands(bands, "smcls_comb_paper", "smcls_comb_old", 200, 600);
        pasteBands(bands, "smcls_comb_new", "smcls_comb_old", "smcls_comb_patch");
        cutBands(bands, "cls_comb",       "cls_comb_new", 110, 200);
        cutBands(bands, "cls_comb_paper", "cls_comb_old", 200, 600);
        pasteBands(bands, "cls_comb_new", "cls_comb_old", "cls_comb_patch");
    }

}

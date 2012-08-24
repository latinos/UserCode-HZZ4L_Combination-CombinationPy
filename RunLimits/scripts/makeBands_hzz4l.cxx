void fixHGG_FREQ(TFile *bands, TString name) {
    TGraphAsymmErrors *graph = bands->Get(name);
    const int numsmpoints = 61;
    float smmasses[numsmpoints] = {110,110.5,111,111.5,112,112.5,113,113.5,114,114.5,115,115.5,116,116.5,117,117.5,118,118.5,119,119.5,120,120.5,121,121.5,122,122.5,123,123.5,124,124.5,125,125.5,126,126.5,127,127.5,128,128.5,129,129.5,130,130.5,131,131.5,132,132.5,133,133.5,134,134.5,135,135.5,136,136.5,137,137.5,138,138.5,139,139.5,140};
    float smxsecs[numsmpoints]  = {22.7112,22.5139,22.3165,22.1192,21.9219,21.7246,21.5272,21.3299,21.1326,20.9352,20.7379,20.5654,20.3928,20.2203,20.0477,19.8752,19.7027,19.5301,19.3576,19.185,19.0125,18.8608,18.7092,18.5576,18.4059,18.2542,18.1026,17.951,17.7993,17.6476,17.496,17.3593,17.2226,17.086,16.9493,16.8126,16.6759,16.5392,16.4026,16.2659,16.1292,16.0095,15.8898,15.7702,15.6505,15.5308,15.4111,15.2914,15.1718,15.0521,14.9324,14.8237,14.715,14.6064,14.4977,14.389,14.2803,14.1716,14.063,13.9543,13.84562};
    float smbrs[numsmpoints]    = {0.00197,0.00198737,0.00200441,0.00202113,0.00203754,0.00205365,0.00206948,0.00208501,0.00210028,0.00211527,0.00213,0.00214331,0.0021563,0.00216899,0.00218138,0.00219348,0.0022053,0.00221685,0.00222815,0.0022392,0.00225,0.00225571,0.00226125,0.00226662,0.00227182,0.00227687,0.00228177,0.00228652,0.00229114,0.00229563,0.0023,0.00229526,0.00229072,0.00228635,0.00228215,0.00227811,0.00227422,0.00227047,0.00226686,0.00226337,0.00226,0.00224526,0.00223124,0.00221791,0.0022052,0.00219308,0.00218151,0.00217044,0.00215986,0.00214972,0.00214,0.00211438,0.00209031,0.00206765,0.00204629,0.00202612,0.00200703,0.00198895,0.0019718,0.0019555,0.00193};
    for (int i = 0; i < graph->GetN(); ++i) {
        double x = graph->GetX()[i];
        double scale = 0.;
        for (int j = 0; j < numsmpoints; ++j) {
            if (fabs(x - smmasses[j]) < 0.001) { 
                scale = 1.0/(smxsecs[i]*smbrs[i]);
            }
        }
        if (scale == 0) std::cout << "ERROR: mass point " << x << " not in tables." << std::endl;
        graph->GetY()[i] *= scale;
        graph->GetEYhigh()[i] *= scale;
        graph->GetEYlow()[i] *= scale;
    }
    
    bands->WriteTObject(graph, name,  "Overwrite");
}
void doIt(TFile *bands, TString name, TString rootfile) {
    if (/*error=*/gSystem->AccessPathName(rootfile)) return;
    makeBands(bands, name, rootfile, 0);
    if (name == "cls_hgg") {
        fixHGG_FREQ(bands, "cls_hgg_obs");
        fixHGG_FREQ(bands, "cls_hgg_median");
        fixHGG_FREQ(bands, "cls_hgg_median_95");
    }
    printBand(bands, name, name+".txt", /*mean=*/false);
}
void makeBands(int hadd=0) {

    TString prefix="./";
    TFile *bands = new TFile("bands.root","RECREATE");

    obs_avg_mode = LogMeanObs;
    doIt(bands, "bayes_comb", "higgsCombineCOMB_BAYES.root");
    obs_avg_mode = MeanObs;
    //doIt(bands, "bayes_mean_comb", "higgsCombineCOMB_BAYES.root");

    obs_avg_mode = MeanObs;
    use_precomputed_quantiles = true;
    //doIt(bands, "cls_hww",      "higgsCombineHWW_FREQ.root");
    doIt(bands, "cls_hzz4l",    "higgsCombineHZZ4L_FREQ.root");
    //doIt(bands, "cls_htt",      "higgsCombineHTT_FREQ.root");
    //doIt(bands, "cls_hgg",      "higgsCombineHGG_FREQ.root");
    //doIt(bands, "cls_hzz2l2q_full",      "higgsCombineHZZ2L2Q_FREQ.root");
    //cutBands(bands, "cls_hzz2l2q_full",  "cls_hzz2l2q", 225., 999.); 
    //importBands(bands, "cls_hzz2l2nu_median", "lands_cls_hww2l2nu.txt", /*hasobs=*/true);
    //bands->WriteTObject(((TGraph*)bands->Get("cls_hzz2l2nu_median_obs"))->Clone("cls_hzz2l2nu_obs"), "cls_hzz2l2nu_obs");

    doIt(bands, "cls_comb",      "higgsCombineCOMBS_FREQ.root");
    zero_is_valid = true;
    doIt(bands, "smcls_comb",   "higgsCombineCOMBS_SMCLS.root");
    zero_is_valid = false;
    use_precomputed_quantiles = false;

    //doIt(bands, "pla_comb",  "higgsCombineCOMB_PLC.root");
    //doIt(bands, "pla_hww",      "higgsCombineHWW_PLC.root");
    //doIt(bands, "pla_hgg",      "higgsCombineHGG_PLC.root");
    //doIt(bands, "pla_htt",      "higgsCombineHTT_PLC.root");
    doIt(bands, "pla_hzz4l",    "higgsCombineHZZ4L_PLC.root");
    //doIt(bands, "pla_hzz2l2q_full",  "higgsCombineHZZ2L2Q_PLC.root");
    //cutBands(bands, "pla_hzz2l2q_full",  "pla_hzz2l2q", 225., 999.); 
    //doIt(bands, "pla_hzz2l2nu",  "higgsCombineHZZ2L2NU_PLC.root");

    /*
    doIt(bands, "ml_hww",  "higgsCombineHWW_ML.root");
    doIt(bands, "ml_hzz2l2nu",  "higgsCombineHZZ2L2NU_ML.root");
    doIt(bands, "ml_hzz4l",  "higgsCombineHZZ4LS2_ML.root");
    doIt(bands, "ml_hzz2l2q",  "higgsCombineHZZ2L2Q_ML.root");
    doIt(bands, "ml_hgg",  "higgsCombineHGG_ML.root");
    doIt(bands, "ml_htt",  "higgsCombineHTT_ML.root");
    doIt(bands, "ml_comb",  "higgsCombineCOMBS2_ML.root");
    */
    //importBands(bands, "ml_comb_obs",     "ml_comb.txt",     false, false);
    //importBands(bands, "ml_hww_obs",      "ml_hww.txt",      false, false);
    //importBands(bands, "ml_htt_obs",      "ml_htt.txt",      false, false);
    //importBands(bands, "ml_hgg_obs",      "ml_hgg.txt",      false, false);
    importBands(bands, "ml_hzz4l_obs",    "ml_hzz4l.txt",    false, false);
    //importBands(bands, "ml_hzz2l2q_obs",  "ml_hzz2l2q.txt",  false, false);
    //importBands(bands, "ml_hzz2l2nu_obs", "ml_hzz2l2nu.txt", false, false);

    //doIt(bands, "pvala_comb",  "higgsCombineCOMB_PLP.root");
    //doIt(bands, "pvala_hww",   "higgsCombineHWW_PLP.root");
    //doIt(bands, "pvala_hgg",   "higgsCombineHGG_PLP.root");
    //doIt(bands, "pvala_htt",      "higgsCombineHTT_PLP.root"); 
    doIt(bands, "pvala_hzz4l",    "higgsCombineHZZ4L_PLP.root");
    //doIt(bands, "pvala_hzz2l2q_full",  "higgsCombineHZZ2L2Q_PLP.root");
    //cutBands(bands, "pvala_hzz2l2q_full",  "pvala_hzz2l2q", 225., 999.); 
    //doIt(bands, "pvala_hzz2l2nu",  "higgsCombineHZZ2L2NU_PLP.root");
}

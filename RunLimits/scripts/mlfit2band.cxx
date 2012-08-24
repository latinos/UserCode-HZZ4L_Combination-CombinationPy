void mlfit2band(const char *name = "COMB") {
    FILE *f = fopen("masses.txt", "r");
    int mass;
    while (fscanf(f, "%d", &mass) == 1) {
        TString fname = TString::Format("%d/mlfit%s.root", mass, name);
        if (/*error=*/gSystem->AccessPathName(fname)) continue;
        gFile = TFile::Open(fname);
        RooFitResult *rfr = (RooFitResult *) gFile->Get("fit_s");
        RooRealVar   *r = (RooRealVar *) rfr->floatParsFinal().find("r");
        gFile->Close();
        printf("%3d %7.3f %7.3f %7.3f\n", mass, r->getVal()+r->getAsymErrorLo(), r->getVal(), r->getVal()+r->getAsymErrorHi());
    }
    fclose(f);
}

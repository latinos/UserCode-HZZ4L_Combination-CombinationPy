TGraph* bestFit(TTree *t, TString x, TString y) {
    t->Draw(y+":"+x, "deltaNLL == 0");
    TGraph *gr0 = (TGraph*) gROOT->FindObject("Graph")->Clone();
    gr0->SetMarkerStyle(34); gr0->SetMarkerSize(2.0);
    if (gr0->GetN() > 1) gr0->Set(1);
    return gr0;
}

TList* contourPlot(TTree *t, TString x, TString y, double pmin, double pmax, TGraph *bestFit) {
    int n = t->Draw(y+":"+x, Form("%f <= quantileExpected && quantileExpected <= %f && quantileExpected != 1",pmin,pmax));
    TGraph *gr = (TGraph*) gROOT->FindObject("Graph")->Clone();
    gr->SetName(Form("grContour_%f",pmax-pmin));
    Double_t x0 = bestFit->GetX()[0], y0 = bestFit->GetY()[0];
    Double_t *xi = gr->GetX(), *yi = gr->GetY();
    int n = gr->GetN();
    for (int i = 0; i < n; ++i) { xi[i] -= x0; yi[i] -= y0; }
    gr->Sort(&TGraph::CompareArg);
    for (int i = 0; i < n; ++i) { xi[i] += x0; yi[i] += y0; }
    TList *ret = new TList();
    ret->Add(gr);
    return ret;
}

TGraph *contourPlot(TTree *t, TString x, TString y, double pmin, double pmax, TGraph *bestFit) {
    int n = t->Draw(y+":"+x, Form("%f <= quantileExpected && quantileExpected <= %f && quantileExpected != 1",pmin,pmax));
    TGraph *gr = (TGraph*) gROOT->FindObject("Graph")->Clone();
    gr->SetName(Form("grContour_%f",pmax-pmin));
    Double_t x0 = bestFit->GetX()[0], y0 = bestFit->GetY()[0];
    Double_t *xi = gr->GetX(), *yi = gr->GetY();
    int n = gr->GetN();
    for (int i = 0; i < n; ++i) { xi[i] -= x0; yi[i] -= y0; }
    gr->Sort(&TGraph::CompareArg);
    for (int i = 0; i < n; ++i) { xi[i] += x0; yi[i] += y0; }
    TList *ret = new TList();
    //   ret->Add(gr);
    return gr;
}

int countGridPointsFromTree(TTree *t, TString x, double xmin = -1, double xmax = -1) {
    if (xmin == xmax) {
        xmin = t->GetMinimum(x);
        xmax = t->GetMaximum(x);
    }
    t->Draw(Form("%s>>h1000(1000,%10g,%10g)", x.Data(),xmin-1e-4,xmax+1e-4), "deltaNLL > 0");
    TH1 *h1000 = (TH1*) gROOT->FindObject("h1000");
    int bins = 0;
    for (int i = 1, n = h1000->GetNbinsX(); i <= n; ++i) {
        if (h1000->GetBinContent(i) != 0) bins++;
    }
    h1000->Delete();
    return bins;
}
TH2 *treeToHist2D(TTree *t, TString x, TString y, TString name, double xmin = -1, double xmax = -1, double ymin = -1, double ymax = -1) {
    if (xmin == xmax) {
        xmin = t->GetMinimum(x);
        xmax = t->GetMaximum(x);
    } 
    if (ymin == ymax) {
        ymin = t->GetMinimum(y);
        ymax = t->GetMaximum(y);
    }
    int xbins = countGridPointsFromTree(t,x,xmin,xmax);
    int ybins = countGridPointsFromTree(t,y,ymin,ymax);
    double dx = (xmax-xmin)/(xbins-1);
    double dy = (ymax-ymin)/(ybins-1);
    xmin -= 0.5*dx; xmax += 0.5*dx;
    ymin -= 0.5*dy; ymax += 0.5*dy;
    if (fabs(xmin) < 1e-5) xmin = 0;
    if (fabs(xmax) < 1e-5) xmax = 0;
    //std::cout << "In making " << name << ", guessed " << xbins << " bins for " << x << " from " << xmin << " to " << xmax << std::endl;
    //std::cout << "In making " << name << ", guessed " << ybins << " bins for " << y << " from " << ymin << " to " << ymax << std::endl;
    t->Draw(Form("2*deltaNLL:%s:%s>>%s_prof(%d,%10g,%10g,%d,%10g,%10g)", y.Data(), x.Data(), name.Data(), xbins, xmin, xmax, ybins, ymin, ymax), "deltaNLL != 0", "PROF");
    TH2 *prof = (TH2*) gROOT->FindObject(name+"_prof");
    TH2D *h2d = new TH2D(name, name, xbins, xmin, xmax, ybins, ymin, ymax);
    for (int ix = 1; ix <= xbins; ++ix) {
        for (int iy = 1; iy <= ybins; ++iy) {
             h2d->SetBinContent(ix, iy, prof->GetBinContent(ix,iy));
        }
    }
    h2d->SetDirectory(0);
    return h2d;
}

TH2D* frameTH2D(TH2D *in){

        double frameValue = 1000;
        if (TString(in->GetName()).Contains("bayes")) frameValue = 0.0;

	Double_t xw = in->GetXaxis()->GetBinWidth(0);
	Double_t yw = in->GetYaxis()->GetBinWidth(0);

	Int_t nx = in->GetNbinsX();
	Int_t ny = in->GetNbinsY();

	Double_t x0 = in->GetXaxis()->GetXmin();
	Double_t x1 = in->GetXaxis()->GetXmax();

	Double_t y0 = in->GetYaxis()->GetXmin();
	Double_t y1 = in->GetYaxis()->GetXmax();

	TH2D *framed = new TH2D(
			Form("%s framed",in->GetName()),
			Form("%s framed",in->GetTitle()),
			nx + 2, x0-xw, x1+xw,
			ny + 2, y0-yw, y1+yw
			);

	//Copy over the contents
	for(int ix = 1; ix <= nx ; ix++){
		for(int iy = 1; iy <= ny ; iy++){
			framed->SetBinContent(1+ix, 1+iy, in->GetBinContent(ix,iy));
		}
	}
	//Frame with huge values
	nx = framed->GetNbinsX();
	ny = framed->GetNbinsY();
	for(int ix = 1; ix <= nx ; ix++){
		framed->SetBinContent(ix,  1, frameValue);
		framed->SetBinContent(ix, ny, frameValue);
	}
	for(int iy = 2; iy <= ny-1 ; iy++){
		framed->SetBinContent( 1, iy, frameValue);
		framed->SetBinContent(nx, iy, frameValue);
	}

	return framed;
}



TList* contourFromTH2(TH2 *h2in, double threshold) {
    std::cout << "Getting contour at threshold " << threshold << " from " << h2in->GetName() << std::endl;
    //http://root.cern.ch/root/html/tutorials/hist/ContourList.C.html
    Double_t contours[1];
    contours[0] = threshold;

    TH2D *h2 = frameTH2D((TH2D*)h2in);

    h2->SetContour(1, contours);

    // Draw contours as filled regions, and Save points
    h2->Draw("CONT Z LIST");
    gPad->Update(); // Needed to force the plotting and retrieve the contours in TGraphs

    // Get Contours
    TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
    TList* contLevel = NULL;

    if (conts == NULL || conts->GetSize() == 0){
        printf("*** No Contours Were Extracted!\n");
        return 0;
    }

    TList *ret = new TList();
    for(int i = 0; i < conts->GetSize(); i++){
        contLevel = (TList*)conts->At(i);
        printf("Contour %d has %d Graphs\n", i, contLevel->GetSize());
        for (int j = 0, n = contLevel->GetSize(); j < n; ++j) {
            TGraph *gr1 = (TGraph*) contLevel->At(j)->Clone();
            ret->Add(gr1);
            //break;
        }
    }
    return ret;
}


void plotScan(){

  TFile *f=new TFile("higgsCombineTest2D.MultiDimFit.mH126.root","READ");
  TTree *t=(TTree*)f->Get("limit");
  TGraph *gr_BestFit=(TGraph*) bestFit(t, "x", "r");
  //  contourPlot(t, "r", "x", double pmin, double pmax,gr_BestFit );
  TH2F* h2D=(TH2F*)treeToHist2D(t, "x", "r", "scan2d", 0.0, 1.0, 0.0, 3.0);
  TGraph *mycont68=(TGraph*)contourFromTH2(h2D, 2.3/2.0 )->At(0);
  mycont68->SetMarkerSize(1.25);
  mycont68->SetLineWidth(1.5);
  TGraph *mycont95=(TGraph*)contourFromTH2(h2D, 6.17/2.0)->At(0);
  mycont95->SetMarkerSize(1.25);
  mycont95->SetLineWidth(1.5);
  mycont95->SetLineStyle(kDashed);

  TGraph *mycont3sigma=(TGraph*)contourFromTH2(h2D, 11.8/2.0)->At(0);
  mycont3sigma->SetMarkerSize(1.25);
  mycont3sigma->SetLineWidth(1.5);
  mycont3sigma->SetLineStyle(kDotted);

  TCanvas *c1=new TCanvas("cc1","CANVAS1",800,800);
  c1->cd();
  gPad->SetRightMargin(0.085);
  h2D->SetXTitle("x (scalar=0, pseudoscalar=1)");
  h2D->SetYTitle("#mu (signal strength)");
  h2D->GetXaxis()->SetLabelSize(0.035);
  h2D->GetYaxis()->SetLabelSize(0.035);
  h2D->Draw("colz");
  mycont68->Draw("L");
  mycont95->Draw("L");
  gr_BestFit->Draw("P");

  /* TCanvas *c2=new TCanvas("cc2","CANVAS2-CONTOUR",800,800);
  c2->cd();
  TFile *fOut=new TFile("test_scan2d.root","RECREATE");
  // TList *lCont=(TList*)contourPlot(t,"x", "r", 0.31, 1.0,gr_BestFit );
  TGraph *grCont=(TGraph*)contourPlot(t,"x", "r", 0.31, 1.0,gr_BestFit );
  grCont->SetMarkerColor(kRed);
  grCont->SetMarkerSize(1.2);
  grCont->Draw("AP");
  */
}

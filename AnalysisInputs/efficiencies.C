#include "TH1F.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <sstream> 
#include "TCutG.h"
#include "TFile.h"
#include "TH2.h"
#include "TPad.h"



//gSystem->AddIncludePath("-I$ROOFITSYS/include/");

//#include "RooCBShape.h"
using namespace RooFit ;
using namespace std;


void efficiencies(){
	
	gSystem->Load("libMathCore");
	
	//--------------------------------------------------------------------
	ifstream fin;
	fin.open("/scratch/hep/ntran/dataFiles/HZZ4L/datasets/datasets_highmass/Nevt_2mu2e.txt");
	
	std::vector< double > masses_V;
	std::vector< double > effs_V;
	
	while (!fin.eof() && fin.good()){

		double mass, val1, val2;
		fin >> mass >> val1 >> val2;
		std::cout << mass << " " << val1 << " " << val2 << std::endl;
		masses_V.push_back( mass );
		effs_V.push_back( val1/val2 );
		
	}
	
	const int sizeV = masses_V.size();
	double masses_A[sizeV];
	double effs_A[sizeV];
	for (int i = 0; i < sizeV; i++){
		
		masses_A[i] = masses_V[i];
		effs_A[i] = effs_V[i];		
		
	}
	
	TGraph* grEff = new TGraph( sizeV, masses_A, effs_A );
	grEff->SetMarkerStyle(20);

	TF1 *polyFunc= new TF1("polyFunc","([0]+[1]*TMath::Erf( (x-[2])/[3] ))*([4]+[5]*x+[6]*x*x)", 110., 600.);
	polyFunc->SetParameters(0.10, 1.42616e-01, 100., 40.,0.6448,0.00226,-1.859e-06);
	//polyFunc1->FixParameter(0,0.10);
	polyFunc->SetLineColor(4);	
	TCanvas *c = new TCanvas("c","c",700,700);
	c->SetGrid();
	TH1F *hr = c->DrawFrame(0.,-0.1,610.,1.);
	grEff->Fit(polyFunc,"Rt");
	grEff->Draw("P");
	polyFunc->Draw("sames");
	c->SaveAs("figs/effs_2mu2e_hi.eps");
	
	
	//TF1 *polyFunc1= new TF1("polyFunc1","[0]+[1]*TMath::Erf( (x-[2])/[3] )", 110., 600.);
	//polyFunc1->SetParameters(0.10, 1.42616e-01, 100., 40.);
	//polyFunc1->FixParameter(0,0.10);
	//polyFunc1->SetLineColor(4);
	
	//TF1 *polyFunc2= new TF1("polyFunc2","[0]+[1]*TMath::Erf( (x-[2])/[3] )", 180., 600.);
	//polyFunc2->SetParameters(0.10, 1.42616e-01, 180, 1.82020e+02);
	//polyFunc2->FixParameter(0,0.10);
	//polyFunc2->SetLineColor(2);
	
	/*
	TCanvas *c = new TCanvas("c","c",700,700);
	c->SetGrid();
	TH1F *hr = c->DrawFrame(0.,0.,610.,1.);
	grEff->Fit(polyFunc1,"Rt");
	//grEff->Fit(polyFunc2,"Rt");
	grEff->Draw("P");
	polyFunc1->Draw("sames");
	//polyFunc2->Draw("sames");
	c->SaveAs("figs/effs_2mu2e.eps");
	*/
}

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
#include "TTree.h"


//gSystem->AddIncludePath("-I$ROOFITSYS/include/");

//#include "RooCBShape.h"
using namespace RooFit ;
using namespace std;


void createEventFiles(){
	
	gSystem->Load("libMathCore");
	
	ifstream fin;
	fin.open("CMSdata/eventList_4.700.txt");
	
	ofstream f_2e2mu;
	f_2e2mu.open("CMSdata/hzz2e2mu_4.700.dat");
	ofstream f_4mu;
	f_4mu.open("CMSdata/hzz4mu_4.700.dat");
	ofstream f_4e;
	f_4e.open("CMSdata/hzz4e_4.700.dat");
	
    long Run,Event,LumiSect;
    int idL1, idL2, idL3, idL4;
    double mass4l, massError;
	TFile *fout = new TFile("massErrors.root","RECREATE");
    TTree *tr = new TTree("passedEvents", "passedEvents");
    tr->Branch("Run",&Run,"Run/l");
    tr->Branch("Event",&Event,"Event/l");
    tr->Branch("LumiSect",&LumiSect,"LumiSect/l");
    tr->Branch("mass4l",&mass4l,"mass4l/D");
    tr->Branch("massError",&massError,"massError/D");
    tr->Branch("idL1",&idL1,"idL1/I");
    tr->Branch("idL2",&idL2,"idL2/I");
    tr->Branch("idL3",&idL3,"idL3/I");
    tr->Branch("idL4",&idL4,"idL4/I");
	
	char output[100];
	int eventCtr = 0;
	if (fin.is_open()) {
		while (!fin.eof()) {
			
			std::string name, channel, period;
			string run, event;
			double mZ1, mZ2, m4l, m4lerr, pT, Y;
			fin >> name >> run >> event >> channel >> mZ1 >> mZ2 >> m4l >> m4lerr >> pT >> Y >> period;
			
			if (channel == "2e2mu"){
				f_2e2mu << m4l << " " << mZ1 << " " << mZ2 << std::endl;
				idL1 = -13; idL2 = 13; idL3 = -11; idL4 = 11;
			}
			else if (channel == "4mu"){
				f_4mu << m4l << " " << mZ1 << " " << mZ2 << std::endl;
				idL1 = -13; idL2 = 13; idL3 = -13; idL4 = 13;				
			}
			else if (channel == "4e"){
				f_4e << m4l << " " << mZ1 << " " << mZ2 << std::endl;
				idL1 = -11; idL2 = 11; idL3 = -11; idL4 = 11;				
			}
			else{
				std::cout << "Invalid channel!" << std::endl;
			}
			
			Run = run;
			Event = event;
			LumiSect = 1;
			mass4l = m4l;
			massError = m4lerr;
			
			tr->Fill();
			
			eventCtr++;
		}
	}
	
	std::cout << "Events: " << eventCtr << std::endl;
	
	fin.close();
	
	f_2e2mu.close();
	f_4mu.close();
	f_4e.close();
	
	fout->Write();
	fout->Close();
}

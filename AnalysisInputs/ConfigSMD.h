//----------> SET INPUT VARIABLES HERE

// Input trees


// With high mass reweights + new MELA (analytical background) + superMELA
TString filePath7TeV = "/afs/cern.ch/user/b/bonato/work/PhysAnalysis/HZZ4L/Trees_191012/PRODFSR_7TeV/";
TString filePath8TeV = "/afs/cern.ch/user/b/bonato/work/PhysAnalysis/HZZ4L/Trees_191012/PRODFSR_8TeV/";
TString filePath7TeVPS = "/afs/cern.ch/user/b/bonato/work/PhysAnalysis/HZZ4L/Trees_191012/JHU_7TeV/";
TString filePath8TeVPS = "/afs/cern.ch/user/b/bonato/work/PhysAnalysis/HZZ4L/Trees_191012/JHU_8TeV/";


// Luminosity, as float and as string to be used in file names, etc.
double lumi7TeV = 5.051;
double lumi8TeV = 12.21;
TString lumistr7TeV = "5.051";
TString lumistr8TeV = "12.21";


// Location of output root files containing data events
TString DataRootFilePath = "../CreateDatacards/CMSdata/"; 
//<----------

const int nPoints7TeV = 22;
int masses7TeV[nPoints7TeV]   = {120,/*125,*/130,140,150,160,170,180,200,210,220,250,300,325,350,400,425,450,475,525,550,575,600}; //FIXME: weights for 125 are wrong.
double mHVal7TeV[nPoints7TeV] = {120,/*125,*/130,140,150,160,170,180,200,210,220,250,300,325,350,400,425,450,475,525,550,575,600};

const int nPoints8TeV = 29;
int masses8TeV[nPoints8TeV]   = {115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,145,150,180,200,250,300,325,350,400,450,500,550,600};
double mHVal8TeV[nPoints8TeV] = {115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,145,150,180,200,250,300,325,350,400,450,500,550,600};


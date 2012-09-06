//----------> SET INPUT VARIABLES HERE

// Input trees
TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/310812/PRODFSR/";
TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/310812/PRODFSR_8TeV/";

// Luminosity, as float and as string to be used in file names, etc.
double lumi7TeV = 5.051;
double lumi8TeV = 5.261;
TString lumistr7TeV = "5.051";
TString lumistr8TeV = "5.261";


// Location of output root files containing data events
TString DataRootFilePath = "../CreateDatacards/CMSdata/"; 
//<----------

//--------------------
// The number and values of mass points for which you have the trees, for 7 and 8 TeV
const int nPoints7TeV = 22;
int masses7TeV[nPoints7TeV]   = {120,125,130,140,150,160,170,180,200,210,220,250,300,325,350,400,425,450,475,/*525,*/550,575,600};
double mHVal7TeV[nPoints7TeV] = {120,125,130,140,150,160,170,180,200,210,220,250,300,325,350,400,425,450,475,/*525,*/550,575,600};

const int nPoints8TeV = 29;
int masses8TeV[nPoints8TeV]   = {115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,145,150,180,200,250,300,325,350,400,450,500,550,600};
double mHVal8TeV[nPoints8TeV] = {115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,145,150,180,200,250,300,325,350,400,450,500,550,600};

//--------------------

// THIS IS TO BE USED TO ADD THE HIGH MASS POINTS
//--------------------
// The number and values of mass points for which you have the trees, for 7 and 8 TeV
//const int nPoints7TeV = 23;
//int masses7TeV[nPoints7TeV]   = {120,125,130,140,150,160,170,180,200,210,220,250,300,325,350,400,425,450,475,525,550,575,600}; //FIXME
//double mHVal7TeV[nPoints7TeV] = {120,125,130,140,150,160,170,180,200,210,220,250,300,325,350,400,425,450,475,525,550,575,600}; //FIXME
//
//const int nPoints8TeV = 37;
//int masses8TeV[nPoints8TeV]   = {115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,145,150,180,200,250,300,325,350,400,450,500,550,600,650,700,750,800,850,900,950,1000};
//double mHVal8TeV[nPoints8TeV] = {115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,145,150,180,200,250,300,325,350,400,450,500,550,600,650,700,750,800,850,900,950,1000};

//--------------------

//----------> SET INPUT VARIABLES HERE

// Input trees

// ICHEP inputs
// TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/240612/PRODFSR";
// TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/240612/PRODFSR_8TeV";

// With LD fixes (phi1, Z1/Z2)
// TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/310812/PRODFSR";
// TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/310812/PRODFSR_8TeV";

// With LD fixes (phi1, Z1/Z2), but use rescaled 7-TeV gZZ for 8 TeV
//TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/310812_ggRescaled/PRODFSR/";
//TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/310812_ggRescaled/PRODFSR_8TeV/";

// With high mass reweights
//TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/310812_HighMassReweight/PRODFSR/";
//TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/310812_HighMassReweight/PRODFSR_8TeV/";

// Experimental
//TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/310812_ggRescaled_HighMassReweight/PRODFSR/";
//TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/310812_ggRescaled_HighMassReweight/PRODFSR_8TeV/";

//171012, HCP unblinding
//TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/171012/PRODFSR/";
//TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/171012/PRODFSR_8TeV/";

//181012, HCP unblinding
// TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/191012/PRODFSR/";
// TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/191012/PRODFSR_8TeV/";

//Final set, with ICHEP data/MC SFs
//TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/261012_ichep/PRODFSR/";
//TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/261012_ichep/PRODFSR_8TeV/";

//Final set, HCP
//TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/261012/PRODFSR/";
//TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/261012/PRODFSR_8TeV/";

//TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/";
//TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/";

//TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/130120/v2/PRODFSR/";
//TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/130120/v2/PRODFSR_8TeV/";

// Legacy paper
//TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/130613b/PRODFSR/";
//TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/130613b/PRODFSR_8TeV/";

TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/130702/PRODFSR/";	 
TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/130702/PRODFSR_8TeV/";


// Luminosity, as float and as string to be used in file names, etc.
double lumi7TeV = 5.051;
double lumi8TeV = 19.79;
TString lumistr7TeV = "5.051";
TString lumistr8TeV = "19.79";


// Location of output root files containing data events
TString DataRootFilePath = "../CreateDatacards/CMSdata/"; 

// THIS IS TO BE USED TO ADD THE HIGH MASS POINTS
//--------------------
// The number and values of mass points for which you have ggH trees, for 7 and 8 TeV
const int nPoints7TeV = 34;
int masses7TeV[nPoints7TeV]   = {120,124,125,126,130,140,150,160,170,180,190,200,210,220,250,275,300,325,350,400,425,450,475,525,550,575,600,650,700,750,800,900,950,1000};
double mHVal7TeV[nPoints7TeV] = {120,124,125,126,130,140,150,160,170,180,190,200,210,220,250,275,300,325,350,400,425,450,475,525,550,575,600,650,700,750,800,900,950,1000};

const int nPoints8TeV = 49;
int masses8TeV[nPoints8TeV]   = {115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,135,140,145,150,160,170,180,190,200,220,250,275,300,325,350,375,400,425,450,475,500,525,550,575,600,650,700,750,800,850,900,950,1000};
double mHVal8TeV[nPoints8TeV] = {115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,135,140,145,150,160,170,180,190,200,220,250,275,300,325,350,375,400,425,450,475,500,525,550,575,600,650,700,750,800,850,900,950,1000};

const int nPoints7TeV_p15 = 12;
int masses7TeV_p15[nPoints7TeV_p15]   = {125,126,400,450,500,550,600,650,700,800,900,1000};
double mHVal7TeV_p15[nPoints7TeV_p15] = {125,126,400,450,500,550,600,650,700,800,900,1000};

const int nPoints8TeV_p15 = 12;
int masses8TeV_p15[nPoints8TeV_p15]   = {125,126,400,450,500,550,600,650,700,800,900,1000};
double mHVal8TeV_p15[nPoints8TeV_p15] = {125,126,400,450,500,550,600,650,700,800,900,1000};

//The number and values of mass points for which you have qqH trees, for 7 and 8 TeV
const int nVBFPoints7TeV = 33; 
int VBFmasses7TeV[nVBFPoints7TeV]   = {115,120,125,130,140,150,160,170,180,190,200,210,220,230,250,275,300,325,350,375,400,425,450,475,500,575,600,650,700,800,900,950,1000};
double mHVBFVal7TeV[nVBFPoints7TeV] = {115,120,125,130,140,150,160,170,180,190,200,210,220,230,250,275,300,325,350,375,400,425,450,475,500,575,600,650,700,800,900,950,1000};

const int nVBFPoints8TeV = 48;
int VBFmasses8TeV[nVBFPoints8TeV]   = {116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,135,140,145,150,160,170,180,190,200,220,250,275,300,325,350,375,400,425,450,475,500,525,550,575,600,650,700,750,800,850,900,950,1000};
double mHVBFVal8TeV[nVBFPoints8TeV] = {116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,135,140,145,150,160,170,180,190,200,220,250,275,300,325,350,375,400,425,450,475,500,525,550,575,600,650,700,750,800,850,900,950,1000};

//The number and values of mass points for which you have VH trees, for 7 and 8 TeV
const int nVHPoints7TeV = 10;
int VHmasses7TeV[nVHPoints7TeV]   = {110,115,120,125,126,130,140,150,160,180}; // FIXME 200 is missing at the moment
double mHVHVal7TeV[nVHPoints7TeV] = {110,115,120,125,126,130,140,150,160,180}; // 

const int nVHPoints8TeV = 11;
int VHmasses8TeV[nVHPoints8TeV]   = {110,115,120,125,126,130,140,150,160,180,200};
double mHVHVal8TeV[nVHPoints8TeV] = {110,115,120,125,126,130,140,150,160,180,200};

double filter_eff_WH_8TeV = 0.0105;
double filter_eff_ZH_8TeV = 0.028;
double filter_eff_ttH_8TeV = 0.0296;

//--------------------


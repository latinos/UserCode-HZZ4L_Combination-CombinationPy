//----------> SET INPUT VARIABLES HERE

// Input trees


// With high mass reweights + new MELA (analytical background) + superMELA
//TString filePath7TeV = "/afs/cern.ch/user/b/bonato/work/PhysAnalysis/HZZ4L/Trees_191012_M126/PRODFSR_7TeV/";
TString filePath7TeV = "root://lxcms02//data/Higgs/rootuplesOut/130702b/PRODFSR/";
TString filePath8TeV = "root://lxcms02//data/Higgs/rootuplesOut/130702b/PRODFSR_8TeV/";
TString filePath7TeVPS = "root://lxcms02//data/Higgs/rootuplesOut/130702b/JHU/";
TString filePath8TeVPS = "root://lxcms02//data/Higgs/rootuplesOut/130702b/JHU_8TeV/";


//--- Flags to control re-computation of KD
bool usePowhegTemplate=false;  // false use analytic bg
bool withPt_ = false;          // Include pT in KD
bool withY_  = false;          //    "    Y  "  "
int sqrts    = 7;              // sqrts, used only for withPt_/withY_


////////////////////////////////////
//--- Really important params --- //
const int mH=126;
const float mzzCutLow=106;
const float mzzCutHigh=141;
int useSqrts=1;              //0=use 7+8TeV; 1=use 7TeV only, 2 use 8TeV only
int altSignal =3; //anything <=2 -> don't do any alternative model templates; 3=0- , 4=2m+, 5=0h+, 6=1+, 7=1-, 8= qq->2+m
TString destDirBase = "../../CreateDatacards/templates2D";
TString destDir; //it must already exist !
const float kdCut=-1.0; //if negative, it is deactivated
//-----


bool extendToHighMass = false; // Include signal samples above 600 GeV
float highMzz=(extendToHighMass?1000:800);
float mBinSize=2.;
string str_mh="8TeV";


//-- binning of 2D template

const int nbinsX=21;
float binsX[nbinsX+1]={0.000, 0.030, 0.060, 0.100, 0.200, 0.300, 0.400, 0.500, 0.550, 0.600, 
		       0.633, 0.666, 0.700, 0.733, 0.766, 0.800, 0.833, 0.866, 0.900, 0.933,
		       0.966, 1.000};
const int nbinsYps=25;
float binsYps[nbinsYps+1]={0.000, 0.100, 0.150, 0.200, 0.233, 0.266, 0.300, 0.333, 0.366, 0.400, 
		       0.433, 0.466, 0.500, 0.533, 0.566, 0.600, 0.633, 0.666, 0.700, 0.733, 
		       0.766, 0.800, 0.850, 0.900, 0.950, 1.000};

const int nbinsYgrav=29;
float binsYgrav[nbinsYgrav+1]={0.000, 0.100, 0.150, 0.175 , 0.200, 0.225, 0.250, 0.275, 0.300, 0.325, 
		       0.350, 0.375, 0.400, 0.425 , 0.450, 0.475, 0.500, 0.525, 0.575, 0.600, 
		       0.633, 0.666, 0.700, 0.733 , 0.766, 0.800, 0.850, 0.900, 0.950, 1.000};




// Luminosity, as float and as string to be used in file names, etc.
double lumi7TeV = 5.051;
double lumi8TeV = 19.79;
TString lumistr7TeV = "5.051";
TString lumistr8TeV = "19.79";


// Location of output root files containing data events
TString DataRootFilePath = "../CreateDatacards/CMSdata/"; 
//<----------

const int nPoints7TeV = 22;
int masses7TeV[nPoints7TeV]   = {120,/*125,*/130,140,150,160,170,180,200,210,220,250,300,325,350,400,425,450,475,525,550,575,600}; //FIXME: weights for 125 are wrong.
double mHVal7TeV[nPoints7TeV] = {120,/*125,*/130,140,150,160,170,180,200,210,220,250,300,325,350,400,425,450,475,525,550,575,600};

const int nPoints8TeV = 29;
int masses8TeV[nPoints8TeV]   = {115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,145,150,180,200,250,300,325,350,400,450,500,550,600};
double mHVal8TeV[nPoints8TeV] = {115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,145,150,180,200,250,300,325,350,400,450,500,550,600};


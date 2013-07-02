//
// pre-load HiggsCSandWidth for compiled macros, and set CMSSW include path.
//
{
  gSystem->Load("libHiggsHiggs_CS_and_Width.so");
  gROOT->LoadMacro("$CMSSW_BASE/src/Higgs/Higgs_CS_and_Width/include/HiggsCSandWidth.h+");
//   gSystem->Load("libZZMatrixElementMELA.so");
//   gROOT->LoadMacro("$CMSSW_BASE/src/ZZMatrixElement/MELA/interface/Mela.h+");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/ ");  
}

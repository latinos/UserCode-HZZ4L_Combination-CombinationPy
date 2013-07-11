void createAllTemplatesPtMWithSlices(TString channel = "all", TString type = "0")
{
  gROOT->ProcessLine(".L makeTemplatesPtSyst.C+");  

  gROOT->ProcessLine("makeTemplatesPtSyst(\"" + channel + "\",0," + type + ",true)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"" + channel + "\",-2," + type + ",true)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"" + channel + "\",-3," + type + ",true)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"" + channel + "\",-5," + type + ",true)");  
  gROOT->ProcessLine("makeTemplatesPtSyst(\"" + channel + "\",-4," + type + ",true)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"" + channel + "\",-6," + type + ",true)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"" + channel + "\",-7," + type + ",true)");
  // gROOT->ProcessLine("makeTemplatesPtSyst(\"" + channel + "\",1," + type + ",true)");
  // gROOT->ProcessLine("makeTemplatesPtSyst(\"" + channel + "\",3," + type + ",true)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"" + channel + "\",4," + type + ",true)");
  // gROOT->ProcessLine("makeTemplatesPtSyst(\"" + channel + "\",5," + type + ",true)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"" + channel + "\",6," + type + ",true)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"" + channel + "\",7," + type + ",true)");

  gROOT->ProcessLine("makeTemplatesPtSyst(\"" + channel + "\",0," + type + ",false)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"" + channel + "\",-2," + type + ",false)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"" + channel + "\",-3," + type + ",false)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"" + channel + "\",-5," + type + ",false)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"" + channel + "\",-4," + type + ",false)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"" + channel + "\",-6," + type + ",false)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"" + channel + "\",-7," + type + ",false)");
  // gROOT->ProcessLine("makeTemplatesPtSyst(\"" + channel + "\",1," + type + ",false)");
  // gROOT->ProcessLine("makeTemplatesPtSyst(\"" + channel + "\",3," + type + ",false)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"" + channel + "\",4," + type + ",false)");
  // gROOT->ProcessLine("makeTemplatesPtSyst(\"" + channel + "\",5," + type + ",false)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"" + channel + "\",6," + type + ",false)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"" + channel + "\",7," + type + ",false)"); 

}

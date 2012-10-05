/*
This script checks the if .root files produced by the combination tool are OK or corrupted. 
Usage:
root -l -q "checkTrees.C+(\"TYPE\", \"inMASSFILE\", \"outMASSFILE\")"
where TYPE can be ASCLS, PLP, PLPE etc.
      inMASSFILE is the .txt file with the mass points
      outMASSFILE is the .txt file which is produced with the mass points to be reprocessed because of file corruption
      (note that the script assumes you are working with a directory structure as the default created by the tool: MASS/higgsCombineHZZL_TYPE.mHMASS.root)
*/


#include <iostream>
#include <fstream>
#include <string>

#include "TFile.h"
#include "TString.h"

using namespace std;

void checkTrees(TString type = "ASCLS", string inMassFile = "masses.txt", string outMassFile = "massesFailed.txt")
{
  TString typeExtend;
  if (type == "ASCLS") typeExtend = "ASCLS.Asymptotic";
  else if (type == "PLP") typeExtend = "PLP.ProfileLikelihood"; 
  else if (type == "PLPE") typeExtend = "PLPE.ProfileLikelihood";

  int nOK(0), nCorrupt(0);
  char c;
  int nLines = 0;
  ifstream myfile(inMassFile.c_str());
  ofstream outFile(outMassFile.c_str());
  int b = 0;

  string str;
  vector<string> str_vec;

  if(!myfile)
    {
      cout<<"Error opening output file"<<endl;
    }
  while(!myfile.eof())
    {
      getline(myfile,str,'\n');
      if (str == "") break;
      //      cout<<b<<" / "<<str<<endl;
      str_vec.push_back(str);
      ++b;
    }
  
  const int nMassPoints(str_vec.size());

  for (int j=0; j<nMassPoints; ++j){
    TString mass(str_vec[j]);
    TFile f(mass+"/higgsCombineHZZ4L_"+typeExtend+".mH"+mass+".root");
    cout<<"Checking mass point "<<mass<<endl;
    if (f.Get("limit")) {
      cout<<"  Tree is ok!"<<endl;
      nOK++;
    }
    else {
      cout<<"  Tree is corrupted..."<<endl;
      outFile << mass << endl;
      nCorrupt++;
    }
  }
  outFile.close();
  
  
  std::cout<<endl<<"Number of mass points: "<<nMassPoints<<std::endl;
  std::cout<<"Number of files corrupted: "<<nCorrupt<<" (see "+(TString)outMassFile+")"<<std::endl;
  
}

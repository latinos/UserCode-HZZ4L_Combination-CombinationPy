//
// Fill fragments with numbers from the wiki.
//
// usage:
// root -q -b zjetRates.cc
//

void zjetRates(int sqrts);


void zjetRates(){
  //  zjetRates(7);
  zjetRates(8);
}



void zjetRates(int sqrts){

  float dijetFraction = 0.2;

  if (sqrts==7) {
    cout << "7TeV  not implemented" << endl;
    return;
  } else {
    //  From wiki:
    // 4e 	6.09 +/- 0.61
    // 4mu 	3.09 +/- 0.45
    // 2mu2e 	6.91 +/- 0.62
    // 2e2mu 	2.34 +/- 0.36

    TString fs[] = {"4e","4mu","2e2mu"};
    float y[]  = {6.09, 3.09,  6.91+2.34};
    //  float ey[] = {0.53, 0.23, sqrt(0.61*0.61+0.17*0.17)};  

    for (int ifs=0; ifs<3; ++ifs) {

      TString outCardName = "CardFragments/zjetRate_";
      outCardName = outCardName + (long) sqrts + "TeV_" + fs [ifs];

      ofstream ofsCard;
      ofsCard.open((outCardName+".txt").Data(),fstream::out);    
      ofsCard << "rate zjets " << y[ifs] << endl;
    
      ofstream ofsCard0;
      ofsCard0.open((outCardName+"_0.txt").Data(),fstream::out);    
      ofsCard0 << "rate zjets " << y[ifs]*(1-dijetFraction) << endl;
    
      ofstream ofsCard1;
      ofsCard1.open((outCardName+"_1.txt").Data(),fstream::out);    
      ofsCard1 << "rate zjets " << y[ifs]*dijetFraction << endl;

    }
  }
}

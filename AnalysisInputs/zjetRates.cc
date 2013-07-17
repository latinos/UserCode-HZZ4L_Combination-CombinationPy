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
    // 4e 	6.15 +/- 0.53
    // 4mu 	3.02 +/- 0.23
    // 2mu2e 	7.06 +/- 0.61
    // 2e2mu 	2.00 +/- 0.17

    TString fs[] = {"4e","4mu","2e2mu"};
    float y[]  = {6.15, 3.02, 7.06+2.00};
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

{
	gSystem->AddIncludePath("-I$ROOFITSYS/include");
	gSystem->Load("HZZ4LRooPdfs_cc.so");
	//gSystem->CompileMacro("signalFitsbatch2.C");
	gROOT->LoadMacro("signalFitsCB.C++");		
	signalFitsCB(3);
	signalFitsCB(1);
	signalFitsCB(2);

}

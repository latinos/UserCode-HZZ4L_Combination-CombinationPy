void newTreeWithAbsError(){
	addBranch("CMSdata_bak/hzz2e2mu_19.63.root", "CMSdata/hzz2e2mu_19.63.root");
	addBranch("CMSdata_bak/hzz4mu_19.63.root", "CMSdata/hzz4mu_19.63.root");
	addBranch("CMSdata_bak/hzz4e_19.63.root", "CMSdata/hzz4e_19.63.root");
	addBranch("CMSdata_bak/hzz2e2mu_5.051.root", "CMSdata/hzz2e2mu_5.051.root");
	addBranch("CMSdata_bak/hzz4mu_5.051.root", "CMSdata/hzz4mu_5.051.root");
	addBranch("CMSdata_bak/hzz4e_5.051.root", "CMSdata/hzz4e_5.051.root");
}
void addBranch(TString oldf, TString newf){
	TFile *fold = new TFile(oldf);
	if(fold==NULL) {cout<<oldf<<" not exist "<<endl; return;}
	TTree * oldt =  (TTree*) fold->Get("data_obs");
	//oldt -> SetName("old");
	if(oldt==NULL) {cout<<" tree data_obs is not there in "<<oldf<<endl; return;}
	TFile *fnew = new TFile(newf, "RECREATE");
	TTree * newt = (TTree*)oldt -> CloneTree();

	Double_t ae;
	TBranch *newBranch = newt->Branch("CMS_zz4l_massRelErr", &ae, "CMS_zz4l_massRelErr/D");

	Double_t CMS_zz4l_mass, CMS_zz4l_massRelErr;
	TBranch        *b_CMS_zz4l_mass;   //!
	TBranch        *b_CMS_zz4l_massRelErr;   //!
	oldt->SetBranchAddress("CMS_zz4l_mass", &CMS_zz4l_mass, &b_CMS_zz4l_mass);
	oldt->SetBranchAddress("CMS_zz4l_massErr", &CMS_zz4l_massRelErr, &b_CMS_zz4l_massRelErr);
	

	//read the number of entries in the t3
	Long64_t nentries = oldt->GetEntries();

	for (Long64_t i = 0; i < nentries; i++){
		oldt->GetEntry(i);
		ae = CMS_zz4l_massRelErr / CMS_zz4l_mass ;
		newBranch->Fill();
	}
	// save only the new version of the tree
	fnew->cd();
	newt->Write();
	fnew->Close();
}

void  tree2() 
{
	TFile *ipf = new TFile("tree.root");
	TTree *tree = (TTree*)ipf->Get("tree");

	Double_t x;
	Double_t e;
	int pid;
	Double_t tof, ctof;
	Double_t tu, td;
	Double_t qu, qd;

	tree->SetBranchAddress("x", &x);
	tree->SetBranchAddress("e",&e);
	tree->SetBranchAddress("tof",&tof);
	tree->SetBranchAddress("ctof", &ctof);
	tree->SetBranchAddress("pid", &pid);
	tree->SetBranchAddress("tu", &tu);
	tree->SetBranchAddress("td", &td);
	tree->SetBranchAddress("qu", &qu);
	tree->SetBranchAddress("qd", &qd);

	//TH1D *tdiff = new TH1D("tdiff","td-tu")
	//TH2D *h = new TH2D("hgtofcx","corrected TOF", 100, -120, 120,100, 3, 3.7);
	//TH1D *h = new TH1D("htofc", "htof", 200, 0, 20);

	//calibration parameters
	Double_t a1 = 3.7439;
	Double_t b1 = -56.4583;
	Double_t a2 = 189.7533;
	Double_t b2 = -7.8178;
	
	Double_t ntof,tx, qx, ce;   

	TFile *opf = new TFile("tree2.root","recreate");
	TTree *opt = new TTree("mytree","tree");
	opt->Branch("x",&x,"x/D");
	opt->Branch("e",&e,"e/D");
	opt->Branch("tof",&tof,"tof/D");
	opt->Branch("ctof",&ctof,"cof/D");
	opt->Branch("pid",&pid,"pid/I");
	opt->Branch("tu",&tu,"tu/D");
	opt->Branch("td",&td,"td/D");
	opt->Branch("qu",&qu,"qu/D");
	opt->Branch("qd",&qd,"qd/D");
	opt->Branch("ntof",&ntof,"ntof/D");
	opt->Branch("tx",&tx,"tx/D");
	opt->Branch("tx",&qx,"qx/D");
	opt->Branch("ce",&ce,"ce/D");


	Long64_t nentries = tree->GetEntries();//get all entries
	for (Long64_t jentry = 0; jentry < nentries; jentry++)
	{
		//loop all entry
		tree->GetEntry(jentry); 

		Double_t ctx = 3.7439*( td-tu -15.08);  //cm
		Double_t d = TMath::Sqrt(502.5*502.5+ctx*ctx);
		ntof = (ctof-26.2117)/d*100; //normalized to 100cm situation

		tx = a1*(td-tu)+b1;
		qx = a2*(TMath::Log(qu/qd))+b2;
		if(pid==1) ce = TMath::Sq(72/ntof);//neutron
		else ce=1.; //gamma

		opt->Fill();//fill new parameter to TTree *opt

		if (jentry % 10000 == 0) 
		{
			cout << "process" << jentry << "of" << nentries  << endl;
			//cout<<"==="<<tx<<"====="<<x<<"Pid="<<pid<<"     "<<"ce="<<ce<<endl;
		}

	}
	
	opt->Write();
	opf->Close();
	
	
	/*
	TCanvas *c1= new TCanvas("c1","hgtofcx");
	hgtofcx->GetXaxis()->SetTitle("cm");
	hgtofcx->GetYaxis()->SetTitle("ns");
	hgtofcx->Draw("colz");

	TCanvas *c2 = new TCanvas("c2","htofc");
	c2->SetLogy();
	htofc->GetXaxis()->SetTitle("ns");
	htofc->GetYaxis()->SetTitle("counts");
	htofc->Draw();

	TCanvas *c3 = new TCanvas("c3","projx of hdx");
	TH1D *hdx1 = hdx->ProjectionX("projx");
	hdx1->Draw();
	hdx1->Fit("gaus","","");
	*/
}
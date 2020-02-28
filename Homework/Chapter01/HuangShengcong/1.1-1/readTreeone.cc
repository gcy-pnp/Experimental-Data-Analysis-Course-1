TH1D *tdiff;

void readTreeone()
{
	TFile *ipf = new TFile("tree.root");
	TTree *tree = (TTree*)ipf->Get("tree");

	Double_t x;
	Double_t e;
	int pid;
	Double_t tof, ctof;
	Double_t tu, td;
	Double_t qu, qd;

	tree->SetBranchAddress("ctof", &ctof);
	tree->SetBranchAddress("pid", &pid);
	tree->SetBranchAddress("tu", &tu);
	tree->SetBranchAddress("td", &td);
	tree->SetBranchAddress("qu", &qu);
	tree->SetBranchAddress("qd", &qd);

	
	tdiff = new TH1D("tdiff","td-tu",140,-20,50);
	TH1D *dtd = new TH1D("dtd","dt/dx",140,-20,50);


	//new data in new rootfile
	//// //calibration parameters
	 Double_t a,b;
	
	//// ... ... ...
	//// //new tree parameters
	//// Double_t tx,qx,ce;
	//// ... ... ...
	Double_t tx;

	TFile *opf = new TFile("tree2.root","recreate");
	 TTree *opt = new TTree("tree","tree");
	opt->Branch("tx", &tx,"tx/D");
	opt->Branch("tu", &tu,"tu/D");
	opt->Branch("ctof", &ctof, "ctof/D");
	opt->Branch("pid", &pid, "pid/I");
	opt->Branch("qu", &qu, "qu/D");
	opt->Branch("qd", &qd, "qd/D");

	//// ... ... ...
	//opt->Print();

	Long64_t nentries = tree->GetEntries();//get all entries
	for (Long64_t jentry = 0; jentry < nentries; jentry++)
	{
		//loop all entry
		tree->GetEntry(jentry); 

		tdiff->Fill(td-tu);

		//// //calculate new parameters
		 tx = tu-td;
		//// ... ... ...

		 opt->Fill();//fill new parameter to TTree* opt
		
		if (jentry % 10000 == 0) cout << "process" << jentry << "of" << nentries << endl;

	}
	
	TCanvas *c1= new TCanvas("c1","tdiff");
	tdiff->Draw();

	for(int i=1; i<tdiff->GetNbinsX(); i++)
	{
		Double_t df = tdiff->GetBinContent(i+1) - tdiff->GetBinContent(i);
		dtd->Fill(tdiff->GetBinLowEdge(i+1),df);  //the X axis and Y axis
	}

	TCanvas *c2 = new TCanvas("c2","dtd");
	dtd->Sumw2(0);
	dtd->Draw();

    opt->Print();
	//ipf->Close(); //can't run
	opt->Write();
	opf->Close();
}
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

	
	TH1D *hqx = new TH1D("hqx","log(qu/qd)",200,-5,5);
	TH1D *lqd = new TH1D("lqd","lq/dx",200,-5,5);

	const Double_t Lambda = 380.; //cm

	Long64_t nentries = tree->GetEntries();//get all entries
	for (Long64_t jentry = 0; jentry < nentries; jentry++)
	{
		//loop all entry
		tree->GetEntry(jentry); 

		//Double_t qx = Lambda/2*TMath::Log(qu/qd);
		 Double_t qx = TMath::Log(qu/qd);
		
		hqx->Fill(qx);

		if (jentry % 10000 == 0) cout << "process" << jentry << "of" << nentries << endl;
	}
	
	TCanvas *c1= new TCanvas("c1","hqx");
	//c1->SetLogy();
	hqx->GetXaxis()->SetTitle("MeV");
	hqx->Draw();

	for(int i=1; i<hqx->GetNbinsX(); i++)
	{
		Double_t df = hqx->GetBinContent(i+1) - hqx->GetBinContent(i);
		lqd->Fill(hqx->GetBinLowEdge(i+1),df);  //the X axis and Y axis
	}

	TCanvas *c2 = new TCanvas("c2","lqd");
	lqd->Sumw2(0);
	lqd->GetXaxis()->SetTitle("MeV");
	lqd->Draw();
	
}
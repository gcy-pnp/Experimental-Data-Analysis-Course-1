
void readTreethree() // TOF  absolutely calibration
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
	tree->SetBranchAddress("ctof", &ctof);
	tree->SetBranchAddress("pid", &pid);
	tree->SetBranchAddress("tu", &tu);
	tree->SetBranchAddress("td", &td);
	tree->SetBranchAddress("qu", &qu);
	tree->SetBranchAddress("qd", &qd);

	
	TH2D *hgtofx = new TH2D("hgtofx","hgtofx",100,-120, 120, 100, 39, 48);
	TH1D *hgctof = new TH1D("hgctof","hgctof",100, 39, 48);


	Long64_t nentries = tree->GetEntries();//get all entries
	for (Long64_t jentry = 0; jentry < nentries; jentry++)
	{
		//loop all entry
		tree->GetEntry(jentry); 

		Double_t ctx = 3.7439*( td-tu -15.08);  //cm
		
		//if(ctof> 0&& ctof< 100)
	//	{
			hgtofx->Fill(ctx,ctof);
			if(ctx<5&&ctx>-5) hgctof->Fill(ctof);
	//	}
		
		if (jentry % 10000 == 0) cout << "process" << jentry << "of" << nentries  << endl;

	}
	
	TCanvas *c1= new TCanvas("c1","hgtofx");
	hgtofx->Draw("colz");

	TCanvas *c2 = new TCanvas("c2","hgctof");
	hgctof->Draw();
	hgctof->Fit("gaus");

/*
	TCanvas *c3 = new TCanvas("c3","projx of hdx");
	TH1D *hdx1 = hdx->ProjectionX("projx");
	hdx1->Draw();
	hdx1->Fit("gaus","","");
*/
}
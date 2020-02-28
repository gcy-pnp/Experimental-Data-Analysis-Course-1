
void readTreetwo() // make a residual graph for the relationship of x and td-tu
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

	
	TH1D *htx = new TH1D("htx","htx",500,-120, 120);
	TH2D *hdx = new TH2D("hdx","hx:htx-hx",100, -20, 20, 500, -120, 120);


	Long64_t nentries = tree->GetEntries();//get all entries
	for (Long64_t jentry = 0; jentry < nentries; jentry++)
	{
		//loop all entry
		tree->GetEntry(jentry); 

		Double_t ctx = 3.7439*( td-tu -15.08);  //cm
		
		htx->Fill(ctx);
		hdx->Fill(ctx-x, x);


		if (jentry % 10000 == 0) cout << "process" << jentry << "of" << nentries <<"/n"<< ctx<<"---------------"<< x << endl;

	}
	
	TCanvas *c1= new TCanvas("c1","htx");
	htx->Draw();

	TCanvas *c2 = new TCanvas("c2","hdx");
	hdx->Draw("colz");

	TCanvas *c3 = new TCanvas("c3","projx of hdx");
	TH1D *hdx1 = hdx->ProjectionX("projx");
	hdx1->Draw();
	hdx1->Fit("gaus","","");

}
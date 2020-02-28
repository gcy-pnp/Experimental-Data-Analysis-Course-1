
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

	
	TH1D *hqx = new TH1D("hqx","hqx",500,-150, 150);
	TH2D *hdx = new TH2D("hdx","hx:hqx-hx",100, -120, 120, 500, -120, 120);


	Long64_t nentries = tree->GetEntries();//get all entries
	for (Long64_t jentry = 0; jentry < nentries; jentry++)
	{
		//loop all entry
		tree->GetEntry(jentry); 
		Double_t temp = TMath::Log(qu/qd);
		Double_t cqx = 189.7533*temp - 7.8178;  //cm
		
		hqx->Fill(cqx);
		hdx->Fill(cqx-x, x);


		if (jentry % 10000 == 0) cout << "process" << jentry << "of" << nentries  << endl;

	}
	
	TCanvas *c1= new TCanvas("c1","hqx");
	hqx->GetXaxis()->SetTitle("cm");
	hqx->Draw();

	TCanvas *c2 = new TCanvas("c2","hdx");
	hdx->GetXaxis()->SetTitle("cm");
	hdx->GetYaxis()->SetTitle("cm");
	hdx->Draw("colz");

	TCanvas *c3 = new TCanvas("c3","projx of hdx");
	TH1D *hdx1 = hdx->ProjectionX("projx");
	hdx1->GetXaxis()->SetTitle("cm");
	hdx1->Draw();
	hdx1->Fit("gaus","","");

}
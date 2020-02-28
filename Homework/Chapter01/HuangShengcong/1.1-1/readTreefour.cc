
void readTreefour() // TOF  absolutely calibration
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

	
TH2D *hgtofcx = new TH2D("hgtofcx","corrected TOF", 100, -120, 120,100, 15, 19);
TH1D *htofc = new TH1D("htofc", "htof", 200, 0, 100);

TH2D *rg = new TH2D("rg","tofc-tof:tof",500,0,120,100,-0.000002,0.000002);//residual graph

	Long64_t nentries = tree->GetEntries();//get all entries
	for (Long64_t jentry = 0; jentry < nentries; jentry++)
	{
		//loop all entry
		tree->GetEntry(jentry); 

		Double_t ctx = 3.7439*( td-tu -15.08);  //cm
		Double_t d = TMath::Sqrt(502.5*502.5+ctx*ctx);
		Double_t tofc = (ctof-26.2117)/d*500; //normalized to 500cm situation

		hgtofcx->Fill(ctx,tofc);//the gamma hits the center of the detector
		htofc->Fill(tofc);

		rg->Fill(tofc-tof,tof);
		
		if (jentry % 10000 == 0) cout << "process" << jentry << "of" << nentries  << endl;

	}
	
	TCanvas *c1= new TCanvas("c1","hgtofcx");
	hgtofcx->GetXaxis()->SetTitle("cm");
	hgtofcx->GetYaxis()->SetTitle("ns");
	hgtofcx->Draw("colz");

	TCanvas *c2 = new TCanvas("c2","htofc");
	c2->SetLogy();
	htofc->GetXaxis()->SetTitle("ns");
	htofc->GetYaxis()->SetTitle("counts");
	htofc->Draw();

	TCanvas *c4 = new TCanvas("c4","rg");
	rg->Draw();
	

	TCanvas *c3 = new TCanvas("c3","projy of rg");
	TH1D *hdx1 = rg->ProjectionY("projy");
	hdx1->Draw();
	hdx1->Fit("gaus","","");

}
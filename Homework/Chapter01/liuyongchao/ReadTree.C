

void ReadTree()
{

	const Double_t D=500.;//cm
	const Double_t L=100.;//cm
	const Double_t dD=0.5;//cm
	const Double_t TRes=1;//ns
	const Double_t Lambda=380.;//cm
	const Double_t QRes=0.1;//
	const Double_t Vsc=7.5;//cm/ns
	const Double_t En0=50.;//MeV
	const Double_t EnRes=50.;//MeV
	const Double_t Eg0=1.;//MeV
	const Double_t Rg=0.3;//ratio of gamma

	TFile *ipf=new TFile("tree.root");
	TTree *tree=(TTree*)ipf->Get("tree");

	Double_t x;
	Double_t e;
	int ng;
	Double_t tof,ctof;
	Double_t tu,td;
	Double_t qu,qd;

	tree->SetBranchAddress("ctof",&ctof);
	tree->SetBranchAddress("ng",&ng);
	tree->SetBranchAddress("tu",&tu);
	tree->SetBranchAddress("td",&td);
	tree->SetBranchAddress("qu",&qu);
	tree->SetBranchAddress("qd",&qd);
	tree->SetBranchAddress("x",&x);
	tree->SetBranchAddress("e",&e);
	//Histogram
	TCanvas *c1=new TCanvas("c1","c1",500,500);//TOF Draw
	TH1D *hTOF=new TH1D("hTOF","Time of flight", 1000,0,150);
	TCanvas *c2=new TCanvas("c2","c2",500,500);//td-tu Draw
	TH1D *tdiff=new TH1D("tdiff","td-tu",140,-20,50);
	TCanvas *c3=new TCanvas("c3","c3",500,500);
	TH1D *dtd=new TH1D("dtd","dt/dx",140,-20,50);

	//新数据写进新的root
	Double_t a,b;
	Double_t tx,qx,ce,te;
	TFile *opf=new TFile("tree2.root","recreate");
	TTree *opt=new TTree("tree2","new tree");
	opt->Branch("tx",&tx,"tx/D");
	opt->Branch("qx",&qx,"qx/D");
	opt->Branch("ce",&ce,"ce/D");
	opt->Branch("e",&e,"e/D");
	opt->Branch("x",&x,"x/D");
	opt->Branch("ng",&ng,"ng/I");
	opt->Branch("te",&te,"te/D");
	opt->Branch("tof",&tof,"tof/D");

	//逐事件读入tree的Branch数据
	Long64_t nentries=tree->GetEntries();
	for (Long64_t jentry = 0; jentry < nentries; jentry++)
	{
		tree->GetEntry(jentry);
		hTOF->Fill(ctof);
		tdiff->Fill(td-tu);
		//calculate new parameters
		//tx=tu-td;
        //''''''''''''''
        if(jentry%100000==0) cout<<"process"<<jentry<<" of "<<nentries<<endl;
	}

	for(int i=1;i<tdiff->GetNbinsX();i++)
	{
		Double_t df=tdiff->GetBinContent(i+1)-tdiff->GetBinContent(i);
		dtd->Fill(tdiff->GetBinCenter(i),df);
	}
	
	c1->cd();
	hTOF->Draw();
	c1->Draw();
	c1->SaveAs("Time of Flight.pdf");

	c2->cd();
	tdiff->Draw();
	c2->Draw();
	c2->SaveAs("td-tu.pdf");

	c3->cd();
	dtd->Sumw2(0);
	dtd->Fit("gaus","","",-14,-9);
	dtd->Draw();
	TF1 *f1=new TF1("f1","[0]*TMath::Exp(-0.5*((x-[1])/[2])^2)",39.5,43);//手动进行高斯参数拟合
	f1->SetParameter(0,-350);
	f1->SetParameter(1,41.5);
	f1->SetParameter(2,0.5);
	dtd->Fit("f1","R+");
	dtd->Draw();
	c3->Draw();
	c3->SaveAs("dtd.pdf");

	TCanvas *c4=new TCanvas("c4","c4",500,500);
	TCanvas *c5=new TCanvas("c5","c5",500,500);
	TCanvas *c6=new TCanvas("c6","c6",500,500);
	TCanvas *c7=new TCanvas("c7","c7",500,500);
	TCanvas *c8=new TCanvas("c8","c8",500,500);
	TCanvas *c9=new TCanvas("c9","c9",500,500);


	TH1D *htx=new TH1D("htx","htx",500,-120,120);
	TH2D *hdx=new TH2D("hdx","hdx-hx:hx",100,-20,20,500,-120,120);
	TH1D *hgctof=new TH1D("hgctof","hgctof",100,39,45);
	TH2D *hgtofx=new TH2D("hgtofx","hgtofx",100,-120,120,100,40,48);
	TH2D *hgtofcx=new TH2D("hgtofcx","corrected TOF",100,-120,120,100,14,19);
	TH1D *htofc=new TH1D("htofc","htof",200,12,100);
	for (Long64_t jentry=0; jentry < nentries; jentry++)
	{
		tree->GetEntry(jentry);
		tx=3.745*(td-tu-(41.4123-12.0510)/2);
		htx->Fill(tx);
		hdx->Fill(tx-x,x);
		if(ctof>40 && ctof<46)
		{
			hgtofx->Fill(tx,ctof);
			if(abs(tx)<5) hgctof->Fill(ctof);
		}
		Double_t d=TMath::Sqrt(502.5*502.5+tx*tx);
		Double_t tofc=(ctof-26.203)/d*500.;
		hgtofcx->Fill(tx,tofc);
		htofc->Fill(tofc);
		tof=(tu+td)/2-L/Vsc-(41.4123-12.0510)/2;
		qx=Lambda/2.*log(qu/qd);
		ce=TMath::Exp(2.*L/Lambda)*qu*qd;
		if(ng==1)
		{
			te=(72.*d/tof/100)*(72.*d/tof/100);
		}	
		else{te=1;}
		opt->Fill();
	}


	c4->cd();
	htx->Draw();
	c4->Draw();
	c4->SaveAs("htx.pdf");

	c5->cd();
	hdx->Draw("colz");
	c5->Draw();
	c5->SaveAs("hdx.pdf");

	c6->cd();
	hgtofx->Draw("colz");
	c6->Draw();
	c6->SaveAs("hgtofx.pdf");

	c7->cd();
	hgctof->Draw();
	hgctof->Fit("gaus");
	c7->Draw();
	c7->SaveAs("hgctof.pdf");

	c8->cd();
	hgtofcx->Draw("colz");
	c8->Draw();
	c8->SaveAs("hgtofcx.pdf");

	c9->cd();
	c9->SetLogy();
	htofc->Draw();
	c9->Draw();
	c9->SaveAs("htofc.pdf");

	ipf->Close();
	opt->Write();
	opf->Close();

	//
	TFile *ipf2=new TFile("tree2.root");
	TTree *tree2=(TTree*)ipf2->Get("tree2");
	tree2->Print();
	TCanvas *c10=new TCanvas("c10","c10",500,500);
	c10->cd();
	tree2->Draw("ce:e>>(2000,0,200,2000,0,200)","","colz");
	c10->Draw();
	c10->SaveAs("readTree_ce_e.pdf");

	TCanvas *c11=new TCanvas("c11","c11",500,500);
	c11->cd();
	tree2->Draw("tx:x>>(1000,-120,120,1000,-120,120)","","colz");
	c11->Draw();
	c11->SaveAs("readTree_tx_x.pdf");

	TCanvas *c12=new TCanvas("c12","c12",500,500);
	c12->cd();
	tree2->Draw("qx:x>>(1000,-120,120,1000,-120,120)","","colz");
	c12->Draw();
	c12->SaveAs("readTree_qx_x.pdf");

	TCanvas *c13=new TCanvas("c13","c13",500,500);
	c13->cd();
	c13->SetLogy();
	tree2->Draw("e>>(2000,-10,200)","","");
	c13->Draw();
	c13->SaveAs("readTree_e.pdf");
	
	TCanvas *c14=new TCanvas("c14","c14",500,500);
	c14->cd();
	c14->SetLogy();
	tree2->Draw("tof>>(2000,-10,100)","","");
	c14->Draw();
	c14->SaveAs("readTree_tof.pdf");

	TCanvas *c15=new TCanvas("c15","c15",500,500);
	c15->cd();
	c15->SetLogy();
	tree2->Draw("te>>(2000,-10,200)","","");
	c15->Draw();
	c15->SaveAs("readTree_te.pdf");

	TCanvas *c16=new TCanvas("c16","c16",500,500);
	c16->cd();
	tree2->Draw("te:e>>(2000,-10,200,2000,-10,200)","","colz");
	c16->Draw();
	c16->SaveAs("readTree_te_e.pdf");


	ipf2->Close();


}

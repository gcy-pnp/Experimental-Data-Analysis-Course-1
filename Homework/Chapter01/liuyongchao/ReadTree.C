

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
	TH1D *hTOF=new TH1D("hTOF","Time of flight", 70,30,100);//飞行时间图
	TCanvas *c2=new TCanvas("c2","c2",500,500);//td-tu Draw
	TH1D *tdiff=new TH1D("tdiff","td-tu",120,-20,50);
	TCanvas *c3=new TCanvas("c3","c3",500,500);
	TH1D *dtd=new TH1D("dtd","dt/dx",120,-20,50);//时间offset修正图
	TCanvas *c17=new TCanvas("c17","c17",500,500);
	TH1D *qdiff=new TH1D("qdiff","log(qu/qd)",60,-0.8,0.8);
	TCanvas *c18=new TCanvas("c18","c18",500,500);
	TH1D *dqd=new TH1D("dqd","dlog(qu/qd)/dx",60,-0.8,0.8);


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
		qdiff->Fill(log(qu/qd));
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

	for(int i=1;i<qdiff->GetNbinsX();i++)
	{
		Double_t dq=qdiff->GetBinContent(i+1)-qdiff->GetBinContent(i);
		dqd->Fill(qdiff->GetBinCenter(i),dq);
	}
	
	c1->cd();
	c1->SetLogy();
	hTOF->Draw();
	// c1->Draw();
	c1->SaveAs("Time of Flight.pdf");

	c2->cd();
	tdiff->Draw();
	// c2->Draw();
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
	// c3->Draw();
	c3->SaveAs("dtd.pdf");

	c17->cd();
	qdiff->Draw();
	c17->SaveAs("qdiff.pdf");

	c18->cd();
	dqd->Sumw2(0);
	dqd->Fit("gaus","","",-0.8,-0.4);
	dqd->Draw();
	TF1 *f2=new TF1("f2","[0]*TMath::Exp(-0.5*((x-[1])/[2])^2)",0.4,0.7);
	f2->SetParameter(0,-430);
	f2->SetParameter(1,0.5);
	f2->SetParameter(2,0.2);
	dqd->Fit("f2","R+");
	dqd->Draw();
	// c3->Draw();
	c18->SaveAs("dqd.pdf");	




	TCanvas *c4=new TCanvas("c4","c4",500,500);
	TCanvas *c5=new TCanvas("c5","c5",500,500);
	TCanvas *c6=new TCanvas("c6","c6",500,500);
	TCanvas *c7=new TCanvas("c7","c7",500,500);
	TCanvas *c8=new TCanvas("c8","c8",500,500);
	TCanvas *c9=new TCanvas("c9","c9",500,500);
	TCanvas *c10=new TCanvas("c10","c10",500,500);
	TCanvas *c11=new TCanvas("c11","c11",500,500);
	TCanvas *c12=new TCanvas("c12","c12",500,500);
	TCanvas *c13=new TCanvas("c13","c13",500,500);


	TH1D *htx=new TH1D("htx","htx",240,-120,120);
	TH2D *hdx=new TH2D("hdx","htx-hx:hx",100,-20,20,500,-120,120);
	TH1D *hgctof=new TH1D("hgctof","hgctof",100,41,45);
	TH2D *hgtofx=new TH2D("hgtofx","hgtofx",100,-120,120,100,40,48);
	TH2D *hgtofcx=new TH2D("hgtofcx","corrected TOF",100,-120,120,100,14,19);
	TH1D *htofc=new TH1D("htofc","htof",200,12,100);
	TH2D *hqdx=new TH2D("hqdx","hqx-hx:hx",100,-100,100,100,-120,120);
	TH1D *hqx=new TH1D("hqx","hqx",150,-150,150);
	TH2D *htede=new TH2D("htede","hte-he:he",220,-20,200,80,-20,20);
	TH1D *hte=new TH1D("hte","hte",110,-20,200);


	for (Long64_t jentry=0; jentry < nentries; jentry++)
	{
		tree->GetEntry(jentry);
		tx=2.*L/(41.4123+12.0510)*(td-tu-(41.4123-12.0510)/2);//////计算tx
		htx->Fill(tx);	//绘制tx一维图
		hdx->Fill(tx-x,x);	//绘制tx-x二维图
		if(ctof>40 && ctof<46)
		{
			hgtofx->Fill(tx,ctof);
			if(abs(tx)<5) hgctof->Fill(ctof);
		}
		tof=(tu+td)/2-26.203;
		Double_t d=TMath::Sqrt(502.5*502.5+tx*tx);
		Double_t tofc=(ctof-26.203)/d*500.;
		hgtofcx->Fill(tx,tofc);
		htofc->Fill(tofc);
		qx=2.*L/(0.512125+0.548872)*(log(qu/qd)-(0.512125-0.548872)/2);//////////计算qx
		hqx->Fill(qx);//qx一位图
		hqdx->Fill(qx-x,x);//qx-x:x二维图
		// ce=TMath::Exp(2.*L/Lambda)*qu*qd;
		//////计算能量
		if(tofc>18.5)
		{
			te=(72.*d/tof/100)*(72.*d/tof/100);
		}	
		else{te=1;}
		htede->Fill(e,te-e);
		hte->Fill(te);
		opt->Fill();
	}


	c4->cd();
	htx->Draw();
	// c4->Draw();
	c4->SaveAs("htx.pdf");

	c5->cd();
	hdx->Draw("colz");
	// c5->Draw();
	c5->SaveAs("htxdx.pdf");

	c6->cd();
	hgtofx->Draw("colz");
	// c6->Draw();
	c6->SaveAs("hgtofx.pdf");

	c7->cd();
	hgctof->Draw();
	hgctof->Fit("gaus");
	// c7->Draw();
	c7->SaveAs("hgctof.pdf");

	c8->cd();
	hgtofcx->Draw("colz");
	// c8->Draw();
	c8->SaveAs("hgtofcx.pdf");

	c9->cd();
	c9->SetLogy();
	htofc->Draw();
	// c9->Draw();
	c9->SaveAs("htofc.pdf");

	c10->cd();
	hqx->Draw();
	c10->SaveAs("hqx.pdf");

	c11->cd();
	hqdx->Draw("colz");
	c11->SaveAs("hqdx.pdf");


	c12->cd();
	htede->Draw("colz");
	c12->SaveAs("htede.pdf");

	c13->cd();
	c13->SetLogy();
	hte->Draw();
	c13->SaveAs("hte.pdf");

	ipf->Close();
	opt->Write();
	opf->Close();



}

void Tree(){
	const Double_t D=500.;//cm
	const Double_t L=100.;//cm
	const Double_t dD=0.5;//cm
	const Double_t TRes=1;//ns
	const Double_t Lambda=380.;//cm
	const Double_t QRes=0.1;//
	const Double_t Vsc=7.5;//cm/ns
	const Double_t En0=100.;//MeV
	const Double_t EnRes=50.;//MeV
	const Double_t Eg0=1.;//MeV
	const Double_t Rg=0.3;//ratio of gamma

	//定义新root文件，声明新tree
	TFile *opf=new TFile("tree.root","recreate");
	TTree *opt=new TTree("tree","tree structure");
	//声明tree结构中要定义的分支
	Double_t x;//入射位置
	Double_t e;//能量
	int ng;//粒子种类,n=1,g=0
	Double_t tof,ctof;
	Double_t qu,qd;
	Double_t tu,td;

	Double_t tu_off=5.5;//time offset
	Double_t td_off=20.4;//time offset

	//将变量添加到tree当中
	opt->Branch("x",&x,"x/D");
	opt->Branch("e",&e,"e/D");
	opt->Branch("tof",&tof,"tof/D");
	opt->Branch("ctof",&ctof,"ctof/D");
	opt->Branch("ng",&ng);
	opt->Branch("tu",&tu,"tu/D");
	opt->Branch("td",&td,"td/D");
	opt->Branch("qu",&qu,"qu/D");
	opt->Branch("qd",&qd,"qd/D");

	TH1D *hctof=new TH1D("hctof","neutron time of flight",1000,0,100);
	TRandom3 *gr=new TRandom3(0);
	for(int i=0;i<100000;i++){
		x=gr->Uniform(-L,L);
		Double_t Dr=D+gr->Uniform(-0.5,0.5)*dD;
		Double_t d=TMath::Sqrt(Dr*Dr+x*x);
		if (gr->Uniform()<Rg)
		{
			ng=0;
			e=Eg0;
			tof=3.333*(d*0.01);
		}
		else{
			ng=1;
			e=gr->Gaus(En0,EnRes/2.35);
			tof=72./TMath::Sqrt(e)*(d*0.01);
		}
		tu=tof+(L-x)/Vsc+gr->Gaus(0,TRes/2.35)+tu_off;
		td=tof+(L+x)/Vsc+gr->Gaus(0,TRes/2.35)+td_off;
		ctof=(tu+td)/2.;
		hctof->Fill(ctof);
		Double_t q0=e*gr->Uniform();
		qu=q0*TMath::Exp(-(L-x)/Lambda);
		qu=gr->Gaus(qu,qu*QRes/2.35);
		qd=q0*TMath::Exp(-(L+x)/Lambda);
		qd=gr->Gaus(qd,qd*QRes/2.35);
		opt->Fill();
	} 

	hctof->Write();
	opt->Write();
	opf->Close();

	//读取tree.root文件并绘制各种图
	//c1 and c2 is the same way
	// TFile *ipf=new TFile("tree.root");
	// // ipf->ls();
	// TTree *tree=(TTree*)ipf->Get("tree");
	// tree->Print();
	// tree->Show(1);
	// tree->Show(100);
	// TCanvas *c1=new TCanvas("c1","c1",500,500);
	// c1->cd();
	// tree->Draw("td-tu>>htx(500,-20,50)");
	// c1->Draw();
	// c1->SaveAs("tree_htx_td_diff_tu.pdf");

	// TCanvas *c2=new TCanvas("c2","c2",500,500);
	// c2->cd();
	// tree->Draw("tof>>hh(1000,0,100)");
	// c2->SetLogy();
	// c2->Draw();
	// c2->SaveAs("tree_hh_tof.pdf");

	// //c3 shows a different one
	// TCanvas *c3=new TCanvas("c3","c3",500,500);
	// TH1D *hh=(TH1D*) ipf->Get("hctof");
	// c3->cd();
	// hh->Draw();
	// c3->Draw();
	// c3->SaveAs("tree_hh_hctof.pdf");

	// TCanvas *c4=new TCanvas("c4","c4",500,500);
	// c4->cd();
	// c4->SetLogy(0);
	// gStyle->SetPalette(1);
	// tree->Draw("ctof:td-tu>>(1000,-20,50,1000,0,200)","","colz");
	// c4->Draw();
	// c4->SaveAs("tree_ctof_td_diff_tu.pdf");

	// TCanvas *c5=new TCanvas("c5","c5",500,500);
	// c5->cd();
	// tree->Draw("ctof:td-tu>>(2000,-20,50,50,40,46)","","colz");
	// c5->Draw();
	// c5->SaveAs("tree_ctof_td_diff_tu_gamma.pdf");

	// TCanvas *c6=new TCanvas("c6","c6",500,500);
	// c6->cd();
	// tree->Draw("td-tu:x","","colz");
	// c6->Draw();
	// c6->SaveAs("td_diff_tu_x.pdf");

	// TCanvas *c7=new TCanvas("c7","c7",500,500);
	// c7->cd();
	// tree->Draw("log(qu/qd):x","","colz");
	// c7->Draw();
	// c7->SaveAs("tree_logqu_qd_x.pdf");
	// ipf->Close();
}

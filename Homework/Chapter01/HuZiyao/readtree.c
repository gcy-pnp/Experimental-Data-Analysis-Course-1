#include<iostream>
#include<cmath>
using namespace std;
////////////////常量
const Double_t D = 500.;//cm, distance between target and the scin.(Center)
const Double_t L = 100.;//cm, half length of the scin.
const Double_t dD = 5.;//cm, thickness of the scin.
const Double_t TRes = 1.;//ns, time resolution(FWHM) of the scintillator.
const Double_t Lambda = 380.;//cm, attenuation lenght of the scin.
const Double_t QRes = 0.1;//relative energy resolution(FWHM) of the scin. 
const Double_t Vsc = 7.5;//ns/cm, speed of light in the scin.
const Double_t En0 = 100;//MeV, average neutron energy
const Double_t EnRes = 50.;//MeV, energy spread of neutron(FWHM)
const Double_t Eg0 = 1;//MeV, gamma energy  
const Double_t Rg = 0.3;//ratio of gamma,ratio of neutron 1-Rg 
TH1D* hTOF;
TH1D* htx;
TH1D* Dhtx;
TH1D* nhtx;
TH1D* hqx;
TH1D* Dhqx;
TH1D* nhqx;

TH2D* Dhtx2;
TH1D* Dhtx1;///投影用
TH2D* Dhqx2;
TH1D* Dhqx1;
void readTree()
{
	TFile* ipf = new TFile("tree.root");	//文件写入内存
	TTree* tree = (TTree*)ipf->Get("tree");	// 得到tree的指针
	
	TFile* opf = new TFile("tree2.root", "recreate");
	TTree* opt = new TTree("tree", "tree");
	//////声明tree的branch变量
	Double_t x, tx, qx;
	Double_t e, ce;
	int pid;
	Double_t tof, ctof, ntof;
	Double_t tu, td;
	Double_t qu, qd;

	//////变量指向branch地址
	tree->SetBranchAddress("ctof", &ctof);
	tree->SetBranchAddress("pid", &pid);
	tree->SetBranchAddress("tu", &tu);
	tree->SetBranchAddress("td", &td);
	tree->SetBranchAddress("qu", &qu);
	tree->SetBranchAddress("qd", &qd);
	tree->SetBranchAddress("tof", &tof);
	tree->SetBranchAddress("x", &x);
	tree->SetBranchAddress("e", &e);

	opt->Branch("x", &x, "x/D");
	opt->Branch("tx", &tx, "tx/D");
	opt->Branch("qx", &qx, "qx/D");
	opt->Branch("e", &e, "e/D");
	opt->Branch("tof", &tof, "tof/D");
	opt->Branch("ctof", &ctof, "ctof/D");
	opt->Branch("pid", &pid, "pid/I");
	opt->Branch("tu", &tu, "tu/D");
	opt->Branch("td", &td, "td/D");
	opt->Branch("qu", &qu, "qu/D");
	opt->Branch("qd", &qd, "qd/D");
	opt->Branch("ce", &ce, "ce/D");
	opt->Branch("ntof", &ntof, "ntof/D");

	//////histogram定义////
	hTOF = new TH1D("hTOF", "time of flight(without corrected)", 1000, 0, 100);
	htx = new TH1D("htx", "tx(without calibration) =tu-td", 150, -45, 15);
	Dhtx = new TH1D("Dhtx", "dn/dtx", 150, -45. - .01, 15. - .01);
	nhtx = new TH1D("nhtx", "tx", 1000, -130, 130);
	hqx = new TH1D("hqx", "log(qu/qd)", 100, -0.8, 0.8);
	Dhqx = new TH1D("Dhqx", "dlog/dx", 100, -.8 - .001, .8 - .001);
	nhqx = new TH1D("nhqx", "qx", 1000, -130, 130);

	Dhqx2 = new TH2D("Dhqx2", "qx-x:x", 500, -100, 100, 500, -130, 130);
	Dhtx2 = new TH2D("Dhtx2", "tx-x:x", 500, -20, 20, 500, -130, 130);
	///////历遍////
	Long64_t nentries = tree->GetEntries();//得到事件总数
	for (Long64_t j = 0; j < nentries; j++)
	{
		tree->GetEntry(j);//将j号的事件branch的值填入变量中;
		htx->Fill(tu - td);
		hqx->Fill(log(qu / qd));
		if (j % 99999 == 0)cout << "Finished percent:" << j * 100. / nentries << "%" << endl;
	}
	/////////对直方图操作处理数据得到a,b/////
	
	for (int i = 1; i < htx->GetNbinsX(); i++)
	{
		Double_t df = htx->GetBinContent(i + 1) - htx->GetBinContent(i);
		Dhtx->Fill(htx->GetBinLowEdge(i + 1), df);
	}
	Dhtx->Sumw2(0);
	//Dhtx->Draw();
	Dhtx->Write();
	//htx->Draw();
	htx->Write();

	//Dhtx->Fit("gaus", "", "", -44, -39);	//			得到meant1=-4.13556e+01	cm
	
	TF1* f1 = new TF1("f1", "[0]*exp(-0.5*((x-[1])/[2])^2)", 10.5, 13.5);
	f1->SetParameter(0, -350);
	f1->SetParameter(1, 12);
	f1->SetParameter(2, 0.5);
	//Dhtx->Fit("f1", "R");						//得到meant2=1.20017e+01		cm
	
	Double_t meant1 = -4.13556e+01;
	Double_t meant2 = 1.20017e+01;
	///////计算a,b
	Double_t a = 2. * L / (meant1 - meant2);		//			= -3.7483156	cm/ns
	Double_t b = (meant2 + meant1) * L / (meant2 - meant1);//	 = -55.013841;		cm
	//////////////////////////////tx=a*(tu-td)+b
	
	///////////根据log(qu/qd)两边界获得cLambda以此确定qx
	for (int i = 1; i < hqx->GetNbinsX(); i++)
	{
		Double_t df = hqx->GetBinContent(i + 1) - hqx->GetBinContent(i);
		Dhqx->Fill(hqx->GetBinLowEdge(i + 1), df);
	}
	Dhqx->Sumw2(0);
	//Dhqx->Draw();
	Dhqx->Write();
	//hqx->Draw();
	hqx->Write();
	//Dhqx->Fit("gaus", "", "", -.75, -.35);			//得到meanq1 = -5.19224e-01

	TF1* f2 = new TF1("f2", "[0]*exp(-0.5*((x-[1])/[2])^2)", 0.37, 0.75);
	f2->SetParameter(0, -300);
	f2->SetParameter(1, 0.5);
	f2->SetParameter(2, 0.2);
	//Dhqx->Fit("f2", "R");								//得到meanq2 = 5.41140e-01

	///////计算cLambda
	Double_t meanq1 = -5.19224e-01;
	Double_t meanq2 = 5.41140e-01;
	Double_t cLambda = 4. * L / (meanq2 - meanq1);		//cm
	////////////////////////////////qx=cLambda*.5*log(qu/qd)


	/////////再次历遍得到刻度后的tx,qx
	for (Long64_t j = 0; j < nentries; j++)
	{
		tree->GetEntry(j);			//将j号的事件branch的值填入变量中;
		tx = a * (tu - td) + b;
		qx = cLambda * .5 * log(qu / qd);
		nhtx->Fill(tx);
		Dhtx2->Fill(tx - x, x);
		nhqx->Fill(qx);
		Dhqx2->Fill(qx - x, x);
		if (j % 99999 == 0)cout << "Finished percent:" << j * 100. / nentries << "%" << endl;
	}
	Dhtx1 = Dhtx2->ProjectionX("projxofDhtx2");
	Dhqx1 = Dhqx2->ProjectionX("projxofDhqx2");
	//nhtx->Draw();
	nhtx->Write();
	//Dhtx2->Draw("colz");
	Dhtx2->Write();
	//Dhtx1->Draw();
	Dhtx1->Write();
	//Dhtx1->Fit("gaus", "", "", -10, 10);			//mean=8.40724e-01 cm		sigma=2.24502e+00	cm

	//nhqx->Draw();
	nhqx->Write();;
	//Dhqx2->Draw("colz");
	Dhqx2->Write();
	//Dhqx1->Draw();
	Dhqx1->Write();
	//Dhqx1->Fit("gaus", "", "", -50, 50);			//mean=-3.32857e-02 cm		sigma=1.13365e+01	cm
	
	////////////////TOF计算ntof和中子能量ce,(只做了tx的)
	TH2D* hgctoftx;//histogram gammar caculated tof Vs tx;
	TH2D* hgctofqx;
	
	TH2D* hgctofCtx = new TH2D("hgctofCtx", "gammar ctof corrected buy tx Vs tx", 500, -130, 130, 100, 3, 3.7);
	
	TH1D* hgcctoftx;//histogram gammar centor caculated tof Divided use tx
	TH1D* hgcctofqx;
	
	TH1D* hTOFCtx = new TH1D("hTOFCtx", "TOF corrected by tx", 1000, 2, 12);
	TH1D* hceCtx = new TH1D("hceCtx", "ce corrected by tx", 1000, 0, 200);

	hgctoftx = new TH2D("hgctoftx", "gammar caculated tof Vs tx", 500, -130, 130, 100, 41, 45);
	hgcctoftx = new TH1D("hgcctoftx", "centor gammar |tx|<5cm cTOF", 100, 41, 45);

	for (Long64_t j = 0; j < nentries; j++)
	{
		tree->GetEntry(j);			//将j号的事件branch的值填入变量中;
		hTOF->Fill(ctof);
		tx = a * (tu - td) + b;
		if (ctof > 41.5 && ctof < 44.5)
		{
			hgctoftx->Fill(tx, ctof);
			if (abs(tx) <= 5) 
			{
				hgcctoftx->Fill(ctof);//中心粒子的ctof
			}
		}
		if (j % 99999 == 0)cout << "Finished percent:" << j * 100. / nentries << "%" << endl;
	}
	//hTOF->Draw();
	hTOF->Write();
	//hgctoftx->Draw("colz");
	hgctoftx->Write();
	//hgcctoftx->Draw();
	hgcctoftx->Write();
	//hgcctoftx->Fit("gaus", "", "", 42, 44);		//mean = 4.29476e+01 ns
	Double_t cTOF0 = 4.29476e+01;
	for (Long64_t j = 0; j < nentries; j++)
	{
		tree->GetEntry(j);			//将j号的事件branch的值填入变量中;
		tx = a * (tu - td) + b;
		Double_t d = sqrt(502.5 * 502.5 + tx * tx);
		ntof = (ctof - cTOF0 + 16.748) * 100. / d;
		if (ntof > 3 && ntof < 3.7)
		{
			hgctofCtx->Fill(tx, ntof);
		}
		hTOFCtx->Fill(ntof);
		if (pid == 1) {
			ce = pow(72. / ntof, 2);			//MeV
			hceCtx->Fill(ce);
		}
		else {
			ce = -1;
		}
		if (j % 99999 == 0)cout << "Finished percent:" << j * 100. / nentries << "%" << endl;
	}
	//hgctofCtx->Draw("colz");
	hgctofCtx->Write();
	//hTOFCtx->Draw();
	hTOFCtx->Write();
	//hceCtx->Draw();
	hceCtx->Write();
	//hceCtx->Fit("gaus");

	///////////////写入新文件
	for (Long64_t j = 0; j < nentries; j++)
	{
		tree->GetEntry(j);
		tx = a * (tu - td) + b;
		qx = cLambda * .5 * log(qu / qd);
		Double_t d = sqrt(502.5 * 502.5 + tx * tx);
		ntof = (ctof - cTOF0 + 16.748) * 100. / d;
		if (pid == 1) {
			ce = pow(72. / ntof, 2);			//MeV
			hceCtx->Fill(ce);
		}
		else {
			ce = -1;
		}
		opt->Fill();

	}
	ipf->Close();
	opt->Write();
	opf->Close();
}
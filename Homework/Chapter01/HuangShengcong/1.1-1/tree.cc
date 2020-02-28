void tree()
{
	const Double_t D = 500; //cm
	const Double_t L = 100; //cm
	const Double_t dD = 5; //cm ���ߵĺ��
	const Double_t TRes = 1; //��˸���ʱ��ֱ���(FWHM)
	const Double_t Lambda = 380; //cm ��˸��˥������
	const Double_t QRes = 0.1; //��˸�����������ֱ��ʣ�FWHM��
	const Double_t Vsc = 7.5; //ns/cm,��˸���еĹ���
	const Double_t En0 = 100; //MeV,ƽ����������
	const Double_t EnRes = 50; //MeV,���ӵ�������ɢ��FWHM��
	const Double_t Eg0 = 1; // MeV,٤������
	const Double_t Rg = 0.3; //�ñ��ʣ�����1-Rg����

	TFile *opf = new TFile("tree.root", "recreate");
	TTree *opt = new TTree("tree", "tree structure");

	Double_t x; //����λ��
	Double_t e; //����
	int pid; //�������࣬n:pid=1,g:pid=0
	Double_t tof, ctof; //����ʵ�ʷ���ʱ�������õ������ӷ���ʱ��
	Double_t tu, td;
	Double_t qu, qd;

	Double_t tu_off = 5.5; //ʱ��ƫ�ƣ�///PMT�Ķ�Խʱ��+�����ϵĴ���ʱ��
	Double_t td_off = 20.4; //ʱ��ƫ��

	opt->Branch("x", &x, "x/D");
	opt->Branch("e", &e, "e/D");
	opt->Branch("tof", &tof, "tof/D");
	opt->Branch("ctof", &ctof, "ctof/D");
	opt->Branch("pid", &pid, "pid/I");
	opt->Branch("tu", &tu, "tu/D");
	opt->Branch("td", &td, "td/D");
	opt->Branch("qu", &qu, "qu/D");
	opt->Branch("qd", &qd, "qd/D");

	TH1D *hctof = new TH1D("hctof", "neutron time of flight", 1000, 0, 100);
	TRandom3 *gr = new TRandom3(0);

	for (int i = 0; i < 100000; i++)
	{
		x = gr->Uniform(-L, L); //�����������Ӿ��ȷֲ�������̽��������
		Double_t Dr = D + gr->Uniform(-0.5, 0.5)*dD; //cm
		Double_t d = TMath::Sqrt(Dr*Dr + x * x);//cm ,flight path
		if (gr->Uniform() < Rg)
		{
			// is gamma
			pid = 0;
			e = Eg0;
			tof = 3.333*(d*0.01);
		}
		else
		{
			//is neutron
			pid = 1;
			e = gr->Gaus(En0, EnRes / 2.35);
			tof = 72 * d*0.01 / (TMath::Sqrt(e));
		}
		tu = tof + (L - x) / Vsc + gr->Gaus(0, TRes / 2.35) + tu_off;
		td = tof + (L + x) / Vsc + gr->Gaus(0, TRes / 2.35) + td_off;
		ctof = (tu + td) / 2; //simplified calculation
		hctof->Fill(ctof);
		
		Double_t q0 = e * gr->Uniform(); //energy of recoil proton in plas
		qu = q0 * TMath::Exp(-(L - x) / Lambda);
		qu = gr->Gaus(qu, qu*QRes / 2.35);
		qd = q0 * TMath::Exp(-(L+x)/Lambda);
		qd = gr->Gaus(qd, qd * QRes / 2.35);
		opt->Fill(); //������õı���ֵ����Tree��

	}
	hctof->Write();
	opt->Write();
	opf->Close();

}
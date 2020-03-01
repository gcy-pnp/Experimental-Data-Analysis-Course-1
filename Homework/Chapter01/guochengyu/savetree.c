TH1D *hce;
TH2D *hgntofx;
TH1D *hntof;

void savetree()
{
    const Double_t D = 500.; //cm, distance between target and the scin.(Center)
    const Double_t L = 100.; //cm, half length of the scin.
    const Double_t dD = 5.;  //cm, thickness of the scin.
    const Double_t tof0 = 3.333 * ((D + dD / 2) * 0.01);
    const Double_t p1 = -11.6999;
    const Double_t p2 = 41.5097;
    const Double_t p3 = -0.541643;
    const Double_t p4 = 0.532719;
    const Double_t Tof0 = 42.9488; //由readtree.c拟合给出

    // 1.打开文件，得到TTree指针
    TFile *ipf = new TFile("tree.root");     //打开ROOT文件
    TTree *tree = (TTree *)ipf->Get("tree"); //得到名字为“tree”的TTree指针

    //2. 声明tree的Branch变量
    Double_t x;
    Double_t e;
    int pid;
    Double_t tof, ctof;
    Double_t tu, td;
    Double_t qu, qd;
    Double_t tx, d, ntof, qx, ce;

    //3. 将变量指向对应Branch的地址
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("e", &e);
    tree->SetBranchAddress("tof", &tof);
    tree->SetBranchAddress("ctof", &ctof);
    tree->SetBranchAddress("pid", &pid);
    tree->SetBranchAddress("tu", &tu);
    tree->SetBranchAddress("td", &td);
    tree->SetBranchAddress("qu", &qu);
    tree->SetBranchAddress("qd", &qd);

    TFile *opf = new TFile("tree2.root", "recreate");
    TTree *opt = new TTree("tree2", "newtree");
    hce = new TH1D("hce", "hce", 1000, -150, 350);
    hgntofx = new TH2D("hgntofx", "corrected TOF", 200, -120, 120, 200, 0, 30);
    hntof = new TH1D("hntof", "htof", 200, 0, 40);

    //4.设置新tree的Branch
    opt->Branch("x", &x, "x/D");
    opt->Branch("e", &e, "e/D");
    opt->Branch("tof", &tof, "tof/D");
    opt->Branch("ctof", &ctof, "ctof/D");
    opt->Branch("pid", &pid, "pid/I");
    opt->Branch("tu", &tu, "tu/D");
    opt->Branch("td", &td, "td/D");
    opt->Branch("qu", &qu, "qu/D");
    opt->Branch("qd", &qd, "qd/D");
    opt->Branch("tx", &tx, "tx/D");
    opt->Branch("qx", &qx, "qx/D");
    opt->Branch("ce", &ce, "ce/D");
    opt->Branch("ntof", &ntof, "ntof/D");

    Long64_t nentries = tree->GetEntries();
    for (Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        //calculate new parameters
        tree->GetEntry(jentry);
        tx = 2 * L / (p2 - p1) * (td - tu - (p1 + p2) / 2);
        qx = 2 * L / (p4 - p3) * (TMath::Log(qu / qd) - (p4 + p3));
        d = TMath::Sqrt((D + dD / 2) * (D + dD / 2) + tx * tx);
        ntof = (ctof + tof0 - Tof0) / d * 100.;
        hgntofx->Fill(tx, ntof); //gamma hits the center of the det.
        hntof->Fill(ntof);
        if (ctof > 44.5)
        {
            ce = (72 / ntof) * (72 / ntof);
            hce->Fill(ce);
        }
        else
        {
            ce = -100;
            hce->Fill(ce);
        }
        opt->Fill(); //fill new parameter to TTree* opt
        if (jentry % 100000 == 0)
            cout << "process " << jentry << " of " << nentries << endl;
    }

    ipf->Close();
    hgntofx->Write();
    hntof->Write();
    hce->Write();
    opt->Write();
    opf->Close();
}
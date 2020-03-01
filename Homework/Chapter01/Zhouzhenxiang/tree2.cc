void tree2()
{
  TFile *ipf=new TFile("tree.root");//打开ROOT文件
  TTree *tree=(TTree*)ipf->Get("tree");//得到tree的指针
  TFile *opf=new TFile("tree2.root","recreate");
  TTree *opt=new TTree("tree2","tree structure");
  Double_t x;
  Double_t e,ce;
  int pid;
  Double_t tof, ctof,ntof ;
  Double_t tu, td, tx;
  Double_t qu, qd, qx;
  const Double_t D=500;
  const Double_t dD=5;


  tree->SetBranchAddress("x",&x);
  tree->SetBranchAddress("e",&e);
  tree->SetBranchAddress("tof",&tof);
  tree->SetBranchAddress("ctof",&ctof);
  tree->SetBranchAddress("pid",&pid);
  tree->SetBranchAddress("tu",&tu);   
  tree->SetBranchAddress("td",&td);
  tree->SetBranchAddress("qu",&qu);   
  tree->SetBranchAddress("qd",&qd);

  opt->Branch("x",&x,"x/D");
  opt->Branch("e",&e,"e/D");
  opt->Branch("ce",&ce,"ce/D");
  opt->Branch("tof",&tof,"tof/D");
  opt->Branch("ctof",&ctof,"ctof/D");
  opt->Branch("ntof",&ntof,"ntof/D");
  opt->Branch("tu",&tu,"tu/D");
  opt->Branch("td",&td,"td/D");
  opt->Branch("tx",&tx,"tx/D");
  opt->Branch("qu",&qu,"qu/D");
  opt->Branch("qd",&qd,"qd/D");
  opt->Branch("qx",&qx,"qx/D");
  opt->Branch("pid",&pid,"pid/I");


  //Histogram

  //逐事件读取tree的branch数据
  Long64_t nentries=tree->GetEntries();//得到事件总数
  Double_t lambda=381.9;

  for(Long64_t jentry=0; jentry<nentries; jentry++) {//对每个事件进行遍历
    tree->GetEntry(jentry);//将第j个事件数据填入对应变量，每次变量值会变成当前事件对应的数据。
    tx=3.737*(td-tu-14.872);
    qx=(lambda/2)*TMath::Log(qu/qd);
    ntof=(ctof-26.20)/(TMath::Sqrt(502.5*502.5+tx*tx))*100;
    ce=72.*72./(ntof*ntof);
    opt->Fill();
    //if(jentry%100000==0) cout<<"process "<<jentry<<" of "<<nentries<<endl;
  }
  opt->Write();
  opf->Close();
}

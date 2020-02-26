#include <TROOT.h> 
#include <TFile.h>
#include <TTree.h> 
#include <TH1D.h> 
#include <TH2D.h> 
#include <TF1.h> 
#include <TFitResult.h> 

#include <iostream> 
using namespace std;

void ana2()
{
  // 常量声明
  const Double_t L = 100.;    //cm, half length of the scin.

  // 1.打开文件，得到TTree指针
  TFile *ipf = new TFile("sim.root"); //打开ROOT文件
  TTree *ipt = (TTree*)ipf->Get("tree"); //得到名字为"tree"的TTree指针

  // 2. 声明tree的Branch变量, 将变量指向对应Branch的地址
  Double_t x;
  Double_t e;
  int pid;
  Double_t tof, ctof;
  Double_t tu, td;
  Double_t qu, qd;
  ipt->SetBranchAddress("x",&x);
  ipt->SetBranchAddress("e",&e);
  ipt->SetBranchAddress("pid",&pid);
  ipt->SetBranchAddress("tof",&tof);
  ipt->SetBranchAddress("ctof",&ctof);
  ipt->SetBranchAddress("tu",&tu);   
  ipt->SetBranchAddress("td",&td);
  ipt->SetBranchAddress("qu",&qu);   
  ipt->SetBranchAddress("qd",&qd);

  TFile *ipf2 = new TFile("ana1.root");
  TFitResultPtr r1 = (TFitResult*)ipf2->Get("r1");
  TFitResultPtr r2 = (TFitResult*)ipf2->Get("r2");
  Double_t txl = r1->Parameter(1);
  Double_t txr = r2->Parameter(1);
  Double_t txoff = (txl+txr)/2;
  Double_t tslope = 2*L/(txr-txl);
  cout << "txoff=" << txoff << "  tslope=" << tslope << endl; 

  //将新数据写入新的ROOT文件 -对应的代码用 ////标出
  TFile *opf = new TFile("ana2.root","recreate");
  TTree *opt=new TTree("tree","tree");

  Double_t tx,qx,ce;
  opt->Branch("tx",&tx,"tx/D");

  TH1D *htx = new TH1D("htx","htx",500,-120,120);
  TH2D *hdx = new TH2D("hdx","htx-hx:hx",100,-20,20,500,-120,120);

  // 4. 逐事件读取tree的branch数据
  Long64_t nentries = ipt->GetEntries(); //得到事件总数
  for(Long64_t jentry=0;jentry<nentries;jentry++) { //对每个事件进行遍历
    ipt->GetEntry(jentry); //将第jentry个事件数据填入对应变量，每次变量值会变成当前事件对应的数据。
    
    tx = tslope*(td-tu-txoff);
    htx->Fill(tx);
    hdx->Fill(tx-x,x);

    opt->Fill();//fill new parameter to TTree* opt

    if(jentry%10000==0) cout<<"process "<<jentry<<" of "<<nentries<<endl;
  }

  TH1D *hdx1 = hdx->ProjectionX("projx of hdx");
  TFitResultPtr r = hdx1->Fit("gaus","S");
  r->SetName("r");
  r->Write();

  opf->Write();
  ipf->Close();
  opf->Close();
}

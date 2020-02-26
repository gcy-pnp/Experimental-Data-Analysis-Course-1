#include <TROOT.h> 
#include <TFile.h>
#include <TTree.h> 
#include <TH1D.h> 
#include <TF1.h> 
#include <TFitResult.h> 

#include <iostream> 
using namespace std;

void ana1()
{
  // 1.打开文件，得到TTree指针
  TFile *ipf = new TFile("sim.root"); //打开ROOT文件
  TTree *ipt = (TTree*)ipf->Get("tree"); //得到名字为"tree"的TTree指针

  // 2. 声明tree的Branch变量, 将变量指向对应Branch的地址
  Double_t tu, td;
  ipt->SetBranchAddress("tu",&tu);   
  ipt->SetBranchAddress("td",&td);

  // 3. 输出文件
  TFile *opf = new TFile("ana1.root","recreate");
  TH1D *tdiff = new TH1D("tdiff","td-tu",140,-20,50);

  // 4. 逐事件读取tree的branch数据
  Long64_t nentries = ipt->GetEntries(); //得到事件总数
  for(Long64_t jentry=0;jentry<nentries;jentry++) { //对每个事件进行遍历
    ipt->GetEntry(jentry); //将第jentry个事件数据填入对应变量，每次变量值会变成当前事件对应的数据。

    tdiff->Fill(td-tu);

    if(jentry%10000==0) cout<<"process "<<jentry<<" of "<<nentries<<endl;
  }

  // 5. 确定均匀分布边界的方法 - 微分法
  TH1D *dtd = new TH1D("dtd","dt/dx",141,-20.25,50.25);
  for(int i=1;i<tdiff->GetNbinsX();i++) {
    Double_t df = tdiff->GetBinContent(i+1)-tdiff->GetBinContent(i);
    dtd->Fill(tdiff->GetBinLowEdge(i+1),df);
  }
  dtd->Sumw2(0);
  TFitResultPtr r1 = dtd->Fit("gaus","S","",-14,-9); //txl
  r1->SetName("r1");

  TF1 *f1 = new TF1("f1","[0]*TMath::Exp(-0.5*((x-[1])/[2])^2)",39.5,43);
  f1->SetParameter(0,-350);
  f1->SetParameter(1,41.5);
  f1->SetParameter(2,0.5);
  TFitResultPtr r2 = dtd->Fit("f1","RS+");
  r2->SetName("r2");

  // 6. 存出文件
  tdiff->Write();
  dtd->Write();
  r1->Write();
  r2->Write();

  ipf->Close();
  opf->Close();
}

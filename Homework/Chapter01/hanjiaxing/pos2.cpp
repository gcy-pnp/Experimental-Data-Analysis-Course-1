// 位置刻度
#include <TROOT.h> 
#include <TFile.h>
#include <TTree.h> 
#include <TH1D.h> 
#include <TH2D.h> 
#include <TF1.h> 
#include <TFitResult.h> 

#include <iostream> 
using namespace std;

void pos2()
{
  // 常量声明
  const Double_t L = 100.; //cm, half length of the scin.

  // Input file 1
  TFile *ipf = new TFile("sim.root");
  TTree *ipt = (TTree*)ipf->Get("tree");

  Double_t x;
  Double_t tu, td;
  ipt->SetBranchAddress("x",&x);
  ipt->SetBranchAddress("tu",&tu);   
  ipt->SetBranchAddress("td",&td);

  // Input file 2
  TFile *ipf2 = new TFile("pos1.root");
  TFitResultPtr r1 = (TFitResult*)ipf2->Get("r1");
  TFitResultPtr r2 = (TFitResult*)ipf2->Get("r2");
  Double_t txl = r1->Parameter(1);
  Double_t txr = r2->Parameter(1);
  Double_t txoff = (txl+txr)/2;
  Double_t tslope = 2*L/(txr-txl);
  cout << "txoff=" << txoff << "  tslope=" << tslope << endl; 

  // Output file
  TFile *opf = new TFile("pos2.root","recreate");
  TTree *opt = new TTree("tree","tree");

  Double_t tx;
  opt->Branch("tx",&tx,"tx/D");

  TH1D *htx = new TH1D("htx","htx",500,-120,120);
  TH2D *hdx = new TH2D("hdx","htx-hx:hx",100,-20,20,500,-120,120);

  // Loop
  Long64_t nentries = ipt->GetEntries();
  for(Long64_t jentry=0;jentry<nentries;jentry++) {
    ipt->GetEntry(jentry);
    if(jentry%10000==0) cout<<"process "<<jentry<<" of "<<nentries<<endl;
    
    tx = tslope*(td-tu-txoff);
    htx->Fill(tx);
    hdx->Fill(tx-x,x);

    opt->Fill();
  }

  // 评估位置刻度结果
  TH1D *hdx1 = hdx->ProjectionX("projx of hdx");
  TFitResultPtr r = hdx1->Fit("gaus","S");
  r->SetName("r");
  r->Write();

  // Write and close
  opf->Write();
  ipf->Close();
  opf->Close();
}

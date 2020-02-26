// 微分法确定td-tu均匀分布边界
#include <TROOT.h> 
#include <TFile.h>
#include <TTree.h> 
#include <TH1D.h> 
#include <TF1.h> 
#include <TFitResult.h> 

#include <iostream> 
using namespace std;

void pos1()
{
  // Input file
  TFile *ipf = new TFile("sim.root");
  TTree *ipt = (TTree*)ipf->Get("tree");

  Double_t tu, td;
  ipt->SetBranchAddress("tu",&tu);   
  ipt->SetBranchAddress("td",&td);

  // Output file
  TFile *opf = new TFile("pos1.root","recreate");
  TH1D *tdiff = new TH1D("tdiff","td-tu",140,-20,50);

  // Loop
  Long64_t nentries = ipt->GetEntries();
  for(Long64_t jentry=0;jentry<nentries;jentry++) {
    ipt->GetEntry(jentry);
    if(jentry%10000==0) cout<<"process "<<jentry<<" of "<<nentries<<endl;

    tdiff->Fill(td-tu);
  }

  // 微分法确定td-tu均匀分布边界
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

  // Write and close
  tdiff->Write();
  dtd->Write();
  r1->Write();
  r2->Write();
  ipf->Close();
  opf->Close();
}

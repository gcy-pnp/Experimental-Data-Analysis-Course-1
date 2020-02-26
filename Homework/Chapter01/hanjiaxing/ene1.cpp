// TOF 绝对刻度 1
#include <TROOT.h> 
#include <TFile.h>
#include <TTree.h> 
#include <TH1D.h> 
#include <TH2D.h> 
#include <TF1.h> 
#include <TFitResult.h> 
#include <TMath.h> 

#include <iostream> 
using namespace std;

void ene1()
{
  // Input file
  TFile *ipf = new TFile("sim.root");
  TTree *ipt = (TTree*)ipf->Get("tree");
  ipt->AddFriend("tree2 = tree","pos2.root");

  Double_t tu, td, tx;
  ipt->SetBranchAddress("tu",&tu);   
  ipt->SetBranchAddress("td",&td);
  ipt->SetBranchAddress("tree2.tx",&tx);

  // Output file
  TFile *opf = new TFile("ene1.root","recreate");

  TH2D *htofx = new TH2D("htofx","htofx",100,-120,120,200,0,200);
  TH2D *hgtofx = new TH2D("hgtofx","hgtofx",100,-120,120,100,39,48);
  TH1D *hgctof = new TH1D("hgctof","hgctof",100,39,48);

  // Loop
  Long64_t nentries = ipt->GetEntries();
  for(Long64_t jentry=0;jentry<nentries;jentry++) {
    ipt->GetEntry(jentry);
    if(jentry%10000==0) cout<<"process "<<jentry<<" of "<<nentries<<endl;

    Double_t ctof = (td + tu)/2;
    htofx->Fill(tx,ctof);
    if(ctof>42&& ctof<44.5) { 
      hgtofx->Fill(tx,ctof);
      if(abs(tx)<5) hgctof->Fill(ctof);//gamma hits the center of the det.
    }   
  }

  // 确定 C
  TFitResultPtr r = hgctof->Fit("gaus","S");
  r->SetName("r");
  r->Write();

  // Write and close
  opf->Write();
  ipf->Close();
  opf->Close();
}

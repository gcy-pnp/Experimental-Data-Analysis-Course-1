// TOF 绝对刻度 2
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

void ene2()
{
  // 常量声明
  const Double_t D = 500.;    //cm, distance between target and the scin.(Center)
  const Double_t dD = 5.;     //cm, thickness of the scin.
  Double_t rD = (D+dD/2);
  Double_t TOF0 = rD*0.01*3.333; //[ns]
  cout << "TOF0=" << TOF0 << endl;

  // Input file 1
  TFile *ipf = new TFile("sim.root");
  TTree *ipt = (TTree*)ipf->Get("tree");
  ipt->AddFriend("tree2 = tree","pos2.root");

  Double_t tu, td, tx;
  ipt->SetBranchAddress("tu",&tu);   
  ipt->SetBranchAddress("td",&td);
  ipt->SetBranchAddress("tree2.tx",&tx);

  // Input file 2
  TFile *ipf2 = new TFile("ene1.root");
  TFitResultPtr r = (TFitResult*)ipf2->Get("r");
  Double_t Toff = TOF0 - r->Parameter(1);
  cout << "Toff=" << Toff << endl;

  // Output file
  TFile *opf = new TFile("ene2.root","recreate");
  TTree *opt = new TTree("tree","tree");

  Double_t tofc, gne;
  opt->Branch("tofc",&tofc,"tofc/D");
  opt->Branch("gne",&tofc,"gne/D");

  TH2D *hgtofcx = new TH2D("hgtofcx","corrected TOF",100,-120,120,100,15,19);
  TH1D *htofc = new TH1D("htofc","htof",200,0,100);
  TH1D *hgne = new TH1D("hgne","hgne",350,-110,240);

  // Loop
  Long64_t nentries = ipt->GetEntries();
  for(Long64_t jentry=0;jentry<nentries;jentry++) {
    ipt->GetEntry(jentry);
    if(jentry%10000==0) cout<<"process "<<jentry<<" of "<<nentries<<endl;

    Double_t ctof = (tu+td)/2;
    Double_t d = TMath::Sqrt(rD*rD+tx*tx);
    Double_t tofc = (ctof+Toff)/d*500.; //normalized to 500cm

    if(ctof>42&& ctof<44.5) { 
      gne = -100;
    }   
    else {
      gne = TMath::Power((72*500*0.01)/tofc,2);
    }
    hgne->Fill(gne);

    hgtofcx->Fill(tx,tofc);//gamma hits the center of the det.
    htofc->Fill(tofc);

    opt->Fill();
  }

  // Write and close
  opf->Write();
  ipf->Close();
  ipf2->Close();
  opf->Close();
}

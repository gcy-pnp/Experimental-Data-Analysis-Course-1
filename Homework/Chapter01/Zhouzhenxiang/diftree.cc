void diftree()
{
  TFile *ipf=new TFile("tree.root");
  TTree *tree=(TTree*)ipf->Get("tree");
  Double_t x;
  Double_t e;
  int pid;
  Double_t tof, ctof;
  Double_t tu, td, tdu;
  Double_t qu, qd, lgqdu;

  tree->SetBranchAddress("ctof",&ctof);
  tree->SetBranchAddress("pid",&pid);
  tree->SetBranchAddress("tu",&tu);   
  tree->SetBranchAddress("td",&td);
  tree->SetBranchAddress("qu",&qu);   
  tree->SetBranchAddress("qd",&qd);

  TH1D *hctof=new TH1D("hctof","Time of flight", 1000,0,100);
  TH1D *htdu=new TH1D("htud","td-tu", 140,-20,50);
  TH1D *hqdu=new TH1D("hqdu","log(qd/qu)", 100,-2,2);

  Long64_t nentries=tree->GetEntries();
  for(Long64_t jentry=0; jentry<nentries; jentry++) {
    tree->GetEntry(jentry);
    tdu = td - tu;
    lgqdu = TMath::Log(qd/qu);
    hctof->Fill(ctof);
    htdu->Fill(tdu);
    hqdu->Fill(lgqdu);
    //opt->Fill()
    //if(jentry%100000==0) cout<<"process "<<jentry<<" of "<<nentries<<endl;
  }
  //htud->Draw();

  
  Int_t i,j,k,l;
  Int_t w1, w2;
  TH1D *dhtdu=new TH1D("dhtdu","differential", 141,-20.25,50.25);
  TH1D *dhqdu=new TH1D("dhqdu","differential", 101,-2.002,2.002);
  const double d1 = 0.5;
  const double d2 = 0.04;
  for(i=0; i<htdu->GetNbinsX(); i++) {
    w1=(htdu->GetBinContent(i+1)-htdu->GetBinContent(i))/d1;
    dhtdu->Fill(htdu->GetBinLowEdge(i+1),w1);
    //if(jentry%100000==0) cout<<"process "<<jentry<<" of "<<nentries<<endl;
  }
  for(i=0; i<hqdu->GetNbinsX(); i++) {
    w2=(hqdu->GetBinContent(i+1)-hqdu->GetBinContent(i))/d2;
    dhqdu->Fill(hqdu->GetBinLowEdge(i+1),w2);
    //if(jentry%100000==0) cout<<"process "<<jentry<<" of "<<nentries<<endl;
  }

  


  //dhtdu->Sumw2(0);
  //dhtdu->Draw();
  dhqdu->Sumw2(0);
  dhqdu->Draw();
}

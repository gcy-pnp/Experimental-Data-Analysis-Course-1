void normaltree()
{
  TFile *ipf=new TFile("tree.root");
  TTree *tree=(TTree*)ipf->Get("tree");
  Double_t x;
  Double_t e;
  int pid;
  Double_t tof, ctof;
  Double_t tu, td;
  Double_t qu, qd;

  tree->SetBranchAddress("ctof",&ctof);
  tree->SetBranchAddress("pid",&pid);
  tree->SetBranchAddress("tu",&tu);   
  tree->SetBranchAddress("td",&td);
  tree->SetBranchAddress("qu",&qu);   
  tree->SetBranchAddress("qd",&qd);

  TH2D *hgtofx=new TH2D("hgtofx","hgtofx",100,-120,120,100,39,48);
  TH1D *hgctof=new TH1D("hgctof","hgctof",100,39,48);

  Long64_t nentries=tree->GetEntries();
  for(Long64_t jentry=0; jentry<nentries; jentry++) {
    tree->GetEntry(jentry);
    Double_t tx=3.737*(td-tu-14.872);
    if(ctof>42&&ctof<44.5){
      hgtofx->Fill(tx,ctof);
      if(abs(tx)<5) hgctof->Fill(ctof);
    }
    //opt->Fill()
    //if(jentry%100000==0) cout<<"process "<<jentry<<" of "<<nentries<<endl;
  }
  hgctof->Draw();

  
 
}

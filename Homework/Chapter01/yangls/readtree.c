#include "TH1.h"
#include "TF1.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <TH2F.h>
#include <TText.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TString.h>
#include <math.h>

using namespace std;
TH1D *hdtd;
TH1D *hdqd;
TH1D *htx;
TH1D *hqx;
TH2D *hdtx;
TH2D *hdqx;
TH1D *htdiff;
TH1D *hqdiff;
TH2D *htxx;
TH2D *hqxx;

TH1D *hTOF;
TH2D *hgtofx;
TH1D *hgctof; 
TH2D *hgtofcx;
TH2D *hgtofcnx;
TH1D *htofc; 
TH1D *hce;
TH2D *hcee;
//scin:detectorT
void readtree(){
    gStyle->SetPalette(55);

    //const define
    const Double_t D = 500.;//cm;target to scin
    const Double_t L = 100.;//cm;half lepidth of scin
    const Double_t dD = 5.;//cm;thickness of scin
    const Double_t TRes = 1.;//ns;time resolution
    const Double_t Lambda = 380.;//cm;attenuation lepidth of scin
    const Double_t QRes = 0.1;//ER
    const Double_t EnRes = 50.;//MeV;FWHM
    const Double_t Vsc = 7.5;//ns/cm;speed of g in scin
    const Double_t En0 = 100.;//MeV;neutron e
    const Double_t Eg0 = 1;//MeV;gamma e
    const Double_t Rg = 0.3;//ratio of gamma;ratio of n =0.7;

    //root file and tree define
    TFile *ipf = new TFile("tree.root");
    TTree *tree = (TTree*)ipf->Get("tree");
    //variation define
    Double_t x;
    Double_t e;
    Double_t tof,ctof;//tof:true time of ftly;c:caculation
    Double_t tu,td;//up and down;time
    Double_t qu,qd;//energy
    int pid;//1:n;0:g

    //tree branch
    tree->SetBranchAddress("x",&x);
    tree->SetBranchAddress("e",&e);
    tree->SetBranchAddress("tof",&tof);
    tree->SetBranchAddress("ctof",&ctof);
    tree->SetBranchAddress("pid",&pid);
    tree->SetBranchAddress("tu",&tu);
    tree->SetBranchAddress("td",&td);
    tree->SetBranchAddress("qu",&qu);
    tree->SetBranchAddress("qd",&qd);

    TFile *opf = new TFile("tree2.root","recreate");
    TTree *opt = new TTree("tree2","tree2");
    Double_t tx,qx,ntof,ce;
    opt->Branch("tx",&tx,"tx/D");
    opt->Branch("qx",&qx,"qx/D");
    opt->Branch("ce",&ce,"ce/D");
    opt->Branch("ntof",&ntof,"ntof/D");
    opt->Branch("x",&x,"x/D");
    opt->Branch("e",&e,"e/D");
    opt->Branch("tof",&tof,"tof/D");
    opt->Branch("ctof",&ctof,"ctof/D");
    opt->Branch("pid",&pid,"pid/I");
    opt->Branch("tu",&tu,"tu/D");
    opt->Branch("td",&td,"td/D");
    opt->Branch("qu",&qu,"qu/D");
    opt->Branch("qd",&qd,"qd/D");

    //tx and qx caculate
    hTOF = new TH1D("hTOF","time of ftlight",1000,0,100);
    hgctof = new TH1D("hgctof","hgctof",100,39,45);
    hgtofx = new TH2D("hgtofx","hgtofx",100,-120,120,100,39,45);

    htdiff = new TH1D("htdiff","td-tu",140,-20,50);
    hdtd = new TH1D("hdtd","dt/dx",140,-20,50);
    htx = new TH1D("htx","htx",500,-120,120);
    htxx = new TH2D("htxx","tx:x",500,-120,120,500,-120,120);
    hdtx = new TH2D("hdtx","htx-hx:hx",100,-25,25,500,-120,120);

    hqdiff = new TH1D("hqdiff","log(qu/qd)",80,-1,1);
    hdqd = new TH1D("hdtd","dq/dx",80,-1,1);
    hqx = new TH1D("hqx","hqx",500,-120,120);
    hqxx = new TH2D("hqxx","qx:x",500,-120,120,500,-120,120);
    hdqx = new TH2D("hdqx","hqx-hx:hx",100,-25,25,500,-120,120);
    Long64_t nentries = tree->GetEntries();
    for(Long64_t jentry = 0; jentry<nentries; jentry++){
        tree->GetEntry(jentry);
	hTOF->Fill(ctof);
	htdiff->Fill(td-tu);
	hqdiff->Fill(log(qu/qd));
        for(int i=1; i<htdiff->GetNbinsX(); i++){
            Double_t dtf = htdiff->GetBinContent(i+1)-htdiff->GetBinContent(i);
            hdtd->Fill(htdiff->GetBinLowEdge(i+1),dtf);
        }
        for(int i=1; i<hqdiff->GetNbinsX(); i++){
            Double_t dqf = hqdiff->GetBinContent(i+1)-hqdiff->GetBinContent(i);
            hdqd->Fill(hqdiff->GetBinLowEdge(i+1),dqf);
        }
    }
    //fit to caculate tr and tl
    TF1 *ftl = new TF1("ftl","gaus",-14,-10);//left
    TF1 *ftr = new TF1("ftr","[0]*TMath::Exp(-0.5*((x-[1])/[2])^2)",40.2,43.5);//right
    ftr->SetParameter(0,-350);
    ftr->SetParameter(1,41.5);
    ftr->SetParameter(2,0.5);
    hdtd->Fit("ftl","RQ");
    hdtd->Fit("ftr","RQ+");

    //fit to caculate qr and ql
    TF1 *fql = new TF1("fql","gaus",-0.75,-0.45);//left
    TF1 *fqr = new TF1("fqr","[0]*TMath::Exp(-0.5*((x-[1])/[2])^2)", 0.49,0.77);//right
    fqr->SetParameter(0,-350);
    fqr->SetParameter(1,0.6);
    fqr->SetParameter(2,0.1);
    hdqd->Fit("fql","RQ");
    hdqd->Fit("fqr","RQ+");

    Double_t tl,tr,toffs,a,b;
    tl=ftl->GetParameter(1); 
    tr=ftr->GetParameter(1); 
    toffs=(tr+tl)/2;
    a=2*L/(tr-tl);
    b=-a*toffs;
    cout<<"tl = "<<tl<<"\t"<<"tr = "<<tr<<"\t"<<"toffset = "<<toffs<<endl;
    cout<<"tx=a*(td-tu)+b:  a = "<< a <<"\t"<<"b = "<< b <<endl;
     
    Double_t ql,qr,qoffs;
    ql=fql->GetParameter(1); 
    qr=fqr->GetParameter(1); 
    qoffs=(qr+ql)/2;
    cout<<"ql = "<<ql<<"\t"<<"qr = "<<qr<<"\t"<<"qoffset = "<<qoffs<<endl;
 
    for(Long64_t jentry = 0; jentry<nentries; jentry++)   {
        tree->GetEntry(jentry);
        tx=a*(td-tu)+b;
        htx->Fill(tx);
        hdtx->Fill(tx-x,x);
        htxx->Fill(x,tx);

        qx=2*L/(qr-ql)*(log(qu/qd)-qoffs);
        hqx->Fill(qx);
        hdqx->Fill(qx-x,x);
        hqxx->Fill(x,qx);
        if(ctof>40&&ctof<50)    {
            hgtofx->Fill(tx,ctof);
            if(abs(tx)<5) hgctof->Fill(ctof);
        }
    }
    TF1 *p1 = new TF1("p1","pol1",-100,100);
    htxx->Fit("p1","RQ");
    Double_t k1 = p1->GetParameter(1);
    cout<<"slope of tx:x = "<<k1<<endl;
    TH1D *hdtx1 = hdtx->ProjectionX("projx of hdtx");
    hdtx1->Fit("gaus","Q");

    TF1 *p2 = new TF1("p2","pol1",-100,100);
    hqxx->Fit("p2","RQ");
    Double_t k2 = p2->GetParameter(1);
    cout<<"slope of qx:x = "<<k2<<endl;
    TH1D *hdqx1 = hdqx->ProjectionX("projx of hdqx");
    hdqx1->Fit("gaus","Q");
    
    //tx draw
    TCanvas *c1 = new TCanvas("c1","tx",0,0,900,600);
    c1->Divide(3,2);
    c1->cd(1);
    htdiff->Draw();
    c1->cd(2);
    hdtd->Draw();
    hdtd->Sumw2(0);
    c1->cd(3);
    htxx->Draw("colz");
    c1->cd(4);
    htx->Draw();
    c1->cd(5);
    hdtx->Draw("colz");
    c1->cd(6);
    hdtx1->Draw();
//    c1->SaveAs("tx.pdf");
//    c1->SaveAs("tx.eps");

    //qx draw
    TCanvas *c2 = new TCanvas("c2","qx",0,0,900,600);
    c2->Divide(3,2);
    c2->cd(1);
    hqdiff->Draw();
    c2->cd(2);
    hdqd->Draw();
    hdqd->Sumw2(0);
    c2->cd(3);
    hqxx->Draw("colz");
    c2->cd(4);
    hqx->Draw();
    c2->cd(5);
    hdqx->Draw("colz");
    c2->cd(6);
    hdqx1->Draw();
//    c2->SaveAs("qx.pdf");
//    c2->SaveAs("qx.eps");

    //ce caculate
    htofc = new TH1D("hgtofc","htof",200,12,100);
    hgtofcx = new TH2D("hgtofx","corrected tof",100,-120,120,100,14,20);
    hgtofcnx = new TH2D("hgtofnx","corrected tof",100,-120,120,100,12,80);

    hce = new TH1D("hce","ce",200,-20,200);
    hcee = new TH2D("hcee","ce:e",200,0,200,200,0,200);

    Double_t TOF0,tof0,CTOF,d;
    TF1 *gt = new TF1("gt","gaus",40,50);
    hgctof->Fit("gt","RQ");
    TOF0=gt->GetParameter(1);
    tof0 = 10*5.025/3.;
    CTOF = tof0-TOF0;
    cout<<"TOF(0) = "<<TOF0<<"ns"<<endl;
    cout<<"TOF0 = "<< tof0<<"ns"<<endl;
    cout<<"C = TOF0 - TOF(0) = "<<CTOF<<"ns"<<endl;
    for(Long64_t jentry = 0; jentry<nentries; jentry++)   {
        tree->GetEntry(jentry);
        tx=a*(td-tu)+b;
        d=TMath::Sqrt(502.5*502.5+tx*tx);
        ntof=(ctof+CTOF)/d*500;
        hgtofcx->Fill(tx,ntof);
        hgtofcnx->Fill(tx,ntof);
        htofc->Fill(ntof);
        if(ntof>20)    ce = (72*d/(ntof*100))*(72*d/(ntof*100));
            else    ce = -10;
        hce->Fill(ce);
        hcee->Fill(e,ce);
        if(jentry%20000==0) cout<<"process"<<jentry<<"of"<<nentries<<endl;
        opt->Fill();
    }
    TF1 *gce = new TF1("gce","gaus",0,200);
    hce->Fit("gce","RQ");
    Double_t meance = gce->GetParameter(1);
    cout << "mean of ce = "<<meance<<endl;
    //calibration tof draw
    TCanvas *c3 = new TCanvas("c3","tofcalibration",0,0,900,600);
    c3->Divide(3,2);
    c3->cd(1);
    c3->cd(1)->SetLogy();
    hTOF->Draw();
    c3->cd(2);
    hgtofx->Draw("colz");
    c3->cd(3);
    hgctof->Draw();
    c3->cd(4);
    c3->cd(4)->SetLogy();
    htofc->Draw();
    c3->cd(5);
    hgtofcx->Draw("colz");
    c3->cd(6);
    hgtofcnx->Draw("colz");
//    c3->SaveAs("tofcalibration.pdf");
//    c3->SaveAs("tofcalibration.eps");
  
    //ce draw
    TCanvas *c4 = new TCanvas("c4","ce",0,0,1200,600);
    c4->Divide(2,1);
    c4->cd(1);
    c4->cd(1)->SetLogy();
    hce->Draw();
    c4->cd(2);
    hcee->Draw("colz");
//    c4->SaveAs("ce.pdf");
//    c4->SaveAs("ce.eps");
//    ipf->Close();
    opt->Write();
    opf->Close();
}

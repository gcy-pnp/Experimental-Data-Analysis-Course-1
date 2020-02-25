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

//scin:detector
void tree(){
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
    TFile *opf = new TFile("tree.root","recreate");
    TTree *opt = new TTree("tree","tree structure");
    //variation define
    Double_t x;
    Double_t e;
    Double_t tof,ctof;//tof:true time of fly;c:caculation
    Double_t tu,td;//up and down;time
    Double_t qu,qd;//energy
    int pid;//1:n;0:g
    //offset
    Double_t tuoff= 5.5;
    Double_t tdoff= 20.4;
    
    //tree branch
    opt->Branch("x",&x,"x/D");
    opt->Branch("e",&e,"e/D");
    opt->Branch("tof",&tof,"tof/D");
    opt->Branch("ctof",&ctof,"ctof/D");
    opt->Branch("pid",&pid,"pid/I");
    opt->Branch("tu",&tu,"tu/D");
    opt->Branch("td",&td,"td/D");
    opt->Branch("qu",&qu,"qu/D");
    opt->Branch("qd",&qd,"qd/D");
    //histogram
    TH1D *hctof = new TH1D("hctof","neutron time of flight",1000,0,100);
    //Random
    TRandom3 *gr = new TRandom3(0);

    for(int i=0;i<100000;i++){
        x = gr->Uniform(-L,L);
        Double_t Dr = D + gr->Uniform(-0.5,0.5)*dD;
        Double_t d = TMath::Sqrt(Dr*Dr+x*x);//m,flight path
        if(gr->Uniform()<Rg){//gamma
            pid = 0;
            e = Eg0;
            tof = 10*(d*0.01)/3;
        }
        else {//n
            pid = 1;
            e = gr->Gaus(En0,EnRes/2.35);
            tof = 72./TMath::Sqrt(e)*(d*0.01);
        }
        tu = tof + (L-x)/Vsc + gr->Gaus(0,TRes/2.35)+tuoff;
        td = tof + (L+x)/Vsc + gr->Gaus(0,TRes/2.35)+tdoff;
        ctof = (tu+td)/2;
        hctof->Fill(ctof);
        Double_t q0 = e*gr->Uniform();
        qu = q0*TMath::Exp(-(L-x)/Lambda);
        qu = gr->Gaus(qu,qu*QRes/2.35);
        qd = q0*TMath::Exp(-(L+x)/Lambda);
        qd = gr->Gaus(qd,qd*QRes/2.35);
        opt->Fill();
    }
    
    hctof->Write();
    opt->Write();
    opf->Close();
}

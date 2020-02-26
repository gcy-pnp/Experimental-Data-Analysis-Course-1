// Simulation
// 在ROOT命令行内解释运行此脚本时无需include
// ROOT和C++的头文件，此处头文件只用于vim的ycm.
#include <TROOT.h> 
#include <TFile.h> 
#include <TTree.h> 
#include <TH1D.h> 
#include <TH2D.h> 
#include <TRandom3.h> 
#include <TMath.h> 

void sim(){
  // 常量声明
  const Double_t D = 500.;    //cm, distance between target and the scin.(Center)
  const Double_t L = 100.;    //cm, half length of the scin.
  const Double_t dD = 5.;     //cm, thickness of the scin.
  const Double_t TRes = 1.;   //ns, time resolution(FWHM) of the scintillator.
  const Double_t Lambda = 380.; //cm, attenuation lenght of the scin.
  const Double_t QRes = 0.1;  //relative energy resolution(FWHM) of the scin. 
  const Double_t Vsc = 7.5;   //ns/cm, speed of light in the scin.
  const Double_t En0 = 100;   //MeV, average neutron energy
  const Double_t EnRes = 50.; //MeV, energy spread of neutron(FWHM)
  const Double_t Eg0 = 1;     //MeV, gamma energy  
  const Double_t Rg = 0.3;    //ratio of gamma,ratio of neutron 1-Rg 

  // Output file
  TFile *opf = new TFile("sim.root","recreate");
  TTree *opt = new TTree("tree","tree structure");
  
  int pid;      //粒子种类，n:pid=1,g:pid=0
  Double_t x;   //入射位置
  Double_t e;   //能量
  Double_t tof; //TOF:粒子实际飞行时间
  Double_t tu, td;
  Double_t qu, qd;
  Double_t tu_off = 5.5;  //time offset
  Double_t td_off = 20.4; //time offset

  opt->Branch("x", &x, "x/D");
  opt->Branch("e", &e, "e/D");
  opt->Branch("tof", &tof, "tof/D");
  opt->Branch("pid", &pid, "pid/I");
  opt->Branch("tu", &tu, "tu/D");
  opt->Branch("td", &td, "td/D");
  opt->Branch("qu", &qu, "qu/D"); 
  opt->Branch("qd", &qd, "qd/D");  

  TH1D *hx = new TH1D("hx","x position in [cm]",1000,0,100);
  TH1D *htof = new TH1D("htof","time of flight in [ns]",1000,0,100);

  // Loop
  TRandom3 *gr = new TRandom3(0);
  for(int i=0;i<100000;i++){
    x = gr->Uniform(-L, L); //均匀入射在中子探测器表面.
    Double_t Dr = D+gr->Uniform(-0.5,0.5)*dD; //粒子在探测器厚度范围内均匀产生光信号
    Double_t d = TMath::Sqrt(Dr*Dr+x*x);  //flight path
    hx->Fill(x);

    if(gr->Uniform()<Rg) { 
      pid = 0;              //gamma
      e = Eg0;              //[MeV]
      tof = 3.333*(d*0.01); //[ns]
    }
    else {
      pid = 1;                         //neutron
      e = gr->Gaus(En0,EnRes/2.35);    //[MeV]
      tof=72./TMath::Sqrt(e)*(d*0.01); //[ns]
    }
    htof->Fill(tof);

    tu = tof+(L-x)/Vsc + gr->Gaus(0,TRes/2.35) + tu_off;
    td = tof+(L+x)/Vsc + gr->Gaus(0,TRes/2.35) + td_off;

    Double_t q0 = e*gr->Uniform(); //energy of recoil proton in plas. 0-En
    qu = q0*TMath::Exp(-(L-x)/Lambda);
    qu = gr->Gaus(qu,qu*QRes/2.35);
    qd = q0*TMath::Exp(-(L+x)/Lambda);
    qd = gr->Gaus(qd,qd*QRes/2.35);

    opt->Fill();
  }

  // Write and close
  hx->Write();
  htof->Write();
  opt->Write();
  opf->Close();
}

void tree(){
  const Double_t D=500;//cm, distance between target and the scin.(Center)
  const Double_t L=100;//cm, half length of the scin.
  const Double_t dD=5;//cm, thickness of the scin.
  const Double_t TRes=1.;//ns, sigma of time of the scin.
  const Double_t Lambda=380.;//cm, attenuation lenght of the scin.
  const Double_t QRes=0.1;//relative energy resolution of the scin. 
  const Double_t Vsc=7.5;//ns/cm, speed of light in the scin.
  const Double_t En0=100;//MeV, average neutron energy
  const Double_t EnRes=50.;//MeV, sigma of En
  const Double_t Eg0=1.;//MeV, gamma energy  
  const Double_t Rg=0.3;//ratio of gamma 

  TFile *opf=new TFile("tree.root","recreate");//新文件tree.root的指针 *opf
  TTree *opt=new TTree("tree","tree structure");//新tree的指针 *opt
  // 在tree结构中定义需要的变量分支
  Double_t x;
  Double_t e;
  int pid;
  Double_t tof,ctof;
  Double_t tu, td;
  Double_t qu, qd;

  Double_t tu_off=5.5;//time offset
  Double_t td_off=20.4;//time offset

  // 将变量分支添加到tree结构中,第一个参数为变量名称，第二个为上面定义的变量地址，第三个为变量的类型说明，D表示Double_t。
  opt->Branch("x", &x, "x/D");//position of neutron
  opt->Branch("e", &e, "e/D");//energy of gamma or neutron
  opt->Branch("tof", &tof, "tof/D");//time of flight
  opt->Branch("ctof",&ctof,"ctof/D");//TOF from exp. data
  opt->Branch("pid", &pid, "pid/I");//1/0 : neutron/gamma  
  opt->Branch("tu", &tu, "tu/D");//time of upper side
  opt->Branch("td", &td, "td/D");//time of bottom side
  opt->Branch("qu", &qu, "qu/D"); 
  opt->Branch("qd", &qd, "qd/D");  
  // histogram
   TH1D *hctof=new TH1D("hctof","neutron time of flight",1000,0,100);//文件中还可存储histogram，graph等
  TRandom3 *gr=new TRandom3(0);
  // 循环，逐事件往tree结构里添加对应分支信息。
  for(int i=0;i<100000;i++){
    x=gr->Uniform(-L, L);
    Double_t Dr=D+gr->Uniform(-0.5,0.5)*dD;
    Double_t d=TMath::Sqrt(Dr*Dr+x*x);//m, flight path
    if(gr->Uniform() < Rg) { //gamma
       pid=0;
       e=Eg0;
       tof=3.333*(d*0.01);
    }
    else {  //neutron
        pid=1;
        e=gr->Gaus(En0, EnRes/2.35); // neutron
        tof=72./TMath::Sqrt(e)*(d*0.01);
    }
    tu=tof+(L-x)/Vsc+gr->Gaus(0,TRes/2.35)+tu_off;
    td=tof+(L+x)/Vsc+gr->Gaus(0,TRes/2.35)+td_off;
    ctof=(tu+td)/2.;
    hctof->Fill(ctof);
    Double_t q0=e*gr->Uniform();
    qu=q0*TMath::Exp(-(L-x)/Lambda);
    qu=gr->Gaus(qu,qu*QRes/2.35);
    qd=q0*TMath::Exp(-(L+x)/Lambda);
    qd=gr->Gaus(qd,qd*QRes/2.35);
    opt->Fill();
  }
  // 将数据写入root文件中
  hctof->Write();
  opt->Write();
  opf->Close();
}

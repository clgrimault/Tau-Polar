
void mix(){
  double pol=-0.156789123;

  //  TString fileName="A1Decay1M_default.root";
  TString fileName="Combined100M_Default.root";
  TString HPlus = "omega_a1p_plus";
  TString HMins = "omega_a1p_minus";


  double wplus = 0.5*(1+pol);
  double wplus = 0.5*(1-pol);
  TFile *_file0 = TFile::Open(fileName);
  TFile *out = new TFile("MixedTemplates.root","RECREATE");
  TH1F  *hp = (TH1F*)_file0->Get(HPlus);
  TH1F  *hm = (TH1F*)_file0->Get(HMins);
  hp->SetName(HPlus);
  hm->SetName(HMins);
  double Integral = hp->Integral() + hm->Integral();
   
  hp->Sumw2();
  hm->Sumw2();

  hp->Scale(1./hp->Integral());
  hm->Scale(1./hm->Integral());

  TH1F *mixed = (TH1F*)hp->Clone("mixed");
  mixed->Add(hm,(1-pol)/(1+pol));
  mixed->Scale(1./mixed->Integral());
  mixed->Scale(Integral);
  mixed->Rebin(2);
  hp->Rebin(2);
  hm->Rebin(2);
  mixed->Write();
  hp->Write();
  hm->Write();
  _file0->Close();
}

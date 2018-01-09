
void mix(){
  double pol=0.5123456789;

  TString fileName="Draw.root";
  TString HPlus = "pi_plus";
  TString HMins = "pi_minus";


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
  mixed->Write();
  hp->Write();
  hm->Write();
  _file0->Close();
}

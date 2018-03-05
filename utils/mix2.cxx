
void mix2(){
  double pol=0.1567891;

  TString fileup="upSyst.root";
  TString filedown="downSyst.root";
  TString filedefault="defSyst.root";
  TString HPlus = "omega_a1p_plus";
  TString HMins = "omega_a1p_minus";

  TString HPlusup = "omega_a1p_plus_up";
  TString HMinsup = "omega_a1p_minus_up";
  TString HPlusdown = "omega_a1p_plus_down";
  TString HMinsdown = "omega_a1p_minus_down";



  double wplus = 0.5*(1+pol);
  double wplus = 0.5*(1-pol);
  TFile *_filedef = TFile::Open(filedefault);
  TFile *_fileup = TFile::Open(fileup);
  TFile *_filedown = TFile::Open(filedown);
  TFile *out = new TFile("MixedTemplates.root","RECREATE");
  TH1F  *hp = (TH1F*)_filedef->Get(HPlus);
  TH1F  *hm = (TH1F*)_filedef->Get(HMins);


  TH1F  *hp_up = (TH1F*)_fileup->Get(HPlus);
  TH1F  *hm_up = (TH1F*)_fileup->Get(HMins);

  TH1F  *hp_down = (TH1F*)_filedown->Get(HPlus);
  TH1F  *hm_down = (TH1F*)_filedown->Get(HMins);


  hp->SetName(HPlus);
  hm->SetName(HMins);

  hp_up->SetName(HPlusup);
  hm_up->SetName(HMinsup);

  hp_down->SetName(HPlusdown);
  hm_down->SetName(HMinsdown);

  double Integral = hp->Integral() + hm->Integral();
  double Integralup = hp_up->Integral() + hm_up->Integral();
  double Integraldown = hp_down->Integral() + hm_down->Integral();
 
  hp->Sumw2();
  hm->Sumw2();

  hp->Scale(1./hp->Integral());
  hm->Scale(1./hm->Integral());

  hp_up->Sumw2();
  hm_up->Sumw2();

  hp_up->Scale(1./hp_up->Integral());
  hm_up->Scale(1./hm_up->Integral());



  hp_down->Sumw2();
  hm_down->Sumw2();

  hp_down->Scale(1./hp_down->Integral());
  hm_down->Scale(1./hm_down->Integral());





  TH1F *mixed = (TH1F*)hp->Clone("mixed");
  mixed->Add(hm,(1-pol)/(1+pol));
  mixed->Scale(1./mixed->Integral());
  mixed->Scale(Integral);
  mixed->Rebin(2);
  hp->Rebin(2);
  hm->Rebin(2);

  hp_up->Rebin(2);
  hm_up->Rebin(2);

  hp_down->Rebin(2);
  hm_down->Rebin(2);




  mixed->Write();
  hp->Write();
  hm->Write();

  hp_up->Write();
  hm_up->Write();

  hp_down->Write();
  hm_down->Write();



  //  filedefault->Close();
}

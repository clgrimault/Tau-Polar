
void ratio(){
  TString fileName1="A1Decay1M_down.root";
  TString fileName2="A1Decay1M_up.root";
  TString HPlus = "omega_a1_plus";
  TString HPlus2 = "omega_a1p_plus";
  TString HMins = "omega_a1_minus";
  TString HMins2 = "omega_a1p_minus";

  TFile *_file1 = TFile::Open(fileName1);
  TFile *_file2 = TFile::Open(fileName2);

  TH1F  *hp1 = (TH1F*)_file1->Get(HPlus);
  TH1F  *hp2 = (TH1F*)_file2->Get(HPlus);
  hp1->Rebin(2);
  hp2->Rebin(2);
  hp1->Sumw2();
  hp2->Sumw2();

  TH1F *h = new TH1F("h","h",25,-1.1,1.1);
  h->Sumw2();
  h->Divide(hp1,hp2);
  h->Sumw2();
  h->Draw();
}

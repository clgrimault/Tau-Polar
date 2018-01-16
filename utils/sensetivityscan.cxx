
#include "TH1.h"
#include "TF1.h"

  TFile *_file0 = TFile::Open("Combined.root");
//TFile *_file0 = TFile::Open("HelicityVals.root");


double Sensetivityrho(Double_t *x, Double_t *par)
{

  double pol=par[0]*x[0];
  TH1F *plus = (TH1F *)_file0->Get("omega_rho_plus");
  TH1F *minus = (TH1F *)_file0->Get("omega_rho_minus");
  plus->Scale(1./plus->Integral());
  minus->Scale(1./minus->Integral());
  TH1F *hf = (TH1F*)plus->Clone("hf");
  TH1F *hg = (TH1F*)plus->Clone("hg");
  hf->Add(minus,1);
  hg->Add(minus,-1); 
  hf->Add(hg,pol);
  hg->Multiply(hg);
 
  hg->Divide(hf);
  return hg->Integral();
}



double Sensetivitymupi(Double_t *x, Double_t *par)
{

  double pol=par[0]*x[0];
  TH1F *plus = (TH1F *)_file0->Get("OmegaMuPi_plus");
  TH1F *minus = (TH1F *)_file0->Get("OmegaMuPi_minus");
  plus->Scale(1./plus->Integral());
  minus->Scale(1./minus->Integral());
  TH1F *hf = (TH1F*)plus->Clone("hf");
  TH1F *hg = (TH1F*)plus->Clone("hg");
  hf->Add(minus,1);
  hg->Add(minus,-1);
  hf->Scale(0.5);  
  hg->Scale(0.5);
  hf->Add(hg,pol);
  hg->Multiply(hg);
  hg->Divide(hf);
 
  return sqrt(hg->Integral());

}


double Sensetivitya1old(Double_t *x, Double_t *par)
{

  double pol=par[0]*x[0];
  TH1F *plus = (TH1F *)_file0->Get("omega_a1_plus");
  TH1F *minus = (TH1F *)_file0->Get("omega_a1_minus");
  plus->Scale(1./plus->Integral());
  minus->Scale(1./minus->Integral());
  TH1F *hf = (TH1F*)plus->Clone("hf");
  TH1F *hg = (TH1F*)plus->Clone("hg");
  hf->Add(minus,1);
  hg->Add(minus,-1); 
 hf->Add(hg,pol);
  hg->Multiply(hg);
 
  hg->Divide(hf);
  return hg->Integral();
}

double Sensetivitya1new(Double_t *x, Double_t *par)
{

  double pol=par[0]*x[0];
  TH1F *plus = (TH1F *)_file0->Get("TRFomegabar_a1_plus");
  TH1F *minus = (TH1F *)_file0->Get("TRFomegabar_a1_minus");
  plus->Scale(1./plus->Integral());
  minus->Scale(1./minus->Integral());
  TH1F *hf = (TH1F*)plus->Clone("hf");
  TH1F *hg = (TH1F*)plus->Clone("hg");
  hf->Add(minus,1);
  hg->Add(minus,-1); 
 hf->Add(hg,pol);
  hg->Multiply(hg);
 
  hg->Divide(hf);
  return hg->Integral();
}



double Sensetivitypi(Double_t *x, Double_t *par)
{
  double pol=par[0]*x[0];
  TH1F *plus = (TH1F *)_file0->Get("pi_plus");
  TH1F *minus = (TH1F *)_file0->Get("pi_minus");
  plus->Scale(1./plus->Integral());
  minus->Scale(1./minus->Integral());
  TH1F *hf = (TH1F*)plus->Clone("hf");
  TH1F *hg = (TH1F*)plus->Clone("hg");
  hf->Add(minus,1);
  hg->Add(minus,-1);
  hf->Scale(0.5);  
  hg->Scale(0.5);
  hf->Add(hg,pol);
  hg->Multiply(hg);
  hg->Divide(hf);
 
  return sqrt(hg->Integral());
}


double Sensetivitypipi(Double_t *x, Double_t *par)
{

  double pol=par[0]*x[0];
  TH1F *plus = (TH1F *)_file0->Get("Omegapipi_plus");
  TH1F *minus = (TH1F *)_file0->Get("Omegapipi_minus");
  plus->Scale(1./plus->Integral());
  minus->Scale(1./minus->Integral());
  TH1F *hf = (TH1F*)plus->Clone("hf");
  TH1F *hg = (TH1F*)plus->Clone("hg");
  hf->Add(minus,1);
  hg->Add(minus,-1);
  hf->Scale(0.5);  
  hg->Scale(0.5);
  hf->Add(hg,pol);
  hg->Multiply(hg);
  hg->Divide(hf);

  return sqrt(hg->Integral());


}



double Sensetivitypirho(Double_t *x, Double_t *par)
{

  double pol=par[0]*x[0];
  TH1F *plus = (TH1F *)_file0->Get("omega_pirho_plus");
  TH1F *minus = (TH1F *)_file0->Get("omega_pirho_minus");



  plus->Scale(1./plus->Integral());
  minus->Scale(1./minus->Integral());
  TH1F *hf = (TH1F*)plus->Clone("hf");
  TH1F *hg = (TH1F*)plus->Clone("hg");
  hf->Add(minus,1);
  hg->Add(minus,-1);
  hf->Add(hg,pol);
  hg->Multiply(hg);
 
  hg->Divide(hf);
  return hg->Integral();
}





double Sensetivitymu(Double_t *x, Double_t *par)
{

  double pol=par[0]*x[0];
   TH1F *plus = (TH1F *)_file0->Get("mu_plus");
   TH1F *minus = (TH1F *)_file0->Get("mu_minus");
  plus->Scale(1./plus->Integral());
  minus->Scale(1./minus->Integral());
  TH1F *hf = (TH1F*)plus->Clone("hf");
  TH1F *hg = (TH1F*)plus->Clone("hg");
  hf->Add(minus,1);
  hg->Add(minus,-1);
  hf->Scale(0.5);  
  hg->Scale(0.5);
  hf->Add(hg,pol);
  hg->Multiply(hg);
  hg->Divide(hf);
 
  return sqrt(hg->Integral());


}



double SensetivitymuTest(Double_t *x, Double_t *par)
{

  double pol=par[0]*x[0];
  TH1F *plus = (TH1F *)_file0->Get("mu_plus");
  TH1F *minus = (TH1F *)_file0->Get("mu_minus");
  plus->Scale(1./plus->Integral());
  minus->Scale(1./minus->Integral());
  TH1F *hf = (TH1F*)plus->Clone("hf");
  TH1F *hg = (TH1F*)plus->Clone("hg");
  hf->Add(minus,0.5);
  hf->Draw();
  hg->Add(minus,-1);
  hf->Add(hg,pol);
  hg->Multiply(hg);
 
  hg->Divide(hf);
  return hg->Integral();
}



void sensetivityscan(){

  TFile *_file0 = TFile::Open("Combined.root");
  double pol=-0.0;
  //   TH1F *plus = (TH1F *)_file0->Get("omega_a1_plus");
  //   TH1F *minus = (TH1F *)_file0->Get("omega_a1_minus");

  //   TH1F *plus = (TH1F *)_file0->Get("omega_a1p_plus");
  //   TH1F *minus = (TH1F *)_file0->Get("omega_a1p_minus");


    // TH1F *plus = (TH1F *)_file0->Get("mu_plus");
    // TH1F *minus = (TH1F *)_file0->Get("mu_minus");

    // TH1F *plus = (TH1F *)_file0->Get("rhobeta_plus");
    // TH1F *minus = (TH1F *)_file0->Get("rhobeta_minus");

    // TH1F *plus = (TH1F *)_file0->Get("pi_plus");
    // TH1F *minus = (TH1F *)_file0->Get("pi_minus");

    TH1F *plus = (TH1F *)_file0->Get("Omegapipi_plus");
    TH1F *minus = (TH1F *)_file0->Get("Omegapipi_minus");

 

// TH1F *plus = (TH1F *)_file0->Get("mass_pirho_plus");
// TH1F *minus = (TH1F *)_file0->Get("mass_pirho_minus");

// TH1F *plus = (TH1F *)_file0->Get("mass_murho_plus");
// TH1F *minus = (TH1F *)_file0->Get("mass_murho_minus");

// TH1F *plus = (TH1F *)_file0->Get("mass_rhorho_plus");
// TH1F *minus = (TH1F *)_file0->Get("mass_rhorho_minus");
 
// TH1F *plus = (TH1F *)_file0->Get("OmegaMuPi_plus");
// TH1F *minus = (TH1F *)_file0->Get("OmegaMuPi_minus");

 
// TH1F *plus = (TH1F *)_file0->Get("pipi_mass_plus");
// TH1F *minus = (TH1F *)_file0->Get("pipi_mass_minus");

// TH1F *plus = (TH1F *)_file0->Get("omega_a1pi_plus");
// TH1F *minus = (TH1F *)_file0->Get("omega_a1pi_minus");

 // TH1F *plus = (TH1F *)_file0->Get("omega_pirho_plus");
 // TH1F *minus = (TH1F *)_file0->Get("omega_pirho_minus");

 // TH1F *plus = (TH1F *)_file0->Get("omega_rho_plus");
 // TH1F *minus = (TH1F *)_file0->Get("omega_rho_minus");

// TH1F *plus = (TH1F *)_file0->Get("omegabar_rho_plus");
// TH1F *minus = (TH1F *)_file0->Get("omegabar_rho_minus");


  plus->Scale(1./plus->Integral());
  minus->Scale(1./minus->Integral());
  TH1F *hf = (TH1F*)plus->Clone("hf");
  TH1F *hg = (TH1F*)plus->Clone("hg");

  hf->Add(minus,1);
  hg->Add(minus,-1);

  hf->Scale(0.5);  
  hg->Scale(0.5);
  // hf->Draw();
   // double hfinteg=hf->Integral();
   // hf->Scale(1./hfinteg);
 

  hf->Add(hg,pol);
  hg->Multiply(hg);

  hg->Divide(hf);
   
  hf->Draw();
  hg->Draw("same");

  std::cout<<"hg  Integral " << sqrt(hg->Integral()) << std::endl;


  TF1 *funcpi = new TF1("Sensetivitypi",Sensetivitypi,-0.95,0.95,1);
  funcpi->SetParameter(0,1);
 
  funcpi->Draw();

  
  TF1 *funcpipi = new TF1("Sensetivitypipi",Sensetivitypipi,-0.95,0.95,1);
  funcpipi->SetParameter(0,1);
  funcpipi->SetLineColor(4);
  funcpipi->Draw("same");

  TF1 *funcmu = new TF1("Sensetivitymu",Sensetivitymu,-0.95,0.95,1);
  funcmu->SetParameter(0,1);
  funcmu->SetLineColor(13);
  funcmu->Draw("same");

  TF1 *funcmupi = new TF1("Sensetivitymupi",Sensetivitymupi,-0.95,0.95,1);
  funcmupi->SetParameter(0,1);
  funcmupi->SetLineColor(3);
  funcmupi->Draw("same");

}

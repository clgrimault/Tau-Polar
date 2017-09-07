
#include "TH1.h"
#include "TF1.h"

  TFile *_file0 = TFile::Open("1.root");
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



double Sensetivitymurho(Double_t *x, Double_t *par)
{

  double pol=par[0]*x[0];
  TH1F *plus = (TH1F *)_file0->Get("omega_murho_plus");
  TH1F *minus = (TH1F *)_file0->Get("omega_murho_minus");
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
 hf->Add(hg,pol);
  hg->Multiply(hg);
 
  hg->Divide(hf);
  return hg->Integral();
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
 hf->Add(hg,pol);
  hg->Multiply(hg);
 
  hg->Divide(hf);
  return hg->Integral();
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
 hf->Add(hg,pol);
  hg->Multiply(hg);
 
  hg->Divide(hf);
  return hg->Integral();
}




void SensetivityScan(){

  TFile *_file0 = TFile::Open("1.root");
  double pol=0.;

  TH1F *plus = (TH1F *)_file0->Get("pi_plus");
  TH1F *minus = (TH1F *)_file0->Get("mu_minus");
  plus->Scale(1./plus->Integral());
  minus->Scale(1./minus->Integral());
  TH1F *hf = (TH1F*)plus->Clone("hf");
  TH1F *hg = (TH1F*)plus->Clone("hg");
   hf->Add(minus,1);
   hg->Add(minus,-1);
   hf->Add(hg,0);
   hg->Multiply(hg);
  
   hg->Divide(hf);

     double integ(0);
  //     for(unsigned int i =1; i < hf->GetNbinsX(); i++){
  //         integ +=      hg->GetBinContent(i)* hg->GetBinWidth(i);
  //         }


  //   std::cout<< "sens  "<< sqrt(integ) <<"  "  <<sqrt(hg->Integral())<<std::endl;

   // hf->Draw();
   //  hg->Draw("same");
 
    //-----------------------------------------------------------------------------------------------------------------
          TF1 *funcpi = new TF1("Sensetivitypi",Sensetivitypi,-1,1,1);
          funcpi->SetParameter(0,1);
          funcpi->Draw();


          TF1 *funcpipi = new TF1("Sensetivitypipi",Sensetivitypipi,-1,1,1);
          funcpipi->SetParameter(0,1);
          funcpipi->SetLineColor(4);
          funcpipi->Draw("same");

          TF1 *funcpirho = new TF1("Sensetivitypirho",Sensetivitypirho,-1,1,1);
          funcpirho->SetParameter(0,1);
          funcpirho->SetLineColor(5);
          funcpirho->Draw("same");

    //      TF1 *funcmu = new TF1("Sensetivitymu",Sensetivitymu,-1,1,1);
    //  funcmu->SetParameter(0,1);
    //   funcmu->Draw("same");

          TF1 *funcrho = new TF1("Sensetivityrho",Sensetivityrho,-0.9,0.9,1);
	  funcrho->SetParameter(0,1);
	  funcrho->Draw("same");


          TF1 *funcmurho = new TF1("Sensetivitymurho",Sensetivitymurho,-0.9,0.9,1);
	  funcmurho->SetParameter(0,1);
	  funcmurho->SetLineColor(3);
	  funcmurho->Draw("same");


     //    TF1 *funca1old = new TF1("Sensetivitya1old",Sensetivitya1old,-1,1,1);
     // funca1old->SetParameter(0,1);
     // funca1old->SetLineColor(3);
     //  funca1old->Draw("same");


     //    TF1 *funca1new = new TF1("Sensetivitya1new",Sensetivitya1new,-1,1,1);
     // funca1new->SetParameter(0,1);
     // funca1new->SetLineColor(5);
     //  funca1new->Draw("same");



    // TH1 * hfmbetaM = fmbetaM->GetHistogram();
    // hfmbetaM->SetName("MomentbetaM");
   


 
}

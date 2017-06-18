
#include "TH1.h"
#include "TF1.h"

  TFile *_file0 = TFile::Open("Helic.root");
double Sensetivity(Double_t *x, Double_t *par)
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
  hg->Multiply(hg);
  hf->Add(hg,pol);
  hg->Divide(hf);
 

 return hg->Integral();
}

void utils(){

  TFile *_file0 = TFile::Open("Helic.root");
  double pol=-0.;

  TH1F *plus = (TH1F *)_file0->Get("rho_plus");
  TH1F *minus = (TH1F *)_file0->Get("rho_minus");
  plus->Scale(1./plus->Integral());
  minus->Scale(1./minus->Integral());
  TH1F *hf = (TH1F*)plus->Clone("hf");
  TH1F *hg = (TH1F*)plus->Clone("hg");
  hf->Add(minus,1);
  hg->Add(minus,-1);
   hg->Multiply(hg);
    hf->Add(hg,pol);
    hg->Divide(hf);

     double integ(0);
    //  for(unsigned int i =1; i < hf->GetNbinsX(); i++){
    //       integ +=      hg->GetBinContent(i)* hg->GetBinWidth(i);
    //      }


     std::cout<< "sens  "<< sqrt(integ) <<"  "  <<sqrt(hg->Integral())<<std::endl;

  hf->Draw();
  hg->Draw("same");
 
     //-----------------------------------------------------------------------------------------------------------------
    // TF1 *func = new TF1("Sensetivity",Sensetivity,-1,1,1);
    // func->SetParameter(0,1);
    //unc->Draw();
    // TH1 * hfmbetaM = fmbetaM->GetHistogram();
    // hfmbetaM->SetName("MomentbetaM");
   


 
}


#include "TH1.h"
#include "TF1.h"

  TFile *_file0 = TFile::Open("HelicityVals1V.root");
double Sensetivity(Double_t *x, Double_t *par)
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
  hg->Multiply(hg);
  hf->Add(hg,pol);
  hg->Divide(hf);
 

 return hg->Integral();
}

void utils(){

  TFile *_file0 = TFile::Open("HelicityVals1V.root");
  double pol=-0.;

  TH1F *plus = (TH1F *)_file0->Get("TRFomegabar_a1_plus");
  TH1F *minus = (TH1F *)_file0->Get("TRFomegabar_a1_minus");
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
     for(unsigned int i =1; i < hf->GetNbinsX(); i++){
       //      integ +=      hg->GetBinContent(i)*hg->GetBinContent(i)*hg->GetBinWidth(i)/ (hf->GetBinContent(i) + pol*hg->GetBinContent(i));
       integ +=      hg->GetBinContent(i)* hg->GetBinWidth(i);

       std::cout<< "hg  "<< hg->GetBinContent(i) << "        hf  "<< hf->GetBinContent(i)   <<"   " <<integ <<std::endl;
     }
     std::cout<< "sens  "<< sqrt(integ) <<"  "  <<sqrt(hg->Integral())<<std::endl;


     //-----------------------------------------------------------------------------------------------------------------
    // TF1 *func = new TF1("Sensetivity",Sensetivity,-1,1,1);
    // func->SetParameter(0,1);
    //unc->Draw();
    // TH1 * hfmbetaM = fmbetaM->GetHistogram();
    // hfmbetaM->SetName("MomentbetaM");
   

 hf->Draw();
  hg->Draw("same");
 
}

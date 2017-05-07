#include <stdio.h>
#include <iostream.h>

void MakeSimplePlots(){

  TFile *_file0 = TFile::Open("output.root");

  TH1F *rho_plus = (TH1F *)_file0->Get("rho_plus");
  TH1F *rho_minus = (TH1F *)_file0->Get("rho_minus");
  rho_plus->Scale(1./rho_plus->Integral());
  rho_minus->Scale(1./rho_minus->Integral());


  TH1F *pi_plus = (TH1F *)_file0->Get("pi_plus");
  TH1F *pi_minus = (TH1F *)_file0->Get("pi_minus");
  pi_plus->Scale(1./pi_plus->Integral());
  pi_minus->Scale(1./pi_minus->Integral());


  TH1F *mu_plus = (TH1F *)_file0->Get("mu_plus");
  TH1F *mu_minus = (TH1F *)_file0->Get("mu_minus");
  mu_plus->Scale(1./mu_plus->Integral());
  mu_minus->Scale(1./mu_minus->Integral());


  TH1F *piom_plus = (TH1F *)_file0->Get("piom_plus");
  TH1F *piom_minus = (TH1F *)_file0->Get("piom_minus");
  piom_plus->Scale(1./piom_plus->Integral());
  piom_minus->Scale(1./piom_minus->Integral());

  TH1F *om_plus = (TH1F *)_file0->Get("om_plus");
  TH1F *om_minus = (TH1F *)_file0->Get("om_minus");
  om_plus->Scale(1./om_plus->Integral());
  om_minus->Scale(1./om_minus->Integral());


  TH1F *ommurho_plus = (TH1F *)_file0->Get("ommurho_plus");
  TH1F *ommurho_minus = (TH1F *)_file0->Get("ommurho_minus");
  ommurho_plus->Scale(1./ommurho_plus->Integral());
  ommurho_minus->Scale(1./ommurho_minus->Integral());


  TCanvas *c1 = new TCanvas("c1","c1",100,100,600,600); 
  rho_plus->SetLineColor(2);
  rho_minus->Draw();
  rho_plus->Draw("same");



  TCanvas *c2 = new TCanvas("c2","c2",100,100,600,600); 
  pi_plus->SetLineColor(2);
  pi_minus->Draw();
  pi_plus->Draw("same");


  TCanvas *c3 = new TCanvas("c3","c3",100,100,600,600); 
  mu_plus->SetLineColor(2);
  mu_minus->Draw();
  mu_plus->Draw("same");


  TCanvas *c4 = new TCanvas("c4","c4",100,100,600,600); 
  piom_plus->SetLineColor(2);
  piom_minus->Draw();
  piom_plus->Draw("same");



  TCanvas *c5 = new TCanvas("c5","c5",100,100,600,600); 
  om_plus->SetLineColor(2);
  om_minus->Draw();
  om_plus->Draw("same");


  TCanvas *c6 = new TCanvas("c6","c6",100,100,600,600); 
  ommurho_plus->SetLineColor(2);
  ommurho_minus->Draw();
  ommurho_plus->Draw("same");

  cout<<" KS:   pi  " << pi_plus->KolmogorovTest(pi_minus)<<endl; 
  cout<<" KS:   mu  " << mu_plus->KolmogorovTest(mu_minus)<<endl; 
  cout<<" KS:   piom  " << piom_plus->KolmogorovTest(piom_minus)<<endl; 
  cout<<" KS:   rho  " << rho_plus->KolmogorovTest(rho_minus)<<endl; 
  cout<<" KS:   ommurho  " << ommurho_plus->KolmogorovTest(ommurho_minus)<<endl; 
  cout<<" KS:   om  " << om_plus->KolmogorovTest(om_minus)<<endl; 

  cout<<" Chi2:   pi  " << pi_plus->Chi2Test(pi_minus,"NORM")<<endl; 



  TCanvas *c7 = new TCanvas("c7","c7",100,100,600,600); 
  TH1F *div1=new TH1F("div1","div1",50,-1,1);
  div1->Divide(rho_plus,rho_minus);

  TH1F *div2=new TH1F("div2","div2",50,-1,1);
  div2->Divide(ommurho_plus,ommurho_minus);
  div1->Draw();
  div2->Draw("same");


}

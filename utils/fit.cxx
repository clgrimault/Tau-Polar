#include "stdlib.h"
#include "TH1.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH2.h"
#include "TTree.h"
#include "TComplex.h"
#include "TMath.h"
#include "TMinuit.h"

static int NUMBER_OF_BINS=40;
static int FirstBin =  1;
static int LastBin  =  NUMBER_OF_BINS ;
static int Min =-1.1;
static int Max = 1.1;
double NMC;


  TString HPlusN = "omega_a1_plus";
  TString HMinsN = "omega_a1_minus";

TString File = "MixedTemplates.root";
TFile *datafile = TFile::Open(File);
TH1F *data  = (TH1F *)datafile->Get("mixed");
TH1F *mins = (TH1F *)datafile->Get(HMinsN);
TH1F *plus= (TH1F *)datafile->Get(HPlusN);

  
double HPlus(0);
double HMins(0);
double sigmap(0);

void PerformHisto(){
  for(unsigned int i=FirstBin; i< LastBin; i++){HPlus+=plus->GetBinContent(i);}
  for(unsigned int i=FirstBin; i< LastBin; i++){HMins+=mins->GetBinContent(i);}
  NMC = data->Integral(FirstBin,LastBin);
}


double N_mc(int iBin, double *par){
  double er=1;//0.612;
  double p = (er-(1-par[0])/(1+par[0]))/(er+(1-par[0])/(1+par[0]));
   return NMC*(plus->GetBinContent(iBin[0])*(1+p)/2/HPlus + mins->GetBinContent(iBin[0])*(1-p)/2/HMins)  ;
}

double N_sig(int iBin, double *par){
     double er=1;
     double p = (er-(1-par[0])/(1+par[0]))/(er+(1-par[0])/(1+par[0]));
     return NMC*(plus->GetBinContent(iBin[0])*(1+p)/2/HPlus + mins->GetBinContent(iBin[0])*(1-p)/2/HMins);
}

int N_data(int iBin){
  return data->GetBinContent(iBin[0]);
}


void fit(string type, string ScanOrFit = "f"){
  PerformHisto();

  gStyle->SetOptStat(0000);
  Double_t par[1],dpar[1];
  TMinuit *fitter = new TMinuit(1);

  if(type == 'c')   fitter->SetFCN(Chi2); fitter->SetErrorDef(1);
  //if(type == 'l')   fitter->SetFCN(Likelihood); fitter->SetErrorDef(0.5);
  if(type == 'l')   fitter->SetFCN(LikelihoodAdv); fitter->SetErrorDef(0.5);
 
  Double_t arglist[2];
  Int_t ierflg = 1;

  fitter->mnparm( 0, "polarization",  -0.15,   0.005,   -0.99, 0.99, ierflg);
  fitter->SetPrintLevel(2);

  arglist[0]=500;
  arglist[1]=500;
 
  fitter->mnexcm("SET STR", arglist ,1,ierflg);


 fitter->Migrad();
 //fitter->Simplex();  
 fitter->mnhess();
 fitter->mnhess();
 fitter->mnhess();
 fitter->mnmnos();
  
  get_parameters(fitter,par,dpar);
 
  TH1F *Sum = new TH1F("Sum","Sum",NUMBER_OF_BINS,-1.1,1.1);
  Sum->Sumw2();
  plus->SetLineColor(3);
  mins->SetLineColor(2);
  data->SetLineWidth(2);
  
  plus->Clone("plusClone");
  mins->Clone("minsClone");
  
  
  
  plusClone->SetName("plusClone");
  minsClone->SetName("minsClone");
  double er=1;
  double p = (er-(1-par[0])/(1+par[0]))/(er+(1-par[0])/(1+par[0]));

   plusClone->Scale(NMC*(1+p)/2/HPlus);
   minsClone->Scale(NMC*(1-p)/2/HMins);

  double  llcheck =0;
  for(int ib = FirstBin; ib <= LastBin; ib++ ){
    
        Sum->SetBinContent(ib, NMC*(plus->GetBinContent(ib)*(1+p)/2/HPlus + mins->GetBinContent(ib)*(1-p)/2/HMins)   );
    
    double binerr1 = NMC*sqrt(plus->GetBinContent(ib)*(1+par[0])/2 + mins->GetBinContent(ib)*(1-par[0])/2 )/(HMins+HPlus);
 


    Sum->SetBinError(ib, sqrt(Sum->GetBinContent(ib)) );
     if(Sum->GetBinContent(ib)!=0)  llcheck+= (data->GetBinContent(ib) - Sum->GetBinContent(ib))*(data->GetBinContent(ib) - Sum->GetBinContent(ib))/Sum->GetBinContent(ib);
  }
  
  
   data->SetStats(0);
   data->SetLineWidth(2);
   data->SetMarkerStyle(20);
   data->SetMarkerSize(0.8);
   
   data->SetTitle("");
   data->SetYTitle("");
   data->SetXTitle("E_{#mu}/E_{#tau}");
   
   Sum->SetLineWidth(2);
   Sum->SetLineColor(4);
   // TCanvas *c = new TCanvas("c","PolFit",100,100,600,600);

   minsClone->SetLineColor(2);
   minsClone->SetLineStyle(7);
   minsClone->SetLineWidth(2);
   minsClone->SetMarkerStyle(20);
   minsClone->SetMarkerSize(1.2);
   minsClone->GetXaxis()->SetTitle("#omega_{a1}");
   minsClone->GetXaxis()->SetLabelFont(42);
   minsClone->GetXaxis()->SetLabelSize(0.05);
   minsClone->GetXaxis()->SetTitleSize(0.05);
   minsClone->GetXaxis()->SetTitleOffset(1.4);
   minsClone->GetXaxis()->SetTitleFont(42);
   minsClone->GetYaxis()->SetLabelFont(42);
   minsClone->GetYaxis()->SetLabelSize(0.05);
   minsClone->GetYaxis()->SetTitleSize(0.05);
   minsClone->GetYaxis()->SetTitleOffset(1.4);
   minsClone->GetYaxis()->SetTitleFont(42);
   minsClone->GetZaxis()->SetLabelFont(42);
   minsClone->GetZaxis()->SetLabelSize(0.05);
   minsClone->GetZaxis()->SetTitleSize(0.05);
   minsClone->GetZaxis()->SetTitleFont(42);
   plusClone->SetLineColor(3);
   plusClone->SetLineStyle(7);
   plusClone->SetLineWidth(2);
   plusClone->SetMarkerStyle(20);
   plusClone->SetMarkerSize(1.2);
   plusClone->GetXaxis()->SetTitle("#omega_{a1}");
   plusClone->GetXaxis()->SetLabelFont(42);
   plusClone->GetXaxis()->SetLabelSize(0.05);
   plusClone->GetXaxis()->SetTitleSize(0.05);
   plusClone->GetXaxis()->SetTitleOffset(1.4);
   plusClone->GetXaxis()->SetTitleFont(42);
   plusClone->GetYaxis()->SetLabelFont(42);
   plusClone->GetYaxis()->SetLabelSize(0.05);
   plusClone->GetYaxis()->SetTitleSize(0.05);
   plusClone->GetYaxis()->SetTitleOffset(1.4);
   plusClone->GetYaxis()->SetTitleFont(42);
   plusClone->GetZaxis()->SetLabelFont(42);
   plusClone->GetZaxis()->SetLabelSize(0.05);
   plusClone->GetZaxis()->SetTitleSize(0.05);
   plusClone->GetZaxis()->SetTitleFont(42);
   Sum->SetStats(0);
   Sum->SetFillColor(5);
   Sum->SetLineWidth(2);
   Sum->GetXaxis()->SetTitle("#omega_{a_{1}}");
   Sum->GetXaxis()->SetTitleSize(0.05);
   Sum->GetXaxis()->SetTitleOffset(0.8);
   Sum->GetYaxis()->SetRange(1,1);

   Sum->Draw("E2");
   data->Draw("SAME");
   minsClone->Draw("HISTSAME");
   plusClone->Draw("HISTSAME");
      // TFile *out = new TFile("outhistos.root","RECREATE");
      // data->SetName("data");
      // TH1F *data1  = (TH1F*) data->Clone("data1");
      // data1->Write();
      // Sum->Write();
      // TH1F *minsClone1=(TH1F*)minsClone->Clone("minsClone1");
      // minsClone1->Write();
      // background->Write();
      // TH1F *plusClone1=(TH1F*)plusClone->Clone("plusClone1");
      // plusClone1->Write();
      // out->Write();
      // out->Close();

   
    
   TLegend *leg = new TLegend(0.5401149,0.4957627,0.9465517,0.8432203,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(82);
   leg->SetLineColor(0);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("NULL","Data","lpf");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Sum","Fit","f");
   entry->SetFillStyle(1001);
   entry->SetLineColor(4);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("plusClone","h_{#tau}  = +1","f");
   entry->SetFillStyle(1001);
   entry->SetLineColor(3);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("minsClone","h_{#tau}  = -1","f");
   entry->SetFillStyle(1001);
   entry->SetLineColor(2);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   
     // leg->Draw();
   
   
   TCanvas *c2 = new TCanvas("c2","Scan POlarization",100,100,500,500);
   // if(ScanOrFit == "s")
   {
     fitter->Command("SCAn 1");
     TGraph *gr = (TGraph*)fitter->GetPlot();
     gr->GetXaxis()->SetTitle("P_{#tau}");
     gr->GetYaxis()->SetTitle("#chi^{2}");
     
     gr->Draw(""); 
   }
   
   //       TCanvas *c3 = new TCanvas("c3","Scan p",100,100,500,500);
   
   //       //     if(ScanOrFit == "s")
   // {
   //   //  fitter->mnimpr() ;
   //    fitter->Command("SCAn 2");
   
   //    //   fitter->Comman("IMProve");
   //    TGraph *gr = (TGraph*)fitter->GetPlot();
   //    gr->GetXaxis()->SetTitle("backgound p");
   //    gr->GetYaxis()->SetTitle("#chi^{2}");
   
   //    gr->Draw("alp"); 
   //  }
   
   // // if(type == 'f')
   // //   {
   // 	 TCanvas *c4 = new TCanvas("c4","contour",100,100,500,500);
   // 	 // fitter->Command("SCAn 1");
   // 	 fitter->SetErrorDef(0.5);
   // 	 TGraph *gr = (TGraph*)fitter->Contour(10,0,1);
   // 	 gr->GetXaxis()->SetTitle("P_{#tau}");
   // 	 gr->GetYaxis()->SetTitle("p_bkg");
   
   // 	 gr->Draw("alp"); 
   // 	 //}
}

void get_parameters(TMinuit * minuit,double *par, double *par_err ){
  for( int i =0 ; i <1; i++){
    minuit->GetParameter(i,par[i], par_err[i]);
  }
}
//  Minimization function definition

void Likelihood(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
   double LL=0;
   double LLtot=0;
   int bin[1]={0};
   double mu = ScaleBackground;
   double B_i(0);
   double N_i(0);
   double S_i(0);
   double Sum(0);
   //  printf("---------------- \n");
  for ( int i = 1; i< 20 + 1; i++){

     bin[0] = i;
     //     B_i = N_bkgns(i);
     N_i = N_data(i);
     S_i = N_sig(i, par); 

     double  mc = N_mc(i, par);
  
     int scale = 1;
     if(mc==0){ scale =0; mc =1;}
     else scale = 1;
     LL+= log(TMath::Poisson(N_i, mc))* TMath::Poisson(B_i,N_bkg(i)*ScaleBackground );// (N_data(i)*log(mc)  - scale*mc) + N_bkgns(i)*log(N_bkg(i)*ScaleBackground) - N_bkg(i)*ScaleBackground;
     //  cout<<" =========================================================================================================================  i " << i << "  " << LL <<endl;
   }

  LL =( LL) - (par[1]-1)*(par[1]-1)/2/0.1/0.1;
  f= -LL;
 }

void LikelihoodAdv(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  double LL=0;
  double LLtot=0;
  int bin[1]={0};
  double mu = 1.;//ScaleBackground;
  double mu = ScaleBackground;
  //double B_i(0);
  double N_i(0);
  double S_i(0);
  
  //printf("----------------  %f \n",mu);
  
  for ( int i = 1; i< 20 + 1; i++){
    bin[0] = i;
    B_i = 0.99*N_bkgns(i);
    N_i = N_data(i);
    S_i = N_sig(i, par);
    double s(0);
    double Sum(0);
    for(int k=0; k< N_i; k++){
      double s(0);
      double nom1 = TMath::Poisson(N_i -k, S_i);
      double nom2 = TMath::Poisson(k, par[1]);
      double nom3 =TMath::Poisson(B_i, mu);
      double denom = TMath::Poisson(B_i+k, mu+par[1])/(mu+par[1]);
	       
      // if(nom2<pow(1,-100))nom2=nom2*pow(10,100);
      // if(denom<pow(1,-100))denom=denom*pow(10,100);


         // cout<<"---  K"<<k<< "  N_i " << N_i <<endl<<endl;
         // cout<<"nom1  "<<nom1<<endl; 
         // cout<<"nom2  "<<nom2<<endl; 
         // cout<<"nom3  "<<nom3<<endl; 
         // cout<<"denom  "<<denom<<endl; 


	 if(denom!=0){ s= (nom1*nom2/denom/(mu+par[1]))*nom3;}//TMath::Poisson(N_i -k, S_i)*TMath::Poisson(k, par[1])*TMath::Poisson(B_i, mu)/TMath::Poisson(B_i+k, mu+par[1])/(mu+par[1]);//  pois(s,N-k)*pois(p,k)*pois(mu,B)/pois(mu+p, B+k);
       else{s==0;}
      if(s < pow(10.,-50)) s==0;
      if(s==s && Sum==Sum) Sum+=s;
      //      if(s==s && Sum==Sum)  
      //    cout<<" i  "<< i  << " k   "<< k  <<"  s  "<<s  << "  Sum "<< Sum  << "  log(Sum)  " <<log(Sum)<<endl;
    }
    //    LL+= log(Sum)*TMath::Gaus(par[1],0.57,sigmap,true)); // if mu  =1 
    // cout<<" i -----" << i <<endl;
    //   cout<<" =========================================================================================================================i  " <<  i << "  " << (Sum) <<endl;
    LL+= log(Sum);//*TMath::Gaus(par[1],1.,0.2,true)); // if mu = scaleBackgound
   
  }
  LL = (LL) - (par[1]-1.)*(par[1]-1.)/2/0.1/0.1;
  f= -LL;
}


void Chi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  double LL=0;
  double LLtot=0;
  int count =0;
  for ( int i = FirstBin; i<=LastBin; i++){
    double  mc = N_mc(i, par);
    if(mc ==0 || N_data(i)==0) continue;
    LL+= (mc -N_data(i))*(mc - N_data(i))/N_data(i);
  }
  //f= LL;    
  LLtot = LL;
  f= LLtot/1.0;
}



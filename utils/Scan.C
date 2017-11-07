
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
#include "TRandom.h"

static int NUMBER_OF_BINS=50;
//static int NUMBER_OF_BINS=10;

static int Min =-1.;
static int Max = 1;

static float Pol = 0;
static float Lum = 50;


 TFile *_file0 = TFile::Open("1.root");

// TFile *data = TFile::Open("LOCAL_COMBINED_templatesdata_default.root");



//        TH1F *template_minus15 = (TH1F *)templateSignalPlus->Get("polarizationtemplates_default_x2minusSignal");
//        TH1F *template_plus15= (TH1F *)templateSignalPlus->Get("polarizationtemplates_default_x2plusSignal");

TH1F *template_minus153pi = (TH1F *)_file0->Get("pi_minus");
TH1F *template_plus153pi= (TH1F *)_file0->Get("pi_plus");

TH1F *template_minus15Incl = (TH1F *)_file0->Get("pi_minus");
TH1F *template_plus15Incl= (TH1F *)_file0->Get("pi_plus");

TH1F *template_minus154pi = (TH1F *)_file0->Get("pi_minus");
TH1F *template_plus154pi= (TH1F *)_file0->Get("pi_plus");


// TH1F *h_data15  = (TH1F *)data->Get("templatesdata_default_RecOmegaA1pData");

TH1F *h3pi = new TH1F("h3pi","h3pi",50,-1.0,1.0);
TH1F *h4pi = new TH1F("h4pi","h4pi",50,-1.0,1.0);
TH1F *hIncl = new TH1F("hIncl","hincl",50,-1.0,1.0);


TH1F *plus3pi  = new TH1F("plus3pi","plus3pi",50,-1.0,1.0);
TH1F *mins3pi  = new TH1F("mins3pi","mins3pi",50,-1.0,1.0);


TH1F *plus4pi  = new TH1F("plus4pi","plus4pi",50,-1.0,1.0);
TH1F *mins4pi  = new TH1F("mins4pi","mins4pi",50,-1.0,1.0);


TH1F *plusIncl = new TH1F("plusIncl","plusIncl",50,-1.0,1.0);
TH1F *minsIncl = new TH1F("minsIncl","minsIncl",50,-1.0,1.0);


plus3pi ->Sumw2();
mins3pi ->Sumw2();
	      
	      
plus4pi ->Sumw2();
mins4pi ->Sumw2();
	      
	      
plusIncl->Sumw2();
minsIncl->Sumw2();


static float Scale = 1;// Lum*h_data15->Integral();

void PerformHisto(){

 
  for(int ib = 1; ib <NUMBER_OF_BINS + 1; ib++ ){
    plus3pi->SetBinContent(ib,0);
    mins3pi->SetBinContent(ib,0);
    plus3pi->SetBinContent(ib, template_plus153pi->GetBinContent(ib));
    mins3pi->SetBinContent(ib, template_minus153pi->GetBinContent(ib));
    plus3pi->SetBinError(ib, sqrt(template_plus153pi->GetBinContent(ib)));
    mins3pi->SetBinError(ib, sqrt(template_minus153pi->GetBinContent(ib)));

    plusIncl->SetBinContent(ib,0);
    minsIncl->SetBinContent(ib,0);
    plusIncl->SetBinContent(ib, template_plus15Incl->GetBinContent(ib));
    minsIncl->SetBinContent(ib, template_minus15Incl->GetBinContent(ib));
    plusIncl->SetBinError(ib, sqrt(template_plus15Incl->GetBinContent(ib)));
    minsIncl->SetBinError(ib, sqrt(template_minus15Incl->GetBinContent(ib)));



    plus4pi->SetBinContent(ib,0);
    mins4pi->SetBinContent(ib,0);
    plus4pi->SetBinContent(ib, template_plus154pi->GetBinContent(ib));
    mins4pi->SetBinContent(ib, template_minus154pi->GetBinContent(ib));
    plus4pi->SetBinError(ib, sqrt(template_plus154pi->GetBinContent(ib)));
    mins4pi->SetBinError(ib, sqrt(template_minus154pi->GetBinContent(ib)));












  }

  plus3pi ->Scale(1/plus3pi->Integral()); 
  mins3pi ->Scale(1/mins3pi->Integral());

  plus4pi ->Scale(1/plus4pi->Integral()); 
  mins4pi ->Scale(1/mins4pi->Integral());

  plusIncl ->Scale(1/plusIncl->Integral()); 
  minsIncl ->Scale(1/minsIncl->Integral());


  h3pi->Divide(mins3pi,plus3pi);
  h4pi->Divide(mins4pi,plus4pi);
  hIncl->Divide(minsIncl,plusIncl);


  h4pi->SetLineColor(2);

  hIncl->SetLineWidth(2);
  hIncl->SetMarkerStyle(21);


}

 


double sens3pi[21];
double esens3pi[21];
double sensincl[21];
double esensincl[21];
double polX[21];
double epolX[21];


void Scan(){
  double polari=-1.0;
  double poldel = 0.1;
  for(int ipol =0; ipol < 21; ipol++){

    Pol = polari +  poldel*ipol;
    
    polX[ipol] = Pol;
    epolX[ipol] = 0;
    PerformHisto();
    TH1F *FakeData3pi = new TH1F("FakeData3pi","FakeData3pi",50,-1,1.);
    TH1F *FakeData4pi = new TH1F("FakeData4pi","FakeData4pi",50,-1,1.);
    TH1F *FakeDataIncl = new TH1F("FakeDataIncl","FakeDataIncl",50,-1,1.);
    TH1F *FakeData2 = new TH1F("FakeData2","FakeData2",50,-1.0,1.);
    //   FakeData4pi->Sumw2();
    //   FakeDataIncl->Sumw2();
    //   FakeData2->Sumw2();
    //   FakeData3pi->Sumw2();
    FakeData3pi->Add(plus3pi,0.5*(1+Pol));
    FakeData3pi->Add(mins3pi,0.5*(1-Pol));
    FakeData3pi->Scale(Scale);
    

    FakeData4pi->Add(plus4pi,0.5*(1+Pol));
    FakeData4pi->Add(mins4pi,0.5*(1-Pol));
    FakeData4pi->Scale(Scale);
    
    FakeDataIncl->Add(plusIncl,0.5*(1+Pol));
    FakeDataIncl->Add(minsIncl,0.5*(1-Pol));
    FakeDataIncl->Scale(Scale);



  //  FakeDataClone =   FakeData->Clone();
  //  FakeDataClone->SetName("FakeDataClone");
  TRandom *rd = new TRandom();
  for(int ibin=1; ibin < FakeData3pi->GetNbinsX(); ibin++ ){
    // FakeData2->SetBinContent(ibin, rd->Poisson(FakeData->GetBinContent(ibin)));
    //    std::cout<<" nevents "<< rd->Poisson(FakeData->GetBinContent(ibin))  << " FakeData->GetBinContent(ibin)   " <<FakeData->GetBinContent(ibin) <<std::endl;
  }
  //  std::cout<<" nevents "<< FakeData2->GetEntries() << std::endl;
  TFile *out= new TFile("sensoutput.root","RECREATE");


//cout<<rd->Poisson(5)<<endl;
  FakeData3pi->Scale(1/( FakeData3pi->GetBinContent(1)+FakeData3pi->GetBinContent(2)+FakeData3pi->GetBinContent(3)+FakeData3pi->GetBinContent(4)+FakeData3pi->GetBinContent(5)+FakeData3pi->GetBinContent(6)+FakeData3pi->GetBinContent(7)+FakeData3pi->GetBinContent(8)+FakeData3pi->GetBinContent(9)+FakeData3pi->GetBinContent(10) + FakeData3pi->GetBinContent(11)+FakeData3pi->GetBinContent(12)+FakeData3pi->GetBinContent(13)+FakeData3pi->GetBinContent(14)+FakeData3pi->GetBinContent(15)+FakeData3pi->GetBinContent(16)+FakeData3pi->GetBinContent(17)+FakeData3pi->GetBinContent(18)+FakeData3pi->GetBinContent(19)+FakeData3pi->GetBinContent(20))/0.1 );

  FakeDataIncl->Scale(1/( FakeDataIncl->GetBinContent(1)+FakeDataIncl->GetBinContent(2)+FakeDataIncl->GetBinContent(3)+FakeDataIncl->GetBinContent(4)+FakeDataIncl->GetBinContent(5)+FakeDataIncl->GetBinContent(6)+FakeDataIncl->GetBinContent(7)+FakeDataIncl->GetBinContent(8)+FakeDataIncl->GetBinContent(9)+FakeDataIncl->GetBinContent(10) + FakeDataIncl->GetBinContent(11)+FakeDataIncl->GetBinContent(12)+FakeDataIncl->GetBinContent(13)+FakeDataIncl->GetBinContent(14)+FakeDataIncl->GetBinContent(15)+FakeDataIncl->GetBinContent(16)+FakeDataIncl->GetBinContent(17)+FakeDataIncl->GetBinContent(18)+FakeDataIncl->GetBinContent(19)+FakeDataIncl->GetBinContent(20))/0.1 );

  FakeData4pi->Scale(1/( FakeData4pi->GetBinContent(1)+FakeData4pi->GetBinContent(2)+FakeData4pi->GetBinContent(3)+FakeData4pi->GetBinContent(4)+FakeData4pi->GetBinContent(5)+FakeData4pi->GetBinContent(6)+FakeData4pi->GetBinContent(7)+FakeData4pi->GetBinContent(8)+FakeData4pi->GetBinContent(9)+FakeData4pi->GetBinContent(10) + FakeData4pi->GetBinContent(11)+FakeData4pi->GetBinContent(12)+FakeData4pi->GetBinContent(13)+FakeData4pi->GetBinContent(14)+FakeData4pi->GetBinContent(15)+FakeData4pi->GetBinContent(16)+FakeData4pi->GetBinContent(17)+FakeData4pi->GetBinContent(18)+FakeData4pi->GetBinContent(19)+FakeData4pi->GetBinContent(20))/0.1 );



//   mins3pi->Scale(Scale);
//   plus3pi->Scale(Scale);
//   plus3pi->Write();
//   mins3pi->Write();

//   mins4pi->Scale(Scale);
//   plus4pi->Scale(Scale);
//   plus4pi->Write();
//   mins4pi->Write();

//   minsIncl->Scale(Scale);
//   plusIncl->Scale(Scale);
//   plusIncl->Write();
//   minsIncl->Write();



//   FakeData3pi->Write();
//   FakeData4pi->Write();
//   FakeDataIncl->Write();
//   h3pi->Write();
//   h4pi->Write();
//   hIncl->Write();

//   out->Write();
//   out->Close();
//   FakeDataIncl->SetLineWidth(2);
//   FakeDataIncl->SetLineWidth(2);
//   FakeDataIncl->SetMarkerStyle(21);

//   FakeData4pi->SetLineColor(2);
//   FakeData3pi->Draw();
//   FakeData4pi->Draw("same");
//   FakeDataIncl->Draw("same");

//   TCanvas *c = new TCanvas("c","c",300,300,600,600);
//   h3pi->Draw();
//   h4pi->Draw("same");
//   hIncl->Draw("same");

  double branch = 0.0911;
  double delta = 0.1;
  double omega;
  double integ=0;
  double err  =0;
  double F;
  for(int ibin=0; ibin < FakeData3pi->GetNbinsX(); ibin++){
    omega = -0.95 + ibin*delta;
    F = FakeData3pi->GetBinContent(ibin);
    eF = FakeData3pi->GetBinError(ibin);
    integ += omega*omega*F*delta/(1+Pol*omega)/(1+Pol*omega);
    err+= pow(omega*omega*eF*delta/(1+Pol*omega)/(1+Pol*omega),2);
    //  std::cout<<" ibin:  " << ibin<<"  integ 3pi  "<< omega*omega*F*delta/(1+Pol*omega)/(1+Pol*omega)<< std::endl;

  }
  sens3pi[ipol]=sqrt(integ);
  esens3pi[ipol]=sqrt(err);
  //  std::cout<<"final 3pi: "<< integ << " pm " <<sqrt(err) <<std::endl;




   delta = 0.1;
   omega =0;
   integ=0;
   err  =0;
   F=0;
   branch = 0.0455  + 0.0911;
  for(int ibin=0; ibin < FakeDataIncl->GetNbinsX(); ibin++){
    omega = -0.95 + ibin*delta;
    F = FakeDataIncl->GetBinContent(ibin);
    eF = FakeDataIncl->GetBinError(ibin);
    integ += omega*omega*F*delta/(1+Pol*omega)/(1+Pol*omega);
    err+= pow(omega*omega*eF*delta/(1+Pol*omega)/(1+Pol*omega),2);
    //     std::cout<<"omega  "<< integ << std::endl;

  }
    
  //std::cout<<"final Incl: "<< integ << " pm " <<sqrt(err) <<std::endl;
  sensincl[ipol]=sqrt(integ);;
  esensincl[ipol]=sqrt(err);
  std::cout<<"sens "<< sensincl[ipol] << " pm " <<sqrt(err) << "    Pol  " << Pol<<std::endl;

   delta = 0.1;
   omega =0;
   integ=0;
   err  =0;
   F=0;
   branch = 0.0455;
  for(int ibin=0; ibin < FakeDataIncl->GetNbinsX(); ibin++){
    omega = -0.95 + ibin*delta;
    F = FakeData4pi->GetBinContent(ibin);
    eF = FakeData4pi->GetBinError(ibin);
    integ += omega*omega*F*delta/(1+Pol*omega)/(1+Pol*omega);
    err+= pow(omega*omega*eF*delta/(1+Pol*omega)/(1+Pol*omega),2);
    //   std::cout<<" ibin:  " << ibin <<"  integ 4pi  "<< omega*omega*F*delta/(1+Pol*omega)/(1+Pol*omega)<< std::endl;
    //  std::cout<<"omega  "<< integ << std::endl;

  }
    
  //  std::cout<<"final 4pi: "<< integ << " pm " <<sqrt(err) <<std::endl;



  }

  TGraphErrors *gr3pi = new TGraphErrors(21,polX,sens3pi,epolX,esens3pi );
  gr3pi->Draw();
  // TGraphErrors *grincl = new TGraphErrors(21,polX,sensincl,epolX,esensincl );
  // grincl->Draw("same");

}

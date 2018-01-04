

void Plotter(){

  TFile *_file0 = TFile::Open("Draw.root");
  TH1F *hm = (TH1F *)_file0->Get("omegabar_rho_minus");
  TH1F *hp = (TH1F *)_file0->Get("omegabar_rho_plus");

  TString Title("#tau #rightarrow #rho #nu");
  TString LabelMinus("#tau_{L}");
  TString LabelPlus("#tau_{R}");


  hm->SetStats(0);
  hp->SetStats(0);
  hm->SetTitle(0);
  hp->SetTitle(0);
  


  Int_t ci;      // for color index setting
  TColor *color; // for color definition with alpha
  ci = TColor::GetColor("#333333");
  hm->SetLineColor(ci);
  hm->SetLineWidth(2);
  hm->GetXaxis()->SetLabelFont(42);
  hm->GetXaxis()->SetLabelSize(0.035);
  hm->GetXaxis()->SetTitleSize(0.035);
  hm->GetXaxis()->SetTitleFont(42);
  hm->GetYaxis()->SetLabelFont(42);
  hm->GetYaxis()->SetTitleSize(0.035);
  hm->GetYaxis()->SetTitleFont(42);
  hm->GetZaxis()->SetLabelFont(42);
  hm->GetZaxis()->SetLabelSize(0.035);
  hm->GetZaxis()->SetTitleSize(0.035);
  hm->GetZaxis()->SetTitleFont(42);


  ci = TColor::GetColor("#6666ff");
  hp->SetLineColor(ci);
  hp->SetLineWidth(2);
  hp->GetXaxis()->SetLabelFont(42);
  hp->GetXaxis()->SetLabelSize(0.035);
  hp->GetXaxis()->SetTitleSize(0.035);
  hp->GetXaxis()->SetTitleFont(42);
  hp->GetYaxis()->SetLabelFont(42);
  hp->GetYaxis()->SetLabelSize(0.035);
  hp->GetYaxis()->SetTitleSize(0.035);
  hp->GetYaxis()->SetTitleFont(42);
  hp->GetZaxis()->SetLabelFont(42);
  hp->GetZaxis()->SetLabelSize(0.035);
  hp->GetZaxis()->SetTitleSize(0.035);
  hp->GetZaxis()->SetTitleFont(42);
   



  hm->DrawNormalized();
  hp->DrawNormalized("same");

  pt = new TPaveText(-0.9494986,0.02809341,-0.6224928,0.03286091,"br");
  pt->SetBorderSize(0);
  pt->SetFillColor(0);
  pt->SetFillStyle(0);

  ci = TColor::GetColor("#333333");
  pt->SetTextColor(ci);
  pt->SetTextFont(42);
  AText = pt->AddText(LabelMinus);
  pt->Draw();


  pt = new TPaveText(0.6540115,0.02815155,0.9810172,0.03291905,"br");
  pt->SetBorderSize(0);
  pt->SetFillColor(0);
  pt->SetFillStyle(0);
  
  ci = TColor::GetColor("#6666ff");
  pt->SetTextColor(ci);
  pt->SetTextFont(42);
  AText = pt->AddText(LabelPlus);
  pt->Draw();

  pt = new TPaveText(-0.2679083,0.04355873,0.2679083,0.04832624,"br");
  pt->SetBorderSize(0);
  pt->SetFillColor(0);
  pt->SetFillStyle(0);

  ci = TColor::GetColor("#660066");
  pt->SetTextColor(ci);
  pt->SetTextFont(42);
  AText = pt->AddText(Title);
  pt->Draw();


}

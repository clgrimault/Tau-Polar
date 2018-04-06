/**
 * Example of use of tauola C++ interfate. Pythia events are
 * generated with a stable tau. Taus are subseuently decay via
 * tauola, plots of polatization observables for tau-> mununu and tau-> pinu
 * are produced.
 *

 */ 
     
#include "Tauola/Log.h"
#include "Tauola/Tauola.h"
#include "Tauola/TauolaHepMCEvent.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h" 
#include "TH3F.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "UserCodes/a1Helper.h"
#include "UserCodes/TauDecaysHelper.h"
#include "UserCodes/TauPolInterface.h"
#include "UserCodes/PolarimetricA1.h"
#include "TLorentzVector.h"
#include "TComplex.h"
#include "TMatrixT.h"
#include "TVectorT.h"
#include "TMatrixTSym.h"
#include "TMath.h"
#include "limits.h"

//pythia header files
#ifdef PYTHIA8180_OR_LATER
#include "Pythia8/Pythia.h" 
#include "Pythia8/Pythia8ToHepMC.h"
#else 
#include "Pythia.h"
#include "HepMCInterface.h"
#endif

//MC-TESTER header files
#include "Generate.h"
#include "HepMCEvent.H"
#include "Setup.H"
 
#include "tauola_print_parameters.h"
using namespace std;
using namespace Pythia8; 
using namespace Tauolapp;


int NumberOfEvents = 20000; 
bool ApplyCut(false);
double pt_cut = 30;  // GeV - approximately correspond to CMS trigger, a cut on the visible decay products
double eta_cut = 10; // - very large value, all events pass - at the moment switched off at all; eta - pseudorapidity;
int EventsToCheck=5;

// elementary test of HepMC typically executed before
// detector simulation based on http://home.fnal.gov/~mrenna/HCPSS/HCPSShepmc.html
// similar test was performed in Fortran
// we perform it before and after Tauola (for the first several events)



void checkMomentumConservationInEvent(HepMC::GenEvent *evt)
{
 	cout<<"List of stable particles: "<<endl;

	double px=0.0,py=0.0,pz=0.0,e=0.0;

	int barcodetau1vertex(0),barcodetau2vertex(0);
	for ( HepMC::GenEvent::particle_const_iterator p = evt->particles_begin();  p != evt->particles_end(); ++p )
	  {
	    if( (*p)->status() == 1 )
	      {
		HepMC::FourVector m = (*p)->momentum();
		HepMC::GenVertex *productProductionvertex = (*p)->production_vertex();
		px+=m.px();
		py+=m.py();
		pz+=m.pz();
		e +=m.e();
	      }
	  }
	cout.precision(6);
	cout.setf(ios_base::floatfield);
}
TMatrixT<double> convertToMatrix(TVectorT<double> V){
  TMatrixT<double> M(V.GetNrows(),1);
  for(int i=0; i < M.GetNrows(); i++){
    M(i,0)=V(i);
  } return M;
}


TLorentzVector 
BoostR(TLorentzVector pB, TLorentzVector frame){
   TMatrixT<double> transform(4,4);
   TMatrixT<double> result(4,1);
   TVectorT<double> vec(4); 
   TVector3 b;
   if(frame.Vect().Mag()==0){ std::cout<<"RH Boost is not set, perfrom calculation in the Lab Frame   "<<std::endl; return pB;}
    if(frame.E()==0){ std::cout<<" Caution: Please check that you perform boost correctly!  " <<std::endl; return pB;} 
   else   b=frame.Vect()*(1/frame.E());
   vec(0)  = pB.E();    vec(1)  = pB.Px();
   vec(2)  = pB.Py();   vec(3)  = pB.Pz();
   double gamma  = 1/sqrt( 1 - b.Mag2());
   transform(0,0)=gamma; transform(0,1) =- gamma*b.X() ;  transform(0,2) =  - gamma*b.Y();  transform(0,3) = - gamma*b.Z(); 
   transform(1,0)=-gamma*b.X(); transform(1,1) =(1+ (gamma-1)*b.X()*b.X()/b.Mag2()) ;  transform(1,2) = ((gamma-1)*b.X()*b.Y()/b.Mag2());  transform(1,3) = ((gamma-1)*b.X()*b.Z()/b.Mag2());
   transform(2,0)=-gamma*b.Y(); transform(2,1) = ((gamma-1)*b.Y()*b.X()/b.Mag2());  transform(2,2) = (1 + (gamma-1)*b.Y()*b.Y()/b.Mag2());  transform(2,3) =  ((gamma-1)*b.Y()*b.Z()/b.Mag2()); 
   transform(3,0)=-gamma*b.Z(); transform(3,1) =((gamma-1)*b.Z()*b.X()/b.Mag2()) ;  transform(3,2) = ((gamma-1)*b.Z()*b.Y()/b.Mag2());  transform(3,3) = (1 + (gamma-1)*b.Z()*b.Z()/b.Mag2()); 
   result=transform*convertToMatrix(vec);
   return TLorentzVector(result(1,0), result(2,0) ,result(3,0), result(0,0));
}

TVector3
Rotate(TVector3 LVec, TVector3 Rot){
  TVector3 vec = LVec;
  vec.RotateZ(0.5*TMath::Pi() - Rot.Phi());  // not 0.5, to avoid warnings about 0 pT
  vec.RotateX(Rot.Theta());
  return vec;
}


void redMinus(TauolaParticle *minus) // this is JAK1
{
  //   
  // this method can be used to redefine branching ratios in decay of tau-
  // either generally, or specific to  tau- with pointer *minus.
  //
  // Pointer *minus can be used to define properties of decays for taus
  // at specific point(s) in the event tree. Example: 
  // vector<TauolaParticle*> x=minus->getMothers();
  // and define special versions depending on x. 
  //
  // Any combination of methods
  // Tauola::setTauBr(int mode, double br);
  //Tauola::setTaukle(double bra1, double brk0, double brk0b,double brks);
   Tauola::setTaukle(1, 0,0,0);
  // can be called here 

   //dec ==2  - muon decay
   //dec ==3  - pion decay
   //dec ==4  - rho decay
   //dec ==5  - a1 decay
   for(unsigned int dec=1; dec <23; dec++){
      double br =0.0; 
      if(dec ==2|| dec ==3|| dec ==4|| dec ==5) br=0.25;
      Tauola::setTauBr(dec, br);
   }

}

void redPlus(TauolaParticle *plus) // this is JAK2
{
  //   
  // this method can be used to redefine branching ratios in decay of tau+
  // either generally, or specific to  tau+ with pointer *plus.
  //
  // Pointer *plus can be used to define properties of decays for tau
  // at specific point(s) in the event tree. Example: 
  // vector<TauolaParticle*> x=plus->getMothers();
  // and define special versions depending on x. 
  //
  // Any combination of methods
  // Tauola::setTauBr(int mode, double br);
  //Tauola::setTaukle(double bra1, double brk0, double brk0b,double brks);
   Tauola::setTaukle(1, 0,0,0);
  // can be called here 
  for(unsigned int dec=1; dec <23; dec++){
     double br =0.0;
      if(dec ==2|| dec ==3|| dec ==4|| dec ==5) br=0.25;
     Tauola::setTauBr(dec, br);
   }
}

void SortPions(std::vector<HepMC::GenParticle > pionsvec)
{

  int npim(0),npip(0);
  int    OSMCPionIndex(0);
  int    SSMCPion1Index(0);
  int    SSMCPion2Index(0);

    HepMC::GenParticle os;
    HepMC::GenParticle ss1;
    HepMC::GenParticle ss2;
    for(int i=0; i<pionsvec.size(); i++){
      if( pionsvec.at(i).pdg_id()== 211) npip++;
      if( pionsvec.at(i).pdg_id()==-211) npim++;
    }
    if(npip == 1 && npim==2){
      int nss=0;
      for(int i=0; i<pionsvec.size(); i++){
	if(pionsvec.at(i).pdg_id()== 211){
	  OSMCPionIndex=i;
	}
	if(pionsvec.at(i).pdg_id()==-211 && nss ==0){
	  nss++;
	  SSMCPion1Index=i;
	}
	if(pionsvec.at(i).pdg_id()==-211 &&  nss == 1){
	  SSMCPion2Index=i;
	}
      }
    }
    if( npip== 2 && npim==1){
      int nss=0;
      for(int i=0; i<pionsvec.size(); i++){
	if(pionsvec.at(i).pdg_id()== -211){
	  
	  OSMCPionIndex=i;
	}
	if(pionsvec.at(i).pdg_id()==211 && nss ==0){
	  nss++;
	  SSMCPion1Index=i;
	}
	if(pionsvec.at(i).pdg_id()==211 &&  nss == 1){
	  SSMCPion2Index=i;
	}
      }
    }
    os=pionsvec.at(OSMCPionIndex);
    ss1=pionsvec.at(SSMCPion1Index);
    ss2=pionsvec.at(SSMCPion2Index);
    
  
    pionsvec.clear();
    pionsvec.push_back(os);
    pionsvec.push_back(ss1);
    pionsvec.push_back(ss2);
}

int main(int argc,char **argv){

  Log::SummaryAtExit();

  // Initialization of pythia
  Pythia pythia;
  Event& event = pythia.event;

  TString path= (TString)std::getenv("PWD") +"/output/";
  TString FileName  = path+"TauolaHelicity_" + TString(argv[1]) + ".root";
  TFile *file = new TFile(FileName,"RECREATE");


  //TFile *file = new TFile("HelicityVals_up.root","RECREATE");

  TH1F *TauMinus_TauPlus_Mass= new TH1F("TauMinus_TauPlus_Mass","M_{#tau^{-}#tau^{+}}",100,0,100);
  TH1F *RhoMinus_RhoPlus_Mass= new TH1F("RhoMinus_RhoPlus_Mass","M_{#rho^{-}#rho^{+}}",100,0,100);
  TH1F *Pions_From_a1_Mass=new TH1F("Pions_From_a1_Mass","M_{#pi#pi#pi}",100,0,3);
  
  
  TH1F *Rho_Pion_From_a1_Mass = new TH1F("Rho_Pion_From_a1_Mass","M_{#rho^{-}#pi^{0}}",100,0,3);
  TH1F *Products_From_a1_Mass = new TH1F("Products_From_a1_Mass","M_{#pi^{-}#pi^{0}#pi^{0}}",100,0,100);
  TH1F *Pions_From_Rho_Mass=new TH1F("Pions_From_Rho_Mass","M_{#pi^{-}#pi^{0}}",100,0,2);
  
  
 
  
  TH2F *ratio_ma1_plus = new TH2F("ratio_ma1_plus", "",100,0,2,100,0,3);
  ratio_ma1_plus->GetXaxis()->SetTitle("m_{a1} truth (GeV/c^{2})"); 
  ratio_ma1_plus->GetYaxis()->SetTitle("m_{a1} plus (GeV/c^{2})");
  ratio_ma1_plus->SetOption("COLZ");
  
  TH2F *ratio_ma1_minus = new TH2F("ratio_ma1_minus", "",100,0,2,100,0,3);
  ratio_ma1_minus->GetXaxis()->SetTitle("m_{a1} truth (GeV/c^{2})"); 
  ratio_ma1_minus->GetYaxis()->SetTitle("m_{a1} minus (GeV/c^{2})");
  ratio_ma1_minus->SetOption("COLZ");
 
  
  //---------------------------- error vs confidence interval ----------------------------------
  TH2F *absolute_error_momentum_tau_plus_vs_x = new TH2F("absolute_error_momentum_tau_plus_vs_x", "|(p^{+}_{#tau}(rec) - p_{#tau}(gen) ) / p_{#tau}(gen)|",100,0,0.5,100,-2,2);
  absolute_error_momentum_tau_plus_vs_x->GetXaxis()->SetTitle("x (GeV/c)"); 
  absolute_error_momentum_tau_plus_vs_x->GetYaxis()->SetTitle("|(p^{+}_{#tau}(rec) - p_{#tau}(gen) ) / p_{#tau}(gen)|");
  absolute_error_momentum_tau_plus_vs_x->SetOption("COLZ");
  
  TH2F *absolute_error_momentum_tau_minus_vs_x = new TH2F("absolute_error_momentum_tau_minus_vs_x", "|(p^{-}_{#tau}(rec) - p_{#tau}(gen) ) / p_{#tau}(gen)|",100,0,0.5,100,-2,2);
  absolute_error_momentum_tau_minus_vs_x->GetXaxis()->SetTitle("x (GeV/c)"); 
  absolute_error_momentum_tau_minus_vs_x->GetYaxis()->SetTitle("|(p^{-}_{#tau}(rec) - p_{#tau}(gen) ) / p_{#tau}(gen)|");
  absolute_error_momentum_tau_minus_vs_x->SetOption("COLZ");
  
  TH2F *absolute_error_momentum_tau_mean_vs_x = new TH2F("absolute_error_momentum_tau_mean_vs_x", "|(p^{mean}_{#tau}(rec) - p_{#tau}(gen) ) / p_{#tau}(gen)|",100,0,0.5,100,-20,50);
  absolute_error_momentum_tau_mean_vs_x->GetXaxis()->SetTitle("x (GeV/c)"); 
  absolute_error_momentum_tau_mean_vs_x->GetYaxis()->SetTitle("|(p^{mean}_{#tau}(rec) - p_{#tau}(gen) ) / p_{#tau}(gen)|");
  absolute_error_momentum_tau_mean_vs_x->SetOption("COLZ");
  
  TH2F *error_vs_x = new TH2F("error_vs_x", "#sqrt{(p^{+}_{#tau}(rec) - p_{#tau}(gen))^{2} + (p^{-}_{#tau}(rec) - p_{#tau}(gen))^{2}}",100,0,0.5,100,0,20);
  error_vs_x->GetXaxis()->SetTitle("x (GeV/c)"); 
  error_vs_x->GetYaxis()->SetTitle("#sqrt{(p^{+}_{#tau}(rec) - p_{#tau}(gen))^{2} + (p^{-}_{#tau}(rec) - p_{#tau}(gen))^{2}}");
  error_vs_x->SetOption("COLZ");
  //---------------------------------------------------------------------------------------------
 
 
  //---------------------- theta CM ---------------------------------
  TH1F *theta_cm = new TH1F("theta_cm","#theta_{GJ}^{CM}",100,0,4);
  theta_cm->GetXaxis()->SetTitle("#theta_{GJ}^{CM} (rad)"); 
  theta_cm->GetYaxis()->SetTitle("#");
  
  TH1F *costheta_cm = new TH1F("costheta_cm","cos(#theta_{GJ}^{CM})",100,-1.1,1.1);
  costheta_cm->GetXaxis()->SetTitle("#cos(theta_{GJ}^{CM}) (rad)"); 
  costheta_cm->GetYaxis()->SetTitle("#");
  
  TH2F *theta_CM_vs_a1_momentum = new TH2F("theta_CM_vs_a1_momentum","#theta^{*} vs p_{a_{1}}",100,0,4,100,0,50);
  theta_CM_vs_a1_momentum->GetYaxis()->SetTitle("p_{a_{1}} (GeV/c)"); 
  theta_CM_vs_a1_momentum->GetXaxis()->SetTitle("#theta^{*} (rad)");
  theta_CM_vs_a1_momentum->SetOption("COLZ");
  
  TH2F *theta_CM_vs_tau_momentum_gen = new TH2F("theta_CM_vs_tau_momentum_gen","#theta^{*} vs p_{#tau}(gen)",100,0,4,100,0,50);
  theta_CM_vs_tau_momentum_gen->GetYaxis()->SetTitle("p_{#tau} (GeV/c)"); 
  theta_CM_vs_tau_momentum_gen->GetXaxis()->SetTitle("#theta^{*} (rad)");
  theta_CM_vs_tau_momentum_gen->SetOption("COLZ");
  //
  TH2F *theta_CM_vs_tau_momentum_rec_pa1_20 = new TH2F("theta_CM_vs_tau_momentum_rec_pa1_20","#theta^{*} vs p_{#tau}(rec)",100,0,4,100,0,50);
  theta_CM_vs_tau_momentum_rec_pa1_20->GetYaxis()->SetTitle("p_{#tau} (GeV/c)"); 
  theta_CM_vs_tau_momentum_rec_pa1_20->GetXaxis()->SetTitle("#theta^{*} (rad)");
  theta_CM_vs_tau_momentum_rec_pa1_20->SetOption("COLZ");
  
  TH2F *theta_CM_vs_tau_momentum_rec_pa1_23 = new TH2F("theta_CM_vs_tau_momentum_rec_pa1_23","#theta^{*} vs p_{#tau}(rec)",100,0,4,100,0,50);
  theta_CM_vs_tau_momentum_rec_pa1_23->GetYaxis()->SetTitle("p_{#tau} (GeV/c)"); 
  theta_CM_vs_tau_momentum_rec_pa1_23->GetXaxis()->SetTitle("#theta^{*} (rad)");
  theta_CM_vs_tau_momentum_rec_pa1_23->SetOption("COLZ");
  
  TH2F *theta_CM_vs_tau_momentum_rec_pa1_25 = new TH2F("theta_CM_vs_tau_momentum_rec_pa1_25","#theta^{*} vs p_{#tau}(rec)",100,0,4,100,0,50);
  theta_CM_vs_tau_momentum_rec_pa1_25->GetYaxis()->SetTitle("p_{#tau} (GeV/c)"); 
  theta_CM_vs_tau_momentum_rec_pa1_25->GetXaxis()->SetTitle("#theta^{*} (rad)");
  theta_CM_vs_tau_momentum_rec_pa1_25->SetOption("COLZ");
  
  TH2F *theta_CM_vs_tau_momentum_rec_pa1_27 = new TH2F("theta_CM_vs_tau_momentum_rec_pa1_27","#theta^{*} vs p_{#tau}(rec)",100,0,4,100,0,50);
  theta_CM_vs_tau_momentum_rec_pa1_27->GetYaxis()->SetTitle("p_{#tau} (GeV/c)"); 
  theta_CM_vs_tau_momentum_rec_pa1_27->GetXaxis()->SetTitle("#theta^{*} (rad)");
  theta_CM_vs_tau_momentum_rec_pa1_27->SetOption("COLZ");
  
  TH2F *theta_CM_vs_tau_momentum_rec_pa1_30 = new TH2F("theta_CM_vs_tau_momentum_rec_pa1_30","#theta^{*} vs p_{#tau}(rec)",100,0,4,100,0,50);
  theta_CM_vs_tau_momentum_rec_pa1_30->GetYaxis()->SetTitle("p_{#tau} (GeV/c)"); 
  theta_CM_vs_tau_momentum_rec_pa1_30->GetXaxis()->SetTitle("#theta^{*} (rad)");
  theta_CM_vs_tau_momentum_rec_pa1_30->SetOption("COLZ");
  
  TH2F *theta_CM_vs_tau_momentum_rec_pa1_35 = new TH2F("theta_CM_vs_tau_momentum_rec_pa1_35","#theta^{*} vs p_{#tau}(rec)",100,0,4,100,0,50);
  theta_CM_vs_tau_momentum_rec_pa1_35->GetYaxis()->SetTitle("p_{#tau} (GeV/c)"); 
  theta_CM_vs_tau_momentum_rec_pa1_35->GetXaxis()->SetTitle("#theta^{*} (rad)");
  theta_CM_vs_tau_momentum_rec_pa1_35->SetOption("COLZ");
  
  TH2F *theta_CM_vs_tau_momentum_rec_pa1_40 = new TH2F("theta_CM_vs_tau_momentum_rec_pa1_40","#theta^{*} vs p_{#tau}(rec)",100,0,4,100,0,50);
  theta_CM_vs_tau_momentum_rec_pa1_40->GetYaxis()->SetTitle("p_{#tau} (GeV/c)"); 
  theta_CM_vs_tau_momentum_rec_pa1_40->GetXaxis()->SetTitle("#theta^{*} (rad)");
  theta_CM_vs_tau_momentum_rec_pa1_40->SetOption("COLZ");
  //
  TH2F *costheta_CM_vs_a1_momentum = new TH2F("costheta_CM_vs_a1_momentum","cos(#theta^{*}) vs p_{a_{1}}",100,-1.1,1.1,100,0,50);
  costheta_CM_vs_a1_momentum->GetYaxis()->SetTitle("p_{a_{1}} (GeV/c)"); 
  costheta_CM_vs_a1_momentum->GetXaxis()->SetTitle("cos(#theta^{*}) (rad)");
  costheta_CM_vs_a1_momentum->SetOption("COLZ");
  
  TH2F *costheta_CM_vs_tau_momentum_gen = new TH2F("costheta_CM_vs_tau_momentum_gen","cos(#theta^{*}) vs p_{#tau}(gen)",100,-1.1,1.1,100,0,50);
  costheta_CM_vs_tau_momentum_gen->GetYaxis()->SetTitle("p_{#tau} (GeV/c)"); 
  costheta_CM_vs_tau_momentum_gen->GetXaxis()->SetTitle("cos(#theta^{*}) (rad)");
  costheta_CM_vs_tau_momentum_gen->SetOption("COLZ");
  //
  TH2F *costheta_CM_vs_tau_momentum_rec_pa1_20 = new TH2F("costheta_CM_vs_tau_momentum_rec_pa1_20","cos(#theta^{*}) vs p_{#tau}(rec)",100,-1.1,1.1,100,0,50);
  costheta_CM_vs_tau_momentum_rec_pa1_20->GetYaxis()->SetTitle("p_{#tau} (GeV/c)"); 
  costheta_CM_vs_tau_momentum_rec_pa1_20->GetXaxis()->SetTitle("cos(#theta^{*}) (rad)");
  costheta_CM_vs_tau_momentum_rec_pa1_20->SetOption("COLZ");
  
  TH2F *costheta_CM_vs_tau_momentum_rec_pa1_23 = new TH2F("costheta_CM_vs_tau_momentum_rec_pa1_23","cos(#theta^{*}) vs p_{#tau}(rec)",100,-1.1,1.1,100,0,50);
  costheta_CM_vs_tau_momentum_rec_pa1_23->GetYaxis()->SetTitle("p_{#tau} (GeV/c)"); 
  costheta_CM_vs_tau_momentum_rec_pa1_23->GetXaxis()->SetTitle("cos(#theta^{*}) (rad)");
  costheta_CM_vs_tau_momentum_rec_pa1_23->SetOption("COLZ");
  
  TH2F *costheta_CM_vs_tau_momentum_rec_pa1_25 = new TH2F("costheta_CM_vs_tau_momentum_rec_pa1_25","cos(#theta^{*}) vs p_{#tau}(rec)",100,-1.1,1.1,100,0,50);
  costheta_CM_vs_tau_momentum_rec_pa1_25->GetYaxis()->SetTitle("p_{#tau} (GeV/c)"); 
  costheta_CM_vs_tau_momentum_rec_pa1_25->GetXaxis()->SetTitle("cos(#theta^{*}) (rad)");
  costheta_CM_vs_tau_momentum_rec_pa1_25->SetOption("COLZ");
  
  TH2F *costheta_CM_vs_tau_momentum_rec_pa1_27 = new TH2F("costheta_CM_vs_tau_momentum_rec_pa1_27","cos(#theta^{*}) vs p_{#tau}(rec)",100,-1.1,1.1,100,0,50);
  costheta_CM_vs_tau_momentum_rec_pa1_27->GetYaxis()->SetTitle("p_{#tau} (GeV/c)"); 
  costheta_CM_vs_tau_momentum_rec_pa1_27->GetXaxis()->SetTitle("cos(#theta^{*}) (rad)");
  costheta_CM_vs_tau_momentum_rec_pa1_27->SetOption("COLZ");
  
  TH2F *costheta_CM_vs_tau_momentum_rec_pa1_30 = new TH2F("costheta_CM_vs_tau_momentum_rec_pa1_30","cos(#theta^{*}) vs p_{#tau}(rec)",100,-1.1,1.1,100,0,50);
  costheta_CM_vs_tau_momentum_rec_pa1_30->GetYaxis()->SetTitle("p_{#tau} (GeV/c)"); 
  costheta_CM_vs_tau_momentum_rec_pa1_30->GetXaxis()->SetTitle("cos(#theta^{*}) (rad)");
  costheta_CM_vs_tau_momentum_rec_pa1_30->SetOption("COLZ");
  
  TH2F *costheta_CM_vs_tau_momentum_rec_pa1_35 = new TH2F("costheta_CM_vs_tau_momentum_rec_pa1_35","cos(#theta^{*}) vs p_{#tau}(rec)",100,-1.1,1.1,100,0,50);
  costheta_CM_vs_tau_momentum_rec_pa1_35->GetYaxis()->SetTitle("p_{#tau} (GeV/c)"); 
  costheta_CM_vs_tau_momentum_rec_pa1_35->GetXaxis()->SetTitle("cos(#theta^{*}) (rad)");
  costheta_CM_vs_tau_momentum_rec_pa1_35->SetOption("COLZ");
  
  TH2F *costheta_CM_vs_tau_momentum_rec_pa1_40 = new TH2F("costheta_CM_vs_tau_momentum_rec_pa1_40","cos(#theta^{*}) vs p_{#tau}(rec)",100,-1.1,1.1,100,0,50);
  costheta_CM_vs_tau_momentum_rec_pa1_40->GetYaxis()->SetTitle("p_{#tau} (GeV/c)"); 
  costheta_CM_vs_tau_momentum_rec_pa1_40->GetXaxis()->SetTitle("cos(#theta^{*}) (rad)");
  costheta_CM_vs_tau_momentum_rec_pa1_40->SetOption("COLZ");
  //-----------------------------------------------------------------
 
  
  //------------------- negative rootsquare ---------------------------
  TH1F *negative_sqrt_a1_mass = new TH1F("negative_sqrt_a1_mass", "m_{a_{1}} => #sqrt{-}",100,0,2);
  negative_sqrt_a1_mass->GetXaxis()->SetTitle("m_{a_{1}} (GeV/c^{2})"); 
  negative_sqrt_a1_mass->GetYaxis()->SetTitle("#");
  
  TH1F *negative_sqrt_a1_momentum = new TH1F("negative_sqrt_a1_momentum", "p_{a_{1}} => #sqrt{-}",100,0,50);
  negative_sqrt_a1_momentum->GetXaxis()->SetTitle("p_{a_{1}} (GeV/c)"); 
  negative_sqrt_a1_momentum->GetYaxis()->SetTitle("#");
  
  TH1F *negative_sqrt_fixed_mass_a1_momentum = new TH1F("negative_sqrt_fixed_mass_a1_momentum", "p_{a_{1}} => #sqrt{-} (fixed mass)",100,0,50);
  negative_sqrt_fixed_mass_a1_momentum->GetXaxis()->SetTitle("p_{a_{1}} (GeV/c)"); 
  negative_sqrt_fixed_mass_a1_momentum->GetYaxis()->SetTitle("#");
  
  
  //------------ errors mean correlations ------------------------
  TH2F *error_a1_mass_vs_a1_momentum = new TH2F("error_a1_mass_vs_a1_momentum", "|p^{mean}_{#tau}(rec) - p_{#tau}(gen)|>20",100,0,2,100,0,50);
  error_a1_mass_vs_a1_momentum->GetXaxis()->SetTitle("m_{a_{1}} (GeV/c^{2})"); 
  error_a1_mass_vs_a1_momentum->GetYaxis()->SetTitle("p_{a_{1}} (GeV/c)");
  error_a1_mass_vs_a1_momentum->SetOption("COLZ");
  
  TH2F *error_a1_mass_vs_theta = new TH2F("error_a1_mass_vs_theta", "|p^{mean}_{#tau}(rec) - p_{#tau}(gen)|>20",100,0,2,100,0,0.02);
  error_a1_mass_vs_theta->GetXaxis()->SetTitle("m_{a_{1}} (GeV/c^{2})"); 
  error_a1_mass_vs_theta->GetYaxis()->SetTitle("#theta (rad)");
  error_a1_mass_vs_theta->SetOption("COLZ");
  
  TH2F *error_theta_vs_a1_momentum = new TH2F("error_theta_vs_a1_momentum", "|p^{mean}_{#tau}(rec) - p_{#tau}(gen)|>20",100,0,0.02,100,0,50);
  error_theta_vs_a1_momentum->GetXaxis()->SetTitle("#theta (rad)"); 
  error_theta_vs_a1_momentum->GetYaxis()->SetTitle("p_{a_{1}} (GeV/c)");
  error_theta_vs_a1_momentum->SetOption("COLZ");
  
  
  
  //----------------------------------- tau informations --------------------------------------------
  TH1F *tau1_gen_Mass = new TH1F("tau1_gen_Mass", "m_{#tau}(gen)",100,0,2);
  tau1_gen_Mass->GetXaxis()->SetTitle("m_{#tau} (GeV/c^{2})"); tau1_gen_Mass->GetYaxis()->SetTitle("#");
  
  TH1F *tau1_gen_Momentum = new TH1F("tau1_gen_Momentum", "p_{#tau}(gen)",100,0,50);
  tau1_gen_Momentum->GetXaxis()->SetTitle("p_{#tau} (GeV/c)"); tau1_gen_Momentum->GetYaxis()->SetTitle("#");
  
  TH1F *tau1_gen_pt = new TH1F("tau1_gen_pt", "p_{T}_{#tau}(gen)",100,0,50);
  tau1_gen_pt->GetXaxis()->SetTitle("p_{T}_{#tau} (GeV/c)"); tau1_gen_pt->GetYaxis()->SetTitle("#");
  //------------------------------------------------------------------------------------------------
  
  
  //----------------------------------- a1 informations -----------------------------------------------------
  TH1F *a1_rec_Mass = new TH1F("a1_rec_Mass", "m_{a_{1}}(rec)",100,0,2);
  a1_rec_Mass->GetXaxis()->SetTitle("m_{a_{1}} (GeV/c^{2})"); a1_rec_Mass->GetYaxis()->SetTitle("#");
  
  TH1F *a1_rec_Momentum = new TH1F("a1_rec_Momentum", "p_{a_{1}}(rec)",100,0,50);
  a1_rec_Momentum->GetXaxis()->SetTitle("p_{a_{1}} (GeV/c)"); a1_rec_Momentum->GetYaxis()->SetTitle("#");
  
  TH1F *a1_rec_pt = new TH1F("a1_rec_pt", "p_{T}_{a1}(rec)",100,0,50);
  a1_rec_pt->GetXaxis()->SetTitle("p_{T}_{a1} (GeV/c)"); a1_rec_pt->GetYaxis()->SetTitle("#");
  //
  TH1F *a1_gen_Mass = new TH1F("a1_gen_Mass", "m_{a_{1}}(gen)",100,0,2);
  a1_gen_Mass->GetXaxis()->SetTitle("m_{a_{1}} (GeV/c^{2})"); a1_gen_Mass->GetYaxis()->SetTitle("#");
  
  TH1F *a1_gen_Momentum = new TH1F("a1_gen_Momentum", "p_{a_{1}}(gen)",100,0,50);
  a1_gen_Momentum->GetXaxis()->SetTitle("p_{a_{1}} (GeV/c)"); a1_gen_Momentum->GetYaxis()->SetTitle("#");
  
  TH1F *a1_gen_pt = new TH1F("a1_gen_pt", "p_{T}_{a1}(gen)",100,0,50);
  a1_gen_pt->GetXaxis()->SetTitle("p_{T}_{a1} (GeV/c)"); a1_gen_pt->GetYaxis()->SetTitle("#");
  
  
  
  //----------------------------------------------  absolute error -------------------------------------------------------
  TH1F *absolute_error_mass_a1 = new TH1F("absolute_error_mass_a1","m_{a1}(rec) - m_{a1}(gen)",100,-1,1);
  absolute_error_mass_a1->GetYaxis()->SetTitle("#");
  
  TH1F *absolute_error_momentum_a1 = new TH1F("absolute_error_momentum_a1","p_{a1}(rec) - p_{a1}(gen)",100,-1,1);
  absolute_error_momentum_a1->GetYaxis()->SetTitle("#");
  //----------------------------------------------------------------------------------------------------------------------
  
  //----------------------------------------------  relative error -------------------------------------------------------
  TH1F *relative_error_mass_a1 = new TH1F("relative_error_mass_a1","(m_{a1}(rec) - m_{a1}(gen) ) / m_{a1}(gen)",100,-1,1);
  relative_error_mass_a1->GetYaxis()->SetTitle("#"); 
  
  TH1F *relative_error_momentum_a1 = new TH1F("relative_error_momentum_a1","(p_{a1}(rec) - p_{a1}(gen) ) / p_{a1}(gen)",100,-1,1);
  relative_error_momentum_a1->GetYaxis()->SetTitle("#"); 
  //----------------------------------------------------------------------------------------------------------------------
  
  //----------------------------------------------  relative error abs---------------------------------------------------
  TH1F *relative_error_abs_mass_a1 = new TH1F("relative_error_abs_mass_a1","|(m_{a1}(rec) - m_{a1}(gen) ) / m_{a1}(gen)|",100,0,1);
  relative_error_abs_mass_a1->GetYaxis()->SetTitle("#");
  
  TH1F *relative_error_abs_momentum_a1 = new TH1F("relative_error_abs_momentum_a1","|(p_{a1}(rec) - p_{a1}(gen) ) / p_{a1}(gen)|",100,0,1);
  relative_error_abs_momentum_a1->GetYaxis()->SetTitle("#"); 
  //----------------------------------------------------------------------------------------------------------------------
  
  
  
  
  //---------------------------------------------------------
  
  //------------ theta vs a1 --------------------------------
  TH2F *theta_GJ_vs_a1_momentum = new TH2F("theta_GJ_vs_a1_momentum","#theta_{GJ} vs p_{a_{1}}",100,0,50,100,0,0.1);
  theta_GJ_vs_a1_momentum->GetXaxis()->SetTitle("p_{a_{1}} (GeV/c)"); 
  theta_GJ_vs_a1_momentum->GetYaxis()->SetTitle("#theta_{GJ} (rad)");
  theta_GJ_vs_a1_momentum->SetOption("COLZ");
  
  TH2F *theta_max_vs_a1_momentum = new TH2F("theta_max_vs_a1_momentum","#theta_{max} vs p_{a_{1}}",100,0,50,100,0,0.1);
  theta_max_vs_a1_momentum->GetXaxis()->SetTitle("p_{a_{1}} (GeV/c)"); 
  theta_max_vs_a1_momentum->GetYaxis()->SetTitle("#theta_{max} (rad)");
  theta_max_vs_a1_momentum->SetOption("COLZ");
  
  TH2F *costheta_GJ_vs_a1_momentum = new TH2F("costheta_GJ_vs_a1_momentum","cos(#theta_{GJ}) vs p_{a_{1}}",100,0,50,100,0.8,1);
  
  TH2F *costheta_max_vs_a1_momentum = new TH2F("costheta_max_vs_a1_momentum","cos(#theta_{max}) vs p_{a_{1}}",100,0,50,100,0.8,1);
  //----------------------------------------------------------
  //------------- theta / theta_max -------------------------
  TH1F *ratio_theta = new TH1F("ratio_theta", "#theta_{GJ} / #theta_{max}",100,0,10);
  ratio_theta->GetYaxis()->SetTitle("#");

  TH2F *theta_max_vs_theta_GJ = new TH2F("theta_max_vs_theta_GJ","#theta_{max} vs #theta_{GJ}",100,0,0.1,100,0,0.1);
  theta_max_vs_theta_GJ->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  theta_max_vs_theta_GJ->GetYaxis()->SetTitle("#theta_{max} (rad)");
  theta_max_vs_theta_GJ->SetOption("COLZ");
  //----------------------------------------------------------

  //-------------- physical meaning of theta --------------------------
  TH1F *Unphysics_theta = new TH1F("Unphysics_theta", "Unphysical #theta" ,100,0,0.2); //rad
  Unphysics_theta->GetXaxis()->SetTitle("#theta (rad)"); Unphysics_theta->GetYaxis()->SetTitle("#");
  
  TH1F *Unphysics_a1_momentum = new TH1F("Unphysics_a1_momentum", "Unphysical p_{a_{1}} (#theta_{GJ} > #theta_{max})" ,100,0,50);
  Unphysics_a1_momentum->GetXaxis()->SetTitle("p_{a_{1}} (GeV/c)"); Unphysics_a1_momentum->GetYaxis()->SetTitle("#");
  
  TH2F *Physical_theta_GJ_vs_a1_momentum = new TH2F("Physical_theta_GJ_vs_a1_momentum","#theta_{GJ} vs p_{a_{1}}",100,0,50,100,0,0.1);
  Physical_theta_GJ_vs_a1_momentum->GetXaxis()->SetTitle("p_{a_{1}} (GeV/c)"); 
  Physical_theta_GJ_vs_a1_momentum->GetYaxis()->SetTitle("#theta_{GJ} (rad)");
  Physical_theta_GJ_vs_a1_momentum->SetOption("COLZ");
  
  TH2F *Physical_theta_max_vs_a1_momentum = new TH2F("Physical_theta_max_vs_a1_momentum","#theta_{max} vs p_{a_{1}}",100,0,50,100,0,0.1);
  Physical_theta_max_vs_a1_momentum->GetXaxis()->SetTitle("p_{a_{1}} (GeV/c)"); 
  Physical_theta_max_vs_a1_momentum->GetYaxis()->SetTitle("#theta_{max} (rad)");
  Physical_theta_max_vs_a1_momentum->SetOption("COLZ");
  
  TH2F *Physical_costheta_GJ_vs_a1_momentum = new TH2F("Physical_costheta_GJ_vs_a1_momentum","cos(#theta_{GJ}) vs p_{a_{1}}",100,0,50,100,0.8,1);
  TH2F *Physical_costheta_max_vs_a1_momentum = new TH2F("Physical_costheta_max_vs_a1_momentum","cos(#theta_{max}) vs p_{a_{1}}",100,0,50,100,0.8,1);
  //---------------------------------------------------------------------------

  //################################# relative errors +/-/average solutions #####################################
  
  //###################################################  PLUS ##################################################################
  //############################################################################################################################
  //----------------------------------------------  absolute error -------------------------------------------------------
  TH1F *absolute_error_momentum_tau_plus = new TH1F("absolute_error_momentum_tau_plus","p^{+}_{#tau}(rec) - p_{#tau}(gen)",100,-45,45);
  absolute_error_momentum_tau_plus->GetYaxis()->SetTitle("#"); 
  //----------------------------------------------------------------------------------------------------------------------
  
  //----------------------------------------------  relative error -------------------------------------------------------
  TH1F *relative_error_momentum_tau_plus = new TH1F("relative_error_momentum_tau_plus","(p^{+}_{#tau}(rec) - p_{#tau}(gen) ) / p_{#tau}(gen)",100,-1,1);
  relative_error_momentum_tau_plus->GetYaxis()->SetTitle("#"); 
  //----------------------------------------------------------------------------------------------------------------------
  
  //----------------------------------------------  relative error abs---------------------------------------------------
  TH1F *relative_error_abs_momentum_tau_plus = new TH1F("relative_error_abs_momentum_tau_plus","|(p^{+}_{#tau}(rec) - p_{#tau}(gen) ) / p_{#tau}(gen)|",100,0,1);
  relative_error_abs_momentum_tau_plus->GetYaxis()->SetTitle("#"); 
  //----------------------------------------------------------------------------------------------------------------------
  
  
  //##################################################  MINUS ##################################################################
  //############################################################################################################################
  //----------------------------------------------  absolute error -------------------------------------------------------
  TH1F *absolute_error_momentum_tau_minus = new TH1F("absolute_error_momentum_tau_minus","p^{-}_{#tau}(rec) - p_{#tau}(gen)",100,-45,45);
  absolute_error_momentum_tau_minus->GetYaxis()->SetTitle("#"); 
  //----------------------------------------------------------------------------------------------------------------------
  
  //----------------------------------------------  relative error -------------------------------------------------------
  TH1F *relative_error_momentum_tau_minus = new TH1F("relative_error_momentum_tau_minus","(p^{-}_{#tau}(rec) - p_{#tau}(gen) ) / p_{#tau}(gen)",100,-1,1);
  relative_error_momentum_tau_minus->GetYaxis()->SetTitle("#"); 
  //----------------------------------------------------------------------------------------------------------------------
  
  //----------------------------------------------  relative error abs---------------------------------------------------
  TH1F *relative_error_abs_momentum_tau_minus = new TH1F("relative_error_abs_momentum_tau_minus","|(p^{-}_{#tau}(rec) - p_{#tau}(gen) ) / p_{#tau}(gen)|",100,0,1);
  relative_error_abs_momentum_tau_minus->GetYaxis()->SetTitle("#"); 
  //----------------------------------------------------------------------------------------------------------------------
  
  //##################################################  MEAN ###################################################################
  //############################################################################################################################
  //----------------------------------------------  absolute error -------------------------------------------------------
  TH1F *absolute_error_momentum_tau_mean = new TH1F("absolute_error_momentum_tau_mean","p^{mean}_{#tau}(rec) - p_{#tau}(gen)",100,-45,45);
  absolute_error_momentum_tau_mean->GetYaxis()->SetTitle("#"); 
  //----------------------------------------------------------------------------------------------------------------------
  
  //----------------------------------------------  relative error -------------------------------------------------------
  TH1F *relative_error_momentum_tau_mean = new TH1F("relative_error_momentum_tau_mean","(p^{mean}_{#tau}(rec) - p_{#tau}(gen) ) / p_{#tau}(gen)",100,-1,1);
  relative_error_momentum_tau_mean->GetYaxis()->SetTitle("#"); 
  //----------------------------------------------------------------------------------------------------------------------
  
  //----------------------------------------------  relative error abs---------------------------------------------------
  TH1F *relative_error_abs_momentum_tau_mean = new TH1F("relative_error_abs_momentum_tau_mean","|(p^{mean}_{#tau}(rec) - p_{#tau}(gen) ) / p_{#tau}(gen)|",100,0,1);
  relative_error_abs_momentum_tau_mean->GetYaxis()->SetTitle("#"); 
  //----------------------------------------------------------------------------------------------------------------------
  
  
  TH2F *error_momentum_tau_mean_vs_momentum_a1 = new TH2F("error_momentum_tau_mean_vs_momentum_a1","|(p_^{mean}{#tau}(rec) - p_{#tau}(gen) ) / p_{#tau}(gen)| vs p_{a1}",100,0,50,100,0,1);
  error_momentum_tau_mean_vs_momentum_a1->GetXaxis()->SetTitle("p_{a1} (GeV/c)");
  error_momentum_tau_mean_vs_momentum_a1->GetYaxis()->SetTitle("|(p^{mean}_{#tau}(rec) - p_{#tau}(gen) ) / p_{#tau}(gen)|");
  //-----------------------------------------------------------------------------

  //--------------- ratio +/-/average solutions ---------------------
  TH1F *ratio_plus = new TH1F("ratio_plus", "p^{+}_{#tau} / p_{#tau}(truth)",100,0,5);
  ratio_plus->GetYaxis()->SetTitle("#");
 
  TH1F *ratio_minus = new TH1F("ratio_minus", "p^{-}_{#tau} / p_{#tau}(truth)",100,0,5);
  ratio_minus->GetYaxis()->SetTitle("#");
  
  TH1F *ratio_mean = new TH1F("ratio_mean", "p^{mean}_{#tau} / p_{#tau}(truth)",100,0,5);
  ratio_mean->GetYaxis()->SetTitle("#");
  //----------------------------------------------------------------

  //---------------------- width ------------------------------------
  TH2F *width_vs_a1_mass = new TH2F("width_vs_a1_mass","(p^{+}_{#tau} - p^{-}_{#tau}) vs m_{a1}",100,0,2,100,0,40);
  width_vs_a1_mass->GetXaxis()->SetTitle("m_{a1} (GeV/c^{2})"); 
  width_vs_a1_mass->GetYaxis()->SetTitle("p^{+}_{#tau} - p^{-}_{#tau} (GeV/c)");
  width_vs_a1_mass->SetOption("COLZ");
  
  TH2F *width_vs_a1_momentum = new TH2F("width_vs_a1_momentum","(p^{+}_{#tau} - p^{-}_{#tau}) vs p_{a1}",100,0,50,100,0,40);
  width_vs_a1_momentum->GetXaxis()->SetTitle("p_{a1} (GeV/c)"); 
  width_vs_a1_momentum->GetYaxis()->SetTitle("p^{+}_{#tau} - p^{-}_{#tau} (GeV/c)");
  width_vs_a1_momentum->SetOption("COLZ");
  
  TH2F *width_vs_a1_fixed_mass = new TH2F("width_vs_a1_fixed_mass","(p^{+}_{#tau} - p^{-}_{#tau}) vs m_{a1} (p_{a1} fixed)",100,0,2,100,0,40);
  width_vs_a1_fixed_mass->GetXaxis()->SetTitle("m_{a1} (GeV/c^{2})"); 
  width_vs_a1_fixed_mass->GetYaxis()->SetTitle("p^{+}_{#tau} - p^{-}_{#tau} (GeV/c)");
  width_vs_a1_fixed_mass->SetOption("COLZ");
  
  TH2F *width_vs_a1_fixed_momentum = new TH2F("width_vs_a1_fixed_momentum","(p^{+}_{#tau} - p^{-}_{#tau}) vs p_{a1} (m_{a1} fixed)",100,0,50,100,0,40);
  width_vs_a1_fixed_momentum->GetXaxis()->SetTitle("p_{a1} (GeV/c)"); 
  width_vs_a1_fixed_momentum->GetYaxis()->SetTitle("p^{+}_{#tau} - p^{-}_{#tau} (GeV/c)");
  width_vs_a1_fixed_momentum->SetOption("COLZ");
  
  
  
  TH2F *width_vs_a1_momentum_theta_1 = new TH2F("width_vs_a1_momentum_theta_1","(p^{+}_{#tau} - p^{-}_{#tau}) vs p_{a1} (0.010 rad < #theta_{GJ} < 0.012 rad)",100,0,50,100,0,40);
  width_vs_a1_momentum_theta_1->GetXaxis()->SetTitle("p_{a1} (GeV/c)"); 
  width_vs_a1_momentum_theta_1->GetYaxis()->SetTitle("p^{+}_{#tau} - p^{-}_{#tau} (GeV/c)");
  width_vs_a1_momentum_theta_1->SetOption("COLZ");
  
  TH2F *width_vs_a1_momentum_theta_2 = new TH2F("width_vs_a1_momentum_theta_2","(p^{+}_{#tau} - p^{-}_{#tau}) vs p_{a1} (0.012 rad < #theta_{GJ} < 0.014 rad)",100,0,50,100,0,40);
  width_vs_a1_momentum_theta_2->GetXaxis()->SetTitle("p_{a1} (GeV/c)"); 
  width_vs_a1_momentum_theta_2->GetYaxis()->SetTitle("p^{+}_{#tau} - p^{-}_{#tau} (GeV/c)");
  width_vs_a1_momentum_theta_2->SetOption("COLZ");
  
  TH2F *width_vs_a1_momentum_theta_3 = new TH2F("width_vs_a1_momentum_theta_3","(p^{+}_{#tau} - p^{-}_{#tau}) vs p_{a1} (0.014 rad < #theta_{GJ} < 0.016 rad)",100,0,50,100,0,40);
  width_vs_a1_momentum_theta_3->GetXaxis()->SetTitle("p_{a1} (GeV/c)"); 
  width_vs_a1_momentum_theta_3->GetYaxis()->SetTitle("p^{+}_{#tau} - p^{-}_{#tau} (GeV/c)");
  width_vs_a1_momentum_theta_3->SetOption("COLZ");
  
  TH2F *width_vs_a1_momentum_theta_4 = new TH2F("width_vs_a1_momentum_theta_4","(p^{+}_{#tau} - p^{-}_{#tau}) vs p_{a1} (0.016 rad < #theta_{GJ} < 0.018 rad)",100,0,50,100,0,40);
  width_vs_a1_momentum_theta_4->GetXaxis()->SetTitle("p_{a1} (GeV/c)"); 
  width_vs_a1_momentum_theta_4->GetYaxis()->SetTitle("p^{+}_{#tau} - p^{-}_{#tau} (GeV/c)");
  width_vs_a1_momentum_theta_4->SetOption("COLZ");
  
  TH2F *width_vs_a1_momentum_theta_5 = new TH2F("width_vs_a1_momentum_theta_5","(p^{+}_{#tau} - p^{-}_{#tau}) vs p_{a1} (0.018 rad < #theta_{GJ} < 0.020 rad)",100,0,50,100,0,40);
  width_vs_a1_momentum_theta_5->GetXaxis()->SetTitle("p_{a1} (GeV/c)"); 
  width_vs_a1_momentum_theta_5->GetYaxis()->SetTitle("p^{+}_{#tau} - p^{-}_{#tau} (GeV/c)");
  width_vs_a1_momentum_theta_5->SetOption("COLZ");

  //####################################################################
  //#######################  tau momentum ##############################
  //####################################################################
  
  //----------------------------------------------  absolute error -------------------------------------------------------
  TH1F *absolute_error_momentum_tau = new TH1F("absolute_error_momentum_tau","p_{#tau}(rec) - p_{#tau}(gen)",100,-45,45);
  absolute_error_momentum_tau->GetYaxis()->SetTitle("#"); 
  //----------------------------------------------------------------------------------------------------------------------
  
  //----------------------------------------------  relative error -------------------------------------------------------
  TH1F *relative_error_momentum_tau = new TH1F("relative_error_momentum_tau","(p_{#tau}(rec) - p_{#tau}(gen) ) / p_{#tau}(gen)",100,-1,1);
  relative_error_momentum_tau->GetYaxis()->SetTitle("#"); 
  //----------------------------------------------------------------------------------------------------------------------
  
  //----------------------------------------------  relative error abs---------------------------------------------------
  TH1F *relative_error_abs_momentum_tau = new TH1F("relative_error_abs_momentum_tau","|(p_{#tau}(rec) - p_{#tau}(gen) ) / p_{#tau}(gen)|",100,0,1);
  relative_error_abs_momentum_tau->GetYaxis()->SetTitle("#"); 
  //----------------------------------------------------------------------------------------------------------------------
  
  
  //---------------------- momentum of tau vs momentum of a1 -----------------------
  TH2F *Momentum_tau_vs_Momentum_a1 = new TH2F("Momentum_tau_vs_Momentum_a1","p_{#tau} vs p_{a1}",100,0,50,100,0,50);
  Momentum_tau_vs_Momentum_a1->GetXaxis()->SetTitle("p_{a1} (GeV/c)"); 
  Momentum_tau_vs_Momentum_a1->GetYaxis()->SetTitle("p_{#tau} (GeV/c)");
  Momentum_tau_vs_Momentum_a1->SetOption("COLZ");
  
  TH2F *Momentum_tau_vs_Momentum_a1_theta_1 = new TH2F("Momentum_tau_vs_Momentum_a1_theta_1","p_{#tau} vs p_{a1} (0.010 rad < #theta_{GJ} < 0.012 rad)",100,0,50,100,0,50);
  Momentum_tau_vs_Momentum_a1_theta_1->GetXaxis()->SetTitle("p_{a1} (GeV/c)"); 
  Momentum_tau_vs_Momentum_a1_theta_1->GetYaxis()->SetTitle("p_{#tau} (GeV/c)");
  Momentum_tau_vs_Momentum_a1_theta_1->SetOption("COLZ");
  
  TH2F *Momentum_tau_vs_Momentum_a1_theta_2 = new TH2F("Momentum_tau_vs_Momentum_a1_theta_2","p_{#tau} vs p_{a1} (0.012 rad < #theta_{GJ} < 0.014 rad)",100,0,50,100,0,50);
  Momentum_tau_vs_Momentum_a1_theta_2->GetXaxis()->SetTitle("p_{a1} (GeV/c)"); 
  Momentum_tau_vs_Momentum_a1_theta_2->GetYaxis()->SetTitle("p_{#tau} (GeV/c)");
  Momentum_tau_vs_Momentum_a1_theta_2->SetOption("COLZ");
  
  TH2F *Momentum_tau_vs_Momentum_a1_theta_3 = new TH2F("Momentum_tau_vs_Momentum_a1_theta_3","p_{#tau} vs p_{a1} (0.014 rad < #theta_{GJ} < 0.016 rad)",100,0,50,100,0,50);
  Momentum_tau_vs_Momentum_a1_theta_3->GetXaxis()->SetTitle("p_{a1} (GeV/c)"); 
  Momentum_tau_vs_Momentum_a1_theta_3->GetYaxis()->SetTitle("p_{#tau} (GeV/c)");
  Momentum_tau_vs_Momentum_a1_theta_3->SetOption("COLZ");
  
  TH2F *Momentum_tau_vs_Momentum_a1_theta_4= new TH2F("Momentum_tau_vs_Momentum_a1_theta_4","p_{#tau} vs p_{a1} (0.016 rad < #theta_{GJ} < 0.018 rad)",100,0,50,100,0,50);
  Momentum_tau_vs_Momentum_a1_theta_4->GetXaxis()->SetTitle("p_{a1} (GeV/c)"); 
  Momentum_tau_vs_Momentum_a1_theta_4->GetYaxis()->SetTitle("p_{#tau} (GeV/c)");
  Momentum_tau_vs_Momentum_a1_theta_4->SetOption("COLZ");
  
  TH2F *Momentum_tau_vs_Momentum_a1_theta_5= new TH2F("Momentum_tau_vs_Momentum_a1_theta_5","p_{#tau} vs p_{a1} (0.018 rad < #theta_{GJ} < 0.020 rad)",100,0,50,100,0,50);
  Momentum_tau_vs_Momentum_a1_theta_5->GetXaxis()->SetTitle("p_{a1} (GeV/c)"); 
  Momentum_tau_vs_Momentum_a1_theta_5->GetYaxis()->SetTitle("p_{#tau} (GeV/c)");
  Momentum_tau_vs_Momentum_a1_theta_5->SetOption("COLZ");
  //---------------------------------------------------------------------

  //---------------------- momentum of tau vs mass of a1 -----------------------
  
  TH2F *Momentum_tau_vs_Mass_a1 = new TH2F("Momentum_tau_vs_Mass_a1","p_{#tau} vs m_{a1}",100,0,2,100,0,50);
  Momentum_tau_vs_Mass_a1->GetXaxis()->SetTitle("m_{a1} (GeV/c^{2})"); 
  Momentum_tau_vs_Mass_a1->GetYaxis()->SetTitle("p_{#tau} (GeV/c)");
  Momentum_tau_vs_Mass_a1->SetOption("COLZ");
  
  TH2F *Momentum_tau_vs_Mass_a1_theta_1 = new TH2F("Momentum_tau_vs_Mass_a1_theta_1","p_{#tau} vs m_{a1} (0.010 rad < #theta_{GJ} < 0.012 rad)",100,0,2,100,0,50);
  Momentum_tau_vs_Mass_a1_theta_1->GetXaxis()->SetTitle("m_{a1} (GeV/c^{2})"); 
  Momentum_tau_vs_Mass_a1_theta_1->GetYaxis()->SetTitle("p_{#tau} (GeV/c)");
  Momentum_tau_vs_Mass_a1_theta_1->SetOption("COLZ");
  
  TH2F *Momentum_tau_vs_Mass_a1_theta_2 = new TH2F("Momentum_tau_vs_Mass_a1_theta_2","p_{#tau} vs m_{a1} (0.012 rad < #theta_{GJ} < 0.014 rad)",100,0,2,100,0,50);
  Momentum_tau_vs_Mass_a1_theta_2->GetXaxis()->SetTitle("m_{a1} (GeV/c^{2})"); 
  Momentum_tau_vs_Mass_a1_theta_2->GetYaxis()->SetTitle("p_{#tau} (GeV/c)");
  Momentum_tau_vs_Mass_a1_theta_2->SetOption("COLZ");
  
  TH2F *Momentum_tau_vs_Mass_a1_theta_3 = new TH2F("Momentum_tau_vs_Mass_a1_theta_3","p_{#tau} vs m_{a1} (0.014 rad < #theta_{GJ} < 0.016 rad)",100,0,2,100,0,50);
  Momentum_tau_vs_Mass_a1_theta_3->GetXaxis()->SetTitle("m_{a1} (GeV/c^{2})"); 
  Momentum_tau_vs_Mass_a1_theta_3->GetYaxis()->SetTitle("p_{#tau} (GeV/c)");
  Momentum_tau_vs_Mass_a1_theta_3->SetOption("COLZ");
  
  TH2F *Momentum_tau_vs_Mass_a1_theta_4= new TH2F("Momentum_tau_vs_Mass_a1_theta_4","p_{#tau} vs m_{a1} (0.016 rad < #theta_{GJ} < 0.018 rad)",100,0,2,100,0,50);
  Momentum_tau_vs_Mass_a1_theta_4->GetXaxis()->SetTitle("m_{a1} (GeV/c^{2})"); 
  Momentum_tau_vs_Mass_a1_theta_4->GetYaxis()->SetTitle("p_{#tau} (GeV/c)");
  Momentum_tau_vs_Mass_a1_theta_4->SetOption("COLZ");
  
  TH2F *Momentum_tau_vs_Mass_a1_theta_5= new TH2F("Momentum_tau_vs_Mass_a1_theta_5","p_{#tau} vs m_{a1} (0.018 rad < #theta_{GJ} < 0.020 rad)",100,0,2,100,0,50);
  Momentum_tau_vs_Mass_a1_theta_5->GetXaxis()->SetTitle("m_{a1} (GeV/c^{2})"); 
  Momentum_tau_vs_Mass_a1_theta_5->GetYaxis()->SetTitle("p_{#tau} (GeV/c)");
  Momentum_tau_vs_Mass_a1_theta_5->SetOption("COLZ");
  //---------------------------------------------------------------------


  //-------------------- momentum of tau (gen) vs theta -----------------------
  TH2F *Momentum_tau_vs_theta = new TH2F("Momentum_tau_vs_theta","p_{#tau} vs #theta",100,0,0.04,100,0,50);
  Momentum_tau_vs_theta->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_vs_theta->GetYaxis()->SetTitle("p_{#tau} (GeV/c)");
  Momentum_tau_vs_theta->SetOption("COLZ");
  
  TH2F *Momentum_tau_vs_theta_pa1_20 = new TH2F("Momentum_tau_vs_theta_pa1_20","p_{#tau} vs #theta (p_{a1} = 20 GeV)",100,0,0.04,100,0,50);
  Momentum_tau_vs_theta_pa1_20->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_vs_theta_pa1_20->GetYaxis()->SetTitle("p_{#tau} (GeV/c)");
  Momentum_tau_vs_theta_pa1_20->SetOption("COLZ");
  
  TH2F *Momentum_tau_vs_theta_pa1_23 = new TH2F("Momentum_tau_vs_theta_pa1_23","p_{#tau} vs #theta (p_{a1} = 23 GeV)",100,0,0.04,100,0,50);
  Momentum_tau_vs_theta_pa1_23->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_vs_theta_pa1_23->GetYaxis()->SetTitle("p_{#tau} (GeV/c)");
  Momentum_tau_vs_theta_pa1_23->SetOption("COLZ");

  TH2F *Momentum_tau_vs_theta_pa1_25 = new TH2F("Momentum_tau_vs_theta_pa1_25","p_{#tau} vs #theta (p_{a1} = 25 GeV)",100,0,0.04,100,0,50);
  Momentum_tau_vs_theta_pa1_25->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_vs_theta_pa1_25->GetYaxis()->SetTitle("p_{#tau} (GeV/c)");
  Momentum_tau_vs_theta_pa1_25->SetOption("COLZ");
  
  TH2F *Momentum_tau_vs_theta_pa1_27 = new TH2F("Momentum_tau_vs_theta_pa1_27","p_{#tau} vs #theta (p_{a1} = 27 GeV)",100,0,0.04,100,0,50);
  Momentum_tau_vs_theta_pa1_27->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_vs_theta_pa1_27->GetYaxis()->SetTitle("p_{#tau} (GeV/c)");
  Momentum_tau_vs_theta_pa1_27->SetOption("COLZ");
  
  TH2F *Momentum_tau_vs_theta_pa1_30 = new TH2F("Momentum_tau_vs_theta_pa1_30","p_{#tau} vs #theta (p_{a1} = 30 GeV)",100,0,0.04,100,0,50);
  Momentum_tau_vs_theta_pa1_30->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_vs_theta_pa1_30->GetYaxis()->SetTitle("p_{#tau} (GeV/c)");
  Momentum_tau_vs_theta_pa1_30->SetOption("COLZ");
  
  TH2F *Momentum_tau_vs_theta_pa1_35 = new TH2F("Momentum_tau_vs_theta_pa1_35","p_{#tau} vs #theta (p_{a1} = 35 GeV)",100,0,0.04,100,0,50);
  Momentum_tau_vs_theta_pa1_35->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_vs_theta_pa1_35->GetYaxis()->SetTitle("p_{#tau} (GeV/c)");
  Momentum_tau_vs_theta_pa1_35->SetOption("COLZ");
  
  TH2F *Momentum_tau_vs_theta_pa1_40 = new TH2F("Momentum_tau_vs_theta_pa1_40","p_{#tau} vs #theta (p_{a1} = 40 GeV)",100,0,0.04,100,0,50);
  Momentum_tau_vs_theta_pa1_40->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_vs_theta_pa1_40->GetYaxis()->SetTitle("p_{#tau} (GeV/c)");
  Momentum_tau_vs_theta_pa1_40->SetOption("COLZ");
  //---------------------------------------------------------------------
  
  //-------------------- momentum of tau (rec) vs theta -----------------------
  /*
  TH2F* Momentum_tau_vs_theta = new TH2F("Momentum_tau_vs_theta","p_{#tau} vs #theta",100,0,0.2,100,0,50);
  Momentum_tau_vs_theta->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_vs_theta->GetYaxis()->SetTitle("p_{#tau} (GeV/c)");
  Momentum_tau_vs_theta->SetOption("COLZ");
  */
  //
  TH2F *Momentum_tau_rec_vs_theta_pa1_20 = new TH2F("Momentum_tau_rec_vs_theta_pa1_20","p_{#tau}(rec) vs #theta (p_{a1} = 20 GeV)",100,0,0.040,100,0,50);
  Momentum_tau_rec_vs_theta_pa1_20->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_vs_theta_pa1_20->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_vs_theta_pa1_20->SetOption("COLZ");
  
  TH2F *Momentum_tau_rec_vs_theta_pa1_23 = new TH2F("Momentum_tau_rec_vs_theta_pa1_23","p_{#tau}(rec) vs #theta (p_{a1} = 23 GeV)",100,0,0.040,100,0,50);
  Momentum_tau_rec_vs_theta_pa1_23->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_vs_theta_pa1_23->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_vs_theta_pa1_23->SetOption("COLZ");
  
  TH2F *Momentum_tau_rec_vs_theta_pa1_25 = new TH2F("Momentum_tau_rec_vs_theta_pa1_25","p_{#tau}(rec) vs #theta (p_{a1} = 25 GeV)",100,0,0.035,100,0,50);
  Momentum_tau_rec_vs_theta_pa1_25->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_vs_theta_pa1_25->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_vs_theta_pa1_25->SetOption("COLZ");
  
  TH2F *Momentum_tau_rec_vs_theta_pa1_27 = new TH2F("Momentum_tau_rec_vs_theta_pa1_27","p_{#tau}(rec) vs #theta (p_{a1} = 27 GeV)",100,0,0.030,100,0,50);
  Momentum_tau_rec_vs_theta_pa1_27->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_vs_theta_pa1_27->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_vs_theta_pa1_27->SetOption("COLZ");
  
  //
  TH2F *Momentum_tau_rec_vs_theta_pa1_30 = new TH2F("Momentum_tau_rec_vs_theta_pa1_30","p_{#tau}(rec) vs #theta (p_{a1} = 30 GeV)",100,0,0.030,100,0,60);
  Momentum_tau_rec_vs_theta_pa1_30->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_vs_theta_pa1_30->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_vs_theta_pa1_30->SetOption("COLZ");
  
  TH2F *Momentum_tau_plus_rec_vs_theta_pa1_30 = new TH2F("Momentum_tau_plus_rec_vs_theta_pa1_30","p^{+}_{#tau}(rec) vs #theta (p_{a1} = 30 GeV)",100,0,0.020,100,0,60);
  Momentum_tau_plus_rec_vs_theta_pa1_30->GetXaxis()->SetTitle("#theta_{GJ} (rad)");
  Momentum_tau_plus_rec_vs_theta_pa1_30->GetYaxis()->SetTitle("p^{+}_{#tau}(rec) (GeV/c)");
  Momentum_tau_plus_rec_vs_theta_pa1_30->SetOption("COLZ");
  
  TH2F *Momentum_tau_minus_rec_vs_theta_pa1_30 = new TH2F("Momentum_tau_minus_rec_vs_theta_pa1_30","p^{-}_{#tau}(rec) vs #theta (p_{a1} = 30 GeV)",100,0,0.020,100,0,60);
  Momentum_tau_minus_rec_vs_theta_pa1_30->GetXaxis()->SetTitle("#theta_{GJ} (rad)");
  Momentum_tau_minus_rec_vs_theta_pa1_30->GetYaxis()->SetTitle("p^{-}_{#tau}(rec) (GeV/c)");
  Momentum_tau_minus_rec_vs_theta_pa1_30->SetOption("COLZ");
  
  TH2F *Mass_a1_above_theta_max_plus_pa1_30 = new TH2F("Mass_a1_above_theta_max_plus_pa1_30","m_{a1}(rec) vs #theta for p^{+}_{#tau} (p_{a1} = 30 GeV)",100,0,0.020,100,0,2);
  Mass_a1_above_theta_max_plus_pa1_30->GetXaxis()->SetTitle("#theta_{GJ} (rad)");
  Mass_a1_above_theta_max_plus_pa1_30->GetYaxis()->SetTitle("m_{a1}(rec) (GeV/c^{2})");
  Mass_a1_above_theta_max_plus_pa1_30->SetOption("COLZ");
  
  TH2F *Mass_a1_above_theta_max_minus_pa1_30 = new TH2F("Mass_a1_above_theta_max_minus_pa1_30","m_{a1}(rec) vs #theta for p^{-}_{#tau} (p_{a1} = 30 GeV)",100,0,0.020,100,0,2);
  Mass_a1_above_theta_max_minus_pa1_30->GetXaxis()->SetTitle("#theta_{GJ} (rad)");
  Mass_a1_above_theta_max_minus_pa1_30->GetYaxis()->SetTitle("m_{a1}(rec) (GeV/c^{2})");
  Mass_a1_above_theta_max_minus_pa1_30->SetOption("COLZ");
  
  //
  TH2F *Momentum_tau_rec_vs_theta_pa1_35 = new TH2F("Momentum_tau_rec_vs_theta_pa1_35","p_{#tau}(rec) vs #theta (p_{a1} = 35 GeV)",100,0,0.025,100,0,60);
  Momentum_tau_rec_vs_theta_pa1_35->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_vs_theta_pa1_35->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_vs_theta_pa1_35->SetOption("COLZ");
  
  TH2F *Momentum_tau_plus_rec_vs_theta_pa1_35 = new TH2F("Momentum_tau_plus_rec_vs_theta_pa1_35","p^{+}_{#tau}(rec) vs #theta (p_{a1} = 35 GeV)",100,0,0.020,100,0,60);
  Momentum_tau_plus_rec_vs_theta_pa1_35->GetXaxis()->SetTitle("#theta_{GJ} (rad)");
  Momentum_tau_plus_rec_vs_theta_pa1_35->GetYaxis()->SetTitle("p^{+}_{#tau}(rec) (GeV/c)");
  Momentum_tau_plus_rec_vs_theta_pa1_35->SetOption("COLZ");
  
  TH2F *Momentum_tau_minus_rec_vs_theta_pa1_35 = new TH2F("Momentum_tau_minus_rec_vs_theta_pa1_35","p^{-}_{#tau}(rec) vs #theta (p_{a1} = 35 GeV)",100,0,0.020,100,0,60);
  Momentum_tau_minus_rec_vs_theta_pa1_35->GetXaxis()->SetTitle("#theta_{GJ} (rad)");
  Momentum_tau_minus_rec_vs_theta_pa1_35->GetYaxis()->SetTitle("p^{-}_{#tau}(rec) (GeV/c)");
  Momentum_tau_minus_rec_vs_theta_pa1_35->SetOption("COLZ");
  
  TH2F *Mass_a1_above_theta_max_plus_pa1_35 = new TH2F("Mass_a1_above_theta_max_plus_pa1_35","m_{a1}(rec) vs #theta for p^{+}_{#tau} (p_{a1} = 35 GeV)",100,0,0.020,100,0,2);
  Mass_a1_above_theta_max_plus_pa1_35->GetXaxis()->SetTitle("#theta_{GJ} (rad)");
  Mass_a1_above_theta_max_plus_pa1_35->GetYaxis()->SetTitle("m_{a1}(rec) (GeV/c^{2})");
  Mass_a1_above_theta_max_plus_pa1_35->SetOption("COLZ");
  
  TH2F *Mass_a1_above_theta_max_minus_pa1_35 = new TH2F("Mass_a1_above_theta_max_minus_pa1_35","m_{a1}(rec) vs #theta for p^{-}_{#tau} (p_{a1} = 35 GeV)",100,0,0.020,100,0,2);
  Mass_a1_above_theta_max_minus_pa1_35->GetXaxis()->SetTitle("#theta_{GJ} (rad)");
  Mass_a1_above_theta_max_minus_pa1_35->GetYaxis()->SetTitle("m_{a1}(rec) (GeV/c^{2})");
  Mass_a1_above_theta_max_minus_pa1_35->SetOption("COLZ");
  
  
  //
  TH2F *Momentum_tau_rec_vs_theta_pa1_40 = new TH2F("Momentum_tau_rec_vs_theta_pa1_40","p_{#tau}(rec) vs #theta (p_{a1} = 40 GeV)",100,0,0.020,100,0,60);
  Momentum_tau_rec_vs_theta_pa1_40->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_vs_theta_pa1_40->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_vs_theta_pa1_40->SetOption("COLZ");
  
  TH2F *Momentum_tau_plus_rec_vs_theta_pa1_40 = new TH2F("Momentum_tau_plus_rec_vs_theta_pa1_40","p^{+}_{#tau}(rec) vs #theta (p_{a1} = 40 GeV)",100,0,0.020,100,0,60);
  Momentum_tau_plus_rec_vs_theta_pa1_40->GetXaxis()->SetTitle("#theta_{GJ} (rad)");
  Momentum_tau_plus_rec_vs_theta_pa1_40->GetYaxis()->SetTitle("p^{+}_{#tau}(rec) (GeV/c)");
  Momentum_tau_plus_rec_vs_theta_pa1_40->SetOption("COLZ");
  
  TH2F *Momentum_tau_minus_rec_vs_theta_pa1_40 = new TH2F("Momentum_tau_minus_rec_vs_theta_pa1_40","p^{-}_{#tau}(rec) vs #theta (p_{a1} = 40 GeV)",100,0,0.020,100,0,60);
  Momentum_tau_minus_rec_vs_theta_pa1_40->GetXaxis()->SetTitle("#theta_{GJ} (rad)");
  Momentum_tau_minus_rec_vs_theta_pa1_40->GetYaxis()->SetTitle("p^{-}_{#tau}(rec) (GeV/c)");
  Momentum_tau_minus_rec_vs_theta_pa1_40->SetOption("COLZ");
  
  TH2F *Mass_a1_above_theta_max_plus_pa1_40 = new TH2F("Mass_a1_above_theta_max_plus_pa1_40","m_{a1}(rec) vs #theta for p^{+}_{#tau} (p_{a1} = 40 GeV)",100,0,0.020,100,0,2);
  Mass_a1_above_theta_max_plus_pa1_40->GetXaxis()->SetTitle("#theta_{GJ} (rad)");
  Mass_a1_above_theta_max_plus_pa1_40->GetYaxis()->SetTitle("m_{a1}(rec) (GeV/c^{2})");
  Mass_a1_above_theta_max_plus_pa1_40->SetOption("COLZ");
  
  TH2F *Mass_a1_above_theta_max_minus_pa1_40 = new TH2F("Mass_a1_above_theta_max_minus_pa1_40","m_{a1}(rec) vs #theta for p^{-}_{#tau} (p_{a1} = 40 GeV)",100,0,0.020,100,0,2);
  Mass_a1_above_theta_max_minus_pa1_40->GetXaxis()->SetTitle("#theta_{GJ} (rad)");
  Mass_a1_above_theta_max_minus_pa1_40->GetYaxis()->SetTitle("m_{a1}(rec) (GeV/c^{2})");
  Mass_a1_above_theta_max_minus_pa1_40->SetOption("COLZ");
 
  
  /*
  TH2F* Momentum_tau_vs_theta_pa1_30 = new TH2F("Momentum_tau_vs_theta_pa1_30","p_{#tau} vs #theta (p_{a1} = 30 GeV)",100,0,0.2,100,0,50);
  Momentum_tau_vs_theta_pa1_30->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_vs_theta_pa1_30->GetYaxis()->SetTitle("p_{#tau} (GeV/c)");
  Momentum_tau_vs_theta_pa1_30->SetOption("COLZ");
  
  TH2F* Momentum_tau_vs_theta_pa1_35 = new TH2F("Momentum_tau_vs_theta_pa1_35","p_{#tau} vs #theta (p_{a1} = 35 GeV)",100,0,0.2,100,0,50);
  Momentum_tau_vs_theta_pa1_35->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_vs_theta_pa1_35->GetYaxis()->SetTitle("p_{#tau} (GeV/c)");
  Momentum_tau_vs_theta_pa1_35->SetOption("COLZ");
  
  TH2F* Momentum_tau_vs_theta_pa1_40 = new TH2F("Momentum_tau_vs_theta_pa1_40","p_{#tau} vs #theta (p_{a1} = 40 GeV)",100,0,0.2,100,0,50);
  Momentum_tau_vs_theta_pa1_40->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_vs_theta_pa1_40->GetYaxis()->SetTitle("p_{#tau} (GeV/c)");
  Momentum_tau_vs_theta_pa1_40->SetOption("COLZ");
  */
  //---------------------------------------------------------------------
  
  //-------------------- momentum of tau (rec) with fixed mass vs theta ----------------------
  TH2F* Momentum_tau_rec_fixed_mass_vs_theta_pa1_20 = new TH2F("Momentum_tau_rec_fixed_mass_vs_theta_pa1_20","p_{#tau}(rec) (fixed m_{a1}) vs #theta (p_{a1} = 20 GeV)",100,0,0.030,100,0,50);
  Momentum_tau_rec_fixed_mass_vs_theta_pa1_20->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_fixed_mass_vs_theta_pa1_20->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_fixed_mass_vs_theta_pa1_20->SetOption("COLZ");
  
  TH2F* Momentum_tau_rec_fixed_mass_vs_theta_pa1_23 = new TH2F("Momentum_tau_rec_fixed_mass_vs_theta_pa1_23","p_{#tau}(rec) (fixed m_{a1}) vs #theta (p_{a1} = 23 GeV)",100,0,0.025,100,0,55);
  Momentum_tau_rec_fixed_mass_vs_theta_pa1_23->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_fixed_mass_vs_theta_pa1_23->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_fixed_mass_vs_theta_pa1_23->SetOption("COLZ");
  
  TH2F* Momentum_tau_rec_fixed_mass_vs_theta_pa1_25 = new TH2F("Momentum_tau_rec_fixed_mass_vs_theta_pa1_25","p_{#tau}(rec) (fixed m_{a1}) vs #theta (p_{a1} = 25 GeV)",100,0,0.025,100,0,60);
  Momentum_tau_rec_fixed_mass_vs_theta_pa1_25->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_fixed_mass_vs_theta_pa1_25->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_fixed_mass_vs_theta_pa1_25->SetOption("COLZ");
  
  TH2F* Momentum_tau_rec_fixed_mass_vs_theta_pa1_27 = new TH2F("Momentum_tau_rec_fixed_mass_vs_theta_pa1_27","p_{#tau}(rec) (fixed m_{a1}) vs #theta (p_{a1} = 27 GeV)",100,0,0.020,100,0,65);
  Momentum_tau_rec_fixed_mass_vs_theta_pa1_27->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_fixed_mass_vs_theta_pa1_27->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_fixed_mass_vs_theta_pa1_27->SetOption("COLZ");
  
  TH2F* Momentum_tau_rec_fixed_mass_vs_theta_pa1_30 = new TH2F("Momentum_tau_rec_fixed_mass_vs_theta_pa1_30","p_{#tau}(rec) (fixed m_{a1}) vs #theta (p_{a1} = 30 GeV)",100,0,0.020,100,0,70);
  Momentum_tau_rec_fixed_mass_vs_theta_pa1_30->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_fixed_mass_vs_theta_pa1_30->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_fixed_mass_vs_theta_pa1_30->SetOption("COLZ");
  
  TH2F* Momentum_tau_rec_fixed_mass_vs_theta_pa1_35 = new TH2F("Momentum_tau_rec_fixed_mass_vs_theta_pa1_35","p_{#tau}(rec) (fixed m_{a1}) vs #theta (p_{a1} = 35 GeV)",100,0,0.020,100,0,80);
  Momentum_tau_rec_fixed_mass_vs_theta_pa1_35->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_fixed_mass_vs_theta_pa1_35->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_fixed_mass_vs_theta_pa1_35->SetOption("COLZ");
  
  TH2F* Momentum_tau_rec_fixed_mass_vs_theta_pa1_40 = new TH2F("Momentum_tau_rec_fixed_mass_vs_theta_pa1_40","p_{#tau}(rec) (fixed m_{a1}) vs #theta (p_{a1} = 40 GeV)",100,0,0.015,100,0,90);
  Momentum_tau_rec_fixed_mass_vs_theta_pa1_40->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_fixed_mass_vs_theta_pa1_40->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_fixed_mass_vs_theta_pa1_40->SetOption("COLZ");
  //---------------------------------------------------------------------

  //-------------------- momentum of tau (rec) vs theta MASS CONDITIONS ----------------------
  
  TH2F *Momentum_tau_rec_vs_theta_pa1_20_small_ma1 = new TH2F("Momentum_tau_rec_vs_theta_pa1_20_small_ma1","p_{#tau}(rec) vs #theta (p_{a1} = 20 GeV and m_{a1} < 1.230 GeV)",100,0,0.040,100,0,50);
  Momentum_tau_rec_vs_theta_pa1_20_small_ma1->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_vs_theta_pa1_20_small_ma1->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_vs_theta_pa1_20_small_ma1->SetOption("COLZ");
  
  TH2F *Momentum_tau_rec_vs_theta_pa1_20_large_ma1 = new TH2F("Momentum_tau_rec_vs_theta_pa1_20_large_ma1","p_{#tau}(rec) vs #theta (p_{a1} = 20 GeV and m_{a1} > 1.230 GeV)",100,0,0.040,100,0,50);
  Momentum_tau_rec_vs_theta_pa1_20_large_ma1->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_vs_theta_pa1_20_large_ma1->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_vs_theta_pa1_20_large_ma1->SetOption("COLZ");
  //
  TH2F *Momentum_tau_rec_vs_theta_pa1_23_small_ma1 = new TH2F("Momentum_tau_rec_vs_theta_pa1_23_small_ma1","p_{#tau}(rec) vs #theta (p_{a1} = 23 GeV and m_{a1} < 1.230 GeV)",100,0,0.040,100,0,50);
  Momentum_tau_rec_vs_theta_pa1_23_small_ma1->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_vs_theta_pa1_23_small_ma1->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_vs_theta_pa1_23_small_ma1->SetOption("COLZ");
  
  TH2F *Momentum_tau_rec_vs_theta_pa1_23_large_ma1 = new TH2F("Momentum_tau_rec_vs_theta_pa1_23_large_ma1","p_{#tau}(rec) vs #theta (p_{a1} = 23 GeV and m_{a1} > 1.230 GeV)",100,0,0.040,100,0,50);
  Momentum_tau_rec_vs_theta_pa1_23_large_ma1->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_vs_theta_pa1_23_large_ma1->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_vs_theta_pa1_23_large_ma1->SetOption("COLZ");
  //
  TH2F *Momentum_tau_rec_vs_theta_pa1_25_small_ma1 = new TH2F("Momentum_tau_rec_vs_theta_pa1_25_small_ma1","p_{#tau}(rec) vs #theta (p_{a1} = 25 GeV and m_{a1} < 1.230 GeV)",100,0,0.035,100,0,50);
  Momentum_tau_rec_vs_theta_pa1_25_small_ma1->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_vs_theta_pa1_25_small_ma1->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_vs_theta_pa1_25_small_ma1->SetOption("COLZ");
  
  TH2F *Momentum_tau_rec_vs_theta_pa1_25_large_ma1 = new TH2F("Momentum_tau_rec_vs_theta_pa1_25_large_ma1","p_{#tau}(rec) vs #theta (p_{a1} = 25 GeV and m_{a1} > 1.230 GeV)",100,0,0.035,100,0,50);
  Momentum_tau_rec_vs_theta_pa1_25_large_ma1->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_vs_theta_pa1_25_large_ma1->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_vs_theta_pa1_25_large_ma1->SetOption("COLZ");
  //
  TH2F *Momentum_tau_rec_vs_theta_pa1_27_small_ma1 = new TH2F("Momentum_tau_rec_vs_theta_pa1_27_small_ma1","p_{#tau}(rec) vs #theta (p_{a1} = 27 GeV and m_{a1} < 1.230 GeV)",100,0,0.030,100,0,50);
  Momentum_tau_rec_vs_theta_pa1_27_small_ma1->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_vs_theta_pa1_27_small_ma1->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_vs_theta_pa1_27_small_ma1->SetOption("COLZ");
  
  TH2F *Momentum_tau_rec_vs_theta_pa1_27_large_ma1 = new TH2F("Momentum_tau_rec_vs_theta_pa1_27_large_ma1","p_{#tau}(rec) vs #theta (p_{a1} = 27 GeV and m_{a1} > 1.230 GeV)",100,0,0.030,100,0,50);
  Momentum_tau_rec_vs_theta_pa1_27_large_ma1->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_vs_theta_pa1_27_large_ma1->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_vs_theta_pa1_27_large_ma1->SetOption("COLZ");
  //
   TH2F *Momentum_tau_rec_vs_theta_pa1_30_small_ma1 = new TH2F("Momentum_tau_rec_vs_theta_pa1_30_small_ma1","p_{#tau}(rec) vs #theta (p_{a1} = 30 GeV and m_{a1} < 1.230 GeV)",100,0,0.035,100,0,50);
  Momentum_tau_rec_vs_theta_pa1_30_small_ma1->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_vs_theta_pa1_30_small_ma1->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_vs_theta_pa1_30_small_ma1->SetOption("COLZ");
  
  TH2F *Momentum_tau_rec_vs_theta_pa1_30_large_ma1 = new TH2F("Momentum_tau_rec_vs_theta_pa1_30_large_ma1","p_{#tau}(rec) vs #theta (p_{a1} = 30 GeV and m_{a1} > 1.230 GeV)",100,0,0.035,100,0,50);
  Momentum_tau_rec_vs_theta_pa1_30_large_ma1->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_vs_theta_pa1_30_large_ma1->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_vs_theta_pa1_30_large_ma1->SetOption("COLZ");
  //
  TH2F *Momentum_tau_rec_vs_theta_pa1_35_small_ma1 = new TH2F("Momentum_tau_rec_vs_theta_pa1_35_small_ma1","p_{#tau}(rec) vs #theta (p_{a1} = 35 GeV and m_{a1} < 1.230 GeV)",100,0,0.030,100,0,50);
  Momentum_tau_rec_vs_theta_pa1_35_small_ma1->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_vs_theta_pa1_35_small_ma1->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_vs_theta_pa1_35_small_ma1->SetOption("COLZ");
  
  TH2F *Momentum_tau_rec_vs_theta_pa1_35_large_ma1 = new TH2F("Momentum_tau_rec_vs_theta_pa1_35_large_ma1","p_{#tau}(rec) vs #theta (p_{a1} = 35 GeV and m_{a1} > 1.230 GeV)",100,0,0.030,100,0,50);
  Momentum_tau_rec_vs_theta_pa1_35_large_ma1->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_vs_theta_pa1_35_large_ma1->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_vs_theta_pa1_35_large_ma1->SetOption("COLZ");
  //
  TH2F *Momentum_tau_rec_vs_theta_pa1_40_small_ma1 = new TH2F("Momentum_tau_rec_vs_theta_pa1_40_small_ma1","p_{#tau}(rec) vs #theta (p_{a1} = 40 GeV and m_{a1} < 1.230 GeV)",100,0,0.030,100,0,50);
  Momentum_tau_rec_vs_theta_pa1_40_small_ma1->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_vs_theta_pa1_40_small_ma1->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_vs_theta_pa1_40_small_ma1->SetOption("COLZ");
  
  TH2F *Momentum_tau_rec_vs_theta_pa1_40_large_ma1 = new TH2F("Momentum_tau_rec_vs_theta_pa1_40_large_ma1","p_{#tau}(rec) vs #theta (p_{a1} = 40 GeV and m_{a1} > 1.230 GeV)",100,0,0.030,100,0,50);
  Momentum_tau_rec_vs_theta_pa1_40_large_ma1->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  Momentum_tau_rec_vs_theta_pa1_40_large_ma1->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)");
  Momentum_tau_rec_vs_theta_pa1_40_large_ma1->SetOption("COLZ");
  //------------------------------------------------------------------------------------
  
  /*
  //-------------------------- Crossing ---------------------------------------------
  //--- pa1 = 30 GeV ---------------
  
  TH2F *crossing_pa1_30_above_theta_mass = new TH2F("crossing_pa1_30_above_theta_mass","m_{#tau} vs m_{a1} (before crossing)",100,0,2,100,0,2);
  crossing_pa1_30_above_theta_mass->GetXaxis()->SetTitle("m_{a1} (GeV/c^{2})"); 
  crossing_pa1_30_above_theta_mass->GetYaxis()->SetTitle("m_{#tau} (GeV/c^{2})");
  crossing_pa1_30_above_theta_mass->SetOption("COLZ");
  
  TH2F *crossing_pa1_30_near_theta_mass = new TH2F("crossing_pa1_30_near_theta_mass","m_{#tau} vs m_{a1} (crossing)",100,0,2,100,0,2);
  crossing_pa1_30_near_theta_mass->GetXaxis()->SetTitle("m_{a1} (GeV/c^{2})"); 
  crossing_pa1_30_near_theta_mass->GetYaxis()->SetTitle("m_{#tau} (GeV/c^{2})");
  crossing_pa1_30_near_theta_mass->SetOption("COLZ");
  
  TH2F *crossing_pa1_30_below_theta_mass = new TH2F("crossing_pa1_30_below_theta_mass","m_{#tau} vs m_{a1} (after crossing)",100,0,2,100,0,2);
  crossing_pa1_30_below_theta_mass->GetXaxis()->SetTitle("m_{a1} (GeV/c^{2})"); 
  crossing_pa1_30_below_theta_mass->GetYaxis()->SetTitle("m_{#tau} (GeV/c^{2})");
  crossing_pa1_30_below_theta_mass->SetOption("COLZ");
  
  //////////
  
  TH2F *crossing_pa1_30_above_theta_momentum = new TH2F("crossing_pa1_30_above_theta_momentum","p_{a1} vs m_{a1} (before crossing)",100,0,2,100,0,50);
  crossing_pa1_30_above_theta_momentum->GetXaxis()->SetTitle("m_{a1} (GeV/c^{2})"); 
  crossing_pa1_30_above_theta_momentum->GetYaxis()->SetTitle("p_{a1} (GeV/c)");
  crossing_pa1_30_above_theta_momentum->SetOption("COLZ");
  
  TH2F *crossing_pa1_30_near_theta_momentum = new TH2F("crossing_pa1_30_near_theta_momentum","p_{a1} vs m_{a1} (crossing)",100,0,2,100,0,50);
  crossing_pa1_30_near_theta_momentum->GetXaxis()->SetTitle("m_{a1} (GeV/c^{2})"); 
  crossing_pa1_30_near_theta_momentum->GetYaxis()->SetTitle("p_{a1} (GeV/c})");
  crossing_pa1_30_near_theta_momentum->SetOption("COLZ");
  
  TH2F *crossing_pa1_30_below_theta_momentum = new TH2F("crossing_pa1_30_below_theta_momentum","p_{a1} vs m_{a1} (after crossing)",100,0,2,100,0,50);
  crossing_pa1_30_below_theta_momentum->GetXaxis()->SetTitle("m_{a1} (GeV/c^{2})"); 
  crossing_pa1_30_below_theta_momentum->GetYaxis()->SetTitle("p_{a1} (GeV/c)");
  crossing_pa1_30_below_theta_momentum->SetOption("COLZ");
  
  //////////////
  TH2F *crossing_pa1_30_above_theta_energy = new TH2F("crossing_pa1_30_above_theta_energy","E_{#tau} vs E_{a1} (before crossing)",100,0,100,100,0,100);
  crossing_pa1_30_above_theta_energy->GetXaxis()->SetTitle("E_{a1} (GeV)"); 
  crossing_pa1_30_above_theta_energy->GetYaxis()->SetTitle("E_{#tau} (GeV)");
  crossing_pa1_30_above_theta_energy->SetOption("COLZ");
  
  TH2F *crossing_pa1_30_near_theta_energy = new TH2F("crossing_pa1_30_near_theta_energy","E_{#tau} vs E_{a1} (crossing)",100,0,100,100,0,100);
  crossing_pa1_30_near_theta_energy->GetXaxis()->SetTitle("E_{a1} (GeV)"); 
  crossing_pa1_30_near_theta_energy->GetYaxis()->SetTitle("E_{#tau} (GeV})");
  crossing_pa1_30_near_theta_energy->SetOption("COLZ");
  
  TH2F *crossing_pa1_30_below_theta_energy = new TH2F("crossing_pa1_30_below_theta_energy","E_{#tau} vs E_{a1} (after crossing)",100,0,100,100,0,100);
  crossing_pa1_30_below_theta_energy->GetXaxis()->SetTitle("E_{a1} (GeV)"); 
  crossing_pa1_30_below_theta_energy->GetYaxis()->SetTitle("E_{#tau} (GeV)");
  crossing_pa1_30_below_theta_energy->SetOption("COLZ");
  
  //--- pa1 = 35 GeV ---------------
  
  TH2F *crossing_pa1_35_above_theta_mass = new TH2F("crossing_pa1_35_above_theta_mass","m_{#tau} vs m_{a1} (before crossing)",100,0,2,100,0,2);
  crossing_pa1_35_above_theta_mass->GetXaxis()->SetTitle("m_{a1} (GeV/c^{2})"); 
  crossing_pa1_35_above_theta_mass->GetYaxis()->SetTitle("m_{#tau} (GeV/c^{2})");
  crossing_pa1_35_above_theta_mass->SetOption("COLZ");
  
  TH2F *crossing_pa1_35_near_theta_mass = new TH2F("crossing_pa1_35_near_theta_mass","m_{#tau} vs m_{a1} (crossing)",100,0,2,100,0,2);
  crossing_pa1_35_near_theta_mass->GetXaxis()->SetTitle("m_{a1} (GeV/c^{2})"); 
  crossing_pa1_35_near_theta_mass->GetYaxis()->SetTitle("m_{#tau} (GeV/c^{2})");
  crossing_pa1_35_near_theta_mass->SetOption("COLZ");
  
  TH2F *crossing_pa1_35_below_theta_mass = new TH2F("crossing_pa1_35_below_theta_mass","m_{#tau} vs m_{a1} (after crossing)",100,0,2,100,0,2);
  crossing_pa1_35_below_theta_mass->GetXaxis()->SetTitle("m_{a1} (GeV/c^{2})"); 
  crossing_pa1_35_below_theta_mass->GetYaxis()->SetTitle("m_{#tau} (GeV/c^{2})");
  crossing_pa1_35_below_theta_mass->SetOption("COLZ");
  
  //////////
  
  TH2F *crossing_pa1_35_above_theta_momentum = new TH2F("crossing_pa1_35_above_theta_momentum","p_{a1} vs m_{a1} (before crossing)",100,0,2,100,0,50);
  crossing_pa1_35_above_theta_momentum->GetXaxis()->SetTitle("m_{a1} (GeV/c^{2})"); 
  crossing_pa1_35_above_theta_momentum->GetYaxis()->SetTitle("p_{a1} (GeV/c)");
  crossing_pa1_35_above_theta_momentum->SetOption("COLZ");
  
  TH2F *crossing_pa1_35_near_theta_momentum = new TH2F("crossing_pa1_35_near_theta_momentum","p_{a1} vs m_{a1} (crossing)",100,0,2,100,0,50);
  crossing_pa1_35_near_theta_momentum->GetXaxis()->SetTitle("m_{a1} (GeV/c^{2})"); 
  crossing_pa1_35_near_theta_momentum->GetYaxis()->SetTitle("p_{a1} (GeV/c})");
  crossing_pa1_35_near_theta_momentum->SetOption("COLZ");
  
  TH2F *crossing_pa1_35_below_theta_momentum = new TH2F("crossing_pa1_35_below_theta_momentum","p_{a1} vs m_{a1} (after crossing)",100,0,2,100,0,50);
  crossing_pa1_35_below_theta_momentum->GetXaxis()->SetTitle("m_{a1} (GeV/c^{2})"); 
  crossing_pa1_35_below_theta_momentum->GetYaxis()->SetTitle("p_{a1} (GeV/c)");
  crossing_pa1_35_below_theta_momentum->SetOption("COLZ");
  
  //////////////
  TH2F *crossing_pa1_35_above_theta_energy = new TH2F("crossing_pa1_35_above_theta_energy","E_{#tau} vs E_{a1} (before crossing)",100,0,100,100,0,100);
  crossing_pa1_35_above_theta_energy->GetXaxis()->SetTitle("E_{a1} (GeV)"); 
  crossing_pa1_35_above_theta_energy->GetYaxis()->SetTitle("E_{#tau} (GeV)");
  crossing_pa1_35_above_theta_energy->SetOption("COLZ");
  
  TH2F *crossing_pa1_35_near_theta_energy = new TH2F("crossing_pa1_35_near_theta_energy","E_{#tau} vs E_{a1} (crossing)",100,0,100,100,0,100);
  crossing_pa1_35_near_theta_energy->GetXaxis()->SetTitle("E_{a1} (GeV)"); 
  crossing_pa1_35_near_theta_energy->GetYaxis()->SetTitle("E_{#tau} (GeV})");
  crossing_pa1_35_near_theta_energy->SetOption("COLZ");
  
  TH2F *crossing_pa1_35_below_theta_energy = new TH2F("crossing_pa1_35_below_theta_energy","E_{#tau} vs E_{a1} (after crossing)",100,0,100,100,0,100);
  crossing_pa1_35_below_theta_energy->GetXaxis()->SetTitle("E_{a1} (GeV)"); 
  crossing_pa1_35_below_theta_energy->GetYaxis()->SetTitle("E_{#tau} (GeV)");
  crossing_pa1_35_below_theta_energy->SetOption("COLZ");
  
  //--- pa1 = 40 GeV ---------------
  
  TH2F *crossing_pa1_40_above_theta_mass = new TH2F("crossing_pa1_40_above_theta_mass","m_{#tau} vs m_{a1} (before crossing)",100,0,2,100,0,2);
  crossing_pa1_40_above_theta_mass->GetXaxis()->SetTitle("m_{a1} (GeV/c^{2})"); 
  crossing_pa1_40_above_theta_mass->GetYaxis()->SetTitle("m_{#tau} (GeV/c^{2})");
  crossing_pa1_40_above_theta_mass->SetOption("COLZ");
  
  TH2F *crossing_pa1_40_near_theta_mass = new TH2F("crossing_pa1_40_near_theta_mass","m_{#tau} vs m_{a1} (crossing)",100,0,2,100,0,2);
  crossing_pa1_40_near_theta_mass->GetXaxis()->SetTitle("m_{a1} (GeV/c^{2})"); 
  crossing_pa1_40_near_theta_mass->GetYaxis()->SetTitle("m_{#tau} (GeV/c^{2})");
  crossing_pa1_40_near_theta_mass->SetOption("COLZ");
  
  TH2F *crossing_pa1_40_below_theta_mass = new TH2F("crossing_pa1_40_below_theta_mass","m_{#tau} vs m_{a1} (after crossing)",100,0,2,100,0,2);
  crossing_pa1_40_below_theta_mass->GetXaxis()->SetTitle("m_{a1} (GeV/c^{2})"); 
  crossing_pa1_40_below_theta_mass->GetYaxis()->SetTitle("m_{#tau} (GeV/c^{2})");
  crossing_pa1_40_below_theta_mass->SetOption("COLZ");
  
  //////////
  
  TH2F *crossing_pa1_40_above_theta_momentum = new TH2F("crossing_pa1_40_above_theta_momentum","p_{a1} vs m_{a1} (before crossing)",100,0,2,100,0,50);
  crossing_pa1_40_above_theta_momentum->GetXaxis()->SetTitle("m_{a1} (GeV/c^{2})"); 
  crossing_pa1_40_above_theta_momentum->GetYaxis()->SetTitle("p_{a1} (GeV/c)");
  crossing_pa1_40_above_theta_momentum->SetOption("COLZ");
  
  TH2F *crossing_pa1_40_near_theta_momentum = new TH2F("crossing_pa1_40_near_theta_momentum","p_{a1} vs m_{a1} (crossing)",100,0,2,100,0,50);
  crossing_pa1_40_near_theta_momentum->GetXaxis()->SetTitle("m_{a1} (GeV/c^{2})"); 
  crossing_pa1_40_near_theta_momentum->GetYaxis()->SetTitle("p_{a1} (GeV/c})");
  crossing_pa1_40_near_theta_momentum->SetOption("COLZ");
  
  TH2F *crossing_pa1_40_below_theta_momentum = new TH2F("crossing_pa1_40_below_theta_momentum","p_{a1} vs m_{a1} (after crossing)",100,0,2,100,0,50);
  crossing_pa1_40_below_theta_momentum->GetXaxis()->SetTitle("m_{a1} (GeV/c^{2})"); 
  crossing_pa1_40_below_theta_momentum->GetYaxis()->SetTitle("p_{a1} (GeV/c)");
  crossing_pa1_40_below_theta_momentum->SetOption("COLZ");
  
  //////////////
  TH2F *crossing_pa1_40_above_theta_energy = new TH2F("crossing_pa1_40_above_theta_energy","E_{#tau} vs E_{a1} (before crossing)",100,0,100,100,0,100);
  crossing_pa1_40_above_theta_energy->GetXaxis()->SetTitle("E_{a1} (GeV)"); 
  crossing_pa1_40_above_theta_energy->GetYaxis()->SetTitle("E_{#tau} (GeV)");
  crossing_pa1_40_above_theta_energy->SetOption("COLZ");
  
  TH2F *crossing_pa1_40_near_theta_energy = new TH2F("crossing_pa1_40_near_theta_energy","E_{#tau} vs E_{a1} (crossing)",100,0,100,100,0,100);
  crossing_pa1_40_near_theta_energy->GetXaxis()->SetTitle("E_{a1} (GeV)"); 
  crossing_pa1_40_near_theta_energy->GetYaxis()->SetTitle("E_{#tau} (GeV})");
  crossing_pa1_40_near_theta_energy->SetOption("COLZ");
  
  TH2F *crossing_pa1_40_below_theta_energy = new TH2F("crossing_pa1_40_below_theta_energy","E_{#tau} vs E_{a1} (after crossing)",100,0,100,100,0,100);
  crossing_pa1_40_below_theta_energy->GetXaxis()->SetTitle("E_{a1} (GeV)"); 
  crossing_pa1_40_below_theta_energy->GetYaxis()->SetTitle("E_{#tau} (GeV)");
  crossing_pa1_40_below_theta_energy->SetOption("COLZ");
  */
  
  //--------------------------- CROSSING ------------------------------------------------------
  TH2F *Theta_vs_momentum_a1_near_crossing = new TH2F("Theta_vs_momentum_a1_near_crossing","p_{a1} vs #theta (crossing)",100,0,50,100,0,0.03);
  Theta_vs_momentum_a1_near_crossing->GetXaxis()->SetTitle("p_{a1} (GeV/c)"); 
  Theta_vs_momentum_a1_near_crossing->GetYaxis()->SetTitle("#theta (rad)");
  Theta_vs_momentum_a1_near_crossing->SetOption("COLZ");
  
  TH2F *Momentum_a1_vs_mass_a1_near_crossing = new TH2F("Momentum_a1_vs_mass_a1_near_crossing","p_{a1} vs m_{a1} (crossing)",100,0,2,100,0,50);
  Momentum_a1_vs_mass_a1_near_crossing->GetXaxis()->SetTitle("m_{a1} (GeV/c^{2})"); 
  Momentum_a1_vs_mass_a1_near_crossing->GetYaxis()->SetTitle("p_{a1} (GeV/c)");
  Momentum_a1_vs_mass_a1_near_crossing->SetOption("COLZ");
  
  //-------------------------------------------------------------------------------------------
  
  TH2F *crossing_pa1_30_theta_mass = new TH2F("crossing_pa1_30_theta_mass","m_{#tau}/m_{a1} vs #theta_{GJ}",100,0,0.03,100,0,4);
  crossing_pa1_30_theta_mass->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  crossing_pa1_30_theta_mass->GetYaxis()->SetTitle("m_{#tau}/m_{a1}");
  crossing_pa1_30_theta_mass->SetOption("COLZ");
  
  TH2F *crossing_pa1_30_theta_momentum = new TH2F("crossing_pa1_30_theta_momentum","p_{a1}/m_{a1} vs #theta_{GJ}",100,0,0.03,100,0,2);
  crossing_pa1_30_theta_momentum->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  crossing_pa1_30_theta_momentum->GetYaxis()->SetTitle("p_{a1}/m_{a1}");
  crossing_pa1_30_theta_momentum->SetOption("COLZ");
  
  TH2F *crossing_pa1_30_theta_energy = new TH2F("crossing_pa1_30_theta_energy","E_{#tau}/E_{a1} vs #theta_{GJ}",100,0,0.03,100,0,1.5);
  crossing_pa1_30_theta_energy->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  crossing_pa1_30_theta_energy->GetYaxis()->SetTitle("E_{#tau}/E_{a1}");
  crossing_pa1_30_theta_energy->SetOption("COLZ");
  //
  TH2F *crossing_pa1_35_theta_mass = new TH2F("crossing_pa1_35_theta_mass","m_{#tau}/m_{a1} vs #theta_{GJ}",100,0,0.02,100,0,4);
  crossing_pa1_35_theta_mass->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  crossing_pa1_35_theta_mass->GetYaxis()->SetTitle("m_{#tau}/m_{a1}");
  crossing_pa1_35_theta_mass->SetOption("COLZ");
  
  TH2F *crossing_pa1_35_theta_momentum = new TH2F("crossing_pa1_35_theta_momentum","p_{a1}/m_{a1} vs #theta_{GJ}",100,0,0.02,100,0,2);
  crossing_pa1_35_theta_momentum->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  crossing_pa1_35_theta_momentum->GetYaxis()->SetTitle("p_{a1}/m_{a1}");
  crossing_pa1_35_theta_momentum->SetOption("COLZ");
  
  TH2F *crossing_pa1_35_theta_energy = new TH2F("crossing_pa1_35_theta_energy","E_{#tau}/E_{a1} vs #theta_{GJ}",100,0,0.02,100,0,1.5);
  crossing_pa1_35_theta_energy->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  crossing_pa1_35_theta_energy->GetYaxis()->SetTitle("E_{#tau}/E_{a1}");
  crossing_pa1_35_theta_energy->SetOption("COLZ");
  //
  TH2F *crossing_pa1_40_theta_mass = new TH2F("crossing_pa1_40_theta_mass","m_{#tau}/m_{a1} vs #theta_{GJ}",100,0,0.015,100,0,4);
  crossing_pa1_40_theta_mass->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  crossing_pa1_40_theta_mass->GetYaxis()->SetTitle("m_{#tau}/m_{a1}");
  crossing_pa1_40_theta_mass->SetOption("COLZ");
  
  TH2F *crossing_pa1_40_theta_momentum = new TH2F("crossing_pa1_40_theta_momentum","p_{a1}/m_{a1} vs #theta_{GJ}",100,0,0.015,100,0,2);
  crossing_pa1_40_theta_momentum->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  crossing_pa1_40_theta_momentum->GetYaxis()->SetTitle("p_{a1}/m_{a1}");
  crossing_pa1_40_theta_momentum->SetOption("COLZ");
  
  TH2F *crossing_pa1_40_theta_energy = new TH2F("crossing_pa1_40_theta_energy","E_{#tau}/E_{a1} vs #theta_{GJ}",100,0,0.015,100,0,1.5);
  crossing_pa1_40_theta_energy->GetXaxis()->SetTitle("#theta_{GJ} (rad)"); 
  crossing_pa1_40_theta_energy->GetYaxis()->SetTitle("E_{#tau}/E_{a1}");
  crossing_pa1_40_theta_energy->SetOption("COLZ");
  
  
  TH2F *Momentum_tau_pa1_30_without_sqrt = new TH2F("Momentum_tau_pa1_30_without_sqrt","p_{#tau}(#sqrt{} = 0) vs #theta_{GJ}",100,0,0.03,100,0,50);
  Momentum_tau_pa1_30_without_sqrt->GetXaxis()->SetTitle("#theta_{GJ} (rad)");
  Momentum_tau_pa1_30_without_sqrt->GetYaxis()->SetTitle("p_{#tau}(#sqrt{} = 0) (GeV/c)");
  Momentum_tau_pa1_30_without_sqrt->SetOption("COLZ");
  
  TH2F *Momentum_tau_pa1_35_without_sqrt = new TH2F("Momentum_tau_pa1_35_without_sqrt","p_{#tau}(#sqrt{} = 0) vs #theta_{GJ}",100,0,0.02,100,0,70);
  Momentum_tau_pa1_35_without_sqrt->GetXaxis()->SetTitle("#theta_{GJ} (rad)");
  Momentum_tau_pa1_35_without_sqrt->GetYaxis()->SetTitle("p_{#tau}(#sqrt{} = 0) (GeV/c)");
  Momentum_tau_pa1_35_without_sqrt->SetOption("COLZ");
  
  TH2F *Momentum_tau_pa1_40_without_sqrt = new TH2F("Momentum_tau_pa1_40_without_sqrt","p_{#tau}(#sqrt{0} = 0) vs #theta_{GJ}",100,0,0.015,100,0,90);
  Momentum_tau_pa1_40_without_sqrt->GetXaxis()->SetTitle("#theta_{GJ} (rad)");
  Momentum_tau_pa1_40_without_sqrt->GetYaxis()->SetTitle("p_{#tau}(#sqrt{} = 0) (GeV/c)");
  Momentum_tau_pa1_40_without_sqrt->SetOption("COLZ");
  
  
  TH2F *Momentum_tau_pa1_30_only_sqrt = new TH2F("Momentum_tau_pa1_30_only_sqrt","p_{#tau}(#sqrt{} = 0) vs #theta_{GJ}",100,0,0.03,100,0,100);
  Momentum_tau_pa1_30_only_sqrt->GetXaxis()->SetTitle("#theta_{GJ} (rad)");
  Momentum_tau_pa1_30_only_sqrt->GetYaxis()->SetTitle("p_{#tau}(#sqrt{}) (GeV/c)");
  Momentum_tau_pa1_30_only_sqrt->SetOption("COLZ");
  
  TH2F *Momentum_tau_pa1_35_only_sqrt = new TH2F("Momentum_tau_pa1_35_only_sqrt","p_{#tau}(#sqrt{} = 0) vs #theta_{GJ}",100,0,0.02,100,0,100);
  Momentum_tau_pa1_35_only_sqrt->GetXaxis()->SetTitle("#theta_{GJ} (rad)");
  Momentum_tau_pa1_35_only_sqrt->GetYaxis()->SetTitle("p_{#tau}(#sqrt{}) (GeV/c)");
  Momentum_tau_pa1_35_only_sqrt->SetOption("COLZ");
  
  TH2F *Momentum_tau_pa1_40_only_sqrt = new TH2F("Momentum_tau_pa1_40_only_sqrt","p_{#tau}(#sqrt{0} = 0) vs #theta_{GJ}",100,0,0.015,100,0,100);
  Momentum_tau_pa1_40_only_sqrt->GetXaxis()->SetTitle("#theta_{GJ} (rad)");
  Momentum_tau_pa1_40_only_sqrt->GetYaxis()->SetTitle("p_{#tau}(#sqrt{}) (GeV/c)");
  Momentum_tau_pa1_40_only_sqrt->SetOption("COLZ");
  
  
  
  //---------------------------- 3D -----------------------------------------
  TH3F *Momentum_tau_vs_theta_vs_Momentum_a1 = new TH3F("Momentum_tau_vs_theta_vs_Momentum_a1","",100,0,0.2,100,0,50,100,0,50);
  Momentum_tau_vs_theta_vs_Momentum_a1->GetXaxis()->SetTitle("#theta_{GJ} (rad)");
  Momentum_tau_vs_theta_vs_Momentum_a1->GetYaxis()->SetTitle("p_{#tau} (GeV/c)");
  Momentum_tau_vs_theta_vs_Momentum_a1->GetZaxis()->SetTitle("p_{a1} (GeV/c)");
  Momentum_tau_vs_theta_vs_Momentum_a1->SetOption("LEGO");
  //-------------------------------------------------------------------------
  
  //-------------------------------------------------- +/- solutions ---------------------------------------------------------
  TH1F *Plus_Minus_solutions = new TH1F("Plus_Minus_solutions","+/- solutions (p^{+}_{#tau}(rec)=0 and p^{-}_{#tau}=1)",2,0,2);
  //--------------------------------------------------------------------------------------------------------------------------
  
  //------------------------------------------- Momentum tau gen vs momentum tau mean ------------------------------------
  TH2F *Momentum_tau_gen_vs_momentum_tau_mean = new TH2F("Momentum_tau_gen_vs_momentum_tau_mean","p_{#tau}(gen) vs p^{mean}_{#tau}(rec) ",150,0,150,100,40,50);
  Momentum_tau_gen_vs_momentum_tau_mean->GetXaxis()->SetTitle("p^{mean}_{#tau}(rec) (GeV/c)"); 
  Momentum_tau_gen_vs_momentum_tau_mean->GetYaxis()->SetTitle("p_{#tau}(gen) (GeV/c)");
  Momentum_tau_gen_vs_momentum_tau_mean->SetOption("COLZ");
  //-----------------------------------------------------------------------------------------------------------------------
  
  //------------------------------------------- Momentum tau gen vs resolution ------------------------------------
  TH2F *Momentum_tau_gen_vs_resolution = new TH2F("Momentum_tau_gen_vs_resolution","p_{#tau}(gen) vs #Delta ",100,-1,3,100,40,50);
  Momentum_tau_gen_vs_resolution->GetYaxis()->SetTitle("p_{#tau}(gen) (GeV/c)"); 
  Momentum_tau_gen_vs_resolution->GetXaxis()->SetTitle("#Delta");
  Momentum_tau_gen_vs_resolution->SetOption("COLZ");
  //-----------------------------------------------------------------------------------------------------------------------
  
  //------------------------------------------- Momentum tau rec vs resolution ------------------------------------
  TH2F *Momentum_tau_rec_vs_resolution = new TH2F("Momentum_tau_rec_vs_resolution","p_{#tau}(rec) vs #Delta",100,-1,3,100,0,50);
  Momentum_tau_rec_vs_resolution->GetYaxis()->SetTitle("p_{#tau}(rec) (GeV/c)"); 
  Momentum_tau_rec_vs_resolution->GetXaxis()->SetTitle("#Delta");
  Momentum_tau_rec_vs_resolution->SetOption("COLZ");
  //-----------------------------------------------------------------------------------------------------------------------
  
  //------------------------------------------- Momentum tau mean rec vs resolution ------------------------------------
  TH2F *Momentum_tau_mean_vs_resolution = new TH2F("Momentum_tau_mean_vs_resolution","p^{mean}_{#tau}(rec) vs #Delta ",100,-1,3,100,0,150);
  Momentum_tau_mean_vs_resolution->GetYaxis()->SetTitle("p^{mean}_{#tau}(rec) (GeV/c)"); 
  Momentum_tau_mean_vs_resolution->GetXaxis()->SetTitle("#Delta");
  Momentum_tau_mean_vs_resolution->SetOption("COLZ");
  //-----------------------------------------------------------------------------------------------------------------------
  
  /*
  //------------------------- omega (helicity)-------------------------------------
  TH1F *Omega_plus_a1 = new TH1F("Omega_plus_a1","#omega^{+}_{a1}",100,-1.1,1.1);
  Omega_plus_a1->GetYaxis()->SetTitle("#");
  
  TH1F *Omega_minus_a1 = new TH1F("Omega_minus_a1","#omega^{-}_{a1}",100,-1.1,1.1);
  Omega_minus_a1->GetYaxis()->SetTitle("#");
  //----------------------------------------------------------------------------------
 */
  
  TH1F *rhobeta_plus= new TH1F("rhobeta_plus","#rho^{+}",50,-1.1,1.1);
  TH1F *rhobeta_minus= new TH1F("rhobeta_minus","#rho^{-}",50,-1.1,1.1);

  TH1F *pi_plus= new TH1F("pi_plus","#pi^{+}",50,-1.1,1.1);
  TH1F *pi_minus= new TH1F("pi_minus","#pi^{-} ",50,-1.1,1.1);

  TH1F *pi_pt= new TH1F("pi_pt","#pi^{-} ",100,0,100);
  TH1F *pi_ptcut= new TH1F("pi_ptcut","#pi^{-} ",100,0,100);
 
  //BT2  (-0.11994077    , 3.76959657E-03)   0.0000000 

  TH1F *pip_plus= new TH1F("pip_plus","#pi^{+}",50,-1.1,1.1);
  TH1F *pip_minus= new TH1F("pip_minus","#pi^{-} ",50,-1.1,1.1);


  TH1F *mu_plus= new TH1F("mu_plus","#mu^{+}",50,-1.1,1.1);
  TH1F *mu_minus= new TH1F("mu_minus","#mu^{-}",50,-1.1,1.1);
  
  TH1F *OmegaMuPi_plus= new TH1F("OmegaMuPi_plus","#omega_{#pi#mu}^{+}",50,-1.1,1.1);
  TH1F *OmegaMuPi_minus= new TH1F("OmegaMuPi_minus","#omega_{#pi#mu}^{-}",50,-1.1,1.1);
 
  TH1F *Omegapipi_plus= new TH1F("Omegapipi_plus","#omega_{#pi#pi}^{+}",50,-1.1,1.1);
  TH1F *Omegapipi_minus= new TH1F("Omegapipi_minus","#omega_{#pi#pi}^{-}",50,-1.1,1.1);
  
  TH1F *pipi_mass_plus= new TH1F("pipi_mass_plus","M_{#pi#pi}^{+}",60,0,120);
  TH1F *pipi_mass_minus= new TH1F("pipi_mass_minus","M_{#pi#pi}^{-}",60,0,120);
  
  TH1F *mupi_mass_plus= new TH1F("mupi_mass_plus","M_{#mu#pi}^{+}",60,0,120);
  TH1F *mupi_mass_minus= new TH1F("mupi_mass_minus","M_{#mu#pi}^{-}",60,0,120);
 
  TH1F *omega_murho_plus= new TH1F("omega_murho_plus","#omega_{#mu#rho}^{+}",50,-1.1,1.1);
  TH1F *omega_murho_minus= new TH1F("omega_murho_minus","#omega_{#mu#rho}^{-}",50,-1.1,1.1);


 TH1F *omega_rho_plus= new TH1F("omega_rho_plus","#omega_{#rho}^{+}",50,-1.1,1.1);
 TH1F *omega_rho_minus= new TH1F("omega_rho_minus","#omega_{#rho}^{-}",50,-1.1,1.1);

 TH1F *omegabar_rho_plus= new TH1F("omegabar_rho_plus","#bar{#omega}_{#rho}^{+}",50,-1.1,1.1);
 TH1F *omegabar_rho_minus= new TH1F("omegabar_rho_minus","#bar{#omega}_{#rho}^{-}",50,-1.1,1.1);


 TH1F *omega_a1_plus= new TH1F("omega_a1_plus","#omega_{a1}^{+}",50,-1.1,1.1);
 TH1F *omega_a1_minus= new TH1F("omega_a1_minus","#omega_{a1}^{-}",50,-1.1,1.1);

 // TH1F *omega_a1_plustest= new TH1F("omega_a1_plustest","#omega_{a1}^{+}",40,-1.1,1.1);
 // TH1F *omega_a1_minustest= new TH1F("omega_a1_minustest","#omega_{a1}^{+}",40,-1.1,1.1);

 TH1F *omega_a1p_plus= new TH1F("omega_a1p_plus","#omega_{a1}^{+}",50,-1.1,1.1);
 TH1F *omega_a1p_minus= new TH1F("omega_a1p_minus","#omega_{a1}^{-}",50,-1.1,1.1);



 TH1F *omegabar_a1_plus= new TH1F("omegabar_a1_plus","#bar{#omega}_{a1}^{+}",50,-1.1,1.1);
 TH1F *omegabar_a1_minus= new TH1F("omegabar_a1_minus","#bar{#omega}_{a1}^{-}",50,-1.1,1.1);
 
 TH1F *TRFomegabar_a1_plus= new TH1F("TRFomegabar_a1_plus","TRF #omega_{a1}^{+}",50,-1.1,1.1);
 TH1F *TRFomegabar_a1_minus= new TH1F("TRFomegabar_a1_minus","TRF #omega_{a1}^{-}",50,-1.1,1.1);

  TH1F *TRFomegabar_a1scalar_plus= new TH1F("TRFomegabar_a1scalar_plus","TRF scalar #omega a1",50,-1.1,1.1);
  TH1F *TRFomegabar_a1scalar_minus= new TH1F("TRFomegabar_a1scalar_minus","TRF scalar  #omega a1",50,-1.1,1.1);
 

 TH2F *cosbetacostheta_plus= new TH2F("cosbetacostheta_plus","cos#beta  cos#theta   {a1}^{+}",50,-1.1,1.1,50,-1.1,1.1);
 TH2F *cosbetacostheta_minus= new TH2F("cosbetacostheta_minus","cos#beta  cos#theta  {a1}^{-}",50,-1.1,1.1,50,-1.1,1.1);
 cosbetacostheta_plus->GetXaxis()->SetTitle("cos#beta"); cosbetacostheta_plus->GetYaxis()->SetTitle("cos#theta");
 cosbetacostheta_minus->GetXaxis()->SetTitle("cos#beta"); cosbetacostheta_minus->GetYaxis()->SetTitle("cos#theta");
 

 TH2F *cosbetacosthetarho_plus= new TH2F("cosbetacosthetarho_plus","cos#beta  cos#theta  #rho^{+}",50,-1.1,1.1,50,-1.1,1.1);
 TH2F *cosbetacosthetarho_minus= new TH2F("cosbetacosthetarho_minus","cos#beta  cos#theta  #rho^{-}" ,50,-1.1,1.1,50,-1.1,1.1);
 cosbetacosthetarho_plus->GetXaxis()->SetTitle("cos#beta"); cosbetacosthetarho_plus->GetYaxis()->SetTitle("cos#theta");
 cosbetacosthetarho_minus->GetXaxis()->SetTitle("cos#beta"); cosbetacosthetarho_minus->GetYaxis()->SetTitle("cos#theta");

 TH2F *TRFcosbetacostheta_plus= new TH2F("TRFcosbetacostheta_plus","cos#beta  cos#theta  {a1}^{+}",50,-1.1,1.1,50,-1.1,1.1);
 TH2F *TRFcosbetacostheta_minus= new TH2F("TRFcosbetacostheta_minus","cos#beta  cos#theta  {a1}^{-}",50,-1.1,1.1,50,-1.1,1.1);
 TRFcosbetacostheta_plus->GetXaxis()->SetTitle("cos#beta"); TRFcosbetacostheta_plus->GetYaxis()->SetTitle("cos#theta");
 TRFcosbetacostheta_minus->GetXaxis()->SetTitle("cos#beta"); TRFcosbetacostheta_minus->GetYaxis()->SetTitle("cos#theta");


 TH1F *omega_a1pi_plus= new TH1F("omega_a1pi_plus","#omega_{a1#pi}^{+}",50,-1.1,1.1);
 TH1F *omega_a1pi_minus= new TH1F("omega_a1pi_minus","#omega_{a1#pi}^{-}",50,-1.1,1.1);

 TH1F *omega_a1mu_plus= new TH1F("omega_a1mu_plus","#omega_{a1#mu}^{+}",50,-1.1,1.1);
 TH1F *omega_a1mu_minus= new TH1F("omega_a1mu_minus","#omega_{a1#mu}^{-}",50,-1.1,1.1);

 TH1F *omega_pirho_plus= new TH1F("omega_pirho_plus","#omega_{#pi#rho}^{+}",50,-1.1,1.1);
 TH1F *omega_pirho_minus= new TH1F("omega_pirho_minus","#omega_{#pi#rho}^{-}",50,-1.1,1.1);
  
 TH1F *omega_rhorho_plus= new TH1F("omega_rhorho_plus","#omega_{#rho#rho}^{+}",50,-1.1,1.1);
 TH1F *omega_rhorho_minus= new TH1F("omega_rhorho_minus","#omega_{#rho#rho}^{-}",50,-1.1,1.1);
  

 TH1F *mass_pirho_plus= new TH1F("mass_pirho_plus","M_{#pi#rho}^{+}",60,0,120);
 TH1F *mass_pirho_minus= new TH1F("mass_pirho_minus","M_{#pi#rho}^{-}",60,0,120);

 TH1F *mass_murho_plus= new TH1F("mass_murho_plus","M_{#mu#rho}^{+}",60,0,120);
 TH1F *mass_murho_minus= new TH1F("mass_murho_minus","M_{#mu#rho}^{-}",60,0,120);

 TH1F *mass_rhorho_plus= new TH1F("mass_rhorho_plus","M_{#rho#rho}^{+}",60,0,120);
 TH1F *mass_rhorho_minus= new TH1F("mass_rhorho_minus","M_{#rho#rho}^{-}",60,0,120);



 TH1F *a1_mcos2gamma_plus= new TH1F("a1_mcos2gamma_plus"," <cos2#gamma>^{+}",50,-0.5,0.5);
 TH1F *a1_mcos2gamma_minus= new TH1F("a1_mcos2gamma_minus"," <cos2#gamma>^{-}",50,-0.5,0.5);



 // TH1F *omega_a1rho_plus= new TH1F("omega_a1rho_plus","#omega_{a1#rho}^{+}",50,-1.1,1.1);
 // TH1F *omega_a1rho_minus= new TH1F("omega_a1rho_minus","#omega_{a1#rho}^{-}",50,-1.1,1.1);


 // TH1F *omega_pirho_plus= new TH1F("omega_pirho_plus","#omega_{#pi#rho}^{+}",50,-1.1,1.1);
 // TH1F *omega_pirho_minus= new TH1F("omega_pirho_minus","#omega_{#pi#rho}^{-}",50,-1.1,1.1);



TH1F *s1= new TH1F("s1","s2",40,0,2);
TH1F *s2= new TH1F("s2","s2",40,0,2);
TH1F *qq= new TH1F("qq","qq",40,0,3);
TH1F *hmag= new TH1F("hmag","hmag",40,0.5,1.5);


 TH2F *s1s2= new TH2F("s1s2","s1s2",20,0,2,20,0,2);




  // Pythia8 HepMC interface depends on Pythia8 version
#ifdef PYTHIA8180_OR_LATER
  HepMC::Pythia8ToHepMC ToHepMC;
#else
  HepMC::I_Pythia8 ToHepMC;
#endif

  //pythia.readString("HadronLevel:all = off");
  pythia.readString("HadronLevel:Hadronize = off");
  pythia.readString("SpaceShower:QEDshowerByL = off");
  pythia.readString("SpaceShower:QEDshowerByQ = off");
  pythia.readString("PartonLevel:ISR = off");
  pythia.readString("PartonLevel:FSR = off");

  // Tauola is currently set to undecay taus. Otherwise, uncomment this line.
  // Uncommenting it will speed up the test significantly
  //  pythia.particleData.readString("15:mayDecay = off");

  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
  pythia.readString("23:onMode = off"); 
  pythia.readString("23:onIfAny = 15");
  //  pythia.readString("23:onIfMatch = 15 -15");

  pythia.init( 11, -11, 92.);          //electron positron collisions

  // Set up Tauola

  // Set Tauola decay mode (if needed)
  //Tauola::setSameParticleDecayMode(2);     //19 and 22 contains K0 
    //   Tauola::setOppositeParticleDecayMode(3); // 20 contains eta

  // Set Higgs scalar-pseudoscalar mixing angle
  //  Tauola::setHiggsScalarPseudoscalarMixingAngle(0.7853);
  //  Tauola::setHiggsScalarPseudoscalarPDG(25);

  Tauola::initialize();
  //  const char* str="hhu1";
  std::cout<<"  "<< TMath::Hash(argv[1])<<std::endl;
 

  //  Tauola::setSeed(time(NULL), 0, 0);
  Tauola::setSeed(TMath::Hash(argv[1]), 0, 0);
  tauola_print_parameters(); // Prints TAUOLA  parameters (residing inside its library): e.g. to test user interface

  // Our default units are GEV and MM, that will be outcome  units after TAUOLA
  // if HepMC unit variables  are correctly set. 
  // with the following coice you can fix the units for final outcome:
  //  Tauola::setUnits(Tauola::GEV,Tauola::MM); 
  //  Tauola::setUnits(Tauola::MEV,Tauola::CM); 

  // Other usefull settings:
  //  Tauola::setEtaK0sPi(0,0,0);  // switches to decay eta K0_S and pi0 1/0 on/off. 
  //  Tauola::setTauLifetime(0.0); //new tau lifetime in mm
    Tauola::spin_correlation.setAll(true);

    Log::LogDebug(true);

      Tauola::setRedefineTauMinus(redMinus);  // activates execution of routine redMinus in TAUOLA interface
      Tauola::setRedefineTauPlus(redPlus);    // activates execution of routine redPlus  in TAUOLA interface
      // Tauola::setRedefineTauMinus(redPlus);  // activates execution of routine redMinus in TAUOLA interface
      // Tauola::setRedefineTauPlus(redMinus);    // activates execution of routine redPlus  in TAUOLA interface



  MC_Initialize();

  // Begin event loop. Generate event.
  for (int iEvent = 0; iEvent < NumberOfEvents; ++iEvent){

    if(iEvent%1000==0) Log::Info()<<"Event: "<<iEvent<<endl;
    if (!pythia.next()) continue;
    //    std::cout<<"-------------- "<<std::endl;
    // Convert event record to HepMC
    HepMC::GenEvent * HepMCEvt = new HepMC::GenEvent();

    // Conversion needed if HepMC use different momentum units
    // than Pythia. However, requires HepMC 2.04 or higher.
    HepMCEvt->use_units(HepMC::Units::GEV,HepMC::Units::MM);

    ToHepMC.fill_next_event(event, HepMCEvt);

    if(iEvent<EventsToCheck)
    {
      cout<<"                                          "<<endl;
      cout<<"Momentum conservation check BEFORE/AFTER Tauola"<<endl;
      checkMomentumConservationInEvent(HepMCEvt);
    }

    // Run TAUOLA on the event
    TauolaHepMCEvent * t_event = new TauolaHepMCEvent(HepMCEvt);

    // Since we let Pythia decay taus, we have to undecay them first.
    t_event->undecayTaus();
    t_event->decayTaus();
    delete t_event; 

    if(iEvent<EventsToCheck)
    {
      checkMomentumConservationInEvent(HepMCEvt);
    }

    int JAK1(0); int SubJAK1(0);
    HepMC::GenParticle *FirstTau;
    std::vector<HepMC::GenParticle > FirstTauProducts;
    int JAK2(0); int SubJAK2(0);
    HepMC::GenParticle *SecondTau;
    std::vector<HepMC::GenParticle > SecondTauProducts;
    std::vector<HepMC::GenParticle > A1Pions;
    std::vector<HepMC::GenParticle > A1Pions1;
    std::vector<HepMC::GenParticle > A1Pions2;
    std::vector<HepMC::GenParticle > SortA1Pions;  //os, ss1, ss2

    HepMC::GenParticle *A1Minus;
    //std::vector<HepMC::GenParticle > A1Minus;
    
    TLorentzVector TauMinus(0,0,0,0);
    TLorentzVector TauPlus(0,0,0,0);
    TLorentzVector a1_gen(0,0,0,0);
     
    for ( HepMC::GenEvent::particle_const_iterator p =HepMCEvt->particles_begin();  p != HepMCEvt->particles_end(); ++p ){  
      if((*p)->pdg_id()==15){
	FirstTau = *p;
	TauMinus += TLorentzVector(FirstTau->momentum().px(), FirstTau->momentum().py(), FirstTau->momentum().pz(), FirstTau->momentum().e());
	
	for ( HepMC::GenEvent::particle_const_iterator d =HepMCEvt->particles_begin();  d != HepMCEvt->particles_end(); ++d ){  
	  if((*d)->pdg_id()!=15){
	    if((*p)->end_vertex() == (*d)->production_vertex()){
	      FirstTauProducts.push_back(**d);

	      if(abs((*d)->pdg_id()) ==  12) {JAK1 =1; 
	      }else if(abs((*d)->pdg_id()) ==  14){ JAK1=2;
	      }else if(abs((*d)->pdg_id())==  211){ JAK1 = 3;/*std::cout<<"pion of negative tau "<< (*d)->pdg_id() <<std::endl;*/}
	      

	      if( abs((*d)->pdg_id())==213 ){
		JAK1 = 4;
		for ( HepMC::GenEvent::particle_const_iterator dd =HepMCEvt->particles_begin();  dd != HepMCEvt->particles_end(); ++dd ){  
		  if( abs((*dd)->pdg_id())!=213  ){
		    if((*d)->end_vertex() == (*dd)->production_vertex()){

		      FirstTauProducts.push_back(**dd);
		    }
		  }
		}
	      }
		if( abs((*d)->pdg_id())==20213 ){
		  JAK1 = 5; int npi(0);
		  A1Minus = *d;
		  a1_gen += TLorentzVector(A1Minus->momentum().px(),A1Minus->momentum().py(), A1Minus->momentum().pz(), A1Minus->momentum().e());
		  for ( HepMC::GenEvent::particle_const_iterator dd =HepMCEvt->particles_begin();  dd != HepMCEvt->particles_end(); ++dd ){  
		    if( abs((*dd)->pdg_id())!=20213  ){
		      if((*d)->end_vertex() == (*dd)->production_vertex()){
		       
			FirstTauProducts.push_back(**dd);
			if(abs((*dd)->pdg_id())==  211)   {
			  A1Pions1.push_back(**dd); npi++;

			}
		      }
		    }
		  }
		   // std::cout<<"a1 pions  of negative tau " << A1Pions.size()<<std::endl;
		   // std::cout<<      A1Pions1.at(2).pdg_id()<< "  " <<A1Pions1.at(2).momentum().px()<<std::endl;
		   // std::cout<<      A1Pions1.at(0).pdg_id()<< "  " <<A1Pions1.at(0).momentum().px()<<std::endl;
		   // std::cout<<      A1Pions1.at(1).pdg_id()<< "  " <<A1Pions1.at(1).momentum().px()<<std::endl;
		   

		  if(npi==3){
		    SubJAK1=51;}
		  else {SubJAK1=52;}
		}
	    }
	  } 
	}
      }
    
      
      if((*p)->pdg_id()==-15){
	SecondTau = *p;
	TauPlus += TLorentzVector(SecondTau->momentum().px(), SecondTau->momentum().py(), SecondTau->momentum().pz(), SecondTau->momentum().e());
	TauMinus_TauPlus_Mass->Fill((TauMinus+TauPlus).M());

	for ( HepMC::GenEvent::particle_const_iterator d =HepMCEvt->particles_begin();  d != HepMCEvt->particles_end(); ++d )
	  {  
	    if((*d)->pdg_id()!=-15){
	     
	      if((*p)->end_vertex() == (*d)->production_vertex()){
		SecondTauProducts.push_back(**d);
		if(abs((*d)->pdg_id()) ==  12) {JAK2 =1; 
		}else if(abs((*d)->pdg_id()) ==  14){ JAK2=2;
		}else if(abs((*d)->pdg_id())==  211){ JAK2=3;/*	std::cout<<"pion of positive tau "<< (*d)->pdg_id() <<std::endl;*/}
	
		
 
		if( abs((*d)->pdg_id())==213 ){
		  JAK2 = 4;
		  for ( HepMC::GenEvent::particle_const_iterator dd =HepMCEvt->particles_begin();  dd != HepMCEvt->particles_end(); ++dd ){  
		    if( abs((*dd)->pdg_id())!=213  ){
		      if((*d)->end_vertex() == (*dd)->production_vertex()){
			SecondTauProducts.push_back(**dd);
		      }
		    }
		  }
		}

		if( abs((*d)->pdg_id())==20213 ){
		  JAK2 = 5; int npi(0);
		  for ( HepMC::GenEvent::particle_const_iterator dd =HepMCEvt->particles_begin();  dd != HepMCEvt->particles_end(); ++dd ){  
		    if( abs((*dd)->pdg_id())!=20213  ){
		      if((*d)->end_vertex() == (*dd)->production_vertex()){
			SecondTauProducts.push_back(**dd);
			if(abs((*dd)->pdg_id())==  211){
			  A1Pions2.push_back(**dd); npi++;
			}
		      }
		    }
		  }

		  if(npi==3) SubJAK2=51; else SubJAK2=52;
		}
		
	      }
	    } 
	  }
      }
    }

    bool HelPlus=false;
    bool HelMinus=false;
    if(Tauola::getHelPlus() == 1 )HelMinus=true;
    if(Tauola::getHelPlus() ==-1)HelPlus=true;

    bool HelPlus2=false;
    bool HelMinus2=false;
    if(Tauola::getHelMinus() == 1 )HelMinus2=true;
    if(Tauola::getHelMinus() ==-1)HelPlus2=true;


    int HelWeightPlus = HelPlus;
    int HelWeightMinus = HelMinus;
    // int HelWeightPlus = 0;
    // int HelWeightMinus = 0;


    // std::cout<<"  two helicities  "<< HelMinus <<"  "<<HelPlus <<std::endl;
    // std::cout<<"  s  "<< HelMinus2 <<"  "<<HelPlus2 <<std::endl;
    int tauHelicity  = Tauola::getHelPlus();

   
    TLorentzVector tau1(0,0,0,0);
    TLorentzVector mu1(0,0,0,0);
    TLorentzVector numu1(0,0,0,0);
    TLorentzVector nutau1(0,0,0,0);
    TLorentzVector nutau2(0,0,0,0);
    TLorentzVector pi2(0,0,0,0);
    TLorentzVector pi1(0,0,0,0);
    TLorentzVector tau2(0,0,0,0);

    TLorentzVector rhopi2(0,0,0,0);
    TLorentzVector rhopi02(0,0,0,0);

    TLorentzVector a1ospi(0,0,0,0);
    TLorentzVector a1ss1pi(0,0,0,0);
    TLorentzVector a1ss2pi(0,0,0,0);
    TLorentzVector a1(0,0,0,0);

  
    
    //-------------------------------------
    TLorentzVector Rho_From_A1(0,0,0,0);
    TLorentzVector Pion_0_From_A1(0,0,0,0);
    TLorentzVector Pion_Minus_From_Rho(0,0,0,0);
    TLorentzVector Pion_0_From_Rho(0,0,0,0);

    
    
   
     //------------------------------------------
    
    

    TauPolInterface Mu1;
    TauPolInterface Mu2;
    TauPolInterface Pi1;
    TauPolInterface Pi2;
    TauPolInterface Rho2;
    TauPolInterface Rho1;
    TauPolInterface A2;
    TauPolInterface A1;

    int taucharge1;
    int taucharge2;

    TauPolInterface TauPolPiPi;
    TauPolInterface TauPolMuPi;
    TauPolInterface TauPolMuRho;
    TauPolInterface TauPolPiRho;
    TauPolInterface TauPolPiA1;
    TauPolInterface TauPolMuA1;
    TauPolInterface TauPolRhoRho;



    a1Helper a1h;
    a1Helper a1hh;
    rhoHelper RhoHelp;
    PolarimetricA1 Polarimetr;
    PolarimetricA1 Polarimetr1;
    TauPolInterface TauPolPi1;
    TauPolInterface TauPolPi2;
   
    TauPolInterface TauPolRho1;

    TauPolInterface TauPolMu1;
    TauPolInterface TauPolRho2;
    TauPolInterface TauPolA1;
 
 
    

    vector<TLorentzVector> tauandprod1,tauandprod2, tauandprodMuon2;
    vector<TLorentzVector> tauandprodMuon1,tauandprodRho,tauandprodRho2,tauandprodA1;
    tau1.SetPxPyPzE(FirstTau->momentum().px(), FirstTau->momentum().py(), FirstTau->momentum().pz(), FirstTau->momentum().e());
    tau2.SetPxPyPzE(SecondTau->momentum().px(), SecondTau->momentum().py(), SecondTau->momentum().pz(), SecondTau->momentum().e());

    bool passed1(true);
    bool passed2(true);
    if(ApplyCut){
      passed1 = false;
      passed2 = false;

      if(JAK1==2){

	for(std::vector<HepMC::GenParticle>::const_iterator a = FirstTauProducts.begin(); a!=FirstTauProducts.end(); ++a){
	  if(abs(a->pdg_id())==13){
	    TLorentzVector  prod(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  );
	    if(prod.Pt() > pt_cut ) passed1=true;
	  }
	}
      }
      if(JAK1==3){
	for(std::vector<HepMC::GenParticle>::const_iterator a = FirstTauProducts.begin(); a!=FirstTauProducts.end(); ++a){
	  if(abs(a->pdg_id())==211){
	    TLorentzVector  prod(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  );
	    if(prod.Pt() > pt_cut ) passed1=true;
	  }
	}
      }

      if(JAK1==4){
	TLorentzVector prod(0,0,0,0);
	for(std::vector<HepMC::GenParticle>::const_iterator a = FirstTauProducts.begin(); a!=FirstTauProducts.end(); ++a){
	  if(abs(a->pdg_id())==211){prod += TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  );}
	  if(abs(a->pdg_id())==111){prod += TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  );}
	}
	if(prod.Pt() > pt_cut ) passed1=true;
      }
      
      if(JAK1==5){
	TLorentzVector prod(0,0,0,0);

	prod+=TLorentzVector(A1Pions1.at(0).momentum().px(), A1Pions1.at(0).momentum().py(), A1Pions1.at(0).momentum().pz(), A1Pions1.at(0).momentum().e());
	prod+=TLorentzVector(A1Pions1.at(1).momentum().px(), A1Pions1.at(1).momentum().py(), A1Pions1.at(1).momentum().pz(), A1Pions1.at(1).momentum().e());
	prod+=TLorentzVector(A1Pions1.at(2).momentum().px(), A1Pions1.at(2).momentum().py(), A1Pions1.at(2).momentum().pz(), A1Pions1.at(2).momentum().e());
	
	if(prod.Pt() > pt_cut ) passed1=true;
      }

      if(JAK2==2){
	for(std::vector<HepMC::GenParticle>::const_iterator a = SecondTauProducts.begin(); a!=SecondTauProducts.end(); ++a){
	  if(abs(a->pdg_id())==13){
	    TLorentzVector  prod(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  );
	    if(prod.Pt() > pt_cut) passed2=true;
	  }
	}
      }

      if(JAK2==3){
	for(std::vector<HepMC::GenParticle>::const_iterator a = SecondTauProducts.begin(); a!=SecondTauProducts.end(); ++a){
	  if(abs(a->pdg_id())==211){
	    TLorentzVector  prod(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  );
	    if(prod.Pt() > pt_cut) passed2=true;
	  }
	}
      }

      if(JAK2==4){
	TLorentzVector prod(0,0,0,0);
	for(std::vector<HepMC::GenParticle>::const_iterator a = SecondTauProducts.begin(); a!=SecondTauProducts.end(); ++a){
	  if(abs(a->pdg_id())==211){prod += TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  );}
	  if(abs(a->pdg_id())==111){prod += TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  );}
	}
	if(prod.Pt() > pt_cut ) passed2=true;
      }

      if(JAK2==5){
	TLorentzVector prod(0,0,0,0);
	prod+=TLorentzVector(A1Pions2.at(0).momentum().px(), A1Pions2.at(0).momentum().py(), A1Pions2.at(0).momentum().pz(), A1Pions2.at(0).momentum().e());
	prod+=TLorentzVector(A1Pions2.at(1).momentum().px(), A1Pions2.at(1).momentum().py(), A1Pions2.at(1).momentum().pz(), A1Pions2.at(1).momentum().e());
	prod+=TLorentzVector(A1Pions2.at(2).momentum().px(), A1Pions2.at(2).momentum().py(), A1Pions2.at(2).momentum().pz(), A1Pions2.at(2).momentum().e());
	if(prod.Pt() > pt_cut ) passed2=true;
      }
    }


    //---------------------------------------------------------------------------
    //--------------------  tau-
    if(JAK1==2){/* vu_muon */
      vector<TLorentzVector> tauandprod;
      tauandprod.push_back(TLorentzVector(FirstTau->momentum().px(), FirstTau->momentum().py(), FirstTau->momentum().pz(), FirstTau->momentum().e()));
      for(std::vector<HepMC::GenParticle>::const_iterator a = FirstTauProducts.begin(); a!=FirstTauProducts.end(); ++a){
	if(abs(a->pdg_id())==13){tauandprod.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );} }
      
      tauandprod1=tauandprod;
      Mu1.Configure(tauandprod1,"lepton");
      
      mu_plus->Fill(Mu1.getOmega(),HelWeightPlus);
      mu_minus->Fill(Mu1.getOmega(),HelWeightMinus);
    }
    
    if(JAK1==3){/* pi^- */
      vector<TLorentzVector> tauandprod;
      tauandprod.push_back(TLorentzVector(FirstTau->momentum().px(), FirstTau->momentum().py(), FirstTau->momentum().pz(), FirstTau->momentum().e()));
       for(std::vector<HepMC::GenParticle>::const_iterator a = FirstTauProducts.begin(); a!=FirstTauProducts.end(); ++a){
	 if(abs(a->pdg_id())==211){tauandprod.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );}}
       tauandprod1 = tauandprod;
       Pi1.Configure(tauandprod1,"pion");
       if(ApplyCut){
	 if(passed1){
	   pip_plus->Fill(Pi1.getOmega(),HelWeightPlus);
	   pip_minus->Fill(Pi1.getOmega(),HelWeightMinus);
	   pi_ptcut->Fill((tauandprod.at(1) ).Pt());
	 }
       }else{
	 pip_plus->Fill(Pi1.getOmega(),HelWeightPlus);
	 pip_minus->Fill(Pi1.getOmega(),HelWeightMinus);
	 pi_pt->Fill((tauandprod.at(1) ) .Pt());
       }


    }
    
  
    

    if(JAK1==4){/* rho^- */
       vector<TLorentzVector> tauandprod;
       tauandprod.push_back(TLorentzVector(FirstTau->momentum().px(), FirstTau->momentum().py(), FirstTau->momentum().pz(), FirstTau->momentum().e()));
       for(std::vector<HepMC::GenParticle>::const_iterator a = FirstTauProducts.begin(); a!=FirstTauProducts.end(); ++a){
     	if(abs(a->pdg_id())==211){
	  tauandprod.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );
	  Pion_Minus_From_Rho += TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e());}
     	
     	if(abs(a->pdg_id())==111){
	  tauandprod.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );
	  Pion_0_From_Rho += TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e());}
	}
       
       Pions_From_Rho_Mass->Fill((Pion_Minus_From_Rho + Pion_0_From_Rho).M());
       
       
       tauandprod1=tauandprod;  
       Rho1.Configure(tauandprod1,"rho");
       TauPolRho1.Configure(tauandprod1,"rho");
       RhoHelp.Configure(tauandprod1);
    	      


       if(ApplyCut){
	 if(passed1){
	   rhobeta_plus->Fill(Rho1.getVisibleOmega(),HelWeightPlus);
	   rhobeta_minus->Fill(Rho1.getVisibleOmega(),HelWeightMinus);
	   
	   omega_rho_plus->Fill(Rho1.getOmega(),HelWeightPlus);
	   omega_rho_minus->Fill(Rho1.getOmega(),HelWeightMinus);
	   
	   omegabar_rho_plus->Fill(Rho1.getOmegabar(),HelWeightPlus);
	   omegabar_rho_minus->Fill(Rho1.getOmegabar(),HelWeightMinus);
	   
	   cosbetacosthetarho_plus->Fill(RhoHelp.getCosbetaRho(),RhoHelp.getCosthetaRho(),HelWeightPlus);
	   cosbetacosthetarho_minus->Fill(RhoHelp.getCosbetaRho(),RhoHelp.getCosthetaRho(),HelWeightMinus);
	 }
       }else{
	 rhobeta_plus->Fill(Rho1.getVisibleOmega(),HelWeightPlus);
	 rhobeta_minus->Fill(Rho1.getVisibleOmega(),HelWeightMinus);
	 
	 omega_rho_plus->Fill(Rho1.getOmega(),HelWeightPlus);
	 omega_rho_minus->Fill(Rho1.getOmega(),HelWeightMinus);
	 
	 omegabar_rho_plus->Fill(Rho1.getOmegabar(),HelWeightPlus);
	 omegabar_rho_minus->Fill(Rho1.getOmegabar(),HelWeightMinus);
	 
	 cosbetacosthetarho_plus->Fill(RhoHelp.getCosbetaRho(),RhoHelp.getCosthetaRho(),HelWeightPlus);
	 cosbetacosthetarho_minus->Fill(RhoHelp.getCosbetaRho(),RhoHelp.getCosthetaRho(),HelWeightMinus);
       }

     }

    
    
   
     if(JAK1==5 /*&& SubJAK1==52*/){
      for(std::vector<HepMC::GenParticle>::const_iterator a = FirstTauProducts.begin(); a!=FirstTauProducts.end(); ++a){
	if(abs(a->pdg_id())==213){
	  Rho_From_A1 += TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e());}
     	if(abs(a->pdg_id())==111){
	  Pion_0_From_A1 += TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e());}
      }
      
      Rho_Pion_From_a1_Mass->Fill((Rho_From_A1 + Pion_0_From_A1).M());
     }
     
     
     
     
    //--------------------------------  a1 
    
   Int_t theta_min_width=0;
   Int_t theta_min_momentum_tau=0;
     /*
    Int_t nb_tau_momentum_plus;
    Int_t nb_tau_momentum_minus;
    */
    
    //Double_t theta_GJ;
    //Double_t theta_max;
    Int_t unphysical_theta = 0;
     if(JAK1==5 && SubJAK1==51){/* a1^- and 3 pions */

       vector<TLorentzVector> particles;
       particles.clear();
       SortPions(A1Pions1);
       int taucharge =  (A1Pions1.at(0).pdg_id()+A1Pions1.at(1).pdg_id()+A1Pions1.at(2).pdg_id() > 0) ? 1 : -1;
       taucharge1=taucharge;

       a1ss1pi.SetPxPyPzE(A1Pions1.at(0).momentum().px(), A1Pions1.at(0).momentum().py(), A1Pions1.at(0).momentum().pz(), A1Pions1.at(0).momentum().e());
       a1ss2pi.SetPxPyPzE(A1Pions1.at(1).momentum().px(), A1Pions1.at(1).momentum().py(), A1Pions1.at(1).momentum().pz(), A1Pions1.at(1).momentum().e());
       a1ospi.SetPxPyPzE(A1Pions1.at(2).momentum().px(), A1Pions1.at(2).momentum().py(), A1Pions1.at(2).momentum().pz(), A1Pions1.at(2).momentum().e());
       
      
       particles.push_back(tau1);
       particles.push_back(a1ospi);
       particles.push_back(a1ss1pi);
       particles.push_back(a1ss2pi);
       
       Polarimetr1.Configure(particles, taucharge);
       A1.Configure(particles,"a1",taucharge);
       a1hh.Configure(particles, a1ospi+a1ss1pi+a1ss2pi);
       TauPolA1.Configure(particles,"a1");
       tauandprod1=particles;
       
       /*
       TLorentzVector Pion1(0,0,0,0);
       TLorentzVector Pion2(0,0,0,0);
       TLorentzVector Pion3(0,0,0,0);
	Pion1+=TLorentzVector(A1Pions1.at(0).momentum().px(), A1Pions1.at(0).momentum().py(), A1Pions1.at(0).momentum().pz(), A1Pions1.at(0).momentum().e());
	Pion2+=TLorentzVector(A1Pions1.at(1).momentum().px(), A1Pions1.at(1).momentum().py(), A1Pions1.at(1).momentum().pz(), A1Pions1.at(1).momentum().e());
	Pion3+=TLorentzVector(A1Pions1.at(2).momentum().px(), A1Pions1.at(2).momentum().py(), A1Pions1.at(2).momentum().pz(), A1Pions1.at(2).momentum().e());
	
	
	
	std::cout<< "Pion1 mass " << a1ospi.M() - Pion1.M() << std::endl;
	std::cout<< "Pion1 mom " << a1ospi.M() - Pion1.M() << std::endl;
	std::cout<< "-------------------------" << std::endl;
	std::cout<< "Pion2 mass " << a1ss1pi.M() - Pion2.M() << std::endl;
	std::cout<< "Pion2 mom " <<  a1ss1pi.P() - Pion2.P() << std::endl;
	std::cout<< "-------------------------" << std::endl;
	std::cout<< "Pion3 mass " << a1ss2pi.M() - Pion3.M() << std::endl;
	std::cout<< "Pion3 mom " << a1ss2pi.P() - Pion3.P() << std::endl;
	std::cout<< "-------------------------" << std::endl;
	std::cout<< "angle a1ss " << tau1.Angle((a1ospi + a1ss1pi + a1ss2pi).Vect())<< std::endl;
	std::cout<< "angle pi " << tau1.Angle((Pion1 + Pion2 + Pion3).Vect())<< std::endl;
	std::cout<< "-------------------------" << std::endl;
	*/
	
	
	//tau1 ---- this is a TLorentzVector of tau1
	TLorentzVector a1LV1 = a1ospi + a1ss1pi + a1ss2pi;
	
	
	TLorentzVector nuLV = tau1 - a1LV1;
	
	Pions_From_a1_Mass->Fill(a1LV1.M());
	
	//Double_t m_a1 = 1.230; 
	
	//std::cout<< "mass of tau " << tau1.M() <<std::endl;
	
	//std::cout<< "mass of a1 (pions) " << a1LV.M() <<std::endl;
	
	//std::cout<< "momentum of a1 (pions) " << a1LV.P() <<std::endl;
	// Gottfried-Jackson angle
	Double_t theta_max1 = asin((tau1.M()*tau1.M() - a1LV1.M()*a1LV1.M())/(2*tau1.M()*a1LV1.P())); 
	
	
	//std::cout<< "theta max " << theta_max <<std::endl;
	
	Double_t theta_GJ1 = tau1.Angle(a1LV1.Vect());
	Double_t theta = a1LV1.Theta();
	 Double_t phi = a1LV1.Phi();
	 
	Double_t theta_CM = asin(sin(theta_GJ1)/sin(theta_max1)); 
	 
	TVector3 n_tau = tau1.Vect()*(1/tau1.Vect().Mag());
	TVector3 n_a1 = a1LV1.Vect()*(1/a1LV1.Vect().Mag());
	
	
	TVector3 v_tau;
	v_tau.SetXYZ(tau1.Px()*(1/tau1.E()),tau1.Py()*(1/tau1.E()),tau1.Pz()*(1/tau1.E()));
	
	Double_t beta_tau = tau1.E()/tau1.M();
	
	TVector3 n_a1CM = a1LV1.Vect()*(1/beta_tau) - v_tau*((tau1.M()*tau1.M() + a1LV1.M()*a1LV1.M())/(2*tau1.M()) );
	
	Double_t scalar_tau_a1 = n_tau.Dot(n_a1CM);
	
	if(scalar_tau_a1<0.){theta_CM = TMath::Pi() - theta_CM;}
	
	//a1ospi+a1ss1pi+a1ss2pi
	
	TVector3 my_nss1 = a1ss1pi.Vect()*(1/a1ss1pi.Vect().Mag());
	TVector3 my_nss2 = a1ss2pi.Vect()*(1/a1ss2pi.Vect().Mag());
	TVector3 my_nPerp = (my_nss1.Cross(my_nss2))*(1/(my_nss1.Cross(my_nss2)).Mag());
	TVector3 my_nL = -1*a1LV1.Vect()*(1/a1LV1.Vect().Mag()); /////////////////// nL = +/- n_a1 ????????????????
	TVector3 my_nT = tau1.Vect()*(1/tau1.Vect().Mag());
	TVector3 my_nPerpCrossnL  = my_nPerp.Cross(my_nL);
	TVector3 my_qvect = a1ospi.Vect()*(1/a1ospi.Vect().Mag());
  
  
	//---------------- Angle beta --------------------
	Double_t cos_beta = my_nT.Dot(my_nPerp);
	Double_t sin_beta = sqrt(1 - cos_beta*cos_beta);
	//------------------------------------------------
  
	//--------------- Angle gamma --------------------
	Double_t sin_gamma = my_nPerpCrossnL.Dot(my_qvect)*(1/my_nPerpCrossnL.Mag());
	Double_t cos_gamma_ana = -1*my_nL.Dot(my_qvect)*(1/my_nPerpCrossnL.Mag());
	Double_t cos_gamma_sin = sqrt(1 - sin_gamma*sin_gamma);
	
	
	
	Double_t my_gamma_ana;
	Double_t my_gamma_sin;
	Double_t a1hh_gamma;
	Double_t a1hh_gammaLF;
	Double_t a1hh_TRF_gamma;
	
	Double_t my_beta;
	Double_t a1hh_beta;
	Double_t a1hh_betaLF;
	Double_t a1hh_TRF_beta;
	
	if(cos_beta > 0. && sin_beta > 0.){my_beta = acos(cos_beta);}
	if(cos_beta < 0. && sin_beta > 0.){my_beta = acos(cos_beta);}
	if(cos_beta > 0. && sin_beta < 0.){my_beta = asin(sin_beta);}
	if(cos_beta < 0. && sin_beta < 0.){my_beta = asin(-1*sin_beta) + TMath::Pi()/2;}
	
	
	if(cos_gamma_ana > 0. && sin_gamma > 0.){my_gamma_ana = acos(cos_gamma_ana);}
	if(cos_gamma_ana < 0. && sin_gamma > 0.){my_gamma_ana = acos(cos_gamma_ana);}
	if(cos_gamma_ana > 0. && sin_gamma < 0.){my_gamma_ana = asin(sin_gamma);}
	if(cos_gamma_ana < 0. && sin_gamma < 0.){my_gamma_ana = asin(-1*sin_gamma) + TMath::Pi()/2;}
	
	if(cos_gamma_sin > 0. && sin_gamma > 0.){my_gamma_sin = acos(cos_gamma_sin);}
	if(cos_gamma_sin < 0. && sin_gamma > 0.){my_gamma_sin = acos(cos_gamma_sin);}
	if(cos_gamma_sin > 0. && sin_gamma < 0.){my_gamma_sin = asin(sin_gamma);}
	if(cos_gamma_sin < 0. && sin_gamma < 0.){my_gamma_sin = asin(-1*sin_gamma) + TMath::Pi()/2;}
	
	//-----------------------------------------------------------
  
	if(a1hh.cosgamma() > 0. && a1hh.singamma() > 0.){a1hh_gamma = acos(a1hh.cosgamma());}
	if(a1hh.cosgamma() < 0. && a1hh.singamma() > 0.){a1hh_gamma = acos(a1hh.cosgamma());}
	if(a1hh.cosgamma() > 0. && a1hh.singamma() < 0.){a1hh_gamma = asin(a1hh.singamma());}
	if(a1hh.cosgamma() < 0. && a1hh.singamma() < 0.){a1hh_gamma = asin(-1*a1hh.singamma()) + TMath::Pi()/2;}
	
	if(a1hh.cosgammaLF() > 0. && a1hh.singammaLF() > 0.){a1hh_gammaLF = acos(a1hh.cosgammaLF());}
	if(a1hh.cosgammaLF() < 0. && a1hh.singammaLF() > 0.){a1hh_gammaLF = acos(a1hh.cosgammaLF());}
	if(a1hh.cosgammaLF() > 0. && a1hh.singammaLF() < 0.){a1hh_gammaLF = asin(a1hh.singammaLF());}
	if(a1hh.cosgammaLF() < 0. && a1hh.singammaLF() < 0.){a1hh_gammaLF = asin(-1*a1hh.singammaLF()) + TMath::Pi()/2;}
	
	if(a1hh.TRF_cosgamma() > 0. && a1hh.TRF_singamma() > 0.){a1hh_TRF_gamma = acos(a1hh.TRF_cosgamma());}
	if(a1hh.TRF_cosgamma() < 0. && a1hh.TRF_singamma() > 0.){a1hh_TRF_gamma = acos(a1hh.TRF_cosgamma());}
	if(a1hh.TRF_cosgamma() > 0. && a1hh.TRF_singamma() < 0.){a1hh_TRF_gamma = asin(a1hh.TRF_singamma());}
	if(a1hh.TRF_cosgamma() < 0. && a1hh.TRF_singamma() < 0.){a1hh_TRF_gamma = asin(-1*a1hh.TRF_singamma()) + TMath::Pi()/2;}
        //
	if(a1hh.cosbeta() > 0. && a1hh.sinbeta() > 0.){a1hh_beta = acos(a1hh.cosbeta());}
	if(a1hh.cosbeta() < 0. && a1hh.sinbeta() > 0.){a1hh_beta = acos(a1hh.cosbeta());}
	if(a1hh.cosbeta() > 0. && a1hh.sinbeta() < 0.){a1hh_beta = asin(a1hh.sinbeta());}
	if(a1hh.cosbeta() < 0. && a1hh.sinbeta() < 0.){a1hh_beta = asin(-1*a1hh.sinbeta()) + TMath::Pi()/2;}
	
	
	Double_t a1hh_sinbetaLF = sqrt(1 - a1hh.cosbetaLF()*a1hh.cosbetaLF());
	
	if(a1hh.cosbetaLF() > 0. && a1hh_sinbetaLF > 0.){a1hh_betaLF = acos(a1hh.cosbetaLF());}
	if(a1hh.cosbetaLF() < 0. && a1hh_sinbetaLF > 0.){a1hh_betaLF = acos(a1hh.cosbetaLF());}
	if(a1hh.cosbetaLF() > 0. && a1hh_sinbetaLF < 0.){a1hh_betaLF = asin(a1hh_sinbetaLF);}
	if(a1hh.cosbetaLF() < 0. && a1hh_sinbetaLF < 0.){a1hh_betaLF = asin(-1*a1hh_sinbetaLF) + TMath::Pi()/2;}
	
	if(a1hh.TRF_cosbeta() > 0. && a1hh.TRF_sinbeta() > 0.){a1hh_TRF_beta = acos(a1hh.TRF_cosbeta());}
	if(a1hh.TRF_cosbeta() < 0. && a1hh.TRF_sinbeta() > 0.){a1hh_TRF_beta = acos(a1hh.TRF_cosbeta());}
	if(a1hh.TRF_cosbeta() > 0. && a1hh.TRF_sinbeta() < 0.){a1hh_TRF_beta = asin(a1hh.TRF_sinbeta());}
	if(a1hh.TRF_cosbeta() < 0. && a1hh.TRF_sinbeta() < 0.){a1hh_TRF_beta = asin(-1*a1hh.TRF_sinbeta()) + TMath::Pi()/2;}
	
	
	
	
	std::cout<< "------------------------------" <<std::endl;
	std::cout<< "a1ss1pi.P() " << a1ss1pi.P() <<std::endl;
	std::cout<< "a1ss2pi.P() " << a1ss2pi.P() <<std::endl;
	
	
	std::cout<< "------------------------------" <<std::endl;
	std::cout<< "my_beta " << my_beta*(180/TMath::Pi()) <<std::endl;
	std::cout<< "a1hh_beta " << a1hh_beta*(180/TMath::Pi())  <<std::endl;
	std::cout<< "a1hh_betaLF " << a1hh_betaLF*(180/TMath::Pi())  <<std::endl;
	std::cout<< "a1hh_TRF_beta " << a1hh_TRF_beta*(180/TMath::Pi()) <<std::endl;
	std::cout<< "------------------------------" <<std::endl;
	std::cout<< "my_gamma_ana " << my_gamma_ana*(180/TMath::Pi())  <<std::endl;
	std::cout<< "my_gamma_sin " << my_gamma_sin*(180/TMath::Pi())  <<std::endl;
	std::cout<< "a1hh_gamma " << a1hh_gamma*(180/TMath::Pi()) <<std::endl;
	std::cout<< "a1hh_gammaLF " <<  a1hh_gammaLF*(180/TMath::Pi()) <<std::endl;
	std::cout<< "a1hh_TRF_gamma " << a1hh_TRF_gamma*(180/TMath::Pi())  <<std::endl;
	
	
	std::cout<< "------------------------------" <<std::endl;
	std::cout<< "my cos_gamma_ana " << cos_gamma_ana  <<std::endl;
	std::cout<< "my cos_gamma_sin " << cos_gamma_sin  <<std::endl;
	std::cout<< "a1hh.cosgamma() " << a1hh.cosgamma()  <<std::endl;
	std::cout<< "a1hh.cosgammaLF() " << a1hh.cosgammaLF()  <<std::endl;
	std::cout<< "a1hh.TRF_cosgamma() " << a1hh.TRF_cosgamma()  <<std::endl;
	std::cout<< "my sin_gamma " << sin_gamma  <<std::endl;
	std::cout<< "a1hh.singamma() " << a1hh.singamma()  <<std::endl;
	std::cout<< "a1hh.singammaLF() " << a1hh.singammaLF() <<std::endl;
	std::cout<< "a1hh.TRF_singamma() " <<  a1hh.TRF_singamma() <<std::endl;
	std::cout<< "------------------------------" <<std::endl;
	std::cout<< "my cos_beta " << cos_beta  <<std::endl;
	std::cout<< "a1hh.cosbeta() " << a1hh.cosbeta() <<std::endl;
	std::cout<< "a1hh.TRF_cosbeta() " << a1hh.TRF_cosbeta()  <<std::endl;
	std::cout<< "a1hh.cosbetaLF() " << a1hh.cosbetaLF()  <<std::endl;
	std::cout<< "my sin_beta " << sin_beta  <<std::endl;
	std::cout<< "a1hh.sinbeta() " << a1hh.sinbeta()  <<std::endl;
	std::cout<< "a1hh.TRF_sinbeta() " << a1hh.TRF_sinbeta() <<std::endl;
	std::cout<< "a1hh_sinbetaLF " << a1hh_sinbetaLF  <<std::endl;
	
	/*
	std::cout<< "------------------------------" <<std::endl;
	std::cout<< "scalar " << scalar_tau_a1 <<std::endl;
	std::cout<< "cos theta CM " <<cos(theta_CM) <<std::endl;
	std::cout<< "theta CM rad " << theta_CM <<std::endl;
	std::cout<< "theta CM deg " << theta_CM*(180/3.14) <<std::endl;
	*/
	
	 /*
	  if(tau1.P()>45){
	  std::cout<< "------------------------------" <<std::endl;
	  std::cout<< "theta GJ " << theta_GJ1 <<std::endl;
	  std::cout<< "p*cos(theta_GJ) " << a1LV1.P()*cos(theta_GJ1) <<std::endl;
	  std::cout<< "p*sin(theta_GJ) " << a1LV1.P()*sin(theta_GJ1) <<std::endl;
	  std::cout<< "px^2 + py^2 " << sqrt(a1LV1.Px()*a1LV1.Px() + a1LV1.Py()*a1LV1.Py()) <<std::endl;
	  std::cout<< "py^2 + pz^2 " << sqrt(a1LV1.Py()*a1LV1.Py() + a1LV1.Pz()*a1LV1.Pz()) <<std::endl;
	  std::cout<< "pz^2 + px^2 " << sqrt(a1LV1.Pz()*a1LV1.Pz() + a1LV1.Px()*a1LV1.Px()) <<std::endl;
	  std::cout<< "pT " << a1LV1.Pt() <<std::endl;
	  std::cout<< "p perp " << a1LV1.Perp() <<std::endl;
	  std::cout<< "pT tau " << tau1.Pt() <<std::endl; 
	  std::cout<< "theta CM " << asin((sin(theta_GJ1))/(sin(theta_max1))) <<std::endl;
	 }
	*/
	
	// General
	
	Double_t a1 = (a1LV1.M()*a1LV1.M() + (a1LV1.P()*sin(theta_GJ1))*(a1LV1.P()*sin(theta_GJ1)));
	Double_t b1 = -1*a1LV1.P()*cos(theta_GJ1)*(a1LV1.M()*a1LV1.M() + tau1.M()*tau1.M());
	Double_t delta1 = (a1LV1.P()*a1LV1.P() + a1LV1.M()*a1LV1.M())*((a1LV1.M()*a1LV1.M() - tau1.M()*tau1.M())*(a1LV1.M()*a1LV1.M() - tau1.M()*tau1.M()) -4*tau1.M()*tau1.M()*a1LV1.P()*a1LV1.P()*sin(theta_GJ1)*sin(theta_GJ1));
	/*
	Double_t a = 4*(a1LV.M()*a1LV.M() + (a1LV.P()*sin(theta_GJ))*(a1LV.P()*sin(theta_GJ)));
	Double_t b = -4*a1LV.P()*cos(theta_GJ)*(a1LV.M()*a1LV.M() + tau1.M()*tau1.M());
	Double_t c = 4*a1LV.E()*a1LV.E()*a1LV.P()*a1LV.P() - (tau1.M()*tau1.M() - a1LV.E()*a1LV.E() -a1LV.P()*a1LV.P())*(tau1.M()*tau1.M() -a1LV.E()*a1LV.E() -a1LV.P()*a1LV.P());
	Double_t delta = b*b - 4*a*c;
	*/
	// Mass a1 constant
	Double_t m_a1 = 1.230; //GeV
	Double_t am1 = (m_a1*m_a1 + (a1LV1.P()*sin(theta_GJ1))*(a1LV1.P()*sin(theta_GJ1)));
	Double_t bm1 = -1*a1LV1.P()*cos(theta_GJ1)*(m_a1*m_a1 + tau1.M()*tau1.M());
	Double_t deltam1 = (a1LV1.P()*a1LV1.P() + m_a1*m_a1)*((m_a1*m_a1 - tau1.M()*tau1.M())*(m_a1*m_a1 - tau1.M()*tau1.M()) -4*tau1.M()*tau1.M()*a1LV1.P()*a1LV1.P()*sin(theta_GJ1)*sin(theta_GJ1));
	
	// Momentum a1 constant
	Double_t p_a1 = 20.; //GeV
	Double_t ap1 = (a1LV1.M()*a1LV1.M() + (p_a1*sin(theta_GJ1))*(p_a1*sin(theta_GJ1)));
	Double_t bp1 = -1*p_a1*cos(theta_GJ1)*(a1LV1.M()*a1LV1.M() + tau1.M()*tau1.M());
	Double_t deltap1 = (p_a1*p_a1 + a1LV1.M()*a1LV1.M())*((a1LV1.M()*a1LV1.M() - tau1.M()*tau1.M())*(a1LV1.M()*a1LV1.M() - tau1.M()*tau1.M()) -4*tau1.M()*tau1.M()*p_a1*p_a1*sin(theta_GJ1)*sin(theta_GJ1));
	  /*
	  std::cout<< "------------------------------------"<<std::endl;
	  std::cout<< "  a1 invariant mass     " << a1LV.M() << "  tau invariant mass " << tau1.M() <<std::endl;
	  std::cout<< "  delta NOT fixed    " << delta << " delta fixed " << deltam <<std::endl;
	  std::cout<< "------------------------------------"<<std::endl;
	  */
	Double_t tau_momentum_minus1 =  (-b1-sqrt(delta1)) / (2*a1);
	
	Double_t tau_momentum_plus1 =  (-b1+sqrt(delta1)) / (2*a1);
	
	Double_t tau_momentum_mean1 = (tau_momentum_minus1 + tau_momentum_plus1) /2;
	
	Double_t tau_momentum_width1 = tau_momentum_plus1 - tau_momentum_minus1;
	
	/*
	if(abs(tau_momentum_plus1 - tau_momentum_minus1)<0.5 && tau1.M()*tau1.M() < 2*tau1.M()*a1LV1.P()*sin(-1*theta_GJ1)) {
	std::cout<< "tau1 mass "<< tau1.M()<<std::endl;
	std::cout<< "a1 mass truth "<< a1LV1.M()<<std::endl;
	std::cout<< "a1 mass + "<< sqrt(tau1.M()*tau1.M() + 2*tau1.M()*a1LV1.P()*sin(theta_GJ1)) <<std::endl;
	std::cout<< "a1 mass - "<< sqrt(tau1.M()*tau1.M() - 2*tau1.M()*a1LV1.P()*sin(theta_GJ1))<<std::endl;
	std::cout<< "----------------" <<std::endl;
	std::cout<< "theta truth "<< theta_GJ1 <<std::endl;
	std::cout<< "theta + "<< asin((tau1.M()*tau1.M() -a1LV1.M()*a1LV1.M())/(2*tau1.M()*a1LV1.P())) <<std::endl;
	std::cout<< "theta - "<< -1*asin((tau1.M()*tau1.M() -a1LV1.M()*a1LV1.M())/(2*tau1.M()*a1LV1.P())) <<std::endl;
	std::cout<< "----------------" <<std::endl;
	std::cout<< "a1 mom truth "<< a1LV1.P()<<std::endl;
	std::cout<< "a1 mom + "<< (a1LV1.M()*a1LV1.M() -tau1.M()*tau1.M())/(2*tau1.M()*sin(theta_GJ1)) <<std::endl;
	std::cout<< "a1 mom - "<< -1*(a1LV1.M()*a1LV1.M() -tau1.M()*tau1.M())/(2*tau1.M()*sin(theta_GJ1)) <<std::endl;
	std::cout<< "----------------" <<std::endl;
	std::cout<< "+      " << tau_momentum_plus1 <<std::endl;
	std::cout<< "-      " << tau_momentum_minus1 <<std::endl;
	std::cout<< "mean   " << tau_momentum_mean1 <<std::endl;
	}
	*/
	
	
	
	
	
	/*
	std::cout<< "a1 mom " << a1LV.P() <<std::endl;
	std::cout<< "a1 mass "<< a1LV.M()<<std::endl;
	std::cout<< "a1 energy "<< a1LV.E()<<std::endl;
	std::cout<< "a1 gamma "<< a1LV.E()/a1LV.M()<<std::endl;
	
	std::cout<< "tau mom "<< tau1.P()<<std::endl;
	std::cout<< "tau mass " << tau1.M() <<std::endl;
	std::cout<< "tau energy "<< tau1.E()<<std::endl;
	std::cout<< "tau gamma "<< tau1.E()/tau1.M()<<std::endl;
	
	std::cout<< "       " <<(a1LV.M()*a1LV.M() - tau1.M()*tau1.M())*(a1LV.M()*a1LV.M() - tau1.M()*tau1.M()) <<std::endl;
	std::cout<< "       " <<-4*tau1.M()*tau1.M()*a1LV.P()*a1LV.P()*sin(theta_GJ)*sin(theta_GJ)<<std::endl;
	std::cout<< "sin    " << sin(theta_GJ) << std::endl;
	std::cout<< "theta   "<< theta_GJ << std::endl;
	std::cout<< "theta max  "<< theta_max << std::endl;
	std::cout<< "----------------------------"<<std::endl;
	*/
	
	Double_t tau_momentum_gen1;
	Double_t tau_momentum_rec1;
	 
	//Double_t nb_tau_momentum[2];
	//nb_tau_momentum[0]=nb_tau_momentum_plus;
	//nb_tau_momentum[1]=nb_tau_momentum_minus;
	
	/*
	std::cout << "a1 gen M " << a1_gen.M() << "  a1 rec M " << a1LV.M() << std::endl;
	std::cout << "a1 gen P " << a1_gen.P() << "  a1 rec P " << a1LV.P() << std::endl;
	*/
	
	//Double_t x;
	/*
	for (Double_t x=0; x<0.5; x += 0.005){
	  Bool_t test = abs(tau1.P() - tau_momentum_plus)<x;
	  if(test == true){
	    tau_momentum_rec = tau_momentum_plus;
	    absolute_error_momentum_tau_plus_vs_x->Fill(x,tau_momentum_rec - tau1.P());
	  }
	  if(test == false){
	    tau_momentum_rec = tau_momentum_minus;
	    absolute_error_momentum_tau_minus_vs_x->Fill(x,tau_momentum_rec - tau1.P());
	  }
	  error_vs_x->Fill(x,sqrt((tau_momentum_plus - tau1.P())*(tau_momentum_plus - tau1.P()) + (tau_momentum_minus - tau1.P())*(tau_momentum_minus - tau1.P())));
	  absolute_error_momentum_tau_mean_vs_x->Fill(x,tau_momentum_mean - tau1.P());
	  std::cout << "x = " << x << std::endl;
	  std::cout << "p+ - p(gen) " << tau_momentum_plus - tau1.P() << std::endl;
	  std::cout << "p- - p(gen) " << tau_momentum_minus - tau1.P()<< std::endl;
	  std::cout << "--------------------------------" << std::endl;
	}
	*/
	Bool_t plus = abs(tau1.P() - tau_momentum_plus1)<0.5;
	
	Int_t nb_tau_momentum_plus1 =0;
	Int_t nb_tau_momentum_minus1 =0;
	
	if (plus == true){
	  nb_tau_momentum_plus1++;
	  tau_momentum_rec1 = tau_momentum_plus1;
	  tau_momentum_gen1 = tau1.P();
	}
	else{
	  nb_tau_momentum_minus1++;
	  tau_momentum_rec1 = tau_momentum_minus1;  
	  tau_momentum_gen1 = tau1.P();
	}
	
	//##############################################################################################################
	//##################################  PAIRS OF TAU #############################################################
	//##############################################################################################################
	if(JAK2==5 && SubJAK2==51){
	vector<TLorentzVector> particles;
	particles.clear();
	SortPions(A1Pions2);
      
	int taucharge =  (A1Pions2.at(0).pdg_id()+A1Pions2.at(1).pdg_id()+A1Pions2.at(2).pdg_id() > 0) ? 1 : -1;
	taucharge2=taucharge;
	a1ss1pi.SetPxPyPzE(A1Pions2.at(0).momentum().px(), A1Pions2.at(0).momentum().py(), A1Pions2.at(0).momentum().pz(), A1Pions2.at(0).momentum().e());
	a1ss2pi.SetPxPyPzE(A1Pions2.at(1).momentum().px(), A1Pions2.at(1).momentum().py(), A1Pions2.at(1).momentum().pz(), A1Pions2.at(1).momentum().e());
	a1ospi.SetPxPyPzE(A1Pions2.at(2).momentum().px(), A1Pions2.at(2).momentum().py(), A1Pions2.at(2).momentum().pz(), A1Pions2.at(2).momentum().e());
      
	particles.push_back(tau2);
	particles.push_back(a1ospi);
	particles.push_back(a1ss1pi);
	particles.push_back(a1ss2pi);
	
	a1h.Configure(particles, a1ospi+a1ss1pi+a1ss2pi);
	TauPolA1.Configure(particles,"a1");
	tauandprod2=particles;
	Polarimetr.Configure(particles, taucharge);
	A2.Configure(particles,"a1", taucharge);
	
	
	
	//tau2 ---- this is a TLorentzVector of tau2
	TLorentzVector a1LV2 = a1ospi + a1ss1pi + a1ss2pi;
	TLorentzVector nuLV2 = tau2 - a1LV2;
	
	//Pions_From_a1_Mass->Fill(a1LV2.M());
	
	//Double_t m_a1 = 1.230; 
	
	//std::cout<< "mass of tau " << tau1.M() <<std::endl;
	
	//std::cout<< "mass of a1 (pions) " << a1LV.M() <<std::endl;
	
	//std::cout<< "momentum of a1 (pions) " << a1LV.P() <<std::endl;
	// Gottfried-Jackson angle
	Double_t theta_max2 = asin((tau2.M()*tau2.M() - a1LV2.M()*a1LV2.M())/(2*tau2.M()*a1LV2.P())); 
	
	
	//std::cout<< "theta max " << theta_max <<std::endl;
	
	Double_t theta_GJ2 = tau2.Angle(a1LV2.Vect());
	
	
	// General
	
	Double_t a2 = (a1LV2.M()*a1LV2.M() + (a1LV2.P()*sin(theta_GJ2))*(a1LV2.P()*sin(theta_GJ2)));
	Double_t b2 = -1*a1LV2.P()*cos(theta_GJ2)*(a1LV2.M()*a1LV2.M() + tau2.M()*tau2.M());
	Double_t delta2 = (a1LV2.P()*a1LV2.P() + a1LV2.M()*a1LV2.M())*((a1LV2.M()*a1LV2.M() - tau2.M()*tau2.M())*(a1LV2.M()*a1LV2.M() - tau2.M()*tau2.M()) -4*tau2.M()*tau2.M()*a1LV2.P()*a1LV2.P()*sin(theta_GJ2)*sin(theta_GJ2));
	
	// Mass a1 constant
	Double_t m_a1 = 1.230; //GeV
	Double_t am2 = (m_a1*m_a1 + (a1LV2.P()*sin(theta_GJ2))*(a1LV2.P()*sin(theta_GJ2)));
	Double_t bm2 = -1*a1LV2.P()*cos(theta_GJ2)*(m_a1*m_a1 + tau2.M()*tau2.M());
	Double_t deltam2 = (a1LV2.P()*a1LV2.P() + m_a1*m_a1)*((m_a1*m_a1 - tau2.M()*tau2.M())*(m_a1*m_a1 - tau2.M()*tau2.M()) -4*tau2.M()*tau2.M()*a1LV2.P()*a1LV2.P()*sin(theta_GJ2)*sin(theta_GJ2));
	
	// Momentum a1 constant
	Double_t p_a1 = 20.; //GeV
	Double_t ap2 = (a1LV2.M()*a1LV2.M() + (p_a1*sin(theta_GJ2))*(p_a1*sin(theta_GJ2)));
	Double_t bp2 = -1*p_a1*cos(theta_GJ2)*(a1LV2.M()*a1LV2.M() + tau2.M()*tau2.M());
	Double_t deltap2 = (p_a1*p_a1 + a1LV2.M()*a1LV2.M())*((a1LV2.M()*a1LV2.M() - tau2.M()*tau2.M())*(a1LV2.M()*a1LV2.M() - tau2.M()*tau2.M()) -4*tau2.M()*tau2.M()*p_a1*p_a1*sin(theta_GJ2)*sin(theta_GJ2));
	  /*
	  std::cout<< "------------------------------------"<<std::endl;
	  std::cout<< "  a1 invariant mass     " << a1LV.M() << "  tau invariant mass " << tau1.M() <<std::endl;
	  std::cout<< "  delta NOT fixed    " << delta << " delta fixed " << deltam <<std::endl;
	  std::cout<< "------------------------------------"<<std::endl;
	  */
	Double_t tau_momentum_minus2 =  (-b2-sqrt(delta2)) / (2*a2);
	
	Double_t tau_momentum_plus2 =  (-b2+sqrt(delta2)) / (2*a2);
	
	Double_t tau_momentum_mean2 = (tau_momentum_minus2 + tau_momentum_plus2) /2;
	
	Double_t tau_momentum_width2 = tau_momentum_plus2 - tau_momentum_minus2;
	
	Double_t tau_momentum_gen2;
	Double_t tau_momentum_rec2;
	
	Bool_t plus = abs(tau2.P() - tau_momentum_plus2)<0.5;
	
	Int_t nb_tau_momentum_plus2 =0;
	Int_t nb_tau_momentum_minus2 =0;
	
	if (plus == true){
	  nb_tau_momentum_plus2++;
	  tau_momentum_rec2 = tau_momentum_plus2;
	  tau_momentum_gen2 = tau2.P();
	}
	else{
	  nb_tau_momentum_minus2++;
	  tau_momentum_rec2 = tau_momentum_minus2;  
	  tau_momentum_gen2 = tau2.P();
	}
/*
	std::cout<< "tau1 gen " << tau1.P() <<std::endl;
	std::cout<< "tau1 rec " << tau_momentum_mean1 <<std::endl;
	std::cout<< "--------------------------" <<std::endl;
	std::cout<< "tau2 gen " << tau2.P() <<std::endl;
	std::cout<< "tau2 rec " << tau_momentum_mean2 <<std::endl;
	std::cout<< "--------------------------" <<std::endl;
	*/
	
	
	
	}
	
	
	
	
	
	//if (delta1 < 1){
	  ratio_ma1_plus->Fill(a1LV1.M(),sqrt(tau1.M()*tau1.M() + 2*tau1.M()*a1LV1.P()*sin(theta_GJ1)));
	  ratio_ma1_minus->Fill(a1LV1.M(),sqrt(tau1.M()*tau1.M() - 2*tau1.M()*a1LV1.P()*sin(theta_GJ1)));
	//}
	
	
	// Fill histogramms
	
	
	//------------------  theta CM ------------------------
	theta_cm->Fill(theta_CM);
	costheta_cm->Fill(cos(theta_CM));
	
	theta_CM_vs_a1_momentum->Fill(theta_CM,a1_gen.P());
	theta_CM_vs_tau_momentum_gen->Fill(theta_CM,tau1.P());
	
	
	costheta_CM_vs_a1_momentum->Fill(cos(theta_CM),a1_gen.P());
	costheta_CM_vs_tau_momentum_gen->Fill(cos(theta_CM),tau1.P());
	
	
	//----------------- negative rootsquare ----------------
	if (delta1<0.){
	  negative_sqrt_a1_mass->Fill(a1LV1.M());
	  negative_sqrt_a1_momentum->Fill(a1LV1.P());}
	if (deltam1<0.){
	  negative_sqrt_fixed_mass_a1_momentum->Fill(a1LV1.P());}
	//-------------------------------------------------------
	
	//------------ errors mean correlations ------------------
	if(abs(tau_momentum_mean1 - tau1.P())>20){
	error_a1_mass_vs_a1_momentum->Fill(a1LV1.M(),a1LV1.P());
	error_a1_mass_vs_theta->Fill(a1LV1.M(),theta_GJ1);
	error_theta_vs_a1_momentum->Fill(theta_GJ1,a1LV1.P());
	}
	//--------------------------------------------------------
	
	//---------------- tau informations ---------------------------
	tau1_gen_Mass->Fill(tau1.M());
	tau1_gen_Momentum->Fill(tau1.P());
	tau1_gen_pt->Fill(tau1.Pt());
	//------------------------------------------------------------
	
	//--------------- a1 informations -------------------------
	a1_rec_Mass->Fill(a1LV1.M());
	a1_rec_Momentum->Fill(a1LV1.P());
	a1_rec_pt->Fill(a1LV1.Pt());
	
	a1_gen_Mass->Fill(a1_gen.M());
	a1_gen_Momentum->Fill(a1_gen.P());
	a1_gen_pt->Fill(a1_gen.Pt());
	//-----------------------------------------------------------------
	
	//----------------------------------------------  absolute error -------------------------------------------------------
	absolute_error_mass_a1->Fill(a1LV1.M() - a1_gen.M());
	absolute_error_momentum_a1->Fill(a1LV1.P() - a1_gen.P());
	//----------------------------------------------------------------------------------------------------------------------
  
	//----------------------------------------------  relative error -------------------------------------------------------
	relative_error_mass_a1->Fill((a1LV1.M() - a1_gen.M())/a1_gen.M());
	relative_error_momentum_a1->Fill((a1LV1.P() - a1_gen.P())/a1_gen.P());
	//----------------------------------------------------------------------------------------------------------------------
  
	//----------------------------------------------  relative error abs---------------------------------------------------
	relative_error_abs_mass_a1->Fill(abs((a1LV1.M() - a1_gen.M())/a1_gen.M()));
	relative_error_abs_momentum_a1->Fill(abs((a1LV1.P() - a1_gen.P())/a1_gen.P()));
	//---------------------------------------------------------------------------------------------------------------------
	
	
	//---------------------------------------------------------
	
	//------------ theta vs a1 --------------------------------
	theta_GJ_vs_a1_momentum->Fill(a1LV1.P(),theta_GJ1);
	theta_max_vs_a1_momentum->Fill(a1LV1.P(),theta_max1);
	costheta_GJ_vs_a1_momentum->Fill(a1LV1.P(),cos(theta_GJ1));
	costheta_max_vs_a1_momentum->Fill(a1LV1.P(),cos(theta_max1));
	//----------------------------------------------------------
	
	//------------- theta / theta_max --------------------------
	ratio_theta->Fill(theta_GJ1 / theta_max1);
	theta_max_vs_theta_GJ->Fill(theta_GJ1,theta_max1);
	//----------------------------------------------------------
	
	//-------------- physical meaning of theta --------------------------
	if (theta_GJ1 > theta_max1){
	  Unphysics_theta->Fill(theta_GJ1);
	  Unphysics_a1_momentum->Fill(a1LV1.P());}
	else{
	  Physical_theta_GJ_vs_a1_momentum->Fill(a1LV1.P(),theta_GJ1);
	  Physical_theta_max_vs_a1_momentum->Fill(a1LV1.P(),theta_max1);
	  Physical_costheta_GJ_vs_a1_momentum->Fill(a1LV1.P(),cos(theta_GJ1));
	  Physical_costheta_max_vs_a1_momentum->Fill(a1LV1.P(),cos(theta_max1));}
	//---------------------------------------------------------------------------
	
	//################################# relative errors +/-/average solutions #####################################
  
	//###################################################  PLUS ##################################################################
	//############################################################################################################################
	//----------------------------------------------  absolute error -------------------------------------------------------
	if(plus == true){absolute_error_momentum_tau_plus->Fill(tau_momentum_plus1 - tau1.P());} 
	//----------------------------------------------------------------------------------------------------------------------
  
	//----------------------------------------------  relative error -------------------------------------------------------
	if(plus == true){relative_error_momentum_tau_plus->Fill((tau_momentum_plus1 - tau1.P())/tau1.P());} 
        //----------------------------------------------------------------------------------------------------------------------
  
	//----------------------------------------------  relative error abs---------------------------------------------------
	if(plus == true){relative_error_abs_momentum_tau_plus->Fill(abs((tau_momentum_plus1 - tau1.P())/tau1.P()));} 
	//----------------------------------------------------------------------------------------------------------------------
  
	//##################################################  MINUS ##################################################################
	//############################################################################################################################
	//----------------------------------------------  absolute error -------------------------------------------------------
	if(plus == false){absolute_error_momentum_tau_minus->Fill(tau_momentum_minus1 - tau1.P());}
	//----------------------------------------------------------------------------------------------------------------------
  
	//----------------------------------------------  relative error -------------------------------------------------------
	if(plus == false){relative_error_momentum_tau_minus->Fill((tau_momentum_minus1 - tau1.P())/tau1.P());}
	//----------------------------------------------------------------------------------------------------------------------
  
	//----------------------------------------------  relative error abs---------------------------------------------------
	if(plus == false){relative_error_abs_momentum_tau_minus->Fill(abs((tau_momentum_minus1 - tau1.P())/tau1.P()));}
	//----------------------------------------------------------------------------------------------------------------------
  
	//##################################################  MEAN ###################################################################
	//############################################################################################################################
	//----------------------------------------------  absolute error -------------------------------------------------------
	absolute_error_momentum_tau_mean->Fill(tau_momentum_mean1 - tau1.P());
	//----------------------------------------------------------------------------------------------------------------------
  
	//----------------------------------------------  relative error -------------------------------------------------------
	relative_error_momentum_tau_mean->Fill((tau_momentum_mean1 - tau1.P())/tau1.P());
	//----------------------------------------------------------------------------------------------------------------------
  
	//----------------------------------------------  relative error abs---------------------------------------------------
	relative_error_abs_momentum_tau_mean->Fill(abs((tau_momentum_mean1 - tau1.P())/tau1.P()));
	//----------------------------------------------------------------------------------------------------------------------
  
  
	error_momentum_tau_mean_vs_momentum_a1->Fill(a1LV1.P(),abs((tau_momentum_mean1 - tau1.P())/tau1.P()));
	//-----------------------------------------------------------------------------
	
	
	//--------------- ratio +/-/average solutions ---------------------
	if (plus == true){ratio_plus->Fill(tau_momentum_plus1 / tau1.P());} 
	if (plus == false){ratio_minus->Fill(tau_momentum_minus1 / tau1.P());}
	ratio_mean->Fill(tau_momentum_mean1 / tau1.P());
	//----------------------------------------------------------------
	
	//---------------------- width ------------------------------------
	width_vs_a1_momentum->Fill(a1LV1.P(),tau_momentum_width1); 
	width_vs_a1_mass->Fill(a1LV1.M(),tau_momentum_width1);
	width_vs_a1_fixed_momentum->Fill(a1LV1.P(),sqrt(deltam1)/am1); // fixed a1 mass
	width_vs_a1_fixed_mass->Fill(a1LV1.M(),sqrt(deltap1)/ap1); // fixed a1 momentum
	
	if(theta_GJ1 >0.010 && theta_GJ1 <0.012){width_vs_a1_momentum_theta_1->Fill(a1LV1.P(),tau_momentum_width1);}
	if(theta_GJ1 >0.012 && theta_GJ1 <0.014){width_vs_a1_momentum_theta_2->Fill(a1LV1.P(),tau_momentum_width1);}
	if(theta_GJ1 >0.014 && theta_GJ1 <0.016){width_vs_a1_momentum_theta_3->Fill(a1LV1.P(),tau_momentum_width1);}
	if(theta_GJ1 >0.016 && theta_GJ1 <0.018){width_vs_a1_momentum_theta_4->Fill(a1LV1.P(),tau_momentum_width1);}
	if(theta_GJ1 >0.018 && theta_GJ1 <0.020){width_vs_a1_momentum_theta_5->Fill(a1LV1.P(),tau_momentum_width1);}
	
	//####################################################################
	//#######################  tau momentum ##############################
	//####################################################################
	
	
	//-------------------------- absolute error --------------------------------
	absolute_error_momentum_tau->Fill(tau_momentum_rec1 - tau1.P());
	//--------------------------------------------------------------------------
	
	//---------------------------  relative error ------------------------------
	relative_error_momentum_tau->Fill((tau_momentum_rec1 - tau1.P())/tau1.P());
	//---------------------------------------------------------------------------
	
	//---------------------------  relative error abs -----------------------------
	relative_error_abs_momentum_tau->Fill(abs((tau_momentum_rec1 - tau1.P())/tau1.P()));
	//-----------------------------------------------------------------------------
	
	//---------------------- momentum of tau vs momentum of a1 -----------------------
	Momentum_tau_vs_Momentum_a1->Fill(a1LV1.P(),tau_momentum_rec1);
	
	if(theta_GJ1 >0.010 && theta_GJ1 <0.012){Momentum_tau_vs_Momentum_a1_theta_1->Fill(a1LV1.P(),tau_momentum_rec1);}
	if(theta_GJ1 >0.012 && theta_GJ1 <0.014){Momentum_tau_vs_Momentum_a1_theta_2->Fill(a1LV1.P(),tau_momentum_rec1);}
	if(theta_GJ1 >0.014 && theta_GJ1 <0.016){Momentum_tau_vs_Momentum_a1_theta_3->Fill(a1LV1.P(),tau_momentum_rec1);}
	if(theta_GJ1 >0.016 && theta_GJ1 <0.018){Momentum_tau_vs_Momentum_a1_theta_4->Fill(a1LV1.P(),tau_momentum_rec1);}
	if(theta_GJ1 >0.018 && theta_GJ1 <0.020){Momentum_tau_vs_Momentum_a1_theta_5->Fill(a1LV1.P(),tau_momentum_rec1);}
	//---------------------------------------------------------------------
	
	//---------------------- momentum of tau vs mass of a1 -----------------------
	Momentum_tau_vs_Mass_a1->Fill(a1LV1.M(),tau_momentum_rec1);
	
	if(theta_GJ1 >0.010 && theta_GJ1 <0.012){Momentum_tau_vs_Mass_a1_theta_1->Fill(a1LV1.M(),tau_momentum_rec1);}
	if(theta_GJ1 >0.012 && theta_GJ1 <0.014){Momentum_tau_vs_Mass_a1_theta_2->Fill(a1LV1.M(),tau_momentum_rec1);}
	if(theta_GJ1 >0.014 && theta_GJ1 <0.016){Momentum_tau_vs_Mass_a1_theta_3->Fill(a1LV1.M(),tau_momentum_rec1);}
	if(theta_GJ1 >0.016 && theta_GJ1 <0.018){Momentum_tau_vs_Mass_a1_theta_4->Fill(a1LV1.M(),tau_momentum_rec1);}
	if(theta_GJ1 >0.018 && theta_GJ1 <0.020){Momentum_tau_vs_Mass_a1_theta_5->Fill(a1LV1.M(),tau_momentum_rec1);}
	//---------------------------------------------------------------------
	
	//-------------------- momentum of tau (gen) vs theta -----------------------
	Momentum_tau_vs_theta->Fill(theta_GJ1,tau_momentum_rec1);
	
	if(a1LV1.P() >19.8 && a1LV1.P() <20.2){Momentum_tau_vs_theta_pa1_20->Fill(theta_GJ1,tau_momentum_rec1);}
	if(a1LV1.P() >22.8 && a1LV1.P() <23.2){Momentum_tau_vs_theta_pa1_23->Fill(theta_GJ1,tau_momentum_rec1);}
	if(a1LV1.P() >24.8 && a1LV1.P() <25.2){Momentum_tau_vs_theta_pa1_25->Fill(theta_GJ1,tau_momentum_rec1);}
	if(a1LV1.P() >26.8 && a1LV1.P() <27.2){Momentum_tau_vs_theta_pa1_27->Fill(theta_GJ1,tau_momentum_rec1);}
	if(a1LV1.P() >29.8 && a1LV1.P() <30.2){Momentum_tau_vs_theta_pa1_30->Fill(theta_GJ1,tau_momentum_rec1);}
	if(a1LV1.P() >34.8 && a1LV1.P() <35.2){Momentum_tau_vs_theta_pa1_35->Fill(theta_GJ1,tau_momentum_rec1);}
	if(a1LV1.P() >39.8 && a1LV1.P() <40.2){Momentum_tau_vs_theta_pa1_40->Fill(theta_GJ1,tau_momentum_rec1);}
	//---------------------------------------------------------------------
	
	//-------------------- momentum of tau (rec) vs theta ----------------------
	if(a1LV1.P() >19.6 && a1LV1.P() <20.4){
	  Momentum_tau_rec_vs_theta_pa1_20->Fill(theta_GJ1,tau_momentum_plus1);
	  Momentum_tau_rec_vs_theta_pa1_20->Fill(theta_GJ1,tau_momentum_minus1);
	  theta_CM_vs_tau_momentum_rec_pa1_20->Fill(theta_GJ1,tau_momentum_plus1);
	  theta_CM_vs_tau_momentum_rec_pa1_20->Fill(theta_GJ1,tau_momentum_minus1);
	  costheta_CM_vs_tau_momentum_rec_pa1_20->Fill(theta_GJ1,tau_momentum_plus1);
	  costheta_CM_vs_tau_momentum_rec_pa1_20->Fill(theta_GJ1,tau_momentum_minus1);}
	if(a1LV1.P() >22.6 && a1LV1.P() <23.4){
	  Momentum_tau_rec_vs_theta_pa1_23->Fill(theta_GJ1,tau_momentum_plus1);
	  Momentum_tau_rec_vs_theta_pa1_23->Fill(theta_GJ1,tau_momentum_minus1);
	  theta_CM_vs_tau_momentum_rec_pa1_23->Fill(theta_GJ1,tau_momentum_plus1);
	  theta_CM_vs_tau_momentum_rec_pa1_23->Fill(theta_GJ1,tau_momentum_minus1);
	  costheta_CM_vs_tau_momentum_rec_pa1_23->Fill(theta_GJ1,tau_momentum_plus1);
	  costheta_CM_vs_tau_momentum_rec_pa1_23->Fill(theta_GJ1,tau_momentum_minus1);}
	if(a1LV1.P() >24.6 && a1LV1.P() <25.4){
	  Momentum_tau_rec_vs_theta_pa1_25->Fill(theta_GJ1,tau_momentum_plus1);
	  Momentum_tau_rec_vs_theta_pa1_25->Fill(theta_GJ1,tau_momentum_minus1);
	  theta_CM_vs_tau_momentum_rec_pa1_25->Fill(theta_GJ1,tau_momentum_plus1);
	  theta_CM_vs_tau_momentum_rec_pa1_25->Fill(theta_GJ1,tau_momentum_minus1);
	  costheta_CM_vs_tau_momentum_rec_pa1_25->Fill(theta_GJ1,tau_momentum_plus1);
	  costheta_CM_vs_tau_momentum_rec_pa1_25->Fill(theta_GJ1,tau_momentum_minus1);}
	if(a1LV1.P() >26.6 && a1LV1.P() <27.4){
	  Momentum_tau_rec_vs_theta_pa1_27->Fill(theta_GJ1,tau_momentum_plus1);
	  Momentum_tau_rec_vs_theta_pa1_27->Fill(theta_GJ1,tau_momentum_minus1);
	  theta_CM_vs_tau_momentum_rec_pa1_27->Fill(theta_GJ1,tau_momentum_plus1);
	  theta_CM_vs_tau_momentum_rec_pa1_27->Fill(theta_GJ1,tau_momentum_minus1);
	  costheta_CM_vs_tau_momentum_rec_pa1_27->Fill(theta_GJ1,tau_momentum_plus1);
	  costheta_CM_vs_tau_momentum_rec_pa1_27->Fill(theta_GJ1,tau_momentum_minus1);}
	if(a1LV1.P() >29.6 && a1LV1.P() <30.4){
	  Momentum_tau_rec_vs_theta_pa1_30->Fill(theta_GJ1,tau_momentum_plus1);
	  Momentum_tau_rec_vs_theta_pa1_30->Fill(theta_GJ1,tau_momentum_minus1);
	  Momentum_tau_plus_rec_vs_theta_pa1_30->Fill(theta_GJ1,tau_momentum_plus1);
	  Momentum_tau_minus_rec_vs_theta_pa1_30->Fill(theta_GJ1,tau_momentum_minus1);
	  theta_CM_vs_tau_momentum_rec_pa1_30->Fill(theta_GJ1,tau_momentum_plus1);
	  theta_CM_vs_tau_momentum_rec_pa1_30->Fill(theta_GJ1,tau_momentum_minus1);
	  costheta_CM_vs_tau_momentum_rec_pa1_30->Fill(theta_GJ1,tau_momentum_plus1);
	  costheta_CM_vs_tau_momentum_rec_pa1_30->Fill(theta_GJ1,tau_momentum_minus1);
	  //std::cout << "theta max 30 " << theta_max << std::endl;
	  if(theta_GJ1 >= theta_max1){Mass_a1_above_theta_max_plus_pa1_30->Fill(theta_GJ1,a1LV1.M());}
	  if(theta_GJ1 < theta_max1){Mass_a1_above_theta_max_minus_pa1_30->Fill(theta_GJ1,a1LV1.M());}}
	if(a1LV1.P() >34.6 && a1LV1.P() <35.4){
	  Momentum_tau_rec_vs_theta_pa1_35->Fill(theta_GJ1,tau_momentum_plus1);
	  Momentum_tau_rec_vs_theta_pa1_35->Fill(theta_GJ1,tau_momentum_minus1);
	  Momentum_tau_plus_rec_vs_theta_pa1_35->Fill(theta_GJ1,tau_momentum_plus1);
	  Momentum_tau_minus_rec_vs_theta_pa1_35->Fill(theta_GJ1,tau_momentum_minus1);
	  theta_CM_vs_tau_momentum_rec_pa1_35->Fill(theta_GJ1,tau_momentum_plus1);
	  theta_CM_vs_tau_momentum_rec_pa1_35->Fill(theta_GJ1,tau_momentum_minus1);
	  costheta_CM_vs_tau_momentum_rec_pa1_35->Fill(theta_GJ1,tau_momentum_plus1);
	  costheta_CM_vs_tau_momentum_rec_pa1_35->Fill(theta_GJ1,tau_momentum_minus1);
	  //std::cout << "theta max 35 " << theta_max << std::endl;
	  if(theta_GJ1 >= theta_max1){Mass_a1_above_theta_max_plus_pa1_35->Fill(theta_GJ1,a1LV1.M());}
	  if(theta_GJ1 < theta_max1){Mass_a1_above_theta_max_minus_pa1_35->Fill(theta_GJ1,a1LV1.M());}}
	if(a1LV1.P() >39.6 && a1LV1.P() <40.4){
	  Momentum_tau_rec_vs_theta_pa1_40->Fill(theta_GJ1,tau_momentum_plus1);
	  Momentum_tau_rec_vs_theta_pa1_40->Fill(theta_GJ1,tau_momentum_minus1);
	  Momentum_tau_plus_rec_vs_theta_pa1_40->Fill(theta_GJ1,tau_momentum_plus1);
	  Momentum_tau_minus_rec_vs_theta_pa1_40->Fill(theta_GJ1,tau_momentum_minus1);
	  theta_CM_vs_tau_momentum_rec_pa1_40->Fill(theta_GJ1,tau_momentum_plus1);
	  theta_CM_vs_tau_momentum_rec_pa1_40->Fill(theta_GJ1,tau_momentum_minus1);
	  costheta_CM_vs_tau_momentum_rec_pa1_40->Fill(theta_GJ1,tau_momentum_plus1);
	  costheta_CM_vs_tau_momentum_rec_pa1_40->Fill(theta_GJ1,tau_momentum_minus1);
	  //std::cout << "theta max 40 " << theta_max << std::endl;
	  if(theta_GJ1 >= theta_max1){Mass_a1_above_theta_max_plus_pa1_40->Fill(theta_GJ1,a1LV1.M());}
	  if(theta_GJ1 < theta_max1){Mass_a1_above_theta_max_minus_pa1_40->Fill(theta_GJ1,a1LV1.M());}}
	//---------------------------------------------------------------------
	
	//-------------------- momentum of tau (rec) with fixed mass vs theta ----------------------
	if(a1LV1.P() >19.6 && a1LV1.P() <20.4){
	  Momentum_tau_rec_fixed_mass_vs_theta_pa1_20->Fill(theta_GJ1,(-bm1+sqrt(deltam1)) / (2*am1));
	  Momentum_tau_rec_fixed_mass_vs_theta_pa1_20->Fill(theta_GJ1,(-bm1-sqrt(deltam1)) / (2*am1));}
	if(a1LV1.P() >22.6 && a1LV1.P() <23.4){
	  Momentum_tau_rec_fixed_mass_vs_theta_pa1_23->Fill(theta_GJ1,(-bm1+sqrt(deltam1)) / (2*am1));
	  Momentum_tau_rec_fixed_mass_vs_theta_pa1_23->Fill(theta_GJ1,(-bm1-sqrt(deltam1)) / (2*am1));}
	if(a1LV1.P() >24.6 && a1LV1.P() <25.4){
	  Momentum_tau_rec_fixed_mass_vs_theta_pa1_25->Fill(theta_GJ1,(-bm1+sqrt(deltam1)) / (2*am1));
	  Momentum_tau_rec_fixed_mass_vs_theta_pa1_25->Fill(theta_GJ1,(-bm1-sqrt(deltam1)) / (2*am1));}
	if(a1LV1.P() >26.6 && a1LV1.P() <27.4){
	  Momentum_tau_rec_fixed_mass_vs_theta_pa1_27->Fill(theta_GJ1,(-bm1+sqrt(deltam1)) / (2*am1));
	  Momentum_tau_rec_fixed_mass_vs_theta_pa1_27->Fill(theta_GJ1,(-bm1-sqrt(deltam1)) / (2*am1));}
	if(a1LV1.P() >29.6 && a1LV1.P() <30.4){
	  Momentum_tau_rec_fixed_mass_vs_theta_pa1_30->Fill(theta_GJ1,(-bm1+sqrt(deltam1)) / (2*am1));
	  Momentum_tau_rec_fixed_mass_vs_theta_pa1_30->Fill(theta_GJ1,(-bm1-sqrt(deltam1)) / (2*am1));}
	if(a1LV1.P() >34.6 && a1LV1.P() <35.4){
	  Momentum_tau_rec_fixed_mass_vs_theta_pa1_35->Fill(theta_GJ1,(-bm1+sqrt(deltam1)) / (2*am1));
	  Momentum_tau_rec_fixed_mass_vs_theta_pa1_35->Fill(theta_GJ1,(-bm1-sqrt(deltam1)) / (2*am1));}
	if(a1LV1.P() >39.6 && a1LV1.P() <40.4){
	  Momentum_tau_rec_fixed_mass_vs_theta_pa1_40->Fill(theta_GJ1,(-bm1+sqrt(deltam1)) / (2*am1));
	  Momentum_tau_rec_fixed_mass_vs_theta_pa1_40->Fill(theta_GJ1,(-bm1-sqrt(deltam1)) / (2*am1));}
	//------------------------------------------------------------------------------------------
	
	//-------------------- momentum of tau (rec) vs theta MASS CONDITIONS ----------------------
	if(a1LV1.P() >19.6 && a1LV1.P() <20.4 && (a1LV1.E()/a1LV1.M())<(tau1.E()/tau1.M())){
	  Momentum_tau_rec_vs_theta_pa1_20_small_ma1->Fill(theta_GJ1,tau_momentum_plus1);
	  Momentum_tau_rec_vs_theta_pa1_20_small_ma1->Fill(theta_GJ1,tau_momentum_minus1);}
	if(a1LV1.P() >19.6 && a1LV1.P() <20.4 && (a1LV1.E()/a1LV1.M())>(tau1.E()/tau1.M())){
	  Momentum_tau_rec_vs_theta_pa1_20_large_ma1->Fill(theta_GJ1,tau_momentum_plus1);
	  Momentum_tau_rec_vs_theta_pa1_20_large_ma1->Fill(theta_GJ1,tau_momentum_minus1);} 
	  
	if(a1LV1.P() >22.6 && a1LV1.P() <23.4 && (a1LV1.E()/a1LV1.M())<(tau1.E()/tau1.M())){
	  Momentum_tau_rec_vs_theta_pa1_23_small_ma1->Fill(theta_GJ1,tau_momentum_plus1);
	  Momentum_tau_rec_vs_theta_pa1_23_small_ma1->Fill(theta_GJ1,tau_momentum_minus1);}
	if(a1LV1.P() >22.6 && a1LV1.P() <23.4 && (a1LV1.E()/a1LV1.M())>(tau1.E()/tau1.M())){
	  Momentum_tau_rec_vs_theta_pa1_23_large_ma1->Fill(theta_GJ1,tau_momentum_plus1);
	  Momentum_tau_rec_vs_theta_pa1_23_large_ma1->Fill(theta_GJ1,tau_momentum_minus1);}
	  
	if(a1LV1.P() >24.6 && a1LV1.P() <25.4 && (a1LV1.E()/a1LV1.M())<(tau1.E()/tau1.M())){
	  Momentum_tau_rec_vs_theta_pa1_25_small_ma1->Fill(theta_GJ1,tau_momentum_plus1);
	  Momentum_tau_rec_vs_theta_pa1_25_small_ma1->Fill(theta_GJ1,tau_momentum_minus1);}
	if(a1LV1.P() >24.6 && a1LV1.P() <25.4 && (a1LV1.E()/a1LV1.M())>(tau1.E()/tau1.M())){
	  Momentum_tau_rec_vs_theta_pa1_25_large_ma1->Fill(theta_GJ1,tau_momentum_plus1);
	  Momentum_tau_rec_vs_theta_pa1_25_large_ma1->Fill(theta_GJ1,tau_momentum_minus1);}
	  
	if(a1LV1.P() >26.6 && a1LV1.P() <27.4 && (a1LV1.E()/a1LV1.M())<(tau1.E()/tau1.M())){
	  Momentum_tau_rec_vs_theta_pa1_27_small_ma1->Fill(theta_GJ1,tau_momentum_plus1);
	  Momentum_tau_rec_vs_theta_pa1_27_small_ma1->Fill(theta_GJ1,tau_momentum_minus1);}
	if(a1LV1.P() >26.6 && a1LV1.P() <27.4 && (a1LV1.E()/a1LV1.M())>(tau1.E()/tau1.M())){
	  Momentum_tau_rec_vs_theta_pa1_27_large_ma1->Fill(theta_GJ1,tau_momentum_plus1);
	  Momentum_tau_rec_vs_theta_pa1_27_large_ma1->Fill(theta_GJ1,tau_momentum_minus1);}
	  
	if(a1LV1.P() >29.6 && a1LV1.P() <30.4 && (a1LV1.E()/a1LV1.M())<(tau1.E()/tau1.M())){
	  Momentum_tau_rec_vs_theta_pa1_30_small_ma1->Fill(theta_GJ1,tau_momentum_plus1);
	  Momentum_tau_rec_vs_theta_pa1_30_small_ma1->Fill(theta_GJ1,tau_momentum_minus1);}
	if(a1LV1.P() >29.6 && a1LV1.P() <30.4 && (a1LV1.E()/a1LV1.M())>(tau1.E()/tau1.M())){
	  Momentum_tau_rec_vs_theta_pa1_30_large_ma1->Fill(theta_GJ1,tau_momentum_plus1);
	  Momentum_tau_rec_vs_theta_pa1_30_large_ma1->Fill(theta_GJ1,tau_momentum_minus1);}  
	  
	  
	if(a1LV1.P() >34.6 && a1LV1.P() <35.4 && (a1LV1.E()/a1LV1.M())<(tau1.E()/tau1.M())){
	  Momentum_tau_rec_vs_theta_pa1_35_small_ma1->Fill(theta_GJ1,tau_momentum_plus1);
	  Momentum_tau_rec_vs_theta_pa1_35_small_ma1->Fill(theta_GJ1,tau_momentum_minus1);}
	if(a1LV1.P() >34.6 && a1LV1.P() <35.4 && (a1LV1.E()/a1LV1.M())>(tau1.E()/tau1.M())){
	  Momentum_tau_rec_vs_theta_pa1_35_large_ma1->Fill(theta_GJ1,tau_momentum_plus1);
	  Momentum_tau_rec_vs_theta_pa1_35_large_ma1->Fill(theta_GJ1,tau_momentum_minus1);}
	  
	if(a1LV1.P() >39.6 && a1LV1.P() <40.4 && (a1LV1.E()/a1LV1.M())<(tau1.E()/tau1.M())){
	  Momentum_tau_rec_vs_theta_pa1_40_small_ma1->Fill(theta_GJ1,tau_momentum_plus1);
	  Momentum_tau_rec_vs_theta_pa1_40_small_ma1->Fill(theta_GJ1,tau_momentum_minus1);}
	if(a1LV1.P() >39.6 && a1LV1.P() <40.4 && (a1LV1.E()/a1LV1.M())>(tau1.E()/tau1.M())){
	  Momentum_tau_rec_vs_theta_pa1_40_large_ma1->Fill(theta_GJ1,tau_momentum_plus1);
	  Momentum_tau_rec_vs_theta_pa1_40_large_ma1->Fill(theta_GJ1,tau_momentum_minus1);}
	//------------------------------------------------------------------------------------
	
	/*
	//-------------------------- Crossing ---------------------------------------------
	//--- pa1 = 30 GeV ---------------
	if(theta_GJ>0 && theta_GJ<0.019 && a1LV.P() >29.6 && a1LV.P() <30.4){crossing_pa1_30_above_theta_mass->Fill(a1LV.M(),tau1.M());}
	if(theta_GJ>0.019 && theta_GJ<0.021 && a1LV.P() >29.6 && a1LV.P() <30.4){crossing_pa1_30_near_theta_mass->Fill(a1LV.M(),tau1.M());}
	if(theta_GJ>0.021 && a1LV.P() >29.6 && a1LV.P() <30.4){crossing_pa1_30_below_theta_mass->Fill(a1LV.M(),tau1.M());}
	
	if(theta_GJ>0 && theta_GJ<0.019 && a1LV.P() >29.6 && a1LV.P() <30.4){crossing_pa1_30_above_theta_momentum->Fill(a1LV.M(),a1LV.P());}
	if(theta_GJ>0.019 && theta_GJ<0.021 && a1LV.P() >29.6 && a1LV.P() <30.4){crossing_pa1_30_near_theta_momentum->Fill(a1LV.M(),a1LV.P());}
	if(theta_GJ>0.021 && a1LV.P() >29.6 && a1LV.P() <30.4){crossing_pa1_30_below_theta_momentum->Fill(a1LV.M(),a1LV.P());}
	
	if(theta_GJ>0 && theta_GJ<0.019 && a1LV.P() >29.6 && a1LV.P() <30.4){crossing_pa1_30_above_theta_energy->Fill(a1LV.E(),tau1.E());}
	if(theta_GJ>0.019 && theta_GJ<0.021 && a1LV.P() >29.6 && a1LV.P() <30.4){crossing_pa1_30_near_theta_energy->Fill(a1LV.E(),tau1.E());}
	if(theta_GJ>0.021 && a1LV.P() >29.6 && a1LV.P() <30.4){crossing_pa1_30_below_theta_energy->Fill(a1LV.E(),tau1.E());}
	
	//------ pa1 = 35 GeV ----------------
	if(theta_GJ>0 && theta_GJ<0.011 && a1LV.P() >34.6 && a1LV.P() <35.4){crossing_pa1_35_above_theta_mass->Fill(a1LV.M(),tau1.M());}
	if(theta_GJ>0.011 && theta_GJ<0.013 && a1LV.P() >34.6 && a1LV.P() <35.4){crossing_pa1_35_near_theta_mass->Fill(a1LV.M(),tau1.M());}
	if(theta_GJ>0.013 && a1LV.P() >34.6 && a1LV.P() <35.4){crossing_pa1_35_below_theta_mass->Fill(a1LV.M(),tau1.M());}
	
	if(theta_GJ>0 && theta_GJ<0.011 && a1LV.P() >34.6 && a1LV.P() <35.4){crossing_pa1_35_above_theta_momentum->Fill(a1LV.M(),a1LV.P());}
	if(theta_GJ>0.011 && theta_GJ<0.013 && a1LV.P() >34.6 && a1LV.P() <35.4){crossing_pa1_35_near_theta_momentum->Fill(a1LV.M(),a1LV.P());}
	if(theta_GJ>0.013 && a1LV.P() >34.6 && a1LV.P() <35.4){crossing_pa1_35_below_theta_momentum->Fill(a1LV.M(),a1LV.P());}
	
	if(theta_GJ>0 && theta_GJ<0.011 && a1LV.P() >34.6 && a1LV.P() <35.4){crossing_pa1_35_above_theta_energy->Fill(a1LV.E(),tau1.E());}
	if(theta_GJ>0.011 && theta_GJ<0.013 && a1LV.P() >34.6 && a1LV.P() <35.4){crossing_pa1_35_near_theta_energy->Fill(a1LV.E(),tau1.E());}
	if(theta_GJ>0.013 && a1LV.P() >34.6 && a1LV.P() <35.4){crossing_pa1_35_below_theta_energy->Fill(a1LV.E(),tau1.E());}
	
	//---------- pa1 = 40 GeV -------------
	if(theta_GJ>0 && theta_GJ<0.0045 && a1LV.P() >39.6 && a1LV.P() <40.4){crossing_pa1_40_above_theta_mass->Fill(a1LV.M(),tau1.M());}
	if(theta_GJ>0.0045 && theta_GJ<0.0055 && a1LV.P() >39.6 && a1LV.P() <40.4){crossing_pa1_40_near_theta_mass->Fill(a1LV.M(),tau1.M());}
	if(theta_GJ>0.0055 && a1LV.P() >39.6 && a1LV.P() <40.4){crossing_pa1_40_below_theta_mass->Fill(a1LV.M(),tau1.M());}
	
	if(theta_GJ>0 && theta_GJ<0.0045 && a1LV.P() >39.6 && a1LV.P() <40.4){crossing_pa1_40_above_theta_momentum->Fill(a1LV.M(),a1LV.P());}
	if(theta_GJ>0.0045 && theta_GJ<0.0055 && a1LV.P() >39.6 && a1LV.P() <40.4){crossing_pa1_40_near_theta_momentum->Fill(a1LV.M(),a1LV.P());}
	if(theta_GJ>0.0055 && a1LV.P() >39.6 && a1LV.P() <40.4){crossing_pa1_40_below_theta_momentum->Fill(a1LV.M(),a1LV.P());}
	
	if(theta_GJ>0 && theta_GJ<0.0045 && a1LV.P() >39.6 && a1LV.P() <40.4){crossing_pa1_40_above_theta_energy->Fill(a1LV.E(),tau1.E());}
	if(theta_GJ>0.0045 && theta_GJ<0.0055 && a1LV.P() >39.6 && a1LV.P() <40.4){crossing_pa1_40_near_theta_energy->Fill(a1LV.E(),tau1.E());}
	if(theta_GJ>0.0055 && a1LV.P() >39.6 && a1LV.P() <40.4){crossing_pa1_40_below_theta_energy->Fill(a1LV.E(),tau1.E());}
	//-------------------------------------------------------------------------------------
	*/
	
	
	
	
	if( abs(delta1 - 0)<0.4){
	  Theta_vs_momentum_a1_near_crossing->Fill(a1LV1.P(),theta_GJ1);
	  Momentum_a1_vs_mass_a1_near_crossing->Fill(a1LV1.M(),a1LV1.P());}
	
	
	
	if(a1LV1.P() >29.6 && a1LV1.P() <30.4){
	  crossing_pa1_30_theta_mass->Fill(theta_GJ1,tau1.M()/a1LV1.M());
	  crossing_pa1_30_theta_momentum->Fill(theta_GJ1,1.230/a1LV1.M());
	  crossing_pa1_30_theta_energy->Fill(theta_GJ1,(tau1.E()/tau1.M())/(a1LV1.E()/a1LV1.M()));
	  Momentum_tau_pa1_30_without_sqrt->Fill(theta_GJ1,-b1/(2*a1));
	  Momentum_tau_pa1_30_only_sqrt->Fill(theta_GJ1,delta1/(2*a1));
	}
	if(a1LV1.P() >34.6 && a1LV1.P() <35.4){
	  crossing_pa1_35_theta_mass->Fill(theta_GJ1,tau1.M()/a1LV1.M());
	  crossing_pa1_35_theta_momentum->Fill(theta_GJ1,1.230/a1LV1.M());
	  crossing_pa1_35_theta_energy->Fill(theta_GJ1,(tau1.E()/tau1.M())/(a1LV1.E()/a1LV1.M()));
	  Momentum_tau_pa1_35_without_sqrt->Fill(theta_GJ1,-b1/(2*a1));
	  Momentum_tau_pa1_35_only_sqrt->Fill(theta_GJ1,delta1/(2*a1));
	}
	if(a1LV1.P() >39.6 && a1LV1.P() <40.4){
	  crossing_pa1_40_theta_mass->Fill(theta_GJ1,tau1.M()/a1LV1.M());
	  crossing_pa1_40_theta_momentum->Fill(theta_GJ1,1.230/a1LV1.M());
	  crossing_pa1_40_theta_energy->Fill(theta_GJ1,(tau1.E()/tau1.M())/(a1LV1.E()/a1LV1.M()));
	  Momentum_tau_pa1_40_without_sqrt->Fill(theta_GJ1,-b1/(2*a1));
	  Momentum_tau_pa1_40_only_sqrt->Fill(theta_GJ1,delta1/(2*a1));
	}
	//------------------------------------------------------------------------------------
	
	//---------------------------- 3D -----------------------------------------
	Momentum_tau_vs_theta_vs_Momentum_a1->Fill(theta_GJ1,tau_momentum_rec1,a1LV1.P());
	//-------------------------------------------------------------------------
	
	//------------------- +/- solutions -----------------------------
	Plus_Minus_solutions->Fill(nb_tau_momentum_minus1);
	//---------------------------------------------------------------
	
	
	//-------------------------------------- Momentum tau gen vs momentum tau mean ------------------------------------
	Momentum_tau_gen_vs_momentum_tau_mean->Fill(tau_momentum_mean1,tau1.P());
	//-----------------------------------------------------------------------------------------------------------------------
	
	
	//------------------------------------------- Momentum tau gen vs resolution ------------------------------------
	Momentum_tau_gen_vs_resolution->Fill((tau_momentum_mean1 - tau1.P())/tau1.P(),tau1.P()); 
	//-----------------------------------------------------------------------------------------------------------------------
  
	//------------------------------------------- Momentum tau rec vs resolution ------------------------------------
	Momentum_tau_rec_vs_resolution->Fill((tau_momentum_mean1 - tau1.P())/tau1.P(),tau_momentum_rec1);
  	//-----------------------------------------------------------------------------------------------------------------------
  
	//------------------------------------------- Momentum tau mean rec vs resolution ------------------------------------
	Momentum_tau_mean_vs_resolution->Fill((tau_momentum_mean1 - tau1.P())/tau1.P(),tau_momentum_mean1);
	//-----------------------------------------------------------------------------------------------------------------------
  
	
	
	
	
	
	/*
	//------------------------- omega (helicity)-------------------------------------
	Omega_plus_a1->Fill(a1LV.getOmega(),HelWeightPlus);
	Omega_minus_a1->Fill(a1LV.getOmega(),HelWeightMinus);
	//----------------------------------------------------------------------------------
	*/
	
	
	//if(theta_GJ > theta_max -0.0001 && theta_GJ <= theta_max){ ambiguity_a1_momentum->Fill(a1LV.P());}

	
	
	
     }
     
    
     /*
     if (theta_GJ > theta_max){ 
       ++unphysical_theta; 
       std::cout<< "unphysical events " << unphysical_theta <<std::endl;
     }
     */
    //--------------------------------  a1
     /*
std::cout << "nb solu + "<< nb_tau_momentum_plus << std::endl;
    std::cout << "nb solu - "<< nb_tau_momentum_minus << std::endl;
*/

 





/*
    if(JAK1==5 && SubJAK1==51){ //a1^- and 3 pions 
       vector<TLorentzVector> particles;
       particles.clear();
       SortPions(A1Pions1);
       int taucharge =  (A1Pions1.at(0).pdg_id()+A1Pions1.at(1).pdg_id()+A1Pions1.at(2).pdg_id() > 0) ? 1 : -1;
       taucharge1=taucharge;

       a1ss1pi.SetPxPyPzE(A1Pions1.at(0).momentum().px(), A1Pions1.at(0).momentum().py(), A1Pions1.at(0).momentum().pz(), A1Pions1.at(0).momentum().e());
       a1ss2pi.SetPxPyPzE(A1Pions1.at(1).momentum().px(), A1Pions1.at(1).momentum().py(), A1Pions1.at(1).momentum().pz(), A1Pions1.at(1).momentum().e());
       a1ospi.SetPxPyPzE(A1Pions1.at(2).momentum().px(), A1Pions1.at(2).momentum().py(), A1Pions1.at(2).momentum().pz(), A1Pions1.at(2).momentum().e());
       particles.push_back(tau1);
       particles.push_back(a1ospi);
       particles.push_back(a1ss1pi);
       particles.push_back(a1ss2pi);
       
       
       TLorentzVector Pion1(0,0,0,0);
       TLorentzVector Pion2(0,0,0,0);
       TLorentzVector Pion3(0,0,0,0);
	Pion1+=TLorentzVector(A1Pions1.at(0).momentum().px(), A1Pions1.at(0).momentum().py(), A1Pions1.at(0).momentum().pz(), A1Pions1.at(0).momentum().e());
	Pion2+=TLorentzVector(A1Pions1.at(1).momentum().px(), A1Pions1.at(1).momentum().py(), A1Pions1.at(1).momentum().pz(), A1Pions1.at(1).momentum().e());
	Pion3+=TLorentzVector(A1Pions1.at(2).momentum().px(), A1Pions1.at(2).momentum().py(), A1Pions1.at(2).momentum().pz(), A1Pions1.at(2).momentum().e());
       
       

       Polarimetr1.Configure(particles, taucharge);
       A1.Configure(particles,"a1",taucharge);

       omega_a1p_minus->Fill(Polarimetr1.getOmegaA1Bar(),HelWeightMinus);  omega_a1p_plus->Fill(Polarimetr1.getOmegaA1Bar(),HelWeightPlus);
       //  std::cout<<"omegabar   "<<Polarimetr1.getOmegaA1Bar() <<std::endl;
       a1hh.Configure(particles, a1ospi+a1ss1pi+a1ss2pi);
       TauPolA1.Configure(particles,"a1");
       tauandprod1=particles;
       //      omega_a1_minus->Fill(a1h.getA1omega(),HelWeightMinus);                                                          omegabar_a1_plus->Fill(a1h.getA1omega(),HelWeightPlus);

       omegabar_a1_minus->Fill(a1hh.getA1omega(),HelWeightMinus);                                                          omegabar_a1_plus->Fill(a1hh.getA1omega(),HelWeightPlus);

      

       
    }
*/




     //--------------------  tau+
     if(JAK2==2){
       vector<TLorentzVector> tauandprodmu2;
       tauandprodmu2.push_back(TLorentzVector(SecondTau->momentum().px(), SecondTau->momentum().py(), SecondTau->momentum().pz(), SecondTau->momentum().e()));
       for(std::vector<HepMC::GenParticle>::const_iterator a = SecondTauProducts.begin(); a!=SecondTauProducts.end(); ++a){
     	if(abs(a->pdg_id())==13){tauandprodmu2.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );} }
       tauandprod2=tauandprodmu2;
       Mu2.Configure(tauandprod2,"lepton");
     }
     
     if(JAK2==3){
       vector<TLorentzVector> tauandprod;
       tauandprod.push_back(TLorentzVector(SecondTau->momentum().px(), SecondTau->momentum().py(), SecondTau->momentum().pz(), SecondTau->momentum().e()));
       for(std::vector<HepMC::GenParticle>::const_iterator a = SecondTauProducts.begin(); a!=SecondTauProducts.end(); ++a){
	 if(abs(a->pdg_id())==211){tauandprod.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );} }
       tauandprod2 = tauandprod;
       Pi2.Configure(tauandprod2,"pion");
       if(ApplyCut){
	 if(passed1){
	   pi_plus->Fill(Pi2.getOmega(),HelWeightPlus);
	   pi_minus->Fill(Pi2.getOmega(),HelWeightMinus);
	 }
       }else{
	 pi_plus->Fill(Pi2.getOmega(),HelWeightPlus);
	 pi_minus->Fill(Pi2.getOmega(),HelWeightMinus);
       }
     }
     
     if(JAK2==4){
       vector<TLorentzVector> tauandprod;
       tauandprod.push_back(TLorentzVector(SecondTau->momentum().px(), SecondTau->momentum().py(), SecondTau->momentum().pz(), SecondTau->momentum().e()));
       //      std::cout<<" ---- "<<std::endl;
       for(std::vector<HepMC::GenParticle>::const_iterator a =SecondTauProducts.begin(); a!=SecondTauProducts.end(); ++a){
	 //	std::cout<<" rho decays:   "<< a->pdg_id() << std::endl;
	 if(abs(a->pdg_id())==211){tauandprod.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );}
	 if(abs(a->pdg_id())==111){tauandprod.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );}
       }
       tauandprod2=tauandprod;  
       Rho2.Configure(tauandprod2,"rho");


     }

    
    if(JAK2==5 && SubJAK2==51){
      vector<TLorentzVector> particles;
      particles.clear();
      SortPions(A1Pions2);
      
      int taucharge =  (A1Pions2.at(0).pdg_id()+A1Pions2.at(1).pdg_id()+A1Pions2.at(2).pdg_id() > 0) ? 1 : -1;
      taucharge2=taucharge;
      a1ss1pi.SetPxPyPzE(A1Pions2.at(0).momentum().px(), A1Pions2.at(0).momentum().py(), A1Pions2.at(0).momentum().pz(), A1Pions2.at(0).momentum().e());
      a1ss2pi.SetPxPyPzE(A1Pions2.at(1).momentum().px(), A1Pions2.at(1).momentum().py(), A1Pions2.at(1).momentum().pz(), A1Pions2.at(1).momentum().e());
      a1ospi.SetPxPyPzE(A1Pions2.at(2).momentum().px(), A1Pions2.at(2).momentum().py(), A1Pions2.at(2).momentum().pz(), A1Pions2.at(2).momentum().e());
      
      particles.push_back(tau2);
      particles.push_back(a1ospi);
      particles.push_back(a1ss1pi);
      particles.push_back(a1ss2pi);
      a1h.Configure(particles, a1ospi+a1ss1pi+a1ss2pi);
      TauPolA1.Configure(particles,"a1");
      tauandprod2=particles;
      Polarimetr.Configure(particles, taucharge);
      A2.Configure(particles,"a1", taucharge);
      
      if(ApplyCut){
	if(passed2){
	  
	  TRFomegabar_a1_minus->Fill(a1h.TRF_vgetA1omega(),HelWeightMinus);                                      TRFomegabar_a1_plus->Fill(a1h.TRF_vgetA1omega(),HelWeightPlus);
	  TRFomegabar_a1scalar_minus->Fill(a1h.TRF_vgetA1omega("scalar"),HelWeightMinus);                  TRFomegabar_a1scalar_plus->Fill(a1h.TRF_vgetA1omega("scalar"),HelWeightPlus);
	  cosbetacostheta_minus->Fill(a1h.cosbeta(),a1h.costhetaLF(),HelWeightMinus);                              cosbetacostheta_plus->Fill(a1h.cosbeta(),a1h.costhetaLF(),HelWeightPlus);
	  TRFcosbetacostheta_minus->Fill(a1h.TRF_cosbeta(),a1h.costhetaLF(),HelWeightMinus);                  TRFcosbetacostheta_plus->Fill(a1h.TRF_cosbeta(),a1h.costhetaLF(),HelWeightPlus);
	  
	  s1->Fill((a1ospi +a1ss1pi).M() );
	  s2->Fill((a1ospi +a1ss2pi).M() );
	  qq->Fill((a1ospi+a1ss1pi+a1ss2pi).M());

	  omega_a1_minus->Fill(Polarimetr.getOmegaA1Bar(),HelWeightMinus);  omega_a1_plus->Fill(Polarimetr.getOmegaA1Bar(),HelWeightPlus);
	  hmag->Fill(Polarimetr.PVC().Vect().Mag());
	  s1s2->Fill((a1ospi+a1ss2pi).M2(),(a1ospi+a1ss1pi).M2());
	  int hel(0);
	  if(HelWeightMinus==1)     hel = -1;
	  if(HelWeightPlus==1)      hel =  1;
	  if(HelWeightPlus==1)      a1_mcos2gamma_plus->Fill(a1h.getMoment(a1h.costhetaLF(),"c2g", 1));                 
	  if(HelWeightMinus==1)     a1_mcos2gamma_minus->Fill(a1h.getMoment(a1h.costhetaLF(),"c2g", -1));                 
	}
      }else{
	TRFomegabar_a1_minus->Fill(a1h.TRF_vgetA1omega(),HelWeightMinus);                                      TRFomegabar_a1_plus->Fill(a1h.TRF_vgetA1omega(),HelWeightPlus);
	TRFomegabar_a1scalar_minus->Fill(a1h.TRF_vgetA1omega("scalar"),HelWeightMinus);                  TRFomegabar_a1scalar_plus->Fill(a1h.TRF_vgetA1omega("scalar"),HelWeightPlus);
	cosbetacostheta_minus->Fill(a1h.cosbeta(),a1h.costhetaLF(),HelWeightMinus);                              cosbetacostheta_plus->Fill(a1h.cosbeta(),a1h.costhetaLF(),HelWeightPlus);
	TRFcosbetacostheta_minus->Fill(a1h.TRF_cosbeta(),a1h.costhetaLF(),HelWeightMinus);                  TRFcosbetacostheta_plus->Fill(a1h.TRF_cosbeta(),a1h.costhetaLF(),HelWeightPlus);
	
	s1->Fill((a1ospi +a1ss1pi).M() );
	s2->Fill((a1ospi +a1ss2pi).M() );
	qq->Fill((a1ospi+a1ss1pi+a1ss2pi).M());
	
	omega_a1_minus->Fill(Polarimetr.getOmegaA1Bar(),HelWeightMinus);  omega_a1_plus->Fill(Polarimetr.getOmegaA1Bar(),HelWeightPlus);
	hmag->Fill(Polarimetr.PVC().Vect().Mag());
	s1s2->Fill((a1ospi+a1ss2pi).M2(),(a1ospi+a1ss1pi).M2());
	int hel(0);
	if(HelWeightMinus==1)     hel = -1;
	if(HelWeightPlus==1)      hel =  1;
	if(HelWeightPlus==1)      a1_mcos2gamma_plus->Fill(a1h.getMoment(a1h.costhetaLF(),"c2g", 1));                 
	if(HelWeightMinus==1)     a1_mcos2gamma_minus->Fill(a1h.getMoment(a1h.costhetaLF(),"c2g", -1));                 
      }
    }
    
    //----------------------------- pairs -----------------

    
    if(JAK1==4 &&JAK2==3 ){
      if(ApplyCut){
	if(passed1 && passed2){
	  
	  
	  mass_pirho_plus->Fill( (tauandprod1.at(1) + tauandprod1.at(2) + tauandprod2.at(1)).M(), HelWeightPlus);
	  mass_pirho_minus->Fill((tauandprod1.at(1) + tauandprod1.at(2) + tauandprod2.at(1)).M(), HelWeightMinus);
	}
      }else{
	mass_pirho_plus->Fill( (tauandprod1.at(1) + tauandprod1.at(2) + tauandprod2.at(1)).M(), HelWeightPlus);
	mass_pirho_minus->Fill((tauandprod1.at(1) + tauandprod1.at(2) + tauandprod2.at(1)).M(), HelWeightMinus);
      }
    }
    
    if(JAK1==4 &&JAK2==2 ){
      if(ApplyCut){
	if(passed1 && passed2){
	  mass_murho_plus->Fill( (tauandprod1.at(1) + tauandprod1.at(2) + tauandprod2.at(1)).M(), HelWeightPlus);
	  mass_murho_minus->Fill((tauandprod1.at(1) + tauandprod1.at(2) + tauandprod2.at(1)).M(),HelWeightMinus);
	}
      }else{
	  mass_murho_plus->Fill( (tauandprod1.at(1) + tauandprod1.at(2) + tauandprod2.at(1)).M(), HelWeightPlus);
	  mass_murho_minus->Fill((tauandprod1.at(1) + tauandprod1.at(2) + tauandprod2.at(1)).M(),HelWeightMinus);
      }
    }
    
    if(JAK1==3 &&JAK2==2 ){
      if(ApplyCut){
	if(passed1 && passed2){
	  mupi_mass_plus->Fill( (tauandprod1.at(1)  + tauandprod2.at(1)).M(), HelWeightPlus);
	  mupi_mass_minus->Fill( (tauandprod1.at(1)  + tauandprod2.at(1)).M(),HelWeightMinus );
	}
      }else{
	  mupi_mass_plus->Fill( (tauandprod1.at(1)  + tauandprod2.at(1)).M(), HelWeightPlus);
	  mupi_mass_minus->Fill( (tauandprod1.at(1)  + tauandprod2.at(1)).M(),HelWeightMinus );
      }
    }
    
    
    if(JAK1 ==3 && JAK2 == 3){
      if(ApplyCut){
	if(passed1 && passed2){
	  if(Pi1.isConfigured() && Pi2.isConfigured()){
	    double OmPiPi=  (Pi1.getOmega() +Pi2.getOmega() )/(1 + Pi1.getOmega()*Pi2.getOmega());
	    TauPolPiPi.ConfigurePair(tauandprod1,"pion",tauandprod2,"pion");
	    Omegapipi_plus->Fill(TauPolPiPi.getCombOmegaBar(),HelWeightPlus);
	    Omegapipi_minus->Fill(TauPolPiPi.getCombOmegaBar(),HelWeightMinus);       
	    pipi_mass_plus->Fill((tauandprod1.at(1)+tauandprod2.at(1)).M(),HelWeightPlus);
	    pipi_mass_minus->Fill((tauandprod1.at(1)+tauandprod2.at(1)).M(),HelWeightMinus);
	  }
	}
      }else{
	double OmPiPi=  (Pi1.getOmega() +Pi2.getOmega() )/(1 + Pi1.getOmega()*Pi2.getOmega());
	TauPolPiPi.ConfigurePair(tauandprod1,"pion",tauandprod2,"pion");
	Omegapipi_plus->Fill(TauPolPiPi.getCombOmegaBar(),HelWeightPlus);
	Omegapipi_minus->Fill(TauPolPiPi.getCombOmegaBar(),HelWeightMinus);       
	pipi_mass_plus->Fill((tauandprod1.at(1)+tauandprod2.at(1)).M(),HelWeightPlus);
	pipi_mass_minus->Fill((tauandprod1.at(1)+tauandprod2.at(1)).M(),HelWeightMinus);
      } 
    }
    
     
     if(JAK1 ==2 && JAK2 == 3){
       if(Mu1.isConfigured() && Pi2.isConfigured()){
        double OmMuPi=  (Mu1.getOmega() +Pi2.getOmega() )/(1 + Mu1.getOmega()*Pi2.getOmega());
        TauPolMuPi.ConfigurePair(tauandprod1,"lepton",tauandprod2,"pion");
	OmegaMuPi_plus->Fill(TauPolMuPi.getCombOmega(),HelWeightPlus);
        OmegaMuPi_minus->Fill(TauPolMuPi.getCombOmega(),HelWeightMinus);       
       }
     }
     
      if(JAK1 ==2 && JAK2 == 4){
        if(Mu1.isConfigured() && Rho2.isConfigured()){
	  double OmMuRho=  (Mu1.getOmega() +Rho2.getOmega() )/(1 + Mu1.getOmega()*Rho2.getOmega());
	  TauPolMuRho.ConfigurePair(tauandprod1,"lepton",tauandprod2,"rho");
 
	  omega_murho_plus->Fill(TauPolMuRho.getCombOmegaBar(),HelWeightPlus);
	  omega_murho_minus->Fill(TauPolMuRho.getCombOmegaBar(),HelWeightMinus);    
        }
      }

      if(JAK1 ==3 && JAK2 == 4){
	if(Pi1.isConfigured() && Rho2.isConfigured()){
	  double OmPiRho=  (Pi1.getOmega() +Rho2.getOmega() )/(1 + Pi1.getOmega()*Rho2.getOmega());
	   TauPolPiRho.ConfigurePair(tauandprod1,"pion",tauandprod2,"rho");
	   
	   omega_pirho_plus->Fill(TauPolPiRho.getCombOmegaBar(),HelWeightPlus);
	   omega_pirho_minus->Fill(TauPolPiRho.getCombOmegaBar(),HelWeightMinus);       
         }
      }
      
      if(JAK1 ==3 && JAK2 == 5 &&   SubJAK2==51){
	if(Pi1.isConfigured() && A2.isConfigured()){
	  TauPolPiA1.ConfigurePair(tauandprod1,"pion",tauandprod2,"a1",taucharge2);
	   omega_a1pi_plus->Fill(TauPolPiA1.getCombOmegaBar(),HelWeightPlus);
	   omega_a1pi_minus->Fill(TauPolPiA1.getCombOmegaBar(),HelWeightMinus);       
	}
       }
      
      
      if(JAK1 ==2 && JAK2 == 5 &&   SubJAK2==51){
	if(Mu1.isConfigured() && A2.isConfigured()){
	  TauPolMuA1.ConfigurePair(tauandprod1,"lepton",tauandprod2,"a1",taucharge2);
	  omega_a1mu_plus->Fill(TauPolMuA1.getCombOmegaBar(),HelWeightPlus);
	   omega_a1mu_minus->Fill(TauPolMuA1.getCombOmegaBar(),HelWeightMinus);       
	}
      }
      

      if(JAK1 ==4 && JAK2 == 4){
	if(Rho1.isConfigured() && Rho2.isConfigured()){
	  TauPolRhoRho.ConfigurePair(tauandprod1,"rho",tauandprod2,"rho");
	  omega_rhorho_plus->Fill(TauPolRhoRho.getCombVisibleOmega(),HelWeightPlus);
	  omega_rhorho_minus->Fill(TauPolRhoRho.getCombVisibleOmega(),HelWeightMinus);       
	  mass_rhorho_plus->Fill((tauandprod1.at(1) + tauandprod1.at(2) + tauandprod2.at(1)+ tauandprod2.at(2)).M(),HelWeightPlus  );
	  mass_rhorho_minus->Fill((tauandprod1.at(1) + tauandprod1.at(2) + tauandprod2.at(1)+ tauandprod2.at(2)).M(),HelWeightMinus );
	}
      }
      
 
      
      // Run MC-TESTER on the event
      HepMCEvent temp_event(*HepMCEvt,false);
      MC_Analyze(&temp_event);
      
      // Print some events at the end of the run
      if(iEvent>=NumberOfEvents-5){  
	// pythia.event.list();
	HepMCEvt->print();
      }

      // Clean up HepMC event
      delete HepMCEvt;  
  }
  
  
  

  
  file->Write();
  file->Close();
  
  pythia.statistics();
  MC_Finalize();
  
  // This is an access to old FORTRAN info on generated tau sample. 
  // That is why it refers to old version number (eg. 2.7) for TAUOLA.
  //Tauola::summary();
}



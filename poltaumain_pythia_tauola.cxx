/**
 * Example of use of tauola C++ interfate. Pythia events are
 * generated with a stable tau. Taus are subseuently decay via
 * tauola, plots of polatization observables for tau-> mununu and tau-> pinu
 * are produced.
 *
 * @author Vladimir Cherepanov
 * @date  01 MAY 2017
 */ 
     
#include "Tauola/Log.h"
#include "Tauola/Tauola.h"
#include "Tauola/TauolaHepMCEvent.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h" 
#include "TF1.h"
#include "TLorentzVector.h"
#include "UserCodes/a1Helper.h"
#include "UserCodes/TauDecaysHelper.h"
#include "UserCodes/TauPolInterface.h"

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

void redMinus(TauolaParticle *minus)
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
  //   Tauola::setTaukle(double bra1, double brk0, double brk0b,double brks);
  Tauola::setTaukle(1, 0,0,0);
  // can be called here 


  for(unsigned int dec=0; dec <23; dec++){
     double br =0.0; 
         if(dec==2 || dec == 3) br=0.5;
	 // if(dec == 3) br=0.99;
     Tauola::setTauBr(dec, br);
   } 

     // Tauola::setTauBr(0, 0.0);     Tauola::setTauBr(1, 0.0);      Tauola::setTauBr(2, 1.0);      Tauola::setTauBr(3, 0.0);      Tauola::setTauBr(4, 0.0);
     // Tauola::setTauBr(5, 0.0);     Tauola::setTauBr(6, 0.0);      Tauola::setTauBr(7, 0.0);      Tauola::setTauBr(8, 0.0);      Tauola::setTauBr(9, 0.0);
     // Tauola::setTauBr(10, 0.0);   Tauola::setTauBr(11, 0.0);    Tauola::setTauBr(12, 0.0);    Tauola::setTauBr(13, 0.0);    Tauola::setTauBr(14, 0.0);
     // Tauola::setTauBr(15, 0.0);   Tauola::setTauBr(16, 0.0);    Tauola::setTauBr(17, 0.0);    Tauola::setTauBr(18, 0.0);    Tauola::setTauBr(19, 0.0);
     // Tauola::setTauBr(20, 0.0);   Tauola::setTauBr(21, 0.0);    Tauola::setTauBr(22, 0.0);


}

void redPlus(TauolaParticle *plus)
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
  // Tauola::setTaukle(double bra1, double brk0, double brk0b,double brks);
  // can be called here 
  for(unsigned int dec=0; dec <23; dec++){
     double br =0.0;
     //      if(dec==3 || dec ==4 || dec ==5) br=0.33;
      if(dec ==4) br=0.99;
     Tauola::setTauBr(dec, br);
   }


  // Tauola::setTauBr(0, 0.0);     Tauola::setTauBr(1, 0.0);      Tauola::setTauBr(2, 0.0);      Tauola::setTauBr(3, 0.5);      Tauola::setTauBr(4, 0.5);
  // Tauola::setTauBr(5, 0.0);     Tauola::setTauBr(6, 0.0);      Tauola::setTauBr(7, 0.0);      Tauola::setTauBr(8, 0.0);      Tauola::setTauBr(9, 0.0);
  // Tauola::setTauBr(10, 0.0);   Tauola::setTauBr(11, 0.0);    Tauola::setTauBr(12, 0.0);    Tauola::setTauBr(13, 0.0);    Tauola::setTauBr(14, 0.0);
  // Tauola::setTauBr(15, 0.0);   Tauola::setTauBr(16, 0.0);    Tauola::setTauBr(17, 0.0);    Tauola::setTauBr(18, 0.0);    Tauola::setTauBr(19, 0.0);
  // Tauola::setTauBr(20, 0.0);   Tauola::setTauBr(21, 0.0);    Tauola::setTauBr(22, 0.0);


}
void SortPions(    std::vector<HepMC::GenParticle > pionsvec)
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
    
    
    // std::cout<<" os.pdgId()  "<<os.pdg_id() <<std::endl;
    // std::cout<<" ss.pdgId()  "<<ss1.pdg_id() <<std::endl;
    // std::cout<<" ss.pdgId()  "<<ss2.pdg_id() <<std::endl;
    
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
  TFile *file = new TFile("HelicityVals.root","RECREATE");


  TH1F *rhobeta_plus= new TH1F("rhobeta_plus","#rho^{+}",50,-1,1);
  TH1F *rhobeta_minus= new TH1F("rhobeta_minus","#rho^{-}",50,-1,1);

  TH1F *pi_plus= new TH1F("pi_plus","#pi^{+}",50,-1,1);
  TH1F *pi_minus= new TH1F("pi_minus","#pi^{-} ",50,-1,1);

  TH1F *mu_plus= new TH1F("mu_plus","#mu^{+}",50,-1,1);
  TH1F *mu_minus= new TH1F("mu_minus","#mu^{-}",50,-1,1);
  
  TH1F *OmegaMuPi_plus= new TH1F("OmegaMuPi_plus","#omega_{#pi#mu}^{+}",50,-2,2);
  TH1F *OmegaMuPi_minus= new TH1F("OmegaMuPi_minus","#omega_{#pi#mu}^{-}",50,-2,2);
 
  TH1F *Omegapipi_plus= new TH1F("Omegapipi_plus","#omega_{#pi#pi}^{+}",50,-1,1);
  TH1F *Omegapipi_minus= new TH1F("Omegapipi_minus","#omega_{#pi#pi}^{-}",50,-1,1);
  
  TH1F *pipi_mass_plus= new TH1F("pipi_mass_plus","M_{#pi#pi}^{+}",60,20,80);
  TH1F *pipi_mass_minus= new TH1F("pipi_mass_minus","M_{#pi#pi}^{-}",60,20,80);
  
 
 TH1F *omega_murho_plus= new TH1F("omega_murho_plus","#omega_{#mu#rho}^{+}",50,-1,1);
 TH1F *omega_murho_minus= new TH1F("omega_murho_minus","#omega_{#mu#rho}^{-}",50,-1,1);


 TH1F *omega_rho_plus= new TH1F("omega_rho_plus","#omega_{#rho}^{+}",50,-1.1,1.1);
 TH1F *omega_rho_minus= new TH1F("omega_rho_minus","#omega_{#rho}^{-}",50,-1.1,1.1);

 TH1F *omegabar_rho_plus= new TH1F("omegabar_rho_plus","#bar{#omega}_{#rho}^{+}",50,-1.1,1.1);
 TH1F *omegabar_rho_minus= new TH1F("omegabar_rho_minus","#bar{#omega}_{#rho}^{-}",50,-1.1,1.1);


 TH1F *omega_a1_plus= new TH1F("omega_a1_plus","#omega_{a1}^{+}",50,-1,1);
 TH1F *omega_a1_minus= new TH1F("omega_a1_minus","#omega_{a1}^{-}",50,-1,1);

 TH1F *omegabar_a1_plus= new TH1F("omegabar_a1_plus","#bar{#omega}_{a1}^{+}",50,-1,1);
 TH1F *omegabar_a1_minus= new TH1F("omegabar_a1_minus","#bar{#omega}_{a1}^{-}",50,-1,1);
 
 TH1F *TRFomegabar_a1_plus= new TH1F("TRFomegabar_a1_plus","TRF #omega_{a1}^{+}",50,-1,1);
 TH1F *TRFomegabar_a1_minus= new TH1F("TRFomegabar_a1_minus","TRF #omega_{a1}^{-}",50,-1,1);

  TH1F *TRFomegabar_a1scalar_plus= new TH1F("TRFomegabar_a1scalar_plus","TRF scalar #omega a1",50,-1,1);
  TH1F *TRFomegabar_a1scalar_minus= new TH1F("TRFomegabar_a1scalar_minus","TRF scalar  #omega a1",50,-1,1);
 

 TH2F *cosbetacostheta_plus= new TH2F("cosbetacostheta_plus","cos#beta  cos#theta   {a1}^{+}",50,-1,1,50,-1,1);
 TH2F *cosbetacostheta_minus= new TH2F("cosbetacostheta_minus","cos#beta  cos#theta  {a1}^{-}",50,-1,1,50,-1,1);
 cosbetacostheta_plus->GetXaxis()->SetTitle("cos#beta"); cosbetacostheta_plus->GetYaxis()->SetTitle("cos#theta");
 cosbetacostheta_minus->GetXaxis()->SetTitle("cos#beta"); cosbetacostheta_minus->GetYaxis()->SetTitle("cos#theta");

 TH2F *cosbetacosthetarho_plus= new TH2F("cosbetacosthetarho_plus","cos#beta  cos#theta  #rho^{+}",50,-1,1,50,-1,1);
 TH2F *cosbetacosthetarho_minus= new TH2F("cosbetacosthetarho_minus","cos#beta  cos#theta  #rho^{-}" ,50,-1,1,50,-1,1);
 cosbetacosthetarho_plus->GetXaxis()->SetTitle("cos#beta"); cosbetacosthetarho_plus->GetYaxis()->SetTitle("cos#theta");
 cosbetacosthetarho_minus->GetXaxis()->SetTitle("cos#beta"); cosbetacosthetarho_minus->GetYaxis()->SetTitle("cos#theta");

 TH2F *TRFcosbetacostheta_plus= new TH2F("TRFcosbetacostheta_plus","cos#beta  cos#theta  {a1}^{+}",50,-1,1,50,-1,1);
 TH2F *TRFcosbetacostheta_minus= new TH2F("TRFcosbetacostheta_minus","cos#beta  cos#theta  {a1}^{-}",50,-1,1,50,-1,1);
 TRFcosbetacostheta_plus->GetXaxis()->SetTitle("cos#beta"); TRFcosbetacostheta_plus->GetYaxis()->SetTitle("cos#theta");
 TRFcosbetacostheta_minus->GetXaxis()->SetTitle("cos#beta"); TRFcosbetacostheta_minus->GetYaxis()->SetTitle("cos#theta");


 TH1F *omega_a1pi_plus= new TH1F("omega_a1pi_plus","#omega_{a1#pi}^{+}",50,-1,1);
 TH1F *omega_a1pi_minus= new TH1F("omega_a1pi_minus","#omega_{a1#pi}^{-}",50,-1,1);

 TH1F *omega_a1mu_plus= new TH1F("omega_a1mu_plus","#omega_{a1#mu}^{+}",50,-1,1);
 TH1F *omega_a1mu_minus= new TH1F("omega_a1mu_minus","#omega_{a1#mu}^{-}",50,-1,1);

 TH1F *omega_pirho_plus= new TH1F("omega_pirho_plus","#omega_{#pi#rho}^{+}",50,-1,1);
 TH1F *omega_pirho_minus= new TH1F("omega_pirho_minus","#omega_{#pi#rho}^{-}",50,-1,1);
  


 // TH1F *omega_a1rho_plus= new TH1F("omega_a1rho_plus","#omega_{a1#rho}^{+}",50,-1,1);
 // TH1F *omega_a1rho_minus= new TH1F("omega_a1rho_minus","#omega_{a1#rho}^{-}",50,-1,1);


 // TH1F *omega_pirho_plus= new TH1F("omega_pirho_plus","#omega_{#pi#rho}^{+}",50,-1,1);
 // TH1F *omega_pirho_minus= new TH1F("omega_pirho_minus","#omega_{#pi#rho}^{-}",50,-1,1);






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

  //  Log::LogDebug(true);

    Tauola::setRedefineTauMinus(redMinus);  // activates execution of routine redMinus in TAUOLA interface
    Tauola::setRedefineTauPlus(redPlus);    // activates execution of routine redPlus  in TAUOLA interface

  MC_Initialize();

  // Begin event loop. Generate event.
  for (int iEvent = 0; iEvent < NumberOfEvents; ++iEvent){

    if(iEvent%1000==0) Log::Info()<<"Event: "<<iEvent<<endl;
    if (!pythia.next()) continue;

    // Convert event record to HepMC
    HepMC::GenEvent * HepMCEvt = new HepMC::GenEvent();

    // Conversion needed if HepMC use different momentum units
    // than Pythia. However, requires HepMC 2.04 or higher.
    HepMCEvt->use_units(HepMC::Units::GEV,HepMC::Units::MM);

    ToHepMC.fill_next_event(event, HepMCEvt);

    if(iEvent<EventsToCheck)
    {
      cout<<"                                          "<<endl;
      cout<<"Momentum conservation chceck BEFORE/AFTER Tauola"<<endl;
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
    std::vector<HepMC::GenParticle > SortA1Pions;  //os, ss1, ss2
    for ( HepMC::GenEvent::particle_const_iterator p =HepMCEvt->particles_begin();  p != HepMCEvt->particles_end(); ++p ){  
      if((*p)->pdg_id()==15){
	FirstTau = *p;
	for ( HepMC::GenEvent::particle_const_iterator d =HepMCEvt->particles_begin();  d != HepMCEvt->particles_end(); ++d ){  
	  if((*d)->pdg_id()!=15){
	    if((*p)->end_vertex() == (*d)->production_vertex()){
	      FirstTauProducts.push_back(**d);
	      if(abs((*d)->pdg_id()) ==  12) {JAK1 =1; 
	      }else if(abs((*d)->pdg_id()) ==  14){ JAK1=2;
	      }else if(abs((*d)->pdg_id())==  211){ JAK1 = 3;}
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
		  for ( HepMC::GenEvent::particle_const_iterator dd =HepMCEvt->particles_begin();  dd != HepMCEvt->particles_end(); ++dd ){  
		    if( abs((*dd)->pdg_id())!=20213  ){
		      if((*d)->end_vertex() == (*dd)->production_vertex()){
		       
			FirstTauProducts.push_back(**dd);
			if(abs((*dd)->pdg_id())==  211)   {
			  A1Pions.push_back(**dd); npi++;
			}
		      }
		    }
		  }
		  if(npi==3) SubJAK1=51; else SubJAK1=52;
		}
	    }
	  } 
	}
      }
    
      
      if((*p)->pdg_id()==-15){
	SecondTau = *p;
	for ( HepMC::GenEvent::particle_const_iterator d =HepMCEvt->particles_begin();  d != HepMCEvt->particles_end(); ++d )
	  {  
	    if((*d)->pdg_id()!=-15){
	      if((*p)->end_vertex() == (*d)->production_vertex()){
		SecondTauProducts.push_back(**d);
		if(abs((*d)->pdg_id()) ==  12) {JAK2 =1; 
		}else if(abs((*d)->pdg_id()) ==  14){ JAK2=2;
		}else if(abs((*d)->pdg_id())==  211){ JAK2 = 3;}

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
			  A1Pions.push_back(**dd); npi++;
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
  
    Log::Debug(5)<<"helicites =  "<<Tauola::getHelPlus()<<" "<<Tauola::getHelMinus()
                 <<" electroweak wt= "<<Tauola::getEWwt()<<endl;

    bool HelPlus=false;
    bool HelMinus=false;
    if(Tauola::getHelPlus() ==1 )HelPlus=true;;
    if(Tauola::getHelPlus() ==-1)HelMinus=true;;
    int HelWeightPlus = HelPlus;
    int HelWeightMinus = HelMinus;


   
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

    TauDecaysHelper Mu1;
    TauDecaysHelper Pi1;
    TauDecaysHelper Pi2;
    TauDecaysHelper Rho2;
    a1Helper a1h;
    rhoHelper RhoHelp;

    TauDecaysHelper RhoTPI;
    TauPolInterface TauPolPi1;
    TauPolInterface TauPolPi2;
    TauPolInterface TauPolPiPi;


    TauPolInterface TauPolMu1;
    TauPolInterface TauPolRho2;
    TauPolInterface TauPolA1;
 
    TauPolInterface TauPolMuPi;
    TauPolInterface TauPolMuRho;

    TauPolInterface TauPolPiRho;
    TauPolInterface TauPolPiA1;
    TauPolInterface TauPolMuA1;

    vector<TLorentzVector> tauandprod1,tauandprod2 ;
    vector<TLorentzVector> tauandprodMuon1,tauandprodRho,tauandprodA1;
    tau1.SetPxPyPzE(FirstTau->momentum().px(), FirstTau->momentum().py(), FirstTau->momentum().pz(), FirstTau->momentum().e());
    tau2.SetPxPyPzE(SecondTau->momentum().px(), SecondTau->momentum().py(), SecondTau->momentum().pz(), SecondTau->momentum().e());
    //    tau1.Print();
    //    tau2.Print();
    if(JAK1==2){
      vector<TLorentzVector> tauandprod;
      tauandprod.push_back(TLorentzVector(FirstTau->momentum().px(), FirstTau->momentum().py(), FirstTau->momentum().pz(), FirstTau->momentum().e()));
      for(std::vector<HepMC::GenParticle>::const_iterator a = FirstTauProducts.begin(); a!=FirstTauProducts.end(); ++a){
	if(abs(a->pdg_id())==13){tauandprod.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );} }
       Mu1.Configure(tauandprod,"lepton");
      tauandprodMuon1=tauandprod;
      TauPolMu1.Configure(tauandprod,"lepton");
      //      mu_plus->Fill(Mu1.getOmega(),HelWeightPlus);
      //      mu_minus->Fill(Mu1.getOmega(),HelWeightMinus);

      mu_plus->Fill(TauPolMu1.getOmega(),HelWeightPlus);
      mu_minus->Fill(TauPolMu1.getOmega(),HelWeightMinus);
    }


    if(JAK1==3){
      vector<TLorentzVector> tauandprod;
      tauandprod.push_back(TLorentzVector(FirstTau->momentum().px(), FirstTau->momentum().py(), FirstTau->momentum().pz(), FirstTau->momentum().e()));
      for(std::vector<HepMC::GenParticle>::const_iterator a = FirstTauProducts.begin(); a!=FirstTauProducts.end(); ++a){
	if(abs(a->pdg_id())==211){tauandprod.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );} }

  

      tauandprod1 = tauandprod;
      Pi1.Configure(tauandprod,"pion");
      TauPolPi1.Configure(tauandprod,"pion");
      pi_plus->Fill(TauPolPi1.getOmega(),HelWeightPlus);
      pi_minus->Fill(TauPolPi1.getOmega(),HelWeightMinus);
    }

    if(JAK2==3){
      vector<TLorentzVector> tauandprod;
      tauandprod.push_back(TLorentzVector(SecondTau->momentum().px(), SecondTau->momentum().py(), SecondTau->momentum().pz(), SecondTau->momentum().e()));
      for(std::vector<HepMC::GenParticle>::const_iterator a = SecondTauProducts.begin(); a!=SecondTauProducts.end(); ++a){
	if(abs(a->pdg_id())==211){tauandprod.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );} }
      tauandprod2 = tauandprod;
      TauPolPi2.Configure(tauandprod,"pion");
      Pi2.Configure(tauandprod,"pion");
    }
    

    
    //    std::cout<<  "Starting a new event loop .... event number  : "<< iEvent << std::endl;
    if(JAK2==4){
      //      std::cout<<"Rho decay is found! ... Jak = "<< JAK2 <<std::endl;

      vector<TLorentzVector> tauandprod;
      tauandprod.push_back(TLorentzVector(SecondTau->momentum().px(), SecondTau->momentum().py(), SecondTau->momentum().pz(), SecondTau->momentum().e()));
      //      std::cout<<" ---- "<<std::endl;
      for(std::vector<HepMC::GenParticle>::const_iterator a = SecondTauProducts.begin(); a!=SecondTauProducts.end(); ++a){
	//	std::cout<<" rho decays:   "<< a->pdg_id() << std::endl;
	if(abs(a->pdg_id())==211){tauandprod.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );}
	if(abs(a->pdg_id())==111){tauandprod.push_back(TLorentzVector(a->momentum().px(), a->momentum().py(), a->momentum().pz(), a->momentum().e()  ) );}
      }
      RhoTPI.Configure(tauandprod,"rho");
      Rho2.Configure(tauandprod,"rho");

      tauandprodRho=tauandprod;  
      TauPolRho2.Configure(tauandprod,"rho");

      RhoHelp.Configure(tauandprod,tauandprod.at(1) + tauandprod.at(2));

      // rhobeta_plus->Fill(TauPolRho2.getOmega(),HelWeightPlus);
      // rhobeta_minus->Fill(TauPolRho2.getOmega(),HelWeightMinus);


      // omega_rho_plus->Fill(RhoTPI.getCosbetaRho(),HelWeightPlus);
      // omega_rho_minus->Fill(RhoTPI.getCosbetaRho(),HelWeightMinus);

      rhobeta_plus->Fill(RhoHelp.getCosBetaTest(),HelWeightPlus);
      rhobeta_minus->Fill(RhoHelp.getCosBetaTest(),HelWeightMinus);
      // std::cout<<"   "<<HelWeightPlus<< "  "<< HelWeightMinus <<std::endl;

      omega_rho_plus->Fill(RhoHelp.getCosbetaRho(),HelWeightPlus);
      omega_rho_minus->Fill(RhoHelp.getCosbetaRho(),HelWeightMinus);



      //   std::cout<<" rho beta 1   "<< RhoTPI.getCosbetaRho()<<std::endl;
     //    std::cout<<" rho TDC    "<< Rho2.getOmega()<<std::endl;
    //     std::cout<<" rho +   "<< TauPolRho2.getOmega()*HelWeightPlus<<std::endl;
     //     std::cout<<" rho -   "<< TauPolRho2.getOmega()*HelWeightMinus<<std::endl;
     // omegabar_rho_plus->Fill(TauPolRho2.getOmegabar(),HelWeightPlus);
     // omegabar_rho_minus->Fill(TauPolRho2.getOmegabar(),HelWeightMinus);
     omegabar_rho_plus->Fill(RhoHelp.getOmegaRhoBar(),HelWeightPlus);
     omegabar_rho_minus->Fill(RhoHelp.getOmegaRhoBar(),HelWeightMinus);

     //     omegabar_rho_plus->Fill(Rho2.getOmegabar(),HelWeightPlus);
     //     omegabar_rho_minus->Fill(Rho2.getOmegabar(),HelWeightMinus);
     //     std::cout<<  " omegabar     "<< Rho2.getOmegabar() <<std::endl;


     cosbetacosthetarho_plus->Fill(Rho2.getCosbetaRho(),Rho2.getCosthetaRho(),HelWeightPlus);
     cosbetacosthetarho_minus->Fill(Rho2.getCosbetaRho(),Rho2.getCosthetaRho(),HelWeightMinus);

    }

    if(JAK2==5 && SubJAK2==51){
      vector<TLorentzVector> particles;
      SortPions(A1Pions);
      a1ospi.SetPxPyPzE(A1Pions.at(0).momentum().px(), A1Pions.at(0).momentum().py(), A1Pions.at(0).momentum().pz(), A1Pions.at(0).momentum().e());
      a1ss1pi.SetPxPyPzE(A1Pions.at(1).momentum().px(), A1Pions.at(1).momentum().py(), A1Pions.at(1).momentum().pz(), A1Pions.at(1).momentum().e());
      a1ss2pi.SetPxPyPzE(A1Pions.at(2).momentum().px(), A1Pions.at(2).momentum().py(), A1Pions.at(2).momentum().pz(), A1Pions.at(2).momentum().e());
      particles.push_back(tau1);
      particles.push_back(a1ospi);
      particles.push_back(a1ss1pi);
      particles.push_back(a1ss2pi);
      a1h.Configure(particles, a1ospi+a1ss1pi+a1ss2pi);
      TauPolA1.Configure(particles,"a1");
      tauandprodA1=particles;
      //      omega_a1_minus->Fill(a1h.getA1omega(),HelWeightMinus);                                                          omega_a1_plus->Fill(a1h.getA1omega(),HelWeightPlus);

      omega_a1_minus->Fill(TauPolA1.getOmega(),HelWeightMinus);                                                          omega_a1_plus->Fill(TauPolA1.getOmega(),HelWeightPlus);
      TRFomegabar_a1_minus->Fill(a1h.TRF_vgetA1omega(),HelWeightMinus);                                      TRFomegabar_a1_plus->Fill(a1h.TRF_vgetA1omega(),HelWeightPlus);
      TRFomegabar_a1scalar_minus->Fill(a1h.TRF_vgetA1omega("scalar"),HelWeightMinus);                  TRFomegabar_a1scalar_plus->Fill(a1h.TRF_vgetA1omega("scalar"),HelWeightPlus);
      cosbetacostheta_minus->Fill(a1h.cosbeta(),a1h.costhetaLF(),HelWeightMinus);                              cosbetacostheta_plus->Fill(a1h.cosbeta(),a1h.costhetaLF(),HelWeightPlus);
      TRFcosbetacostheta_minus->Fill(a1h.TRF_cosbeta(),a1h.costhetaLF(),HelWeightMinus);                  TRFcosbetacostheta_plus->Fill(a1h.TRF_cosbeta(),a1h.costhetaLF(),HelWeightPlus);
      omegabar_a1_minus->Fill(a1h.vgetA1omega("bar")  ,HelWeightMinus);                                          omegabar_a1_plus->Fill(a1h.vgetA1omega("bar"),HelWeightPlus);
      //      std::cout<<"MomentSFunction "<< a1h.MomentSFunction(0.89,"WA") <<std::endl;
      s1s2->Fill((a1ospi+a1ss2pi).M2(),(a1ospi+a1ss1pi).M2());
    }
 


    if(JAK1 ==3 && JAK2 == 3){
      if(Pi1.isConfigured() && Pi2.isConfigured()){
      double OmPiPi=  (Pi1.getOmega() +Pi2.getOmega() )/(1 + Pi1.getOmega()*Pi2.getOmega());
      TauPolPiPi.ConfigurePair(tauandprod1,"pion",tauandprod2,"pion");

 



      //      Omegapipi_plus->Fill(OmPiPi,HelWeightPlus);
      //      Omegapipi_minus->Fill(OmPiPi,HelWeightMinus);       

      //      std::cout<<" ----  "<<std::endl;
      //      tauandprod1.at(1).Print();
      //      tauandprod2.at(1).Print();
      Omegapipi_plus->Fill(TauPolPiPi.getCombOmega(),HelWeightPlus);
      Omegapipi_minus->Fill(TauPolPiPi.getCombOmega(),HelWeightMinus);       
      pipi_mass_plus->Fill(TauPolPiPi.getVisiblePairLV().M(),HelWeightPlus);
      pipi_mass_minus->Fill(TauPolPiPi.getVisiblePairLV().M(),HelWeightMinus);

      }
    }

    if(JAK1 ==2 && JAK2 == 3){
      if(Mu1.isConfigured() && Pi2.isConfigured()){
      double OmMuPi=  (Mu1.getOmega() +Pi2.getOmega() )/(1 + Mu1.getOmega()*Pi2.getOmega());
      TauPolMuPi.ConfigurePair(tauandprodMuon1,"lepton",tauandprod2,"pion");

 



      //      OmegaMuPi_plus->Fill(OmMuPi,HelWeightPlus);
      //      OmegaMuPi_minus->Fill(OmMuPi,HelWeightMinus);       

      OmegaMuPi_plus->Fill(TauPolMuPi.getCombOmega(),HelWeightPlus);
      OmegaMuPi_minus->Fill(TauPolMuPi.getCombOmega(),HelWeightMinus);       
      }
    }

    if(JAK1 ==2 && JAK2 == 4){
      if(Mu1.isConfigured() && Rho2.isConfigured()){
      double OmMuRho=  (Mu1.getOmega() +Rho2.getOmega() )/(1 + Mu1.getOmega()*Rho2.getOmega());
      TauPolMuRho.ConfigurePair(tauandprodMuon1,"lepton",tauandprodRho,"rho");
      //      omega_murho_plus->Fill(OmMuRho,HelWeightPlus);
      //      omega_murho_minus->Fill(OmMuRho,HelWeightMinus);       

    

      omega_murho_plus->Fill(TauPolMuRho.getCombOmegaBar(),HelWeightPlus);
      omega_murho_minus->Fill(TauPolMuRho.getCombOmegaBar(),HelWeightMinus);    

      // if( fabs(TauPolMuRho.getCombOmega())> 1)std::cout<<" mu:    "<< TauPolMuRho.getOmega("first") << "  rho:   "<<  TauPolMuRho.getOmega("second") << " combined    " <<TauPolMuRho.getCombOmega() <<std::endl;
      }
    }
    if(JAK1 ==3 && JAK2 == 4){
      if(Pi1.isConfigured() && Rho2.isConfigured()){
      double OmPiRho=  (Pi1.getOmega() +Rho2.getOmega() )/(1 + Pi1.getOmega()*Rho2.getOmega());
      TauPolPiRho.ConfigurePair(tauandprod1,"pion",tauandprodRho,"rho");


      omega_pirho_plus->Fill(TauPolPiRho.getCombOmega(),HelWeightPlus);
      omega_pirho_minus->Fill(TauPolPiRho.getCombOmega(),HelWeightMinus);       
      }
    }


    if(JAK1 ==3 && JAK2 == 5 &&   SubJAK2==51){
      if(Pi1.isConfigured() && a1h.isConfigured()){
	double OmPiA1=  (Pi1.getOmega() +a1h.TRF_vgetA1omega() )/(1 + Pi1.getOmega()*a1h.TRF_vgetA1omega());
	TauPolPiA1.ConfigurePair(tauandprod1,"pion",tauandprodA1,"a1");
	omega_a1pi_plus->Fill(TauPolPiA1.getCombOmega(),HelWeightPlus);
	omega_a1pi_minus->Fill(TauPolPiA1.getCombOmega(),HelWeightMinus);       




      //      omega_a1pi_plus->Fill(OmPiA1,HelWeightPlus);
      //      omega_a1pi_minus->Fill(OmPiA1,HelWeightMinus);       
      }
    }


    if(JAK1 ==2 && JAK2 == 5 &&   SubJAK2==51){
      if(Mu1.isConfigured() && a1h.isConfigured()){
	double OmMuA1=  (Mu1.getOmega() +a1h.TRF_vgetA1omega() )/(1 + Mu1.getOmega()*a1h.TRF_vgetA1omega());
	TauPolMuA1.ConfigurePair(tauandprodMuon1,"lepton",tauandprodA1,"a1");
	omega_a1mu_plus->Fill(TauPolMuA1.getCombOmegaBar(),HelWeightPlus);
	omega_a1mu_minus->Fill(TauPolMuA1.getCombOmegaBar(),HelWeightMinus);       
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


// Omegapipi_plus->Write();
// Omegapipi_minus->Write(); 
// Omegapirho_plus->Write();
// Omegapirho_minus->Write();
// rho_plus->Write();
// rho_minus->Write();
// Omegamurho_minus->Write();
// Omegamurho_plus->Write();
// mu_plus->Write();
// mu_minus->Write();
// OmegaMuPi_plus->Write();
// OmegaMuPi_minus->Write();
// pi_plus->Write();
// pi_minus->Write();
// omega_a1_minus->Write();
// omega_a1_plus->Write();
// omegabar_a1_minus->Write();
// omegabar_a1_plus->Write();
// cosbetacostheta_minus->Write();
// cosbetacostheta_plus->Write();
// TRFomegabar_a1_plus->Write();
// TRFomegabar_a1_minus->Write();
// TRFomegabar_a1scalar_plus->Write();
// TRFomegabar_a1scalar_minus->Write();
// omega_a1pi_plus->Write();
// omega_a1pi_minus->Write();
// omega_a1mu_plus->Write();
// omega_a1mu_minus->Write();



  file->Write();
  file->Close();

  pythia.statistics();
  MC_Finalize();
 
  // This is an access to old FORTRAN info on generated tau sample. 
  // That is why it refers to old version number (eg. 2.7) for TAUOLA.
  //Tauola::summary();
}



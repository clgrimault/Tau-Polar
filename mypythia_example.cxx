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
#include "TLorentzVector.h"

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

int NumberOfEvents = 2000000; 
int EventsToCheck=10;

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
  // Tauola::setTaukle(double bra1, double brk0, double brk0b,double brks);
  // can be called here 
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
}

int main(int argc,char **argv){

  Log::SummaryAtExit();

  // Initialization of pythia
  Pythia pythia;
  Event& event = pythia.event;
  TFile *file = new TFile("output.root","RECREATE");

  TH1F *pi_plus= new TH1F("pi_plus","#pi  plus",50,0,1);
  TH1F *pi_minus= new TH1F("pi_minus","#pi minus",50,0,1);

  TH1F *piom_plus= new TH1F("piom_plus","#pi  plus",50,-1,1);
  TH1F *piom_minus= new TH1F("piom_minus","#pi minus",50,-1,1);


  TH1F *mu_plus= new TH1F("mu_plus","#mu plus",50,0,1);
  TH1F *mu_minus= new TH1F("mu_minus","#mu  minus",50,0,1);
  
  TH1F *om_plus= new TH1F("om_plus","#omega plus",50,-1,1);
  TH1F *om_minus= new TH1F("om_minus","#omega  minus",50,-1,1);
  

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
    Tauola::setSameParticleDecayMode(2);     //19 and 22 contains K0 
    Tauola::setOppositeParticleDecayMode(3); // 20 contains eta

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

  //  Tauola::setRedefineTauMinus(redMinus);  // activates execution of routine redMinus in TAUOLA interface
  //  Tauola::setRedefineTauPlus(redPlus);    // activates execution of routine redPlus  in TAUOLA interface

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
    TLorentzVector tau1(0,0,0,0);

    TLorentzVector mu1(0,0,0,0);
    TLorentzVector pi2(0,0,0,0);

    TLorentzVector numu1(0,0,0,0);
    TLorentzVector nutau1(0,0,0,0);
    TLorentzVector nutau2(0,0,0,0);


    TLorentzVector tau2(0,0,0,0);
    int barcodetau1vertex(0),barcodetau2vertex(0);
    for ( HepMC::GenEvent::particle_const_iterator p =HepMCEvt->particles_begin();  p != HepMCEvt->particles_end(); ++p )
      {

	if((*p)->pdg_id()==15){
	  HepMC::GenVertex * tau1endvertex =(*p)->end_vertex(); 
	  barcodetau1vertex = tau1endvertex->barcode();  
	  tau1.SetPxPyPzE((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
	}
	if((*p)->pdg_id()==-15){
	  HepMC::GenVertex * tau2endvertex =(*p)->end_vertex(); 
	  tau2.SetPxPyPzE((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
	  barcodetau2vertex = tau2endvertex->barcode();  
	}
	if( (*p)->status() == 1 )
	  {
	    HepMC::FourVector m = (*p)->momentum();
	    HepMC::GenVertex *productProductionvertex = (*p)->production_vertex();
	    if(productProductionvertex->barcode()==barcodetau1vertex ){
	      if(abs((*p)->pdg_id())==13)mu1.SetPxPyPzE((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
	      if(abs((*p)->pdg_id())==14)numu1.SetPxPyPzE((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
	      if(abs((*p)->pdg_id())==16)nutau1.SetPxPyPzE((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
	    }
	    if(productProductionvertex->barcode()==barcodetau2vertex ){
	      if(abs((*p)->pdg_id())==16)nutau2.SetPxPyPzE((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
	      if(abs((*p)->pdg_id())==211)pi2.SetPxPyPzE((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
	    }
	  }
      }
    Log::Debug(5)<<"helicites =  "<<Tauola::getHelPlus()<<" "<<Tauola::getHelMinus()
                 <<" electroweak wt= "<<Tauola::getEWwt()<<endl;
    double x1p(0),x2p(0);
    double x1m(0),x2m(0);
    double omp(0),omm(0),OmegaP(0), OmegaM(0);
    if(Tauola::getHelPlus() ==1){
      if(tau2.E()!=0 && pi2.E()!=0) {x1p=2*pi2.E()/tau2.E() - 1; pi_plus->Fill(pi2.E()/tau2.E()); piom_plus->Fill(x1p);}
      if(tau1.E()!=0 && mu1.E()!=0){x2p =2*mu1.E()/tau1.E()-1;  mu_plus->Fill(mu1.E()/tau1.E());}
      if(x1p!=0 && x2p!=0){OmegaP = (x1p+x2p)/(1+x1p*x2p);om_plus->Fill(OmegaP);}
    }

    if(Tauola::getHelPlus() ==-1 ){
      if(tau2.E()!=0 && pi2.E()!=0) { x1m = 2*pi2.E()/tau2.E()-1; pi_minus->Fill(pi2.E()/tau2.E()); piom_minus->Fill(x1m);}
      if(tau1.E()!=0 && mu1.E()!=0) {x2m = 2*mu1.E()/tau1.E()-1; mu_minus->Fill(mu1.E()/tau1.E());}
      if(x1m!=0 && x2m!=0){OmegaM = (x1m+x2m)/(1+x1m*x2m); om_minus->Fill(OmegaM);}
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
  pi_plus->Write();
  pi_minus->Write();

  mu_plus->Write();
  mu_minus->Write();
  om_plus->Write();
  om_minus->Write();
  piom_plus->Write();
  piom_minus->Write();

  file->Write();
  file->Close();
  pythia.statistics();
  MC_Finalize();

  // This is an access to old FORTRAN info on generated tau sample. 
  // That is why it refers to old version number (eg. 2.7) for TAUOLA.
  //Tauola::summary();
}


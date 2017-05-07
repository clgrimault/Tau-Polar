/**
 * Example of use of tauola C++ interface. Pythia events are
 * generated with a stable tau. Taus are subsequently decay via
 * tauola.
 *
 * @author Nadia Davidson
 * @date 17 June 2008
 */

#include "Tauola/Log.h"
#include "Tauola/Tauola.h"
#include "Tauola/TauolaHepMCEvent.h"
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

int NumberOfEvents = 1000;
int EventsToCheck=20;

// elementary test of HepMC typically executed before
// detector simulation based on http://home.fnal.gov/~mrenna/HCPSS/HCPSShepmc.html
// similar test was performed in Fortran
// we perform it before and after Tauola (for the first several events)
void checkMomentumConservationInEvent(HepMC::GenEvent *evt)
{
	//cout<<"List of stable particles: "<<endl;

	double px=0.0,py=0.0,pz=0.0,e=0.0;
	
	for ( HepMC::GenEvent::particle_const_iterator p = evt->particles_begin();
	      p != evt->particles_end(); ++p )
	{
		if( (*p)->status() == 1 )
		{
			HepMC::FourVector m = (*p)->momentum();
			px+=m.px();
			py+=m.py();
			pz+=m.pz();
			e +=m.e();
			//(*p)->print();
		}
	}
  cout.precision(6);
  cout.setf(ios_base::floatfield);
	cout<<endl<<"Vector Sum: "<<px<<" "<<py<<" "<<pz<<" "<<e<<endl;
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
  //Tauola::setTauBr(3, 0.5);      Tauola::setTauBr(4, 0.5);
     Tauola::setTauBr(0, 0.0);     Tauola::setTauBr(1, 0.0);      Tauola::setTauBr(2, 0.0);      Tauola::setTauBr(3, 0.5);      Tauola::setTauBr(4, 0.0);
     Tauola::setTauBr(5, 0.0);     Tauola::setTauBr(6, 0.0);      Tauola::setTauBr(7, 0.0);      Tauola::setTauBr(8, 0.0);      Tauola::setTauBr(9, 0.0);
     Tauola::setTauBr(10, 0.0);   Tauola::setTauBr(11, 0.0);    Tauola::setTauBr(12, 0.0);    Tauola::setTauBr(13, 0.0);    Tauola::setTauBr(14, 0.0);
     Tauola::setTauBr(15, 0.0);   Tauola::setTauBr(16, 0.0);    Tauola::setTauBr(17, 0.0);    Tauola::setTauBr(18, 0.0);    Tauola::setTauBr(19, 0.0);
     Tauola::setTauBr(20, 0.0);   Tauola::setTauBr(21, 0.0);    Tauola::setTauBr(22, 0.0);

}

int main(int argc,char **argv){

  Log::SummaryAtExit();

  // Initialization of pythia
  Pythia pythia;
  Event& event = pythia.event;

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
  //  Tauola::setOppositeParticleDecayMode(3); // 20 contains eta
  Tauola::setOppositeParticleDecayMode(0); // 20 contains eta

  // Set Higgs scalar-pseudoscalar mixing angle
  //  Tauola::setHiggsScalarPseudoscalarMixingAngle(0.7853);
  //  Tauola::setHiggsScalarPseudoscalarPDG(25);

  Tauola::initialize();  
   // Tauola::setTauBr(0, 0);     Tauola::setTauBr(1, 0);      Tauola::setTauBr(2, 0);      Tauola::setTauBr(3, 0);      Tauola::setTauBr(4, 0);
   // Tauola::setTauBr(5, 0);     Tauola::setTauBr(6, 0);      Tauola::setTauBr(7, 0);      Tauola::setTauBr(8, 0);      Tauola::setTauBr(9, 0);
   // Tauola::setTauBr(10, 0);   Tauola::setTauBr(11, 0);    Tauola::setTauBr(12, 0);    Tauola::setTauBr(13, 0);    Tauola::setTauBr(14, 0);
   // Tauola::setTauBr(15, 0);   Tauola::setTauBr(16, 0);    Tauola::setTauBr(17, 0);    Tauola::setTauBr(18, 0);    Tauola::setTauBr(19, 0);
   // Tauola::setTauBr(20, 0);   Tauola::setTauBr(21, 0);    Tauola::setTauBr(22, 0);
  tauola_print_parameters(); // Prints TAUOLA  parameters (residing inside its library): e.g. to test user interface

  // Our default units are GEV and MM, that will be outcome  units after TAUOLA
  // if HepMC unit variables  are correctly set. 
  // with the following coice you can fix the units for final outcome:
  //  Tauola::setUnits(Tauola::GEV,Tauola::MM); 
  //  Tauola::setUnits(Tauola::MEV,Tauola::CM); 

  // Other usefull settings:
  //  Tauola::setEtaK0sPi(0,0,0);  // switches to decay eta K0_S and pi0 1/0 on/off. 
  //  Tauola::setTauLifetime(0.0); //new tau lifetime in mm
  //  Tauola::spin_correlation.setAll(false);

  //  Log::LogDebug(true);

  //  Tauola::setRedefineTauMinus(redMinus);  // activates execution of routine redMinus in TAUOLA interface
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

    int JAK1(0);
    HepMC::GenParticle *FirstTau;
    std::vector<HepMC::GenParticle > FirstTauProducts;
    int JAK2(0);
    HepMC::GenParticle *SecondTau;
    std::vector<HepMC::GenParticle > SecondTauProducts;
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
	      }
	    } 
	  }
      }
    }
  
    // for(std::vector<HepMC::GenParticle>::const_iterator a = FirstTauProducts.begin(); a!=FirstTauProducts.end(); ++a){
    //   cout<<"  00000000000  Loop Over  FirstTauProducts  " << a->pdg_id()<< "   status " <<a->status()<<  "   JAK1   " << JAK1<<endl;
    // }
    
    

  

    // for(std::vector<HepMC::GenParticle>::const_iterator a = SecondTauProducts.begin(); a!=SecondTauProducts.end(); ++a){

    //   cout<<"  00000000000  Loop Over   SecondTauProducts   " << a->pdg_id()<< "  ststua   "<<a->status()<<   "   px  " <<a->momentum().px() << "   SecondTau  "<< SecondTau->momentum().px()<<  "   JAK2   " << JAK2<<endl;

    // }






    Log::Debug(5)<<"helicites =  "<<Tauola::getHelPlus()<<" "<<Tauola::getHelMinus()
                  <<" electroweak wt= "<<Tauola::getEWwt()<<endl;

    // Run MC-TESTER on the event
    HepMCEvent temp_event(*HepMCEvt,false);
    MC_Analyze(&temp_event);

    // Print some events at the end of the run
    if(iEvent>=NumberOfEvents-5){  
      pythia.event.list();
      HepMCEvt->print();
    }

    // Clean up HepMC event
    delete HepMCEvt;  
  }

  pythia.statistics();
  MC_Finalize();

  // This is an access to old FORTRAN info on generated tau sample. 
  // That is why it refers to old version number (eg. 2.7) for TAUOLA.
  // Tauola::summary();
}


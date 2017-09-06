#ifndef TauPolInterface_h
#define TauPolInterface_h
#include <vector>
#include "TLorentzVector.h"
#include "TComplex.h"
#include "TMatrixT.h"
#include "TVectorT.h"
#include "TMatrixTSym.h"
#include <string.h>
#include <vector> 
#include "rhoHelper.h"
#include "a1Helper.h"

using namespace std;
class TauPolInterface {
 public:
  TauPolInterface();
  TauPolInterface(vector<TLorentzVector> TauAndProd, string type);
  TauPolInterface(vector<TLorentzVector> TauAndProd1, string type1,vector<TLorentzVector> TauAndProd2, string type2);
  ~TauPolInterface();


  void Configure(vector<TLorentzVector> TauAndProd, string type );
  void ConfigurePair(vector<TLorentzVector> TauAndProd1, string type1,vector<TLorentzVector> TauAndProd2, string type2);
  bool  isConfigured();
  bool  isPairConfigured();
  TLorentzVector Boost(TLorentzVector pB, TLorentzVector frame);
  
  double getOmega(string which="");
  double getOmegabar(string which="");
  double getCombOmega();
  double getCombOmegaBar();

 private:

  void SetupLeg(string which="");


  double mrho;
  double mpi;
  double mtau;
  double ma1;
  bool debug;
  TMatrixT<double> convertToMatrix(TVectorT<double> V);

  TLorentzVector TauLV;
  TLorentzVector ProductLV;
  TLorentzVector TauRhoPi;
  TLorentzVector TauRhoPi0;
  TLorentzVector InvisibleLV;
  TLorentzVector DPF_TauLV;
  TLorentzVector DPF_TauRhoPi;
  TLorentzVector DPF_TauRhoPi0;
  TLorentzVector DPF_InvisibleLV;
  string type_;


  TLorentzVector TauLV1,TauLV2;
  TLorentzVector ProductLV1,ProductLV2;
  TLorentzVector TauRhoPi1,TauRhoPi2;
  TLorentzVector TauRhoPi01,TauRhoPi02;
  TLorentzVector InvisibleLV1,InvisibleLV2;
  TLorentzVector DPF_TauLV1,DPF_TauLV2;
  TLorentzVector DPF_TauRhoPi1,DPF_TauRhoPi2;
  TLorentzVector DPF_TauRhoPi01,DPF_TauRhoPi02;
  TLorentzVector DPF_InvisibleLV1,DPF_InvisibleLV2;
  string type1_,type2_;
};
#endif

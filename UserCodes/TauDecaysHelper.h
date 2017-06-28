#ifndef TauDecaysHelper_h
#define TauDecaysHelper_h
#include <vector>
#include "TLorentzVector.h"
#include "TComplex.h"
#include "TMatrixT.h"
#include "TVectorT.h"
#include "TMatrixTSym.h"
#include <string.h>
#include <vector>
#include "TLorentzVector.h"
using namespace std;
class TauDecaysHelper {
 public:
  TauDecaysHelper();
  TauDecaysHelper(vector<TLorentzVector> TauAndProd, string type);
  ~TauDecaysHelper();

  void Configure(vector<TLorentzVector> TauAndProd, string type );
  bool  isConfigured();
  TLorentzVector Boost(TLorentzVector pB, TLorentzVector frame);
  double getOmega();
  double getCosbetaRho();
  double getCosthetaRho();
 private:
  double mrho;
  double mpi;
  double mtau;
  double ma1;

  TMatrixT<double> convertToMatrix(TVectorT<double> V);
  TLorentzVector TauLV;
  TLorentzVector ProductLV;
  TLorentzVector TauRhoPi;
  TLorentzVector TauRhoPi0;
  TLorentzVector InvisibleLV;
  string type_;
};
#endif

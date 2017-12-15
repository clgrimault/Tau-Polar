// -*- C++ -*-
//
// 
/**\class rhoHelper.h rhoHelper.cc
 Description: 
*/
//
// Original Author:  Vladimir Cherepanov 
//         Created:  Mon Sep 4 13:49:02 CET 2017
//
//

#ifndef rhoHelper_h
#define rhoHelper_h


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


class rhoHelper {
 
 public:
  rhoHelper();
  rhoHelper(vector<TLorentzVector> TauRhoandProd, int helicity);
  rhoHelper(vector<TLorentzVector> TauRhoandProd, TLorentzVector RefernceFrame, int helicity);
  ~rhoHelper();
  void Configure(vector<TLorentzVector> TauRhoandProd, int helicity);
  void Configure(vector<TLorentzVector> TauRhoandProd, TLorentzVector RefernceFrame, int helicity);
  bool isConfigured();
  void Setup(vector<TLorentzVector> TauRhoandProd, TLorentzVector ReferenceFrame, int helicity );
  void Initialize(TLorentzVector t, TLorentzVector mu);
  std::vector<TLorentzVector> getBoosted(){return TauRhoandProd_RF;}
  TLorentzVector Boost(TLorentzVector pB, TLorentzVector frame);
  TVector3 Rotate(TVector3 LVec, TVector3 Rot);

  TVector3 RotationVector;
  //====================
  double getCosthetaRho(); 
  double getSinthetaRho();

  double getCosbetaRho();
  double getSinbetaRho();

  double getCosBetaTest();
  double getSinBetaTest();

  double TFK_cosbeta();
  double TFK_sinbeta();

  double TFK_costheta();
  double TFK_sintheta();

  bool toplot;


  double getUltrarel_cospsiLF();
  double getSinpsiLF();
  double DPF_cosalpha(); 
  double DPF_sinalpha(); 

  double DPF_cosalphaTest(); 
  double DPF_sinalphaTest(); 



  TVector3 nL();
  TVector3 nT();
  TVector3 nPerp();
  TVector3 ns();
  TVector3 nTRotated();
  TVector3 nPerpRotated();
  TLorentzVector sLV();


  double getOmegaRho();
  double getOmegaRhoBar();

 private:
  double mpi;
  double mpi0;
  double mtau;
  double coscab;
  double mrho;
  bool debug;
  int TauHelicity;


  vector<TLorentzVector> TauRhoandProd_RF;
  TLorentzVector TauLV;
  TLorentzVector TauRhoPi;
  TLorentzVector TauRhoPi0;
  TLorentzVector ProductLV;
  TLorentzVector DPF_TauRhoPi;
  TLorentzVector DPF_TauRhoPi0;
  TLorentzVector InvisibleLV;
  TLorentzVector DPF_TauLV;
  TLorentzVector DPF_InvisibleLV;
  bool isValid_;
  TMatrixT<double> convertToMatrix(TVectorT<double> V);
};
#endif

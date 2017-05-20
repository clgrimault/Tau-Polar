#ifndef a1Helper_h
#define a1Helper_h


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


class a1Helper {

 public:
  a1Helper(vector<TLorentzVector> TauA1andProd);
  a1Helper(vector<TLorentzVector> TauA1andProd, TLorentzVector RefernceFrame);
  ~a1Helper();
  void Configure(TLorentzVector OSPion, TLorentzVector SSPion1, TLorentzVector SSPion2,TLorentzVector TauA1, TLorentzVector TauMu );
  void Setup(vector<TLorentzVector> TauA1andProd, TLorentzVector ReferenceFrame );

  void Initialize(TLorentzVector t, TLorentzVector mu);
  bool OmegaIsValid(){return isValid_;}
  std::vector<TLorentzVector> getBoosted(){return TauA1andProd_RF;}


  void SetParametersReco(TLorentzVector tau, TLorentzVector mu );
  void SetFrame(TLorentzVector Vec );
  TLorentzVector Boost(TLorentzVector pB, TLorentzVector frame);

  double costheta();
  double costheta1();
  float CosBeta();
  double CosBeta1();
  std::vector<float> Sin2Cos2Gamma(TLorentzVector p1,TLorentzVector p2, TLorentzVector p3);
  float CosPsi();





  double lambda(double x, double y, double z);
  double Scalar(TLorentzVector p1, TLorentzVector p2);



  //--------------------------- Hadronic current ---------------------
  float WA(TLorentzVector s1,TLorentzVector s2,TLorentzVector s3,float QQ);
  float WC(TLorentzVector s1,TLorentzVector s2,TLorentzVector s3,float QQ);
  float WD(TLorentzVector s1,TLorentzVector s2,TLorentzVector s3,float QQ);
  float WE(TLorentzVector s1,TLorentzVector s2,TLorentzVector s3,float QQ);
  TComplex  F(float s1, float s2,float QQ);
  double VV1(double s1 ,double s2,double  s3, double Q);
  double VV2(double s1 ,double s2, double s3,double  Q);
  double V1V2(double s1, double s2, double s3, double Q);
  double h0(double s1, double s2, double s3, double Q);
  double h(double s1 ,double s2, double s3, double Q);


  TComplex  BreitWigner(double Q, string type="rho");
  TComplex  BRho(double Q);
  TComplex F1(double s1, double s2, double Q);
  TComplex F2(double s1, double s2, double Q);
  TComplex F4(double s1, double s2, double s3, double Q);
  TComplex   Conjugate(TComplex a);
  double  Widths(double Q, string type="rho");
  double ppi(double QQ);
  double ga1(double  Q);
  double Mass(string type="rho");


  TComplex  BWa1(float QQ);
  TComplex  BWrho(float QQ);
  TComplex  BWrhoPrime(float QQ);
  float GammaA1(float QQ);
  float gForGammaA1(float QQ);
  float GammaRho(float QQ);
  float  GammaRhoPrime(float QQ);



  double GetOmegaA1();



 private:
  double mpi;
  double mpi0;
  double mtau;
  double coscab;
  double mrho;
  double mrhoprime;
  double ma1;
  double mpiprime;
  double Gamma0rho; 
  double Gamma0rhoprime; 
  double Gamma0a1;
  double Gamma0piprime;
  double fpi;
  double fpiprime;
  double gpiprimerhopi;
  double grhopipi;
  double beta;

  const TLorentzVector a1pos;
  const TLorentzVector a1pss1;
  const TLorentzVector a1pss2;



  vector<TLorentzVector> TauA1andProd_RF;



  TLorentzVector Boost_;

  TLorentzVector KFitTau_;
  TLorentzVector RecoMuon_;
  TLorentzVector Tau1_;
  TLorentzVector Tau2_;

  TLorentzVector TauMu1_;
  TLorentzVector TauMu2_;

  bool Flag_;

  TLorentzVector Z_;
  TLorentzVector TauA1_;
  TLorentzVector TauMu_;


  TLorentzVector A1ZFrame_;
  TLorentzVector OSPionZFrame_;
  TLorentzVector SSPion1ZFrame_;
  TLorentzVector SSPion2ZFrame_;
  bool isValid_;

  TMatrixT<double> convertToMatrix(TVectorT<double> V);


};
#endif

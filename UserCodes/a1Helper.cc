#include "a1Helper.h"
#include <iostream>

a1Helper::a1Helper(vector<TLorentzVector> TauA1andProd){
  if(TauA1andProd.size()!=4){
    std::cout<<" Warning!! Size of input vector != 4 !! "<<std::endl;
  }
  TLorentzVector fakeboost(0,0,0,0);
  Setup(TauA1andProd,fakeboost);
}


a1Helper::a1Helper(vector<TLorentzVector> TauA1andProd, TLorentzVector RefernceFrame){
  if(TauA1andProd.size()!=4){
    std::cout<<" Warning!! Size of input vector != 4 !! "<<std::endl;
  }
  Setup(TauA1andProd,RefernceFrame);
}


void 
a1Helper::Setup(vector<TLorentzVector> TauA1andProd, TLorentzVector ReferenceFrame){
   mpi   = 0.13957018; // GeV 
   mpi0 = 0.1349766;   // GeV
   mtau = 1.776; // GeV
   coscab = 0.975; 
   mrho = 0.773; // GeV
   mrhoprime = 1.370; // GeV
   ma1 = 1.251; // GeV
   mpiprime = 1.300; // GeV
   Gamma0rho  =0.145; // GeV
   Gamma0rhoprime = 0.510; // GeV
   Gamma0a1 = 0.599; // GeV
   Gamma0piprime = 0.3; // GeV
   fpi= 0.093; // GeV
   fpiprime = 0.08; //GeV
   gpiprimerhopi = 5.8; //GeV
   grhopipi = 6.08;  //GeV
   beta = -0.145;

   for(int i=0; i<TauA1andProd.size(); i++){
     TauA1andProd_RF.push_back(Boost(TauA1andProd.at(i),ReferenceFrame));
   }
   LFosPionLV  = TauA1andProd.at(1);
   LFss1pionLV =TauA1andProd.at(2);
   LFss2pionLV =TauA1andProd.at(3);
   LFa1LV = LFosPionLV+LFss1pionLV+LFss2pionLV;
   LFtauLV = TauA1andProd.at(0);
   LFQ= LFa1LV.M();

   _osPionLV   = TauA1andProd_RF.at(1);
   _ss1pionLV =TauA1andProd_RF.at(2);
   _ss2pionLV =TauA1andProd_RF.at(3);
   _a1LV = _osPionLV+_ss1pionLV+_ss2pionLV;
   _tauLV = TauA1andProd_RF.at(0);
   _s12 = _ss1pionLV +_ss2pionLV;
   _s13 = _ss1pionLV + _osPionLV;
   _s23 = _ss2pionLV + _osPionLV;
   _s1  =  _s23.M2(); 
   _s2  =  _s13.M2();
   _s3  =  _s12.M2();
   _Q = _a1LV.M();
}

void 
a1Helper::Configure(TLorentzVector OSPion, TLorentzVector SSPion1, TLorentzVector SSPion2,TLorentzVector TauA1, TLorentzVector TauMu  ){


  Z_ =TauA1 + TauMu;
  TLorentzVector A1LabFrame = OSPion + SSPion1 + SSPion2;

   // OSPionZFrame_ = Boost(OSPion,Z_);
   // SSPion1ZFrame_= Boost(SSPion1,Z_);
   // SSPion2ZFrame_= Boost(SSPion2,Z_);
   // A1ZFrame_     = Boost(A1LabFrame,Z_);

   OSPionZFrame_ = OSPion;
   SSPion1ZFrame_= SSPion1;
   SSPion2ZFrame_= SSPion2;
   A1ZFrame_     = A1LabFrame;

   TauA1_ = TauA1;
  // std::cout<<"A1 Lab Frame  "<< A1LabFrame.Px() << "  " <<A1LabFrame.Py() << " "<< A1LabFrame.Pz() << "  " <<A1LabFrame.E()<<std::endl;
  // std::cout<<"A1 Z frame    "<< A1ZFrame_.Px() << "  "<< A1ZFrame_.Py() << " "<< A1ZFrame_.Pz() << "  " <<A1ZFrame_.E()<<std::endl;
  // std::cout<<"Z   "<< Z_.Px() << "  "<< Z_.Py() << " "<< Z_.Pz() << "  " <<Z_.E()<<std::endl;

}




void 
a1Helper::SetParametersReco(TLorentzVector tau, TLorentzVector mu ){
 Initialize(tau,mu);
}
void 
a1Helper::SetFrame(TLorentzVector vec){
  Boost_ = vec;
}



a1Helper::~a1Helper(){
}



void 
a1Helper::Initialize(TLorentzVector t, TLorentzVector mu){
  RecoMuon_=mu;
  KFitTau_=t;
}





double 
a1Helper::lambda(double x, double y, double z){
    return x*x +y*y +z*z - 2*x*y - 2*x*z - 2*z*y;
}
TLorentzVector 
a1Helper::Boost(TLorentzVector pB, TLorentzVector frame){
   TMatrixT<double> transform(4,4);
   TMatrixT<double> result(4,1);
   TVectorT<double> vec(4); 
   TVector3 b;
   if(frame.Vect().Mag()==0){ std::cout<<" Boost is not set, perfrom calculation in the Lab Frame   "<<std::endl; return pB;}
    if(frame.E()==0){ std::cout<<" Caution: Please check that you perform boost correctly!  " <<std::endl; return pB;} 
   else   b=frame.Vect()*(1/frame.E());
   vec(0)  = pB.E();    vec(1)  = pB.Px();
   vec(2)  = pB.Py();  vec(3)  = pB.Pz();
   double gamma  = 1/sqrt( 1 - b.Mag2());
   transform(0,0)=gamma; transform(0,1) =- gamma*b.X() ;  transform(0,2) =  - gamma*b.Y();  transform(0,3) = - gamma*b.Z(); 
   transform(1,0)=-gamma*b.X(); transform(1,1) =(1+ (gamma-1)*b.X()*b.X()/b.Mag2()) ;  transform(1,2) = ((gamma-1)*b.X()*b.Y()/b.Mag2());  transform(1,3) = ((gamma-1)*b.X()*b.Z()/b.Mag2());
   transform(2,0)=-gamma*b.Y(); transform(2,1) = ((gamma-1)*b.Y()*b.X()/b.Mag2());  transform(2,2) = (1 + (gamma-1)*b.Y()*b.Y()/b.Mag2());  transform(2,3) =  ((gamma-1)*b.Y()*b.Z()/b.Mag2()); 
   transform(3,0)=-gamma*b.Z(); transform(3,1) =((gamma-1)*b.Z()*b.X()/b.Mag2()) ;  transform(3,2) = ((gamma-1)*b.Z()*b.Y()/b.Mag2());  transform(3,3) = (1 + (gamma-1)*b.Z()*b.Z()/b.Mag2()); 
   result=transform*convertToMatrix(vec);
   return TLorentzVector(result(1,0), result(2,0) ,result(3,0), result(0,0));
}
double 
a1Helper::Scalar(TLorentzVector p1, TLorentzVector p2){
    return p1.Vect()*p2.Vect();
}

//---------------------------------------  hadronic current ---------------------------
double 
a1Helper::WA(){

//   float SS1 = (s2+s3).M()*(s2+s3).M();
//   float SS2 = (s1+s3).M()*(s1+s3).M();
//   float SS3 = (s1+s2).M()*(s1+s2).M();
  // std::cout<<  "WA:   VV1()   " << VV1() <<std::endl; 
  // std::cout<<  "WA:    F1().Rho2()  " << F1().Rho2()<<std::endl; 
  // std::cout<<  "WA:   VV2()   " << VV2()<<std::endl; 
  // std::cout<<  "WA:    F2().Rho2()  " << F2().Rho2()<<std::endl; 
  // std::cout<<  "WA:    V1V2()  " << V1V2()<<std::endl; 
  // std::cout<<  "WA:     ( F1()*Conjugate(F2()) ).Re() " << ( F1()*Conjugate(F2()) ).Re()<<std::endl; 
  // std::cout<<  "WA:  ( Conjugate(F1())*F2() ).Re()    " <<  ( Conjugate(F1())*F2() ).Re()  <<std::endl; 

  return  VV1()*F1().Rho2() + VV2()*F2().Rho2()  + 2*V1V2()*( F1()*Conjugate(F2()) ).Re();
 }

 double 
a1Helper::WC(){
   return  -(-VV1() + 2*h() )*F1().Rho2() - (-VV2() + 2*h())*F2().Rho2()   -   (-2*V1V2() - 4*h())*( F1()*Conjugate(F2()) ).Re();
 } 

 double
 a1Helper::WD(){
   double QQ = _Q*_Q;
   double undersqrt1 = VV1()  -h();
   double undersqrt2 = VV2()  -h();

//   if(undersqrt1 < 0) undersqrt1 =0;
//   if(undersqrt2 < 0) undersqrt2 =0;

   return  -sqrt(h()) * ( 2 * sqrt(undersqrt1) * F1().Rho2() - 2*sqrt(undersqrt2)*F2().Rho2()  
				     + (QQ - mpi*mpi + _s3)*(_s1 - _s2 )*( F1()*Conjugate(F2()) ).Re()/QQ/sqrt(h0() ) );
 }

 double
 a1Helper::WE(){
  return  3*sqrt(h()*h0())*( F1()*Conjugate(F2()) ).Im();
 }

double
 a1Helper::WSA(){
  double QQ = _Q*_Q;
  return  QQ*F4().Rho2();
 }
double
 a1Helper::WSB(){
  //  double QQ = _Q*_Q;
   double undersqrt1 = VV1()  -h();
   double undersqrt2 = VV2()  -h();
   return  -2*_Q* (sqrt(undersqrt1) * (F1()*Conjugate(F4())).Re() +   sqrt(undersqrt2)*(F2()*Conjugate(F4())).Re()  );
 }
double
 a1Helper::WSD(){
  double QQ = _Q*_Q;
  return  2*sqrt(QQ*h())* ( (F1()*Conjugate(F4())).Re() - (F2()*Conjugate(F4())).Re()   );
 }
double
 a1Helper::WSC(){
  //  double QQ = _Q*_Q;
   double undersqrt1 = VV1()  -h();
   double undersqrt2 = VV2()  -h();
   return  2*_Q* (sqrt(undersqrt1) * (F1()*Conjugate(F4())).Im() +   sqrt(undersqrt2)*(F2()*Conjugate(F4())).Im()  );
 }
double
 a1Helper::WSE(){
  double QQ = _Q*_Q;
   return  -2*sqrt(QQ*h())* ( (F1()*Conjugate(F4())).Im() - (F2()*Conjugate(F4())).Im()   );
 }



double
a1Helper::cosgammaLF(){
  double QQ=LFQ*LFQ;
  // double B1 = (pow(_ss1pionLV.E()*_tauLV.E()   - _ss1pionLV.Vect().Dot(_a1LV.Vect()),2 ) - QQ*mpi*mpi)/QQ;
  // double B2 = (pow(_ss2pionLV.E()*_tauLV.E()   - _ss2pionLV.Vect().Dot(_a1LV.Vect()),2 ) - QQ*mpi*mpi)/QQ;
  double B3 = (pow(LFosPionLV.E()*LFtauLV.E()   - LFosPionLV.Vect().Dot(LFa1LV.Vect()),2 ) - QQ*mpi*mpi)/QQ;

  // double T = 0.5*sqrt(-lambda(B1,B2,B3));
  // double A1=(_a1LV.E()*_a1LV.Vect().Dot(_ss1pionLV.Vect()) - _ss1pionLV.E()*_a1LV.P()*_a1LV.P())/QQ;
  // double A2=(_a1LV.E()*_a1LV.Vect().Dot(_ss2pionLV.Vect()) - _ss2pionLV.E()*_a1LV.P()*_a1LV.P())/QQ;
  double A3=(LFa1LV.E() * LFa1LV.Vect().Dot(LFosPionLV.Vect()) - LFosPionLV.E()*LFa1LV.P()*LFa1LV.P()) / LFQ;
  std::cout<<"sqrt B3 "<< sqrt(B3)<<std::endl;
  std::cout<<"A3 "<< A3<<std::endl;

  std::cout<< "fuck 1  "   <<LFosPionLV.E()*LFtauLV.E()   - LFosPionLV.Vect().Dot(LFa1LV.Vect())<<std::endl;
  std::cout<< "fuck 2  "   <<LFa1LV.E() * LFa1LV.Vect().Dot(LFosPionLV.Vect()) - LFosPionLV.E()*LFa1LV.P()*LFa1LV.P()<<std::endl;
                                         
  std::cout<< "QQ "   << LFQ*LFQ<< "  _Q_Q  "<< _Q*_Q << std::endl;

  if(B3<=0 || cosbetaLF() >=1){std::cout<<"Warning! In a1Helper::cosgamma square root <=0! return 0"<<std::endl; return 0;}
  return A3/LFa1LV.P()/sqrt(B3)/sqrt(1 - cosbetaLF()*cosbetaLF());
}

double
a1Helper::singammaLF(){
  double QQ=LFQ*LFQ;
   double B1 = (pow(LFss1pionLV.E()*LFa1LV.E()   - LFss1pionLV.Vect().Dot(LFa1LV.Vect()),2 ) - QQ*mpi*mpi)/QQ;
  double B2 = (pow(LFss2pionLV.E()*LFa1LV.E()   - LFss2pionLV.Vect().Dot(LFa1LV.Vect()),2 ) - QQ*mpi*mpi)/QQ;
   double B3 = (pow(LFosPionLV.E()*LFa1LV.E()   - LFosPionLV.Vect().Dot(LFa1LV.Vect()),2 ) - QQ*mpi*mpi)/QQ;

  double T = 0.5*sqrt(-lambda(B1,B2,B3));

  double A1=(LFa1LV.E()*LFa1LV.Vect().Dot(LFss1pionLV.Vect()) - LFss1pionLV.E()*LFa1LV.P()*LFa1LV.P())/QQ;
  //  double A2=(_a1LV.E()*_a1LV.Vect().Dot(_ss2pionLV.Vect()) - _ss2pionLV.E()*_a1LV.P()*_a1LV.P())/QQ;
  double A3=(LFa1LV.E()*LFa1LV.Vect().Dot(LFosPionLV.Vect()) - LFosPionLV.E()*LFa1LV.P()*LFa1LV.P())/QQ;

  if(A3==0 || T==0){std::cout<<"Warning! In a1Helper::singamma denominator ==0! return 0"<<std::endl; return 0;}
  double scale = -(B3*A1/A3 - 0.5*(B2 - B1 - B3))/T;
  std::cout<<"scale  " << scale <<std::endl;
  return cosgammaLF()*scale;
}
double
a1Helper::cos2gamma(){
   return singamma()*singamma()   -     cosgamma()*cosgamma();
}

double
a1Helper::sin2gamma(){
  return 2*singamma()*cosgamma();
}
double 
a1Helper::cospsiLF(){
  double QQ = LFQ*LFQ;
  double s = 4*LFtauLV.E()*LFtauLV.E();
  double x = 2*LFa1LV.E()/sqrt(s);
  if(x*x  - 4*QQ/s <= 0 ){std::cout<<"Warning! In a1Helper::cospsi root square <=0! return 0"<<std::endl; return 0;}
  return    ( x*(mtau*mtau + QQ)  - 2*QQ  )   /   ( mtau*mtau  - QQ   ) / sqrt(x*x  - 4*QQ/s); 
}
double 
a1Helper::sinpsiLF(){
  if(cospsiLF()*cospsiLF() > 1  )std::cout<<"Warning! In a1Helper::sinpsi root square <=0! return nan"<<std::endl;
  return    sqrt(1 - cospsiLF()*cospsiLF());
}

double 
a1Helper::ultrarel_cospsiLF(){
  double QQ = LFQ*LFQ;
  double cos = (costhetaLF()*(mtau*mtau  + QQ)   + (mtau*mtau  - QQ))/(costhetaLF()*(mtau*mtau  - QQ)   + (mtau*mtau  + QQ));
  return cos;
}

double 
a1Helper::costhetaLF(){
  double QQ = LFQ*LFQ;
  double x = LFa1LV.E()/LFtauLV.E();
  double s = 4*LFtauLV.E()*LFtauLV.E();
  if( 1 - 4*mtau*mtau/s  <= 0 ){std::cout<<"Warning! In a1Helper::costheta root square <=0! return 0"<<std::endl; return 0;}
  return (2*x*mtau*mtau - mtau*mtau - QQ)/((mtau*mtau - QQ)*sqrt(1 - 4*mtau*mtau/s));
}
double 
a1Helper::sinthetaLF(){
  if( costhetaLF()*costhetaLF() > 1 ) std::cout<<"Warning! In a1Helper::sin heta root square <=0! return nan"<<std::endl; 
  return sqrt(1- costhetaLF()*costhetaLF());
}

double 
a1Helper::cosbetaLF(){
  double QQ = LFQ*LFQ;
  double B1 = (pow(LFss1pionLV.E()*LFa1LV.E()   - Scalar(LFss1pionLV,LFa1LV),2 ) - QQ*mpi*mpi)/QQ;
  double B2 = (pow(LFss2pionLV.E()*LFa1LV.E()   - Scalar(LFss2pionLV,LFa1LV),2 ) - QQ*mpi*mpi)/QQ;
  double B3 = (pow(LFosPionLV.E()*LFa1LV.E()   -   Scalar(LFosPionLV,LFa1LV),2 ) - QQ*mpi*mpi)/QQ;

  TVector3 ss1pionVect = LFss1pionLV.Vect();
  TVector3 ss2pionVect = LFss2pionLV.Vect();
  TVector3 ospionVect = LFosPionLV.Vect();
  float T = 0.5*sqrt(-lambda(B1,B2,B3));
  // std::cout<<" B1  " << B1 <<std::endl;
  // std::cout<<" B2  " << B2 <<std::endl;
  // std::cout<<" B3  " << B3 <<std::endl;
  // std::cout<<"-lambda(B1,B2,B3)" << -lambda(B1,B2,B3) <<std::endl;
  // std::cout<<"  T  " << T <<std::endl;
  if(T==0 || LFa1LV.P()==0){std::cout<<" Warning!  Can not compute cosbetaLF, denominator =0; return 0; "<<std::endl; return 0;}
  return ospionVect.Dot(ss1pionVect.Cross(ss2pionVect)) /LFa1LV.P()/T;
}

double
a1Helper::VV1(){ //  this is -V1^{2}
  double QQ = _Q*_Q;
  return  _s2 - 4*mpi*mpi + pow(_s3 - _s1,2)/4/QQ;
}

double
a1Helper::VV2(){ //  this is -V2^{2}
  double QQ = _Q*_Q;
  return  _s1 - 4*mpi*mpi + pow(_s3 - _s2,2)/4/QQ;
}

double
a1Helper::V1V2(){  // this is -V1V2
  double QQ = _Q*_Q;
  return  (QQ/2 - _s3 - 0.5*mpi*mpi) + (_s3 - _s1)*(_s3 - _s2)/4/QQ;
}


double
a1Helper::h0(){ // this is -3sqrt{h0}/2
  double QQ = _Q*_Q;
  return -4*mpi*mpi + pow(2*mpi*mpi - _s1 - _s2,2)/QQ;
}

double
a1Helper::h(){
  double QQ = _Q*_Q;
  return (_s1*_s2*_s3 - mpi*mpi*pow(QQ - mpi*mpi,2))/h0()/QQ;  // this is sqrt{h}
}



TComplex 
a1Helper::F1(){
  TComplex scale(0, -2*sqrt(2)/3/fpi);
  TComplex res = scale*BreitWigner(_Q,"a1")*BRho(sqrt(_s2));
  //  std::cout<<"  BreitWigner(_Q,a1)  " << BreitWigner(_Q,"a1") << " BRho(_s2)  " << BRho(sqrt(_s2))<< std::endl;
  return res;
}


TComplex 
a1Helper::F2(){
  TComplex scale(0, -2*sqrt(2)/3/fpi);
  TComplex res = scale*BreitWigner(_Q,"a1")*BRho(sqrt(_s1));
  return res;
}

TComplex 
a1Helper::F4(){
  TComplex scale(0, -gpiprimerhopi*grhopipi*fpiprime/2/pow(mrho,4)/pow(mpiprime,3));
  TComplex res = scale*BreitWigner(_Q,"piprime")*(_s1*(_s2-_s3)*BRho(sqrt(_s1)) + _s2*(_s1-_s3)*BRho(sqrt(_s2)));
  return res;
}


TComplex 
a1Helper::BRho(double Q){
  //  std::cout<<"BRho:      BreitWigner(Q) " << BreitWigner(Q) << " BreitWigner(Q,rhoprime) " << BreitWigner(Q,"rhoprime")<< std::endl;
  return (BreitWigner(Q) + beta*BreitWigner(Q,"rhoprime"))/(1+beta);
}

TComplex 
a1Helper::BreitWigner(double Q, string type){
  double QQ=Q*Q;
  double re,im;
  double m = Mass(type);
  double g  = Widths(Q,type);
  re = (m*m*pow(m*m - QQ,2))/(pow(m*m - QQ,2) + m*m*g*g);
  im = m*m*m*g/(pow(m*m - QQ,2) + m*m*g*g);
  TComplex out(re,im);
  return out;
}

double
a1Helper::Widths(double Q, string type){
  double QQ = Q*Q;
  double Gamma;
  Gamma = Gamma0rho*mrho*pow( ppi(QQ)  / ppi(mrho*mrho), 3) /sqrt(QQ);
  if(type == "rhoprime"){
    Gamma=Gamma0rhoprime*QQ/mrhoprime/mrhoprime;
 }
  if(type == "a1"){
    Gamma=Gamma0a1*ga1(Q)/ga1(ma1);
 }
  if(type == "piprime"){
    Gamma = Gamma0piprime*pow( sqrt(QQ)/mpiprime  ,5)*pow( (1-mrho*mrho/QQ)/(1-mrho*mrho/mpiprime/mpiprime) ,3);
  }
  //  std::cout<< " Widths :   type   " << type << " Gamma  " << Gamma << "  QQ  "<< QQ <<std::endl;
  return Gamma;
}
double a1Helper::ga1(double  Q){
  double QQ = Q*Q;
  return (QQ > pow(mrho + mpi,2)) ?  QQ*(1.623 + 10.38/QQ - 9.32/QQ/QQ   + 0.65/QQ/QQ/QQ)  : 4.1*pow(QQ - 9*mpi*mpi,3)*(  1 - 3.3*(QQ - 9*mpi*mpi)  + 5.8*pow(QQ - 9*mpi*mpi,2)  );
}
double
a1Helper::Mass(string type){
  double m = mrho;
  if(type == "rhoprime") return mrhoprime; 
  if(type == "a1") return ma1;
  if(type == "piprime") return mpiprime;
  //std::cout<< "  type   " << type << " Mass  " << std::endl;
  return m;
}


double a1Helper::ppi(double QQ){  if(QQ < 4*mpi*mpi) std::cout<<"Warning! Can not compute ppi(Q); root square <0 ; return nan  "; return 0.5*sqrt(QQ - 4*mpi*mpi);}


double a1Helper::getf(){
  double QQ=_Q*_Q;
  double l  = 0.5*(mtau*mtau + QQ)/sqrt(QQ);
  //  double l0= 0.5*(mtau*mtau - QQ)/sqrt(QQ);
  //------ separate for simple debugging 
  double line1 =   -2*l   *   ( 2*WA()/3   + 0.5*(3*cospsiLF()*cospsiLF()   -1)  *  ( WA()*(3*cosbeta()*cosbeta() -1 )/6    - 0.5*WC()*sinbeta()*sinbeta()*cos2gamma()   + 0.5*WD()* sinbeta()*sinbeta()* sin2gamma() )   )/sqrt(QQ);
  double line2 = mtau*mtau*WA()/QQ + mtau*mtau  *  (  WSA() +  cospsiLF()*sinbeta()*(   WSB() *cosgamma()      - WSD() * singamma())     )/QQ + WE()*cosbeta()*cospsiLF();
  double res = line1+ line2;

  // std::cout<< "f:  line 1   " << line1 <<std::endl;
  // std::cout<< "f:   line 2   " << line2 <<std::endl;

  return res;
}
double a1Helper::getg(){
  double QQ=_Q*_Q;
  //  double l  = 0.5*(mtau*mtau + QQ)/sqrt(QQ);
  double l0= 0.5*(mtau*mtau - QQ)/sqrt(QQ);
  //------ separate for simple debugging 
  double line1 =   -2*l0   * costhetaLF()*  ( 2*WA()/3   + 0.5*(3*cospsiLF()*cospsiLF()   -1)  *  ( WA()*(3*cosbeta()*cosbeta() -1 )/6    -  0.5*WC()*sinbeta()*sinbeta()*cos2gamma()   + 0.5*WD()* sinbeta()*sinbeta()* sin2gamma() )   )/sqrt(QQ);
  double line2 = mtau*mtau*WA()*costhetaLF()/QQ  +    sqrt(mtau*mtau/QQ )  * sinthetaLF()* ( 0.5*WA()*2* sinbeta()*cosbeta()* cosalpha() -  
                                             WC()*sinbeta()* (sinalpha()* sin2gamma() + cos2gamma()* cosalpha()*    cosbeta() )    -    WD()*sinbeta()*( sinalpha()*cos2gamma() + sin2gamma()* cosalpha()*cosbeta()  )- 2*cospsiLF()*sinpsiLF()  *
				            (WA()*(3*cosbeta() *cosbeta() -1 )/6   -    0.5*WC()*sinbeta()* sinbeta()* cos2gamma()+ 0.5*WD()*sinbeta()* sinbeta()* cos2gamma() + WD()*sinbeta()* sinbeta()* sin2gamma())/3   );

  double line3  =  sqrt(mtau*mtau/QQ ) *sinthetaLF()* (WE()*(cosbeta()*sinpsiLF() + sinbeta()*cosalpha()) +cosbeta()*sinalpha()*(WSC()*cosgamma() - WSE()*singamma()) + cosalpha()*(WSC()*singamma() + WSE()*cosgamma()));
  double line4  =  -WE()*costhetaLF()*cosbeta()*cospsiLF() + mtau*mtau*costhetaLF()*(WSA() + cospsiLF()*sinbeta()  * (WSB()*cosgamma() - WSD()* singamma()  ) )/QQ;
  double line5  =  sqrt(mtau*mtau/QQ)*sinthetaLF() *  ( sinpsiLF()*sinbeta()*( WSB()* cosgamma() - WSD()* singamma()) + cosbeta()*cosalpha()*(WSD()*singamma() - WSB()*cosgamma()  ) + 
							sinalpha()*(WSD()*cosgamma() + WSB()*singamma())          );
  double res = line1+ line2 + line3 + line4 + line5;
  // std::cout<< "g:  line 1   " << line1 <<std::endl;
  // std::cout<< "g:   line 2   " << line2 <<std::endl;
  // std::cout<< "g:   line 3   " << line3 <<std::endl;
  // std::cout<< "g:   line 4   " << line4 <<std::endl;
  // std::cout<< "g:   line 5   " << line5 <<std::endl;
  return res;
}



// double
// a1Helper::GetOmegaA1(){
//         isValid_ = false;
//         double omega(-999.);
// 	TLorentzVector pi1 = SSPion1ZFrame_;
// 	TLorentzVector pi2 = SSPion2ZFrame_;
// 	TLorentzVector pi3 = OSPionZFrame_;
// 	TLorentzVector a = A1ZFrame_;
// 	float mtau = 1.777;
// 	float cospsi  = CosPsi();
// 	float sinpsi  = sqrt(1 - cospsi*cospsi);
// 	float sin2psi = 2*sinpsi*cospsi;
	
// 	float sin2gamma = Sin2Cos2Gamma(pi1,pi2,pi3).at(0);
// 	float cos2gamma = Sin2Cos2Gamma(pi1,pi2,pi3).at(1);
	
	
// 	float cosbeta=CosBeta();
// 	float sinbeta = sqrt(1  - cosbeta*cosbeta);
	
	  
// 	float cstheta=costheta();
// 	float sintheta = sqrt(1 - cstheta*cstheta);

// 	double RR  = mtau*mtau/a.M()/a.M();
// 	float U = 0.5*(3*cospsi*cospsi - 1)*(1 - RR);
// 	float V = 0.5*(3*cospsi*cospsi - 1)*(1 + RR)*cstheta + 0.5*3*sin2psi*sintheta*sqrt(RR);
// 	//	std::cout<< " cospsi  "<< cospsi <<" costheta  "<<costheta <<" sin2psi "<<sin2psi <<" sintheta  "<<sintheta  <<" RR  "<< RR <<std::endl;

// 	float Wa =WA(pi1,pi2,pi3,a.M()*a.M());
// 	float Wc =WC(pi1,pi2,pi3,a.M()*a.M());
// 	float Wd =WD(pi1,pi2,pi3,a.M()*a.M());
// 	float We =WE(pi1,pi2,pi3,a.M()*a.M());

// 	float fa1 = (2  + RR + 0.5*(3*cosbeta*cosbeta - 1)*U)*Wa/3 - 0.5*sinbeta*sinbeta*cos2gamma*U*Wc + 0.5*sinbeta*sinbeta*sin2gamma*U*Wd + cospsi*cosbeta*We;
// 	float ga1 = (cstheta*(RR -2) - 0.5*(3*cosbeta*cosbeta - 1)*V)*Wa/3 + 0.5*sinbeta*sinbeta*cos2gamma*V*Wc - 0.5*sinbeta*sinbeta*sin2gamma*V*Wd -cosbeta*(cstheta*cospsi + sintheta*sinpsi*sqrt(RR))*We;
// 	//	std::cout<< " U  "<< U <<" V  "<<V <<" Wa "<<Wa <<" Wc  "<<Wc  <<" Wd  "<< Wd <<" We  "<< We  <<" f "<< fa1 <<" g "<< ga1 <<std::endl;
// 	omega = ga1/fa1;
// 	if(omega > 0 or omega < 0) isValid_ = true;
// 	return omega;
// }


TVector3
a1Helper::nPerp(){
  if(_ss1pionLV.Vect().Cross(_ss2pionLV.Vect()).Mag()==0){ std::cout<<"  Can not return nPerp, same sign pions seem to be parallel in a1 rest frame, return 0,0,0  "<<std::endl; return TVector3(0,0,0);}

  TVector3 nss1= _ss1pionLV.Vect()*(1/_ss1pionLV.Vect().Mag());
  TVector3 nss2= _ss2pionLV.Vect()*(1/_ss2pionLV.Vect().Mag());
  return   (nss1.Cross(nss2))*(1/(nss1.Cross(nss2)).Mag());
}
TVector3
a1Helper::nL(){
  return   -LFa1LV.Vect()*(1/LFa1LV.Vect().Mag());
}
TVector3
a1Helper::nT(){
  return   _tauLV.Vect()*(1/_tauLV.Vect().Mag());
}

double a1Helper::cosalpha(){
   TVector3 nLCrossnT  = nL().Cross(nT());
  TVector3 nLCrossnPerp  = nL().Cross(nPerp());

  if(nLCrossnPerp.Mag() ==0 || nLCrossnT.Mag() ==0){std::cout<<" Can not compute cos alpha, one denominator is 0, return cos alpha =0  "<< std::endl; return 0;}
  return nLCrossnT.Dot(nLCrossnPerp)/nLCrossnT.Mag()/nLCrossnPerp.Mag();
}
double a1Helper::sinalpha(){
  TVector3 nLCrossnT  = nL().Cross(nT());
  TVector3 nLCrossnPerp  = nL().Cross(nPerp());
  if(nLCrossnPerp.Mag() ==0 || nLCrossnT.Mag() ==0){std::cout<<" Can not compute sin alpha, one denominator is 0, return sin alpha =0  "<< std::endl; return 0;}
  return -nT().Dot(nLCrossnPerp)/nLCrossnT.Mag()/nLCrossnPerp.Mag();
}
double a1Helper::cosbeta(){
  return nL().Dot(nPerp());
}
double a1Helper::sinbeta(){
  if(cosbeta()*cosbeta() > 1 ){std::cout<<"Warning! Can not compute sin beta! return 0"<<std::endl; return 0;}
  return sqrt(1 - cosbeta()*cosbeta());
}

double a1Helper::cosgamma(){
  TVector3 nLCrossnPerp  = nL().Cross(nPerp());

  TVector3 qvect = _osPionLV.Vect()*(1/_osPionLV.Vect().Mag());
  qvect.Print();
  if(nLCrossnPerp.Mag()==0) { std::cout<<"Warning! Can not compute cos gamma, denominator =0, return 0  "<< std::endl; return 0; }
  return -nL()*qvect/nLCrossnPerp.Mag();
}

double a1Helper::singamma(){
  TVector3 nLCrossnPerp  = nL().Cross(nPerp());
  TVector3 qvect = _osPionLV.Vect()*(1/_osPionLV.Vect().Mag());

  if(nLCrossnPerp.Mag()==0) { std::cout<<"Warning! Can not compute cos gamma, denominator =0, return 0  "<< std::endl; return 0; }
  return qvect*nLCrossnPerp/nLCrossnPerp.Mag();
}



TComplex 
a1Helper::Conjugate(TComplex a){
  return TComplex(a.Re(), -a.Im());
}
TMatrixT<double> a1Helper::convertToMatrix(TVectorT<double> V){
  TMatrixT<double> M(V.GetNrows(),1);
  for(int i=0; i < M.GetNrows(); i++){
    M(i,0)=V(i);
  } return M;
}
//  cosb = nL*nPerp
// cosgamma = -nL*q3/|nLCrossnPerp|
// cosgamma = q3*nLCrossnPer/|nLCrossnPerp|
// double 
// a1Helper::CosBeta1(){
//   TLorentzVector p1 = SSPion1ZFrame_;
//   TLorentzVector p2 = SSPion2ZFrame_;
//   TLorentzVector p3 = OSPionZFrame_;


//   //std::cout<<"  cosbeta1  ================================= "<<std::endl;
//   double mpi  = 0.139;
// //   double E = p1.E() +  p2.E() +  p3.E(); 
// //   double P = p1.P() +  p2.P() +  p3.P(); 
// //   double QQ = E*E - P*P;

//   TLorentzVector a1 = p1+p2+p3;
//   //  float P = 
//   TLorentzVector s12 = p1+p2;
//   TLorentzVector s13 = p1+p3;
//   TLorentzVector s23 = p2+p3;

//   double QQ = a1.E()*a1.E() - a1.P()*a1.P();


//   TLorentzVector p1Timesp2(p1.Py()*p2.Pz() - p1.Pz()*p2.Py(),p1.Pz()*p2.Px() - p1.Px()*p2.Pz(),p1.Px()*p2.Py() - p1.Py()*p2.Px(),1);
//   float mm=a1.M()*a1.M();
//   float mm12=s12.M()*s12.M();
//   float mm13=s13.M()*s13.M();
//   float mm23=s23.M()*s23.M();
//   float mmpi=mpi*mpi;

//   float l1  = lambda( mm, mm12 , mmpi);
//   float l2  = lambda( mm, mm13 , mmpi);
//   float l3  = lambda( mm, mm23 , mmpi);

 
//   double cbeta = /*8*a1.M()*a1.M()**/Scalar(p3,p1Timesp2)*a1.P()/sqrt(-lambda(l1,l2,l3));

// //   std::cout<<"  QQ  "<< QQ<<std::endl;
// //   std::cout<<"  mm12  "<< mm12<<std::endl;
// //   std::cout<<"  mm13  "<< mm13<<std::endl;
// //   std::cout<<"  mm23  "<< mm23<<std::endl;
// //   std::cout<<"  mm  "<< mm<<std::endl;

// //   std::cout<<"  lambda1 = "<<l1<<std::endl;
// //   std::cout<<"  lambda2 = "<<l2<<std::endl;
// //   std::cout<<"  lambda3 = "<<l3<<std::endl;
// //   std::cout<<"  lambda  = "<<-lambda(l1,l2,l3)<<std::endl;

  
// //   std::cout<<"  Scalar(p3,p1Timesp2)*a1.P()  "<< Scalar(p3,p1Timesp2)<<std::endl;
// //   std::cout<<"  a1.P()  "<< a1.P() <<std::endl;

// //   std::cout<<"  a1.P()  "<< a1.P() << "   *  "<</*8*a1.M()*a1.M()**/Scalar(p3,p1Timesp2)/sqrt(-lambda(l1,l2,l3)) <<std::endl;



// //   std::cout<<"  cbeta  "<< cbeta  <<std::endl;
//   return cbeta;
// }



// std::vector<float> 
// a1Helper::Sin2Cos2Gamma(TLorentzVector p1,TLorentzVector p2, TLorentzVector p3){

//   std::vector<float> sin2cos2;
//   float mpi  = 0.139;
//   TLorentzVector a1 = p1+p2+p3;
//   float QQ = a1.E()*a1.E() - a1.P()*a1.P();

//   float B1 = (pow(p1.E()*a1.E()   - Scalar(p1,a1),2 ) - QQ*mpi*mpi)/QQ;
//   float B2 = (pow(p2.E()*a1.E()   - Scalar(p2,a1),2 ) - QQ*mpi*mpi)/QQ;
//   float B3 = (pow(p3.E()*a1.E()   - Scalar(p3,a1),2 ) - QQ*mpi*mpi)/QQ;

//   float T = 0.5*sqrt(-lambda(B1,B2,B3));

//   float A1=(a1.E()*Scalar(a1,p1) - p1.E()*a1.P()*a1.P())/QQ;
//   float A2=(a1.E()*Scalar(a1,p2) - p2.E()*a1.P()*a1.P())/QQ;
//   float A3=(a1.E()*Scalar(a1,p3) - p3.E()*a1.P()*a1.P())/QQ;


//   float cosgamma = A3/a1.P()/sqrt(B3)/sqrt(1 - CosBeta()*CosBeta());
//   float singamma = -cosgamma*(B3*A1/A3 - 0.5*(B2 - B1 - B3))/T;

//   sin2cos2.push_back(2*singamma*cosgamma);
//   sin2cos2.push_back(2*cosgamma*cosgamma - 1);
//   return sin2cos2;

// }



// float 
// a1Helper::CosPsi(){

//   TLorentzVector p1 = SSPion1ZFrame_;
//   TLorentzVector p2 = SSPion2ZFrame_;
//   TLorentzVector p3 = OSPionZFrame_;
//   TLorentzVector Z=Z_;

//   float mtau =1.777;
//   TLorentzVector a1 = p1+p2+p3;
//   float QQ = a1.E()*a1.E() - a1.P()*a1.P();
//   float cos = (costheta()*(mtau*mtau  + QQ)   + (mtau*mtau  - QQ))/(costheta()*(mtau*mtau  - QQ)   + (mtau*mtau  + QQ));
//   return cos;

// }


// double 
// a1Helper::costheta(){
  
//   TLorentzVector p1 = SSPion1ZFrame_;
//   TLorentzVector p2 = SSPion2ZFrame_;
//   TLorentzVector p3 = OSPionZFrame_;
//   TLorentzVector Z=Z_;

  
//   double zmass = Z.M();
//   double mtau = 1.777;
//   TLorentzVector a1 = p1+p2+p3;
//   float QQ = a1.E()*a1.E() - a1.P()*a1.P();

//   double x = a1.E()/TauA1_.E();
//   double ctheta = (2*x*mtau*mtau - mtau*mtau - QQ)/((mtau*mtau - QQ)*sqrt(1 - 4*mtau*mtau/zmass/zmass));
//   // std::cout<<"p1 "<< p1.Px() << "  " <<p1.Py() << " "<< p1.Pz() << "  " <<p1.M()<<std::endl;
//   // std::cout<<"p2 "<< p2.Px() << "  "<< p2.Py() << " "<< p2.Pz() << "  " <<p2.M()<<std::endl;
//   // std::cout<<"p3 "<< p3.Px() << "  "<< p3.Py() << " "<< p3.Pz() << "  " <<p3.M()<<std::endl;
//   // std::cout<<"a1 "<< a1.Px() << "  "<< a1.Py() << " "<< a1.Pz() << "  " <<a1.M()<<"  " <<a1.M()*a1.M() <<std::endl;
//   // std::cout<<"QQ "<< QQ << " x "<< x <<" zmass  " <<zmass <<std::endl;

//   //  std::cout<<"2*x*mtau*mtau/mtau*mtau - QQ "<< 2*x*mtau*mtau/(mtau*mtau - QQ)<< "  "<< x <<std::endl;
//   return ctheta;
// }


// double 
// a1Helper::costheta1(){
  
//   TLorentzVector p1 = SSPion1ZFrame_;
//   TLorentzVector p2 = SSPion2ZFrame_;
//   TLorentzVector p3 = OSPionZFrame_;
//   TLorentzVector Z=Z_;

//   double zmass = Z.M();
//   double mt = 1.777;
//   TLorentzVector a1 = p1+p2+p3;
//   double ma  = a1.M();
//   double diffmass = mt*mt - ma*ma;
//   float QQ = a1.E()*a1.E() - a1.P()*a1.P();

//   double x = 2*a1.E()/zmass;
//   double ctheta = 4*mt*mt*a1.E()/zmass/diffmass   - (mt*mt + ma*ma)/diffmass;
//   return ctheta;

// }


// float 
// a1Helper::CosBeta(){
//   TLorentzVector p1 = SSPion1ZFrame_;
//   TLorentzVector p2 = SSPion2ZFrame_;
//   TLorentzVector p3 = OSPionZFrame_;

//   float mpi  = 0.139;
// //   float E = p1.E() +  p2.E() +  p3.E(); 
// //   float P = p1.P() +  p2.P() +  p3.P(); 
// //   float QQ = E*E - P*P;


// //   std::cout<<"  cosbeta --------------------- "<<std::endl;

//   TLorentzVector a1 = p1+p2+p3;
//   float QQ = a1.E()*a1.E() - a1.P()*a1.P();

//   float B1 = (pow(p1.E()*a1.E()   - Scalar(p1,a1),2 ) - QQ*mpi*mpi)/QQ;
//   float B2 = (pow(p2.E()*a1.E()   - Scalar(p2,a1),2 ) - QQ*mpi*mpi)/QQ;
//   float B3 = (pow(p3.E()*a1.E()   - Scalar(p3,a1),2 ) - QQ*mpi*mpi)/QQ;

//   float T = 0.5*sqrt(-lambda(B1,B2,B3));

//   TLorentzVector p1Timesp2(p1.Py()*p2.Pz() - p1.Pz()*p2.Py(),p1.Pz()*p2.Px() - p1.Px()*p2.Pz(),p1.Px()*p2.Py() - p1.Py()*p2.Px(),1);
    
    
//   float cbeta = Scalar(p3,p1Timesp2)/a1.P()/T;

// //   std::cout<<"  T  "<< T<<std::endl;
// //   std::cout<<"  B1  "<< B1<<std::endl;
// //   std::cout<<"  B2  "<< B2<<std::endl;
// //   std::cout<<"  B3  "<< B3<<std::endl;
// //   std::cout<<"  QQ  "<< QQ<<std::endl;
// //   std::cout<<"  Scalar(p3,p1Timesp2)  "<< Scalar(p3,p1Timesp2)<<std::endl;
// //   std::cout<<"  cbeta  "<< cbeta  <<std::endl;


//   return cbeta;

// }

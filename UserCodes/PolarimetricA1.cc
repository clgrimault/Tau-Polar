#include "PolarimetricA1.h"
#include <iostream>

PolarimetricA1::PolarimetricA1(){
}

PolarimetricA1::PolarimetricA1(vector<TLorentzVector> TauA1andProd){
  if(TauA1andProd.size()!=4){
    std::cout<<" Warning!! Size of a1 input vector != 4 !! "<<std::endl;
  }
  TLorentzVector fakeboost(0,0,0,0);
  Setup(TauA1andProd,fakeboost);
}


PolarimetricA1::PolarimetricA1(vector<TLorentzVector> TauA1andProd, TLorentzVector RefernceFrame){
  if(TauA1andProd.size()!=4){
    std::cout<<" Warning!! Size of a1 input vector != 4 !! "<<std::endl;
  }
  Setup(TauA1andProd,RefernceFrame);
}


void 
PolarimetricA1::Setup(vector<TLorentzVector> TauA1andProd, TLorentzVector ReferenceFrame){
   mpi  = 0.13957018; // GeV 
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
   debug  = false;

   TVector3 RotVector = TauA1andProd.at(0).Vect();
   ReferenceFrame.SetVect(Rotate(ReferenceFrame.Vect(),RotVector));
   _tauAlongZLabFrame = TauA1andProd.at(0);
   _tauAlongZLabFrame.SetVect(Rotate(_tauAlongZLabFrame.Vect(),RotVector));
   
   for(int i=0; i<TauA1andProd.size(); i++){
     TLorentzVector Rotated = TauA1andProd.at(i);
     Rotated.SetVect(Rotate(Rotated.Vect(),RotVector));
     TauA1andProd_RF.push_back(Boost(Rotated,ReferenceFrame));
   }
   LFosPionLV  = TauA1andProd.at(1);
   LFss1pionLV = TauA1andProd.at(2);
   LFss2pionLV = TauA1andProd.at(3);
   LFa1LV = LFosPionLV+LFss1pionLV+LFss2pionLV;
   LFtauLV = TauA1andProd.at(0);
   LFQ= LFa1LV.M();

  

   _osPionLV   = TauA1andProd_RF.at(1);
   _ss1pionLV  = TauA1andProd_RF.at(2);
   _ss2pionLV  = TauA1andProd_RF.at(3);
   _a1LV       = _osPionLV+_ss1pionLV+_ss2pionLV;
   _tauLV      = TauA1andProd_RF.at(0);
   _nuLV      = _tauLV - _a1LV;
   _s12 = _ss1pionLV  + _ss2pionLV;
   _s13 = _ss1pionLV  + _osPionLV;
   _s23 = _ss2pionLV  + _osPionLV;
   _s1  =  _s23.M2(); 
   _s2  =  _s13.M2();
   _s3  =  _s12.M2();
   _Q = _a1LV.M();

   // std::cout<<"tau, a1, pi1,pi2,pi3,nu  ";
   // _tauLV.Print();
   // _a1LV.Print();
   // _ss1pionLV.Print();
   // _ss2pionLV.Print();
   // _osPionLV.Print();
   // _nuLV.Print();
}

void 
PolarimetricA1::subSetup(double s1, double s2, double s3, double Q){
   _s1  =   s1;
   _s2  =   s2;
   _s3  =   s3;
   _Q = Q;
}



void 
PolarimetricA1::Configure(vector<TLorentzVector> TauA1andProd){

  if(TauA1andProd.size()!=4){
    std::cout<<" Warning!! Size of input vector != 4 !! "<<std::endl;
  }
  TLorentzVector fakeboost(0,0,0,0);
  Setup(TauA1andProd,fakeboost);

}

void 
PolarimetricA1::Configure(vector<TLorentzVector> TauA1andProd, TLorentzVector RefernceFrame){
  if(TauA1andProd.size()!=4){
    std::cout<<" a1 helper:  Warning!! Size of input vector != 4!   Size = "<< TauA1andProd.size()<<std::endl;
  }
  Setup(TauA1andProd,RefernceFrame);

}
bool
PolarimetricA1::isConfigured(){
  if(TauA1andProd_RF.size()!=4){ std::cout<<"Error:   PolarimetricA1 is not Configured! Check  the size of input vector!  Size =  "<< TauA1andProd_RF.size() <<std::endl; return false;} return true;
}



void 
PolarimetricA1::SetParametersReco(TLorentzVector tau, TLorentzVector mu ){
 Initialize(tau,mu);
}
void 
PolarimetricA1::SetFrame(TLorentzVector vec){
  Boost_ = vec;
}



PolarimetricA1::~PolarimetricA1(){
}



void 
PolarimetricA1::Initialize(TLorentzVector t, TLorentzVector mu){
  RecoMuon_=mu;
  KFitTau_=t;
}





double 
PolarimetricA1::lambda(double x, double y, double z){
    return x*x +y*y +z*z - 2*x*y - 2*x*z - 2*z*y;
}
TLorentzVector 
PolarimetricA1::Boost(TLorentzVector pB, TLorentzVector frame){
   TMatrixT<double> transform(4,4);
   TMatrixT<double> result(4,1);
   TVectorT<double> vec(4); 
   TVector3 b;
   if(frame.Vect().Mag()==0){ std::cout<<"PolarimetricA1  Boost is not set, perfrom calculation in the Lab Frame   "<<std::endl; return pB;}
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
PolarimetricA1::Scalar(TLorentzVector p1, TLorentzVector p2){
    return p1.Vect()*p2.Vect();
}
double 
PolarimetricA1::MomentSFunction(double s, string type){
  int cells(20);
  //  double s = Q*Q;
  double intx(0);
  double m1 = mpi;
  double m2 = mpi;
  double m3 = mpi;

  double m13(0);
  double integral(0);

  double da1(0), db1(0);
  vector<double> set;
  set.push_back(_s1);
  set.push_back(_s2);
  set.push_back(_s3);
  set.push_back(_Q);
  double  stepx  = (pow(sqrt(s)-m2,2) - pow( m1+m3,2) ) / cells;
  for(unsigned int i=1;i<cells + 1;i++){ 
    da1 = pow(m1+m3,2) + stepx*(i-1);
    db1 = pow(m1+m3,2) + stepx*i;
    m13 = 0.5*(da1 + db1);
    double  E3s = (m13 - m1*m1 + m3*m3)/(2*sqrt(m13));  
    double  E2s = (s   - m13  -m2*m2)/(2*sqrt(m13));  
    double  m23max =pow (E2s+E3s,2) - pow( sqrt(E2s*E2s - m2*m2) - sqrt(E3s*E3s - m3*m3),2);
    double  m23min =  pow(E2s+E3s,2) - pow( sqrt(E2s*E2s - m2*m2) + sqrt(E3s*E3s - m3*m3),2);
    double  stepy = (m23max - m23min)/cells;
    double da2(0), db2(0);
    double inty(0);
    double m23(0);
    double m12(0);
    for(unsigned int j=1;j<cells + 1;j++){ 
      da2 = m23min + stepy*(j-1);
      db2 = m23min + stepy*j;
      m23 = 0.5*(da2 + db2);
      m12 = s +m1*m1 + m2*m2 + m3*m3 - m13 - m23;
      subSetup(m23,m13,m12,sqrt(s)); 
      // if(s >1.88 && s < 1.90)    std::cout<<"  WD=    "<< WD() << "        m23 = "<< m23 << "       m13= " << m13 <<    "     m12=   "<< m12 << "    sqrts=  " << sqrt(s) <<std::endl;
      // if(s >1.88 && s < 1.90)    std::cout<<"  F1  "<< F1() << "  F2 "<< F2() << " h0 " << h0()<<    "  VV1()-h()=   "<< VV1()-h() << " VV2() -h()=  " <<VV2()-h() <<std::endl;
      double SFunction(0);
      if(type=="WA")SFunction=WA();
      if(type=="WC")SFunction=WC();
      if(type=="WSA")SFunction=WSA();
      if(type=="WSB")SFunction=WSB();
      if(type=="WD"  ){
	if(m23 > m13)SFunction=WD();
	else SFunction=-WD();
      }
      if(type=="WE"){
	if(m23 > m13)SFunction=WE();
	else SFunction=-WE();
      }
      if(type=="WSD"){
	if(m23 > m13)SFunction=WSD();
	else SFunction=-WSD();
      }
      //      std::cout<<"SFunction  "<< SFunction<< std::endl;
      inty+=stepx*stepy*SFunction;
    }
    intx+=inty;
  }
  integral = intx;
  // std::cout<<"check return value  "<< s << "   integral=  " << integral << "  type  "<<   type  << "  :  "<<  WD()<< std::endl;
  subSetup(set.at(0),set.at(1),set.at(2),set.at(3));

  return integral;
}

double PolarimetricA1::K1(double ct, double QQ, int hel){
  if(debug){if(fabs(ct) > 1) std::cout<<"Warning! K1: |ct| > 1 "<<std::endl;}
  return   1 - hel*ct - mtau*mtau*(1+hel*ct)/QQ/QQ;
}
double PolarimetricA1::K2(double ct, double QQ, int hel){
  if(debug){if(fabs(ct) > 1) std::cout<<"Warning! K1: |ct| > 1 "<<std::endl;}
  return   mtau*mtau*(1+hel*ct)/QQ/QQ;
}
double PolarimetricA1::K3(double ct, double QQ, int hel){
  if(debug){if(fabs(ct) > 1) std::cout<<"Warning! K1: |ct| > 1 "<<std::endl;}
  return   1 - hel*ct;
}
double PolarimetricA1::K1bar(double ct, double QQ, int hel){
  if(debug){if(fabs(ct) > 1) std::cout<<"Warning! K1bar: |ct| > 1 "<<std::endl;}
  double cpsi = (ct*(mtau*mtau  + QQ)   + (mtau*mtau  - QQ))/(ct*(mtau*mtau  - QQ)   + (mtau*mtau  + QQ));
  if(debug){if(fabs(cpsi) > 1) std::cout<<"Warning! K1bar: |cpsi| > 1 "<<std::endl;}
  return  K1(ct,QQ,hel)*0.5*(3*cpsi*cpsi - 1) - 3*sqrt(mtau*mtau/QQ/QQ)*cpsi*sqrt(1-cpsi*cpsi)*sqrt(1-ct*ct)*hel;

}
double PolarimetricA1::K2bar(double ct, double QQ, int hel){
  if(debug){if(fabs(ct) > 1) std::cout<<"Warning! K1bar: |ct| > 1 "<<std::endl;}
  double cpsi = (ct*(mtau*mtau  + QQ)   + (mtau*mtau  - QQ))/(ct*(mtau*mtau  - QQ)   + (mtau*mtau  + QQ));
  if(debug){if(fabs(cpsi) > 1) std::cout<<"Warning! K1bar: |cpsi| > 1 "<<std::endl;}
  return  K2(ct,QQ,hel)*cpsi  + sqrt(mtau*mtau/QQ/QQ)*sqrt(1-cpsi*cpsi)*sqrt(1-ct*ct)*hel;

}
 double PolarimetricA1::K3bar(double ct, double QQ, int hel){
  if(debug){if(fabs(ct) > 1) std::cout<<"Warning! K1bar: |ct| > 1 "<<std::endl;}
  double cpsi = (ct*(mtau*mtau  + QQ)   + (mtau*mtau  - QQ))/(ct*(mtau*mtau  - QQ)   + (mtau*mtau  + QQ));
  if(debug){if(fabs(cpsi) > 1) std::cout<<"Warning! K1bar: |cpsi| > 1 "<<std::endl;}
  return  K3(ct,QQ,hel)*cpsi  - sqrt(mtau*mtau/QQ/QQ)*sqrt(1-cpsi*cpsi)*sqrt(1-ct*ct)*hel;

}
 


double 
PolarimetricA1::getMoment(double ct, string type, int hel){
  int cells(20);
  double qqmin  = 0.4;
  double qqmax = 3.0;
  vector<double> set;
  set.push_back(_s1);
  set.push_back(_s2);
  set.push_back(_s3);
  set.push_back(_Q);
  double  stepqq  = ( qqmax - qqmin) / cells;
  double kern(0);
  double atQQa(0);
  double atQQb(0);
  double atQQ(0);
  double integral(0);
  for(unsigned int i=1;i<cells + 1;i++){ 
    atQQa = qqmin + stepqq*(i-1);
    atQQb = qqmin + stepqq*i;
    atQQ = 0.5*(atQQa + atQQb);

    if(type=="one") kern = (2*K1(ct,atQQ,hel) + 3*K2(ct,atQQ,hel))*MomentSFunction(atQQ,"WA");
    if(type=="beta") kern = 0.2*K1bar(ct,atQQ,hel)*MomentSFunction(atQQ,"WA");
    if(type=="c2g") kern = -0.5*K1bar(ct,atQQ,hel)*MomentSFunction(atQQ,"WC");
    if(type=="s2g") kern = 0.5*K1bar(ct,atQQ,hel)*MomentSFunction(atQQ,"WD");
    if(type=="cb") kern = K3bar(ct,atQQ,hel)*MomentSFunction(atQQ,"WE");
    //    std::cout<<"  kern "<< kern << "   "<< stepqq <<std::endl;
    integral += kern*stepqq;

  }
  //  subSetup(set.at(0),set.at(1),set.at(2),set.at(3));
  return integral;
}
 



//---------------------------------------  hadronic current ---------------------------
double 
PolarimetricA1::WA(){
   return  VV1()*F1().Rho2() + VV2()*F2().Rho2()  + 2*V1V2()*( F1()*Conjugate(F2()) ).Re();

 }

 double 
PolarimetricA1::WC(){
   return  -(-VV1() + 2*h() )*F1().Rho2() - (-VV2() + 2*h())*F2().Rho2()   -   (-2*V1V2() - 4*h())*( F1()*Conjugate(F2()) ).Re();
 } 

 double
 PolarimetricA1::WD(){
   double QQ = _Q*_Q;
   double undersqrt1 = VV1()  -h();
   double undersqrt2 = VV2()  -h();

//   if(undersqrt1 < 0) undersqrt1 =0;
//   if(undersqrt2 < 0) undersqrt2 =0;

//   std::cout<<" WD(1)  "<< - 2 * sqrt(undersqrt1) * F1().Rho2()*sqrt(h())  << "  WD(2)   "<<sqrt(h())*2*sqrt(undersqrt2)*F2().Rho2()  << "   WD(3)   "<<  -sqrt(h())* (QQ - mpi*mpi + _s3)*(_s1 - _s2 )*( F1()*Conjugate(F2()) ).Re()/QQ/sqrt(h0() ) <<endl;
  

   return  -sqrt(h()) * ( 2 * sqrt(undersqrt1) * F1().Rho2() - 2*sqrt(undersqrt2)*F2().Rho2()  
			  + (QQ - mpi*mpi + _s3)*(_s1 - _s2 )*( F1()*Conjugate(F2()) ).Re()/QQ/sqrt(h0() ) );


 }

 double
 PolarimetricA1::WE(){
  return  3*sqrt(h()*h0())*( F1()*Conjugate(F2()) ).Im();
 }

double
 PolarimetricA1::WSA(){
  double QQ = _Q*_Q;
  return  QQ*F4().Rho2();
 }
double
 PolarimetricA1::WSB(){
  //  double QQ = _Q*_Q;
   double undersqrt1 = VV1()  -h();
   double undersqrt2 = VV2()  -h();
   return  -2*_Q* (sqrt(undersqrt1) * (F1()*Conjugate(F4())).Re() +   sqrt(undersqrt2)*(F2()*Conjugate(F4())).Re()  );
 }
double
 PolarimetricA1::WSD(){
  double QQ = _Q*_Q;
  return  2*sqrt(QQ*h())* ( (F1()*Conjugate(F4())).Re() - (F2()*Conjugate(F4())).Re()   );
 }
double
 PolarimetricA1::WSC(){
  //  double QQ = _Q*_Q;
   double undersqrt1 = VV1()  -h();
   double undersqrt2 = VV2()  -h();
   return  2*_Q* (sqrt(undersqrt1) * (F1()*Conjugate(F4())).Im() +   sqrt(undersqrt2)*(F2()*Conjugate(F4())).Im()  );
 }
double
 PolarimetricA1::WSE(){
  double QQ = _Q*_Q;
   return  -2*sqrt(QQ*h())* ( (F1()*Conjugate(F4())).Im() - (F2()*Conjugate(F4())).Im()   );
 }



double
PolarimetricA1::cosgammaLF(){
  double QQ=LFQ*LFQ;
  // double B1 = (pow(_ss1pionLV.E()*_tauLV.E()   - _ss1pionLV.Vect().Dot(_a1LV.Vect()),2 ) - QQ*mpi*mpi)/QQ;
  // double B2 = (pow(_ss2pionLV.E()*_tauLV.E()   - _ss2pionLV.Vect().Dot(_a1LV.Vect()),2 ) - QQ*mpi*mpi)/QQ;
  double B3 = (pow(LFosPionLV.E()*LFtauLV.E()   - LFosPionLV.Vect().Dot(LFa1LV.Vect()),2 ) - QQ*mpi*mpi)/QQ;

  // double T = 0.5*sqrt(-lambda(B1,B2,B3));
  // double A1=(_a1LV.E()*_a1LV.Vect().Dot(_ss1pionLV.Vect()) - _ss1pionLV.E()*_a1LV.P()*_a1LV.P())/QQ;
  // double A2=(_a1LV.E()*_a1LV.Vect().Dot(_ss2pionLV.Vect()) - _ss2pionLV.E()*_a1LV.P()*_a1LV.P())/QQ;
  double A3=(LFa1LV.E() * LFa1LV.Vect().Dot(LFosPionLV.Vect()) - LFosPionLV.E()*LFa1LV.P()*LFa1LV.P()) / LFQ;
  // std::cout<<"sqrt B3 "<< sqrt(B3)<<std::endl;
  // std::cout<<"A3 "<< A3<<std::endl;

  // std::cout<< "fuck 1  "   <<LFosPionLV.E()*LFtauLV.E()   - LFosPionLV.Vect().Dot(LFa1LV.Vect())<<std::endl;
  // std::cout<< "fuck 2  "   <<LFa1LV.E() * LFa1LV.Vect().Dot(LFosPionLV.Vect()) - LFosPionLV.E()*LFa1LV.P()*LFa1LV.P()<<std::endl;
                                         
  // std::cout<< "QQ "   << LFQ*LFQ<< "  _Q_Q  "<< _Q*_Q << std::endl;

  if(B3<=0 || cosbetaLF() >=1){std::cout<<"Warning! In PolarimetricA1::cosgamma square root <=0! return 0"<<std::endl; return 0;}
  return A3/LFa1LV.P()/sqrt(B3)/sqrt(1 - cosbetaLF()*cosbetaLF());
}

double
PolarimetricA1::singammaLF(){
  double QQ=LFQ*LFQ;
   double B1 = (pow(LFss1pionLV.E()*LFa1LV.E()   - LFss1pionLV.Vect().Dot(LFa1LV.Vect()),2 ) - QQ*mpi*mpi)/QQ;
  double B2 = (pow(LFss2pionLV.E()*LFa1LV.E()   - LFss2pionLV.Vect().Dot(LFa1LV.Vect()),2 ) - QQ*mpi*mpi)/QQ;
   double B3 = (pow(LFosPionLV.E()*LFa1LV.E()   - LFosPionLV.Vect().Dot(LFa1LV.Vect()),2 ) - QQ*mpi*mpi)/QQ;

  double T = 0.5*sqrt(-lambda(B1,B2,B3));

  double A1=(LFa1LV.E()*LFa1LV.Vect().Dot(LFss1pionLV.Vect()) - LFss1pionLV.E()*LFa1LV.P()*LFa1LV.P())/QQ;
  //  double A2=(_a1LV.E()*_a1LV.Vect().Dot(_ss2pionLV.Vect()) - _ss2pionLV.E()*_a1LV.P()*_a1LV.P())/QQ;
  double A3=(LFa1LV.E()*LFa1LV.Vect().Dot(LFosPionLV.Vect()) - LFosPionLV.E()*LFa1LV.P()*LFa1LV.P())/QQ;

  if(A3==0 || T==0){std::cout<<"Warning! In PolarimetricA1::singamma denominator ==0! return 0"<<std::endl; return 0;}
  double scale = -(B3*A1/A3 - 0.5*(B2 - B1 - B3))/T;
  //  std::cout<<"scale  " << scale <<std::endl;
  return cosgammaLF()*scale;
}
double
PolarimetricA1::cos2gamma(){
   return singamma()*singamma()   -     cosgamma()*cosgamma();
}

double
PolarimetricA1::sin2gamma(){
  return 2*singamma()*cosgamma();
}
double 
PolarimetricA1::cospsiLF(){
  double QQ = LFQ*LFQ;
  double s = 4*LFtauLV.E()*LFtauLV.E();
  double x = 2*LFa1LV.E()/sqrt(s);
  if(x*x  - 4*QQ/s <= 0 ){if(debug){std::cout<<"Warning! In PolarimetricA1::cospsi root square <=0! return 0"<<std::endl;} return 0;}
  return    ( x*(mtau*mtau + QQ)  - 2*QQ  )   /   ( mtau*mtau  - QQ   ) / sqrt(x*x  - 4*QQ/s); 
}
double 
PolarimetricA1::sinpsiLF(){
  if(cospsiLF()*cospsiLF() > 1  ){if(debug){std::cout<<"Warning! In PolarimetricA1::sinpsi root square <=0! return nan"<<std::endl;}}
  return    sqrt(1 - cospsiLF()*cospsiLF());
}

double 
PolarimetricA1::ultrarel_cospsiLF(){
  double QQ = LFQ*LFQ;
  double cos = (costhetaLF()*(mtau*mtau  + QQ)   + (mtau*mtau  - QQ))/(costhetaLF()*(mtau*mtau  - QQ)   + (mtau*mtau  + QQ));
  return cos;
}

double 
PolarimetricA1::costhetaLF(){
  double QQ = LFQ*LFQ;
  double x = LFa1LV.E()/LFtauLV.E();
  double s = 4*LFtauLV.E()*LFtauLV.E();
  if( 1 - 4*mtau*mtau/s  <= 0 ){if(debug){std::cout<<"Warning! In PolarimetricA1::costheta root square <=0! return 0"<<std::endl;} return 0;}
  return (2*x*mtau*mtau - mtau*mtau - QQ)/( (mtau*mtau - QQ)*sqrt(1 - 4*mtau*mtau/s) );
}
double 
PolarimetricA1::sinthetaLF(){
  if( costhetaLF()*costhetaLF() > 1 ) {if(debug){std::cout<<"Warning! In PolarimetricA1::sin theta root square <=0! return nan;   costheta = "<< costhetaLF()<<std::endl; }}
  return sqrt(1- costhetaLF()*costhetaLF());
}

double 
PolarimetricA1::cosbetaLF(){
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
  if(T==0 || LFa1LV.P()==0){if(debug){std::cout<<" Warning!  Can not compute cosbetaLF, denominator =0; return 0; "<<std::endl;} return 0;}
  return ospionVect.Dot(ss1pionVect.Cross(ss2pionVect)) /LFa1LV.P()/T;
}


double
PolarimetricA1::V1(){ 
  double QQ = _Q*_Q;
  return  4*mpi*mpi -_s2  - pow(_s3  - _s1,2)/4/QQ;
}

double
PolarimetricA1::V2(){ 
  double QQ = _Q*_Q;
  return  4*mpi*mpi -_s1  - pow(_s3  - _s2,2)/4/QQ;
}


double
PolarimetricA1::VV12(){ 
  double QQ = _Q*_Q;
  return  -(QQ/2 - _s3 - 0.5*mpi*mpi) - (_s3 - _s1)*(_s3 - _s2)/4/QQ;
}
double
PolarimetricA1::JJ(){
  double QQ = _Q*_Q;
  return  (V1()*BreitWigner(sqrt(_s2),"rho").Rho2() + V2()*BreitWigner(sqrt(_s1),"rho").Rho2()  + VV12()*( BreitWigner(sqrt(_s1),"rho")*Conjugate(BreitWigner(sqrt(_s2),"rho")) + BreitWigner(sqrt(_s2),"rho")*Conjugate(BreitWigner(sqrt(_s1),"rho"))  ))*f3(sqrt(QQ)).Rho2();
  // std::cout<<" FORM1  "<<f3(sqrt(QQ))* BreitWigner(sqrt(_s2),"rho") <<std::endl;
  // std::cout<<" FORM2  "<<f3(sqrt(QQ))* BreitWigner(sqrt(_s1),"rho") <<std::endl;

}


//  double
// PolarimetricA1::JN(){ //  this is -V1^{2}
//   double QQ = _Q*_Q;
//   return  (V1()*BreitWigner(sqrt(_s2),"rho") + V2()*BreitWigner(sqrt(_s1),"rho")  + VV12()*( BreitWigner(sqrt(_s1),"rho")*Conjugate(BreitWigner(sqrt(_s2),"rho")) + BreitWigner(sqrt(_s2),"rho")*Conjugate(BreitWigner(sqrt(_s1),"rho"))  ))*f3().Rho2();
// }


TLorentzVector
PolarimetricA1::PTenzor5(TLorentzVector aR, TLorentzVector aI, TLorentzVector bR, TLorentzVector bI, TLorentzVector c){ 
  TComplex a4(aR.E(), aI.E());  TComplex a1(aR.Px(),aI.Px());   TComplex a2(aR.Py(),aI.Py());   TComplex a3(aR.Pz(),aI.Pz());
  TComplex b4(bR.E(), bI.E());  TComplex b1(bR.Px(),bI.Px());   TComplex b2(bR.Py(),bI.Py());   TComplex b3(bR.Pz(),bI.Pz());
  // double  c1 = c.Px();   double  c2 = c.Py();   double  c3 = c.Pz();   double  c4 = c.E();
  
  //  TComplex c0(cR.E(), cI.E());  TComplex c1(cR.Px(),cI.Px());   TComplex c2(cR.Py(),cI.Py());   TComplex c3(cR.Pz(),cI.Pz());
  double  c1 = c.Px();   double  c2 = c.Py();   double  c3 = c.Pz();   double  c4 = c.E();

  double d34 = (a3*b4 - a4*b3).Im();
  double d24 = (a2*b4 - a4*b2).Im();  
  double d23 = (a2*b3 - a3*b2).Im();
  double d14 = (a1*b4 - a4*b1).Im();
  double d13 = (a1*b3 - a3*b1).Im();
  double d12 = (a1*b2 - a2*b1).Im();

  double PIAX1 = 2*( c2*d34 - c3*d24 + c4*d23);
  double PIAX2 = 2*(-c1*d34 + c3*d14 - c4*d13);
  double PIAX3 = 2*( c1*d24 - c2*d14 + c4*d12);
  double PIAX4 = 2*(-c1*d23 + c2*d13 - c3*d12);

  //std::cout<<" a0, a1, a2, a3  "<< a0<< a1 <<a2<<a3<<std::endl;
  //0:  (  a1*(b2*c3  - b3*c2) - a2*( b1 * c3   - b3 *c1 )  + a3*( b1* c2   - b2* c1)  )
  //1: -(  a0*(b2*c3  - b3*c2) - a2*( b0 * c3   - b3 *c0 )  + a3*( b0* c2   - b2* c0)  )
  //2:  (  a0*(b1*c3  - b3*c1) - a1*( b0 * c3   - b3 *c0 )  + a3*( b0* c1   - b1* c0)  )
  //3: -(  a0*(b1*c2  - b2*c1) - a1*( b0 * c2   - b2 *c0 )  + a2*( b0* c1   - b1* c0)  )

  // TComplex d0 = (  a1*(b2*c3  - b3*c2) - a2*( b1 * c3   - b3 *c1 )  + a3*( b1* c2   - b2* c1)  );
  // TComplex d1 =-(  a0*(b2*c3  - b3*c2) - a2*( b0 * c3   - b3 *c0 )  + a3*( b0* c2   - b2* c0)  );
  // TComplex d2 = (  a0*(b1*c3  - b3*c1) - a1*( b0 * c3   - b3 *c0 )  + a3*( b0* c1   - b1* c0)  );
  // TComplex d3 =-(  a0*(b1*c2  - b2*c1) - a1*( b0 * c2   - b2 *c0 )  + a2*( b0* c1   - b1* c0)  );


  //0:   (  a.Px()*(b.Py()*c.Pz() - b.Pz()*c.Py()) - a.Py()*( b.Px() * c.Pz()  - b.Pz() *c.Px() ) + a.Pz()*( b.Px()* c.Py()   - b.Py()* c.Px())  )
  //1: - (  a.E()*(b.Py()*c.Pz()  - b.Pz()*c.Py()) - a.Py()*( b.E() * c.Pz()   - b.Pz() *c.E() )  + a.Pz()*( b.E() * c.Py()   - b.Py()* c.E())  )
  //2:   (  a.E()*(b.Px()*c.Pz()  - b.Pz()*c.Px()) - a.Px()*( b.E() * c.Pz()   - b.Pz() *c.E() )  + a.Pz()*( b.E() * c.Px()   - b.Px()* c.E())  )
  //3: - (  a.E()*(b.Px()*c.Py()  - b.Py()*c.Px()) - a.Px()*( b.E() * c.Py()   - b.Py() *c.E() )  + a.Py()*( b.E() * c.Px()   - b.Px()* c.E())  )
  // TLorentzVector d(- (  a.E()*(b.Py()*c.Pz()  - b.Pz()*c.Py()) - a.Py()*( b.E() * c.Pz()   - b.Pz() *c.E() )  + a.Pz()*( b.E() * c.Py()   - b.Py()* c.E())  ),
  // 	   	      (  a.E()*(b.Px()*c.Pz()  - b.Pz()*c.Px()) - a.Px()*( b.E() * c.Pz()   - b.Pz() *c.E() )  + a.Pz()*( b.E() * c.Px()   - b.Px()* c.E())  ),
  // 		    - (  a.E()*(b.Px()*c.Py()  - b.Py()*c.Px()) - a.Px()*( b.E() * c.Py()   - b.Py() *c.E() )  + a.Py()*( b.E() * c.Px()   - b.Px()* c.E())  ),
  // 		      (  a.Px()*(b.Py()*c.Pz() - b.Pz()*c.Py()) - a.Py()*( b.Px() * c.Pz()  - b.Pz() *c.Px() ) + a.Pz()*( b.Px()* c.Py()   - b.Py()* c.Px())  ));
  // 
  TLorentzVector d(PIAX1,PIAX2,PIAX3,PIAX4);
    //TLorentzVector d(2*d1.Im(),2*d2.Im(),2*d3.Im(),2*d0.Im());
  // d.Print();
  return d;
}


TComplex
PolarimetricA1::f3(double Q){ 
  return  (coscab*2*sqrt(2)/3/fpi)*BreitWigner(Q,"a1");
}


TLorentzVector
PolarimetricA1::PolarimetricVector(){ 


   TLorentzVector q1= _ss1pionLV;
   TLorentzVector q2= _ss2pionLV;
   TLorentzVector q3= _osPionLV;
   TLorentzVector a1 = q1+q2+q3;
   TLorentzVector N = _nuLV;
   TLorentzVector P = _tauLV;
   double s1 = (q2+q3).M2();
   double s2 = (q1+q3).M2();


   // std::cout<<"-------------"<<std::endl;
   // P.Print();
   // N.Print();
   // a1.Print();
   // std::cout<<"============="<<std::endl;
   // q1.Print();
   // q2.Print();
   // q3.Print();

   TLorentzVector vec1 = q1 - q3 -  a1* (a1*(q1-q3)/a1.M2());
   TLorentzVector vec2 = q2 - q3 -  a1* (a1*(q2-q3)/a1.M2());
//    std::cout<<" s1,s2,QQ  "<< s1 <<"  "<<s2<<"  "<<a1.M2()<<std::endl;

// 0.742717  0.433661  1.61186


   std::cout<<" F3PI1   "<< F3PI(1,1.61186,0.742717,0.433661) <<std::endl;
   std::cout<<" F3PI2   "<< F3PI(2,1.61186,0.742717,0.433661) <<std::endl;
   std::cout<<" F3PI3   "<< F3PI(3,1.61186,0.742717,0.433661) <<std::endl;

   TComplex BWProd1 = f3(a1.M())*BreitWigner(sqrt(s2),"rho");
   TComplex BWProd2 = f3(a1.M())*BreitWigner(sqrt(s1),"rho");

   double BWProd1Re = BWProd1.Re();   double BWProd1Im = BWProd1.Im();
   double BWProd2Re = BWProd2.Re();   double BWProd2Im = BWProd2.Im();

 
   TLorentzVector PT5 = PTenzor5(JConjRe( q1,  q2,  q3,  a1), JConjIm( q1,  q2,  q3,  a1), JRe( q1,  q2,  q3,  a1), JIm( q1,  q2,  q3,  a1),N);
 
   double omega = P*PTenzor(q1,q2,q3,a1,N) - P*PT5;

   //    std::cout<<"P  and P5  omega  "<< omega<<std::endl;
   // PTenzor(q1,q2,q3,a1,N).Print();
   // PT5.Print();


   TLorentzVector out =  (P.M()*P.M()*  (PT5 - PTenzor(q1,q2,q3,a1,N))  -  P*(  P*PT5  -  P*PTenzor(q1,q2,q3,a1,N)))*(1/omega/P.M());

   // std::cout<<"1st:  "<< P.M()*P.M()*  (PT5 - PTenzor(q1,q2,q3,a1,N))/omega/P.M() <<std::endl;
   // std::cout<<"2nd:  "<< P*(  P*PT5  -  P*PTenzor(q1,q2,q3,a1,N))* (1/omega/P.M()) <<std::endl;
   // ( (P.M()/omega)*(PT5 - PTenzor(q1,q2,q3,a1,N))).Print();
   // (P*(  P*PT5  -  P*PTenzor(q1,q2,q3,a1,N))).Print();
   // out.Vect().Print();
   return out;
}

TLorentzVector
PolarimetricA1::PTenzor(TLorentzVector q1, TLorentzVector q2, TLorentzVector q3, TLorentzVector a1, TLorentzVector N){ 
  double s1 = (q2+q3).M2();
  double s2 = (q1+q3).M2();
  double s3 = (q2+q3).M2();
  
  TLorentzVector vec1 = q1 - q3 -  a1* (a1*(q1-q3)/a1.M2());
  TLorentzVector vec2 = q2 - q3 -  a1* (a1*(q2-q3)/a1.M2());



  TComplex L1 = f3(a1.M())*BreitWigner(sqrt(s2),"rho");
  TComplex L2 = f3(a1.M())*BreitWigner(sqrt(s1),"rho");

  TComplex CL1 = Conjugate(f3(a1.M())*BreitWigner(sqrt(s2),"rho"));
  TComplex CL2 = Conjugate(f3(a1.M())*BreitWigner(sqrt(s1),"rho"));



  TComplex factor1=ConjJN(q1,q2,q3,a1,N)*L1 + JN(q1,q2,q3,a1,N)*CL1;
  TComplex factor2=ConjJN(q1,q2,q3,a1,N)*L2 + JN(q1,q2,q3,a1,N)*CL2;

  TLorentzVector  Ptenz= 2*BreitWigner(sqrt(s1),"rho").Rho2()*(vec2*N)*vec2 + 2*BreitWigner(sqrt(s2),"rho").Rho2()*(vec1*N)*vec1 + 2*(BreitWigner(sqrt(s2),"rho")*Conjugate(BreitWigner(sqrt(s1),"rho"))  ).Re() *((vec1*N)*vec2 + (vec2*N)*vec1) - JJ()*N;

  TLorentzVector out  = 2*(factor1*vec1 + factor2*vec2 - JJ()*N);
  return out;
}


TComplex
PolarimetricA1::JN(TLorentzVector q1, TLorentzVector q2, TLorentzVector q3, TLorentzVector a1, TLorentzVector N){ 
  double s1 = (q2+q3).M2();
  double s2 = (q1+q3).M2();

   
  TLorentzVector vec1 = q1 - q3 -  a1* (a1*(q1-q3)/a1.M2());
  TLorentzVector vec2 = q2 - q3 -  a1* (a1*(q2-q3)/a1.M2());

  double prod1 = vec1*N;
  double prod2 = vec2*N;

  TComplex BWProd1 = f3(a1.M())*BreitWigner(sqrt(s2),"rho");
  TComplex BWProd2 = f3(a1.M())*BreitWigner(sqrt(s1),"rho");

  TComplex out  = BWProd1*prod1 + BWProd2*prod2;
  return out;
}
TComplex
PolarimetricA1::ConjJN(TLorentzVector q1, TLorentzVector q2, TLorentzVector q3, TLorentzVector a1, TLorentzVector N){ 
  double s1 = (q2+q3).M2();
  double s2 = (q1+q3).M2();
  
  TLorentzVector vec1 = q1 - q3 -  a1* (a1*(q1-q3)/a1.M2());
  TLorentzVector vec2 = q2 - q3 -  a1* (a1*(q2-q3)/a1.M2());

  double prod1 = vec1*N;
  double prod2 = vec2*N;
  TComplex BWProd1 = Conjugate(f3(a1.M())*BreitWigner(sqrt(s2),"rho"));
  TComplex BWProd2 = Conjugate(f3(a1.M())*BreitWigner(sqrt(s1),"rho"));
  TComplex out  = BWProd1*prod1 + BWProd2*prod2;
  return out;
}

TLorentzVector PolarimetricA1::JRe(TLorentzVector q1, TLorentzVector q2, TLorentzVector q3, TLorentzVector a1){

  double s1 = (q2+q3).M2();
  double s2 = (q1+q3).M2();
  double s3 = (q2+q3).M2();
  
  TLorentzVector vec1 = q1 - q3 -  a1* (a1*(q1-q3)/a1.M2());
  TLorentzVector vec2 = q2 - q3 -  a1* (a1*(q2-q3)/a1.M2());
    
  TComplex BWProd1 = f3(a1.M())*BreitWigner(sqrt(s2),"rho");
  TComplex BWProd2 = f3(a1.M())*BreitWigner(sqrt(s1),"rho");

  TLorentzVector out = vec1*BWProd1.Re() + vec2*BWProd2.Re();
  return out;
}
TLorentzVector PolarimetricA1::JIm(TLorentzVector q1, TLorentzVector q2, TLorentzVector q3, TLorentzVector a1){

  double s1 = (q2+q3).M2();
  double s2 = (q1+q3).M2();
  double s3 = (q2+q3).M2();
  
  TLorentzVector vec1 = q1 - q3 -  a1* (a1*(q1-q3)/a1.M2());
  TLorentzVector vec2 = q2 - q3 -  a1* (a1*(q2-q3)/a1.M2());
    
  TComplex BWProd1 = f3(a1.M())*BreitWigner(sqrt(s2),"rho");
  TComplex BWProd2 = f3(a1.M())*BreitWigner(sqrt(s1),"rho");

  TLorentzVector out = vec1*BWProd1.Im() + vec2*BWProd2.Im();
  return out;
}

TLorentzVector PolarimetricA1::JConjRe(TLorentzVector q1, TLorentzVector q2, TLorentzVector q3, TLorentzVector a1){

  double s1 = (q2+q3).M2();
  double s2 = (q1+q3).M2();
  double s3 = (q2+q3).M2();
  
  TLorentzVector vec1 = q1 - q3 -  a1* (a1*(q1-q3)/a1.M2());
  TLorentzVector vec2 = q2 - q3 -  a1* (a1*(q2-q3)/a1.M2());
    
  TComplex BWProd1 = Conjugate(f3(a1.M())*BreitWigner(sqrt(s2),"rho"));
  TComplex BWProd2 = Conjugate(f3(a1.M())*BreitWigner(sqrt(s1),"rho"));
  //  std::cout<< BWProd1 << BWProd2<<std::endl;
  TLorentzVector out = vec1*BWProd1.Re() + vec2*BWProd2.Re();
  return out;
}
TLorentzVector PolarimetricA1::JConjIm(TLorentzVector q1, TLorentzVector q2, TLorentzVector q3, TLorentzVector a1){

  double s1 = (q2+q3).M2();
  double s2 = (q1+q3).M2();
  double s3 = (q2+q3).M2();
  
  TLorentzVector vec1 = q1 - q3 -  a1* (a1*(q1-q3)/a1.M2());
  TLorentzVector vec2 = q2 - q3 -  a1* (a1*(q2-q3)/a1.M2());
    
  TComplex BWProd1 = Conjugate(f3(a1.M())*BreitWigner(sqrt(s2),"rho"));
  TComplex BWProd2 = Conjugate(f3(a1.M())*BreitWigner(sqrt(s1),"rho"));

  TLorentzVector out = vec1*BWProd1.Im() + vec2*BWProd2.Im();
  return out;
}


double
PolarimetricA1::VV1(){ //  this is -V1^{2}
  double QQ = _Q*_Q;
  return  _s2 - 4*mpi*mpi + pow(_s3 - _s1,2)/4/QQ;
}



double
PolarimetricA1::VV2(){ //  this is -V2^{2}
  double QQ = _Q*_Q;
  return  _s1 - 4*mpi*mpi + pow(_s3 - _s2,2)/4/QQ;
}

double
PolarimetricA1::V1V2(){  // this is -V1V2
  double QQ = _Q*_Q;
  return  (QQ/2 - _s3 - 0.5*mpi*mpi) + (_s3 - _s1)*(_s3 - _s2)/4/QQ;
}


double
PolarimetricA1::h0(){ // this is -3sqrt{h0}/2
  double QQ = _Q*_Q;
  return -4*mpi*mpi + pow(2*mpi*mpi - _s1 - _s2,2)/QQ;
}

double
PolarimetricA1::h(){
  double QQ = _Q*_Q;
  return (_s1*_s2*_s3 - mpi*mpi*pow(QQ - mpi*mpi,2))/h0()/QQ;  // this is sqrt{h}
}



TComplex 
PolarimetricA1::F1(){
  TComplex scale(0, -2*sqrt(2)/3/fpi);
  TComplex res = scale*BreitWigner(_Q,"a1")*BRho(sqrt(_s2));
  //  std::cout<<"  BreitWigner(_Q,a1)  " << BreitWigner(_Q,"a1") << " BRho(_s2)  " << BRho(sqrt(_s2))<< std::endl;
  return res;
}


TComplex 
PolarimetricA1::F2(){
  TComplex scale(0, -2*sqrt(2)/3/fpi);
  TComplex res = scale*BreitWigner(_Q,"a1")*BRho(sqrt(_s1));
  return res;
}

TComplex 
PolarimetricA1::F4(){
  TComplex scale(0, -gpiprimerhopi*grhopipi*fpiprime/2/pow(mrho,4)/pow(mpiprime,3));
  TComplex res = scale*BreitWigner(_Q,"piprime")*(_s1*(_s2-_s3)*BRho(sqrt(_s1)) + _s2*(_s1-_s3)*BRho(sqrt(_s2)));
  return res;
}


TComplex 
PolarimetricA1::BRho(double Q){
  //  std::cout<<"BRho:      BreitWigner(Q) " << BreitWigner(Q) << " BreitWigner(Q,rhoprime) " << BreitWigner(Q,"rhoprime")<< std::endl;
  return (BreitWigner(Q) + beta*BreitWigner(Q,"rhoprime"))/(1+beta);
}

TComplex 
PolarimetricA1::BreitWigner(double Q, string type){
  double QQ=Q*Q;
  double re,im;
  double m = Mass(type);
  double g  = Widths(Q,type);
  // re = (m*m*pow(m*m - QQ,2))/(pow(m*m - QQ,2) + m*m*g*g); // 
  // im = m*m*m*g/(pow(m*m - QQ,2) + m*m*g*g);
  re = (m*m*(m*m - QQ))/(pow(m*m - QQ,2) + m*m*g*g);
  im = Q*g/(pow(m*m - QQ,2) + m*m*g*g);
 

  TComplex out(re,im);
  return out;
}

double
PolarimetricA1::Widths(double Q, string type){
  double QQ = Q*Q;
  double Gamma;
  Gamma = Gamma0rho*mrho*pow( ppi(QQ)  / ppi(mrho*mrho), 3) /sqrt(QQ);
  if(type == "rhoprime"){
    Gamma=Gamma0rhoprime*QQ/mrhoprime/mrhoprime;
 }
  if(type == "a1"){
    Gamma=Gamma0a1*ga1(Q)/ga1(ma1);
    //    Gamma=Gamma0a1*ma1*ga1(Q)/ga1(ma1)/Q;
 }
  if(type == "piprime"){
    Gamma = Gamma0piprime*pow( sqrt(QQ)/mpiprime  ,5)*pow( (1-mrho*mrho/QQ)/(1-mrho*mrho/mpiprime/mpiprime) ,3);
  }
  //  std::cout<< " Widths :   type   " << type << " Gamma  " << Gamma << "  QQ  "<< QQ <<std::endl;
  return Gamma;
}
double PolarimetricA1::ga1(double  Q){
  double QQ = Q*Q;
  return (QQ > pow(mrho + mpi,2)) ?  QQ*(1.623 + 10.38/QQ - 9.32/QQ/QQ   + 0.65/QQ/QQ/QQ)  : 4.1*pow(QQ - 9*mpi*mpi,3)*(  1 - 3.3*(QQ - 9*mpi*mpi)  + 5.8*pow(QQ - 9*mpi*mpi,2)  );
}
double
PolarimetricA1::Mass(string type){
  double m = mrho;
  if(type == "rhoprime") return mrhoprime; 
  if(type == "a1") return ma1;
  if(type == "piprime") return mpiprime;
  //std::cout<< "  type   " << type << " Mass  " << std::endl;
  return m;
}


double PolarimetricA1::ppi(double QQ){  if(QQ < 4*mpi*mpi) std::cout<<"Warning! Can not compute ppi(Q); root square <0 ; return nan  "; return 0.5*sqrt(QQ - 4*mpi*mpi);}


 double PolarimetricA1::getf(){
   double QQ=_Q*_Q;
   double l  = 0.5*(mtau*mtau + QQ)/sqrt(QQ);
   double line1 =   -2*l   *   ( 2*WA()/3   + 0.5*(3*cospsiLF()*cospsiLF()   -1)  *  ( WA()*(3*cosbeta()*cosbeta() -1 )/6    - 0.5*WC()*sinbeta()*sinbeta()*cos2gamma()   + 0.5*WD()* sinbeta()*sinbeta()* sin2gamma() )   )/sqrt(QQ);
   double line2 = mtau*mtau*WA()/QQ + mtau*mtau  *  (  WSA() +  cospsiLF()*sinbeta()*(   WSB() *cosgamma()      - WSD() * singamma())     )/QQ + WE()*cosbeta()*cospsiLF();
   double res = line1+ line2;

   return res;
 }
 double PolarimetricA1::getg(){
   double QQ=_Q*_Q;
   //  double l  = 0.5*(mtau*mtau + QQ)/sqrt(QQ);
   double l0= 0.5*(mtau*mtau - QQ)/sqrt(QQ);
   double line1 =   -2*l0   * costhetaLF()*  ( 2*WA()/3   + 0.5*(3*cospsiLF()*cospsiLF()   -1)  *  ( WA()*(3*cosbeta()*cosbeta() -1 )/6    -  0.5*WC()*sinbeta()*sinbeta()*cos2gamma()   + 0.5*WD()* sinbeta()*sinbeta()* sin2gamma() )   )/sqrt(QQ);
   double line2 = mtau*mtau*WA()*costhetaLF()/QQ  +    sqrt(mtau*mtau/QQ )  * sinthetaLF()* ( 0.5*WA()*2* sinbeta()*cosbeta()* cosalpha() -  
                                              WC()*sinbeta()* (sinalpha()* sin2gamma() + cos2gamma()* cosalpha()*    cosbeta() )    -    WD()*sinbeta()*( sinalpha()*cos2gamma() + sin2gamma()* cosalpha()*cosbeta()  )- 2*cospsiLF()*sinpsiLF()  *
 				            (WA()*(3*cosbeta() *cosbeta() -1 )/6   -    0.5*WC()*sinbeta()* sinbeta()* cos2gamma()+ 0.5*WD()*sinbeta()* sinbeta()* cos2gamma() + WD()*sinbeta()* sinbeta()* sin2gamma())/3   );

   double line3  =  sqrt(mtau*mtau/QQ ) *sinthetaLF()* (WE()*(cosbeta()*sinpsiLF() + sinbeta()*cosalpha()) +cosbeta()*sinalpha()*(WSC()*cosgamma() - WSE()*singamma()) + cosalpha()*(WSC()*singamma() + WSE()*cosgamma()));
   double line4  =  -WE()*costhetaLF()*cosbeta()*cospsiLF() + mtau*mtau*costhetaLF()*(WSA() + cospsiLF()*sinbeta()  * (WSB()*cosgamma() - WSD()* singamma()  ) )/QQ;
   double line5  =  sqrt(mtau*mtau/QQ)*sinthetaLF() *  ( sinpsiLF()*sinbeta()*( WSB()* cosgamma() - WSD()* singamma()) + cosbeta()*cosalpha()*(WSD()*singamma() - WSB()*cosgamma()  ) + 
 							sinalpha()*(WSD()*cosgamma() + WSB()*singamma())          );
   double res = line1+ line2 + line3 + line4 + line5;
   return res;
 }


 double PolarimetricA1::vgetf(TString type){
   double QQ=_Q*_Q;
   double RR  = mtau*mtau/QQ; double R = sqrt(RR);
   float U = 0.5*(3*cospsiLF()*cospsiLF() - 1)*(1 - RR);
   double B = 0.5*(3*cosbeta()*cosbeta() - 1);

   double fA =  WA()*(2+RR + B*U)/3;
   double fC = -WC()*0.5*U*sinbeta()*sinbeta()* cos2gamma();
   double fD = WD()*0.5*U*sinbeta()*sinbeta()* sin2gamma();
   double fE =  WE()*cospsiLF()*cosbeta();

   double res = fA+fC+fD+fE;

   return res;
 }


 double PolarimetricA1::vgetg(TString type){
   double QQ=_Q*_Q;
   double RR  = mtau*mtau/QQ; double R = sqrt(RR);
   float U = 0.5*(3*cospsiLF()*cospsiLF() - 1)*(1 - RR);
   float V = 0.5*(3*cospsiLF()*cospsiLF() - 1)*(1 + RR)*costhetaLF() + 0.5*3*2*cospsiLF()* sinpsiLF()*sinthetaLF()*R;
   double B = 0.5*(3*cosbeta()*cosbeta() - 1);
   double fact =0;
   if(type == "bar") fact =1;
   double gA =  WA()*(costhetaLF()*(RR - 2)   - B*V)/3                                                               +     fact*WA()*0.5*R*sinthetaLF()*cosalpha()*2*sinbeta()*cosbeta();
   double gC =  WC()*0.5*V*sinbeta()*sinbeta()* cos2gamma()                                                 -      fact*WC()*R*sinthetaLF()*sinbeta()*(sinalpha()*sin2gamma()  -  cos2gamma()*cosalpha()*cosbeta() ) ;
   double gD = -WD()*0.5*V*sinbeta()*sinbeta()* sin2gamma()                                                 -      fact*WD()*R*sinthetaLF()*sinbeta()*(sinalpha()*cos2gamma() + sin2gamma()* cosalpha()*cosbeta()  );
   double gE = - WE()*cosbeta()*( costhetaLF()*cospsiLF() + R*sinthetaLF()*sinpsiLF())             +     fact*WE()*R*sinthetaLF()*sinbeta()*cosalpha();


   double res = gA+gC+gD+gE;
   return res;
 }

 double PolarimetricA1::vgetfscalar(TString type){
   double QQ=_Q*_Q;
   double RR  = mtau*mtau/QQ; double R = sqrt(RR);
   float U = 0.5*(3*cospsiLF()*cospsiLF() - 1)*(1 - RR);
   float V = 0.5*(3*cospsiLF()*cospsiLF() - 1)*(1 + RR)*costhetaLF() + 0.5*3*2*cospsiLF()* sinpsiLF()*sinthetaLF()*R;
   double B = 0.5*(3*cosbeta()*cosbeta() - 1);

   double fA =  WA()*(2+RR + B*U)/3;
   double fC = -WC()*0.5*U*sinbeta()*sinbeta()* cos2gamma();
   double fD = WD()*0.5*U*sinbeta()*sinbeta()* sin2gamma();
   double fE =  WE()*cospsiLF()*cosbeta();
   double fSA = WSA()*RR;
   double fSB = WSB()*RR*cospsiLF()*sinbeta()*cosgamma();
   double fSC = 0;
   double fSD = -WSD()*RR*cospsiLF()*sinbeta()*singamma();
   double fSE =0;

   double res = fA+fC+fD+fE  + fSA + fSB + fSC  + fSD + fSE;

   return res;
 }
 double PolarimetricA1::vgetgscalar(TString type){
   double QQ=_Q*_Q;
   double RR  = mtau*mtau/QQ; double R = sqrt(RR);
   float U = 0.5*(3*cospsiLF()*cospsiLF() - 1)*(1 - RR);
   float V = 0.5*(3*cospsiLF()*cospsiLF() - 1)*(1 + RR)*costhetaLF() + 0.5*3*2*cospsiLF()* sinpsiLF()*sinthetaLF()*R;
   double B = 0.5*(3*cosbeta()*cosbeta() - 1);
   double fact =0;
   if(type == "bar") fact =1;

  
   double gA =  WA()*(costhetaLF()*(RR - 2)   - B*V)/3                                                               +     fact*WA()*0.5*R*sinthetaLF()*cosalpha()*2*sinbeta()*cosbeta();
   double gC =  WC()*0.5*V*sinbeta()*sinbeta()* cos2gamma()                                                 -      fact*WC()*R*sinthetaLF()*sinbeta()*(sinalpha()*sin2gamma()  -  cos2gamma()*cosalpha()*cosbeta() ) ;
   double gD = -WD()*0.5*V*sinbeta()*sinbeta()* sin2gamma()                                                 -      fact*WD()*R*sinthetaLF()*sinbeta()*(sinalpha()*cos2gamma() + sin2gamma()* cosalpha()*cosbeta()  );
   double gE = - WE()*cosbeta()*( costhetaLF()*cospsiLF() + R*sinthetaLF()*sinpsiLF())             +     fact*WE()*R*sinthetaLF()*sinbeta()*cosalpha();
   double gSA =WSA()*RR*costhetaLF();
   double gSB =WSB()*R*(R*cospsiLF()*costhetaLF()*sinbeta()*cosgamma() + sinthetaLF()* ( sinpsiLF()*sinbeta()*cosgamma()  -  cosbeta()* cosalpha()* cosgamma() + sinalpha()*singamma())   );
   double gSC = WSC()*R*sinthetaLF()*(cosbeta()*sinalpha()*cosgamma() + cosalpha()*singamma());
   double gSD = WSD()*R*(sinthetaLF()*(cosbeta()*cosalpha()*singamma() + sinalpha()*cosgamma() - sinpsiLF()*sinbeta()*singamma()  )      - R*costhetaLF()*cospsiLF()*sinbeta()*singamma() );
   double gSE = -WSE()*R*sinthetaLF()*(cosbeta()*sinalpha()*singamma() -  cosalpha()*cosgamma());
   double res = gA+gC+gD+gE;

   return res;
 }

 double PolarimetricA1::TRF_vgetf(TString type){
   double QQ=_Q*_Q;
   double RR  = mtau*mtau/QQ; double R = sqrt(RR);

   double cb = TRF_cosbeta();     double ct = costhetaLF();    double ca = TRF_cosalpha();   double cg = TRF_cosgamma();  
   double sb = TRF_sinbeta();     double st =  sinthetaLF();    double sa = TRF_sinalpha();   double sg = TRF_singamma();  
   double s2g  = 2*sg*cg; double c2g = cg*cg - sg*sg;
   double Bb = 0.5*(cb*cb + 1);
   double fact=0;
   if(type=="scalar") fact=1;

   double fA =  WA()*(Bb*(1 - RR) + RR);
   double fC = -WC()*0.5*sb*sb*c2g*(1- RR);
   double fD = WD()*0.5*(1-RR)*sb*sb*s2g;
   double fE =  WE()*cb;
   double fSA = WSA()*RR;
   double fSB = WSB()*RR*sb*cg;
   double fSC = 0;
   double fSD = -WSD()*RR*sb*sg;
   double fSE = 0;
   double res = fA+fC+fD+fE+fact*(fSA+fSB+ fSC+fSD+fSE);

   return res;
 }
 double PolarimetricA1::TRF_vgetg(TString type){
   double QQ=_Q*_Q;
   double RR  = mtau*mtau/QQ; double R = sqrt(RR);

   double cb = TRF_cosbeta();     double ct = costhetaLF();    double ca = TRF_cosalpha();   double cg = TRF_cosgamma();  
   double sb = TRF_sinbeta();     double st = sinthetaLF();    double sa = TRF_sinalpha();   double sg = TRF_singamma();  
   double s2g  = 2*sg*cg; double c2g = cg*cg - sg*sg;
   double s2b  = 2*sb*cb; double c2b = cb*cb - sb*sb;
   double Bb = 0.5*(cb*cb + 1);
   double fact=0;
   if(type=="scalar") fact=1;

   double gA =  WA()*(RR*ct - Bb*ct*(1+RR)) + WA()*0.5*R*st*s2b*ca ;
   double gC =  WC()*(0.5*ct*(1+RR)  *sb*sb*c2g)  - WC()*R*st*sb*( sa*s2g - c2g*ca*cb );
   double gD = -WD()*(0.5*ct*(1+RR)*sb*sb*s2g)  - WD()*R*st*sb*( sa*c2g + s2g*ca*cb  );
   double gE =  -WE()*ct*cb + WE()*R*st*sb*ca ;

   double gSA = WSA()*RR*ct;
   double gSB = WSB()*(RR*ct*sb*cg - R*st*(cb*ca*cg - sa*sg  ));
   double gSC = WSC()*R*st*(cb*sa*cg + ca*sg);
   double gSD = WSD()*(R*st*(cb*ca*sg + sa*cg) - RR*ct*sb*sg);
   double gSE = -WSE()*R*st*(cb*sa*sg - ca*cg);

   double res = gA+gC+gD+gE+  fact*(gSA + gSB + gSC+ gSD + gSE);

   return res;
 }

void PolarimetricA1::debugger(){
   double QQ=_Q*_Q;
   double RR  = mtau*mtau/QQ; double R = sqrt(RR);

   double cb = TRF_cosbeta();     double ct = costhetaLF();    double ca = TRF_cosalpha();   double cg = TRF_cosgamma();  
   double sb = TRF_sinbeta();     double st =  sinthetaLF();    double sa = TRF_sinalpha();   double sg = TRF_singamma();  
   double s2g  = 2*sg*cg; double c2g = cg*cg - sg*sg;
   double Bb = 0.5*(cb*cb + 1);

   double fA =  WA()*(Bb*(1 - RR) + RR);
   double fC = -WC()*0.5*sb*sb*c2g*(1- RR);
   double fD = WD()*0.5*(1-RR)*sb*sb*s2g;
   double fE =  WE()*cb;
   double fSA = WSA()*RR;
   double fSB = WSB()*RR*sb*cg;
   double fSC = 0;
   double fSD = -WSD()*RR*sb*sg;
   double fSE = 0;
   double s2b  = 2*sb*cb; double c2b = cb*cb - sb*sb;

 
   double gA =  WA()*(RR*ct - Bb*ct*(1+RR)) + WA()*0.5*R*st*s2b*ca ;
   double gC = WC()*(ct*(1+RR)  *s2b*c2g)  - WC()*R*st*sb*( sa*s2g - c2g*ca*cb    );
   double gD = -WD()*(ct*(1+RR)*sb*sb*s2g) +WD()* R*st*sb*( sa*c2g + s2g*ca*cb  );
   double gE =  WE()*(R*st*sb*ca) - WE()*ct*cb;

   double gSA = WSA()*RR*ct;
   double gSB = WSB()*(RR*ct*sb*sg - R*st*(cb*ca*cg - sa*sg  ));
   double gSC = WSC()*R*st*(cb*sa*cg + ca*sg);
   double gSD = WSD()*(R*st*(cb*ca*sg + sa*cg) - RR*ct*sb*sg);
   double gSE = -WSE()*st*(cb*sa*sg - ca*cg);

   std::cout<<"  TRF_f"<<std::endl; 
   std::cout<<" fa + fc + fd + fe   "<< fA+fC+fD+fE;  std::cout<<"   TRF g non alpha   "<< WA()*(RR*ct - Bb*ct*(1+RR))+  WC()*(ct*(1+RR)  *s2b*c2g) -WD()*(ct*(1+RR)*sb*sb*s2g)- WE()*ct*cb <<std::endl;
   std::cout<<" fA   "<<  fA   << "    gA  + gAa  " << WA()*(RR*ct - Bb*ct*(1+RR)) <<"  +  " << WA()*0.5*R*st*s2b*ca<<std::endl;
   std::cout<<" fC  "<<  fC   << "     gC  + gCa  " <<  WC()*(ct*(1+RR)  *s2b*c2g) <<"  +  " << - WC()*R*st*sb*( sa*s2g - c2g*ca*cb    )<<std::endl;
   std::cout<<" fD    "<<  fD   << "    gD  + gDa " << -WD()*(ct*(1+RR)*sb*sb*s2g)<<"  +  " << WD()* R*st*sb*( sa*c2g + s2g*ca*cb  )<<std::endl;
   std::cout<<" fE    "<<  fE   << "    gE  + gEa  " << - WE()*ct*cb<<"  +  " << WE()*(R*st*sb*ca)<<std::endl;

   std::cout<<"  TRF_g"<<std::endl; 
   std::cout<<" ga + gc + gd + ge   "<< gA+gC+gD+gE<<std::endl;
}

double PolarimetricA1::TRF_vgetA1omega(TString type){
  if(TRF_vgetf(type)==0){ if(debug){std::cout<<"Warning!  Can not return vomegascalar; f(0)=0; return -50;  "<<std::endl;} return -50;}
  return TRF_vgetg(type)/TRF_vgetf(type);
}

double PolarimetricA1::result(){
  //  PolarimetricVector().Vect().Print();
   return nTAlongZLabFrame()*PolarimetricVector().Vect();
}


double PolarimetricA1::vgetA1omegascalar(TString type){
  if(vgetfscalar(type)==0){ if(debug){std::cout<<"Warning!  Can not return vomegascalar; f(0)=0; return -5;  "<<std::endl;} return -5;}
  return vgetgscalar(type)/vgetfscalar(type);
}
double PolarimetricA1::vgetA1omega(TString type){
  if(vgetf(type)==0){ if(debug){std::cout<<"Warning!  Can not return vomega; f(0)=0; return -5;  "<<std::endl; }return -5;}
  return vgetg(type)/vgetf(type);
}
double PolarimetricA1::getA1omegaBar(){
  if(getf()==0){ if(debug){std::cout<<"Warning!  Can not return omega; f(0)=0; return -5;  "<<std::endl;} return -5;}
  return getg()/getf();
}
double
PolarimetricA1::getA1omega(){
  double QQ=_Q*_Q;
  double RR  = mtau*mtau/QQ;
  float U = 0.5*(3*cospsiLF()*cospsiLF() - 1)*(1 - RR);
  float V = 0.5*(3*cospsiLF()*cospsiLF() - 1)*(1 + RR)*costhetaLF() + 0.5*3*2*cospsiLF()* sinpsiLF()*sinthetaLF()*sqrt(RR);
  
  float fa1 = (2  + RR + 0.5*(3*cosbeta()*cosbeta()- 1)*U)*WA()/3 - 0.5*sinbeta()*sinbeta()*cos2gamma()*U*WC() + 0.5*sinbeta()*sinbeta()*sin2gamma()*U*WD() + cospsiLF()*cosbeta()*WE();
  float ga1 = (costhetaLF()*(RR -2) - 0.5*(3*cosbeta()*cosbeta() - 1)*V)*WA()/3 + 0.5*sinbeta()*sinbeta()*cos2gamma()*V*WC() - 0.5*sinbeta()*sinbeta()*sin2gamma()*V*WD() -cosbeta()*(costhetaLF()*cospsiLF() + sinthetaLF()*sinpsiLF()*sqrt(RR))*WE();

  double omega = ga1/fa1;
  if(omega > 0 or omega < 0) 	return omega;
  return -999;
}
TLorentzVector
PolarimetricA1::sLV(){
  double QQ = _Q*_Q;
  double l0 = 0.5*(mtau*mtau + QQ)/sqrt(QQ);
  double l   = 0.5*(mtau*mtau  - QQ)/sqrt(QQ);
  return TLorentzVector(sinthetaLF(),0,-l0*costhetaLF()/mtau,-l*costhetaLF()/mtau);
}


TVector3
PolarimetricA1::nPerp(){
  if(_ss1pionLV.Vect().Cross(_ss2pionLV.Vect()).Mag()==0){if(debug){ std::cout<<"  Can not return nPerp, same sign pions seem to be parallel in a1 rest frame, return 0,0,0  "<<std::endl;} return TVector3(0,0,0);}

  TVector3 nss1= _ss1pionLV.Vect()*(1/_ss1pionLV.Vect().Mag());
  TVector3 nss2= _ss2pionLV.Vect()*(1/_ss2pionLV.Vect().Mag());
  return   (nss1.Cross(nss2))*(1/(nss1.Cross(nss2)).Mag());
}

TVector3
PolarimetricA1::ns(){
  return   sLV().Vect()*(1/sLV().Vect().Mag());
}
TVector3
PolarimetricA1::nL(){
  return   -LFa1LV.Vect()*(1/LFa1LV.Vect().Mag());
}
TVector3
PolarimetricA1::nT(){
  return   _tauLV.Vect()*(1/_tauLV.Vect().Mag());
}
TVector3
PolarimetricA1::nTAlongZLabFrame(){
  return _tauAlongZLabFrame.Vect()*(1/_tauAlongZLabFrame.Vect().Mag());
}


double  
PolarimetricA1::TRF_cosalpha(){
   TVector3 nTCrossns  = nT().Cross(ns());
   TVector3 nTCrossnPerp  = nT().Cross(nPerp());

   if(nTCrossns.Mag() ==0 || nTCrossnPerp.Mag() ==0){if(debug){std::cout<<" Can not compute cos alpha, one denominator is 0, return TRF cos alpha =0  "<< std::endl; }return 0;}
  return nTCrossns.Dot(nTCrossnPerp)/nTCrossns.Mag()/nTCrossnPerp.Mag();
}
double  
PolarimetricA1::TRF_sinalpha(){
   TVector3 nTCrossns  = nT().Cross(ns());
   TVector3 nTCrossnPerp  = nT().Cross(nPerp());

   if(nTCrossns.Mag() ==0 || nTCrossnPerp.Mag() ==0){if(debug){std::cout<<" Can not compute sin alpha, one denominator is 0, return TRF sin alpha =0  "<< std::endl; }return 0;}
  return -ns().Dot(nTCrossnPerp)/nTCrossns.Mag()/nTCrossnPerp.Mag();

}


double PolarimetricA1::TRF_cosbeta(){
  return nT().Dot(nPerp());
}
double PolarimetricA1::TRF_sinbeta(){
  if(TRF_cosbeta()*TRF_cosbeta() > 1 ){if(debug){std::cout<<"Warning! Can not compute TRF sin beta! return 0"<<std::endl;} return 0;}
  return sqrt(1 - TRF_cosbeta()*TRF_cosbeta());
}

double PolarimetricA1::TRF_cosgamma(){
  TVector3 nTCrossnPerp  = nT().Cross(nPerp());

  TVector3 qvect = _osPionLV.Vect()*(1/_osPionLV.Vect().Mag());
  //  qvect.Print();
  if(nTCrossnPerp.Mag()==0) { if(debug){std::cout<<"Warning! Can not compute TRF cos gamma, denominator =0, return 0  "<< std::endl;} return 0; }
  return -nT()*qvect/nTCrossnPerp.Mag();
}

double PolarimetricA1::TRF_singamma(){
  TVector3 nTCrossnPerp  = nT().Cross(nPerp());
  TVector3 qvect = _osPionLV.Vect()*(1/_osPionLV.Vect().Mag());

  if(nTCrossnPerp.Mag()==0) { if(debug){std::cout<<"Warning! Can not compute TRF  sin gamma, denominator =0, return 0  "<< std::endl;} return 0; }
  return qvect*nTCrossnPerp/nTCrossnPerp.Mag();
}


 // double  TRF_cosbeta();      double  TRF_cosalpha();   double  TRF_cosgamma();  
 // double TRF_sinbeta();        double TRF_sinalpha();    double  TRF_singamma();  

double PolarimetricA1::cosalpha(){
   TVector3 nLCrossnT  = nL().Cross(nT());
   TVector3 nLCrossnPerp  = nL().Cross(nPerp());

   if(nLCrossnPerp.Mag() ==0 || nLCrossnT.Mag() ==0){if(debug){std::cout<<" Can not compute cos alpha, one denominator is 0, return cos alpha =0  "<< std::endl;} return 0;}
  return nLCrossnT.Dot(nLCrossnPerp)/nLCrossnT.Mag()/nLCrossnPerp.Mag();
}
double PolarimetricA1::sinalpha(){
  TVector3 nLCrossnT  = nL().Cross(nT());
  TVector3 nLCrossnPerp  = nL().Cross(nPerp());
  if(nLCrossnPerp.Mag() ==0 || nLCrossnT.Mag() ==0){if(debug){std::cout<<" Can not compute sin alpha, one denominator is 0, return sin alpha =0  "<< std::endl; }return 0;}
  return -nT().Dot(nLCrossnPerp)/nLCrossnT.Mag()/nLCrossnPerp.Mag();
}
double PolarimetricA1::cosbeta(){
  return nL().Dot(nPerp());
}
double PolarimetricA1::sinbeta(){
  if(cosbeta()*cosbeta() > 1 ){if(debug){std::cout<<"Warning! Can not compute sin beta! return 0"<<std::endl;} return 0;}
  return sqrt(1 - cosbeta()*cosbeta());
}

double PolarimetricA1::cosgamma(){
  TVector3 nLCrossnPerp  = nL().Cross(nPerp());

  TVector3 qvect = _osPionLV.Vect()*(1/_osPionLV.Vect().Mag());
  //  qvect.Print();
  if(nLCrossnPerp.Mag()==0) { if(debug){std::cout<<"Warning! Can not compute cos gamma, denominator =0, return 0  "<< std::endl; }return 0; }
  return -nL()*qvect/nLCrossnPerp.Mag();
}

double PolarimetricA1::singamma(){
  TVector3 nLCrossnPerp  = nL().Cross(nPerp());
  TVector3 qvect = _osPionLV.Vect()*(1/_osPionLV.Vect().Mag());

  if(nLCrossnPerp.Mag()==0) { if(debug){std::cout<<"Warning! Can not compute cos gamma, denominator =0, return 0  "<< std::endl;} return 0; }
  return qvect*nLCrossnPerp/nLCrossnPerp.Mag();
}



TComplex 
PolarimetricA1::Conjugate(TComplex a){
  return TComplex(a.Re(), -a.Im());
}
TMatrixT<double> PolarimetricA1::convertToMatrix(TVectorT<double> V){
  TMatrixT<double> M(V.GetNrows(),1);
  for(int i=0; i < M.GetNrows(); i++){
    M(i,0)=V(i);
  } return M;
}
TVector3
PolarimetricA1::Rotate(TVector3 LVec, TVector3 Rot){
  TVector3 vec = LVec;
  vec.RotateZ(0.5*TMath::Pi() - Rot.Phi());  // not 0.5, to avoid warnings about 0 pT
  vec.RotateX(Rot.Theta());
  return vec;
}



//------- L-wave BreightWigner for rho
TComplex 
PolarimetricA1::BWIGML(double S, double M,  double G, double m1, double m2, int L){
  int IPOW;
  double MP = pow(m1+m2,2);
  double MM = pow(m1-m2,2);
  double MSQ = M*M;
  double W = sqrt(S);
  double WGS =0.0;
  double QS,QM;
  if(W > m1+m2){
    QS = sqrt(fabs( (S  - MP)*(S  - MM)))/W;
    QM = sqrt(fabs( (MSQ - MP)*(MSQ - MM)))/M;
    IPOW = 2*L +1;
    WGS=G*(MSQ/W)*pow(QS/QM, IPOW);
  }

  // std::cout<<" MSQ  "<<MSQ <<std::endl;
  // std::cout<<" WGS:  "<<G*(MSQ/W) <<" *  " <<QS/QM <<std::endl;
 TComplex out;
 out = TComplex(MSQ,0)/TComplex(MSQ - S, -WGS) ;
 return out;
}

TComplex
PolarimetricA1::FPIKM(double W, double XM1, double XM2){
  double ROM  = 0.773;
  double ROG  = 0.145;
  double ROM1 = 1.370;
  double ROG1 = 0.510;
  double BETA1=-0.145;
  double PIM  = 0.140;
  
  double S=W*W;
  int L =1; // P-wave
  TComplex out = (BWIGML(S,ROM,ROG,XM1,XM2,L) + BETA1*BWIGML(S,ROM1,ROG1,XM1,XM2,L))/(1+BETA1);
  return out;
  
} 


TComplex
PolarimetricA1::F3PI(double IFORM,double QQ,double SA,double SB){
  double MRO = 0.7743;
  double GRO = 0.1491;
  double MRP = 1.370 ;
  double GRP = 0.386 ;
  double MF2 = 1.275;
  double GF2 = 0.185;
  double MF0 = 1.186;
  double GF0 = 0.350;
  double MSG = 0.860;
  double GSG = 0.880;
  double MPIZ = mpi0;
  double MPIC = mpi;


  double M1 = mpi;
  double M2 = mpi;
  double M3 = mpi;

  double M1SQ = M1*M1;
  double M2SQ = M2*M2;
  double M3SQ = M3*M3;
  
  
  TComplex  BT1 = TComplex(1.,0.);
  TComplex  BT2 = TComplex(0.12,0.)*TComplex(1, 0.99*TMath::Pi(), true);//  TComplex(1, 0.99*TMath::Pi(), true);   Real part must be equal to one, stupid polar implemenation in root
  TComplex  BT3 = TComplex(0.37,0.)*TComplex(1, -0.15*TMath::Pi(), true);
  TComplex  BT4 = TComplex(0.87,0.)*TComplex(1, 0.53*TMath::Pi(), true);
  TComplex  BT5 = TComplex(0.71,0.)*TComplex(1, 0.56*TMath::Pi(), true);
  TComplex  BT6 = TComplex(2.10,0.)*TComplex(1, 0.23*TMath::Pi(), true);
  TComplex  BT7 = TComplex(0.77,0.)*TComplex(1, -0.54*TMath::Pi(), true);

  TComplex  F3PIFactor(0.,0.); // initialize to zero

  if(IFORM == 1 || IFORM == 2 ){
  double S1 = SA;
  double S2 = SB;
  double S3 = QQ-SA-SB+M1SQ+M2SQ+M3SQ;
  // std::cout<<"S1  "<< BT1  <<std::endl;
  // std::cout<<"S2  "<< BT5  <<std::endl;
  // std::cout<<"S3  "<< BT6  <<std::endl;
  //Lorentz invariants for all the contributions:
  double F134 = -(1./3.)*((S3-M3SQ)-(S1-M1SQ));
  double F15A = -(1./2.)*((S2-M2SQ)-(S3-M3SQ));
  double F15B = -(1./18.)*(QQ-M2SQ+S2)*(2.*M1SQ+2.*M3SQ-S2)/S2;
  double F167 = -(2./3.);

  // Breit Wigners for all the contributions:

 
  TComplex  FRO1 = BWIGML(S1,MRO,GRO,M2,M3,1);
  TComplex  FRP1 = BWIGML(S1,MRP,GRP,M2,M3,1);
  TComplex  FRO2 = BWIGML(S2,MRO,GRO,M3,M1,1);
  TComplex  FRP2 = BWIGML(S2,MRP,GRP,M3,M1,1);
  TComplex  FF21 = BWIGML(S1,MF2,GF2,M2,M3,2);
  TComplex  FF22 = BWIGML(S2,MF2,GF2,M3,M1,2);
  TComplex  FSG2 = BWIGML(S2,MSG,GSG,M3,M1,0);
  TComplex  FF02 = BWIGML(S2,MF0,GF0,M3,M1,0);

  // std::cout<<"FRO1  "<< FRO1  <<std::endl;
  // std::cout<<" FRP1 "<< FRP1  <<std::endl;
  // std::cout<<"FRO2  "<< FRO2  <<std::endl;
  // std::cout<<"FRP2  "<< FRP2  <<std::endl;
  // std::cout<<"FF21  "<< FF21  <<std::endl;
  // std::cout<<"FF22  "<<FF22   <<std::endl;
  // std::cout<<"FSG2  "<<FSG2   <<std::endl;
  // std::cout<<"FF02  "<<FF02   <<std::endl;


  F3PIFactor = BT1*FRO1+BT2*FRP1+
    BT3*TComplex(F134,0.)*FRO2+BT4*TComplex(F134,0.)*FRP2
    -BT5*TComplex(F15A,0.)*FF21-BT5*TComplex(F15B,0.)*FF22
    -BT6*TComplex(F167,0.)*FSG2-BT7*TComplex(F167,0.)*FF02;

  } else if (IFORM == 3 ){

    double S3 = SA;
    double S1 = SB;
    double S2 = QQ-SA-SB+M1SQ+M2SQ+M3SQ;

    double F34A = (1./3.)*((S2-M2SQ)-(S3-M3SQ));
    double F34B = (1./3.)*((S3-M3SQ)-(S1-M1SQ));
    double F35A = -(1./18.)*(QQ-M1SQ+S1)*(2.*M2SQ+2.*M3SQ-S1)/S1;
    double F35B =  (1./18.)*(QQ-M2SQ+S2)*(2.*M3SQ+2.*M1SQ-S2)/S2;
    double F36A = -(2./3.);
    double F36B =  (2./3.);

    //C Breit Wigners for all the contributions:
    TComplex  FRO1 = BWIGML(S1,MRO,GRO,M2,M3,1);
    TComplex  FRP1 = BWIGML(S1,MRP,GRP,M2,M3,1);
    TComplex  FRO2 = BWIGML(S2,MRO,GRO,M3,M1,1);
    TComplex  FRP2 = BWIGML(S2,MRP,GRP,M3,M1,1);
    TComplex  FF21 = BWIGML(S1,MF2,GF2,M2,M3,2);
    TComplex  FF22 = BWIGML(S2,MF2,GF2,M3,M1,2);
    TComplex  FSG1 = BWIGML(S1,MSG,GSG,M2,M3,0);
    TComplex  FSG2 = BWIGML(S2,MSG,GSG,M3,M1,0);
    TComplex  FF01 = BWIGML(S1,MF0,GF0,M2,M3,0);
    TComplex  FF02 = BWIGML(S2,MF0,GF0,M3,M1,0);
    
    F3PIFactor = 
      BT3*(TComplex(F34A,0.)*FRO1+TComplex(F34B,0.)*FRO2)+
      BT4*(TComplex(F34A,0.)*FRP1+TComplex(F34B,0.)*FRP2)
      -BT5*(TComplex(F35A,0.)*FF21+TComplex(F35B,0.)*FF22)
      -BT6*(TComplex(F36A,0.)*FSG1+TComplex(F36B,0.)*FSG2)
      -BT7*(TComplex(F36A,0.)*FF01+TComplex(F36B,0.)*FF02);
    
    // F3PIFactor = TComplex(0.,0.);
 

  }

  TComplex FORMA1 = FA1A1P(QQ);
  TComplex F3PIFactor_ret =  F3PIFactor*FORMA1;

  return F3PIFactor;
} 


TComplex
PolarimetricA1::FA1A1P(double XMSQ){
  double  XM1 = 1.275000;
  double  XG1 =0.700 ; 
  double  XM2 = 1.461000 ;
  double  XG2 = 0.250; 
  TComplex BET = TComplex(0.00,0.);

  double GG1 = XM1*XG1/(1.3281*0.806);
  double GG2 = XM2*XG2/(1.3281*0.806);
  double XM1SQ = XM1*XM1;
  double XM2SQ = XM2*XM2;

  double GF = WGA1(XMSQ);
  double FG1 = GG1*GF;
  double FG2 = GG2*GF;
  TComplex F1 = TComplex(-XM1SQ,0.0)/TComplex(XMSQ-XM1SQ,FG1);
  TComplex F2 = TComplex(-XM2SQ,0.0)/TComplex(XMSQ-XM2SQ,FG2);
  TComplex FA1A1P = F1+BET*F2;

  return FA1A1P;
}


double 
PolarimetricA1::WGA1(double QQ){
// C mass-dependent M*Gamma of a1 through its decays to 
// C.   [(rho-pi S-wave) + (rho-pi D-wave) + 
// C.    (f2 pi D-wave) + (f0pi S-wave)]
// C.  AND simple K*K S-wave
  double  MKST = 0.894;
  double  MK = 0.496;
  double  MK1SQ = (MKST+MK)*(MKST+MK);
  double  MK2SQ = (MKST-MK)*(MKST-MK);
  //C coupling constants squared:
  double   C3PI = 0.2384*0.2384;
  double   CKST = 4.7621*4.7621*C3PI;
// C Parameterization of numerical integral of total width of a1 to 3pi.
// C From M. Schmidtler, CBX-97-64-Update.
  double  S = QQ;
  double  WG3PIC = WGA1C(S);
  double  WG3PIN = WGA1N(S);

  //C Contribution to M*Gamma(m(3pi)^2) from S-wave K*K, if above threshold
  double  GKST = 0.0;
  if(S > MK1SQ) GKST = sqrt((S-MK1SQ)*(S-MK2SQ))/(2.*S);

  double WGA1_ret = C3PI*(WG3PIC+WG3PIN)+CKST*GKST;
  return WGA1_ret;
}


double
PolarimetricA1::WGA1C(double S){
  double STH,Q0,Q1,Q2,P0,P1,P2,P3,P4,G1_IM;

  Q0 =   5.80900; Q1 =  -3.00980; Q2 =   4.57920;
  P0 = -13.91400; P1 =  27.67900; P2 = -13.39300;
  P3 =   3.19240; P4 =  -0.10487; STH=0.1753;


  if(S < STH){
    G1_IM=0.0;
  }else if(S > STH && S < 0.823){
    G1_IM = Q0*   pow(S-STH,3)   *(1. + Q1*(S-STH) + Q2*pow(S-STH,2));
  }
  else{
    G1_IM = P0 + P1*S + P2*S*S+ P3*S*S*S + P4*S*S*S*S;
  }

  double WGA1C_ret = G1_IM;
  return WGA1C_ret;
}

double
PolarimetricA1::WGA1N(double S){
  double STH,Q0,Q1,Q2,P0,P1,P2,P3,P4,G1_IM;
  Q0 =   6.28450;Q1 =  -2.95950;Q2 =   4.33550;
  P0 = -15.41100;P1 =  32.08800;P2 = -17.66600;
  P3 =   4.93550;P4 =  -0.37498;STH   = 0.1676;
  
  if(S < STH){
    G1_IM = 0.0;
  }else if(S > STH && S < 0.823){
    G1_IM = Q0*pow(S-STH,3)*(1. + Q1*(S-STH) + Q2*pow(S-STH,2));
  }
  else{
    G1_IM = P0 + P1*S + P2*S*S+ P3*S*S*S + P4*S*S*S*S;
  }
  double WGA1N_ret = G1_IM;
  return WGA1N_ret;
}


// C=======================================================================
//       COMPLEX FUNCTION FA1A1P(XMSQ)
// C     ==================================================================
// C     complex form-factor for a1+a1prime.                       AJW 1/98
// C     ==================================================================

//       REAL XMSQ
//       REAL PKORB,WGA1
//       REAL XM1,XG1,XM2,XG2,XM1SQ,XM2SQ,GG1,GG2,GF,FG1,FG2
//       COMPLEX BET,F1,F2
//       INTEGER IFIRST/0/

//       IF (IFIRST.EQ.0) THEN
//         IFIRST = 1

// C The user may choose masses and widths that differ from nominal:


// C scale factors relative to nominal:
//         GG1 = XM1*XG1/(1.3281*0.806)
//         GG2 = XM2*XG2/(1.3281*0.806)

//         XM1SQ = XM1*XM1
//         XM2SQ = XM2*XM2
//       END IF

//       GF = WGA1(XMSQ)
//       FG1 = GG1*GF
//       FG2 = GG2*GF
//       F1 = CMPLX(-XM1SQ,0.0)/CMPLX(XMSQ-XM1SQ,FG1)
//       F2 = CMPLX(-XM2SQ,0.0)/CMPLX(XMSQ-XM2SQ,FG2)
//       FA1A1P = F1+BET*F2

//       RETURN
//       END
// C=======================================================================




//       COMPLEX FUNCTION FPIKM(W,XM1,XM2)
// C **********************************************************
// C     PION FORM FACTOR
// C **********************************************************
//       COMPLEX BWIGM
//       REAL ROM,ROG,ROM1,ROG1,BETA1,PI,PIM,S,W
//       EXTERNAL BWIG
//       DATA  INIT /0/
// C
// C ------------ PARAMETERS --------------------
//       IF (INIT.EQ.0 ) THEN
//       INIT=1
//       PI=3.141592654
//       PIM=.140
//       ROM=0.773
//       ROG=0.145
//       ROM1=1.370
//       ROG1=0.510
// C      BETA1=-0.145
//       BETA1=0.
//       ENDIF
// C -----------------------------------------------
//       S=W**2
//       FPIKM=(BWIGM(S,ROM,ROG,XM1,XM2)+BETA1*BWIGM(S,ROM1,ROG1,XM1,XM2))
//      & /(1+BETA1)
// C      FPIKM=(BWIGM(S,ROM1,ROG1,XM1,XM2))
// C     & /(1+BETA1)
//       RETURN
//       END



//       COMPLEX FUNCTION BWIGML(S,M,G,M1,M2,L)
// C **********************************************************
// C     L-WAVE BREIT-WIGNER
// C **********************************************************
//       REAL S,M,G,M1,M2
//       INTEGER L,IPOW
//       REAL MSQ,W,WGS,MP,MM,QS,QM

//       MP = (M1+M2)**2
//       MM = (M1-M2)**2
//       MSQ = M*M
//       W = SQRT(S)
//       WGS = 0.0
//       IF (W.GT.(M1+M2)) THEN
//         QS=SQRT(ABS((S   -MP)*(S   -MM)))/W
//         QM=SQRT(ABS((MSQ -MP)*(MSQ -MM)))/M
//         IPOW = 2*L+1
//         WGS=G*(MSQ/W)*(QS/QM)**IPOW
//       ENDIF

//       BWIGML=CMPLX(MSQ,0.)/CMPLX(MSQ-S,-WGS)

//       RETURN
//       END










// C=======================================================================
//       COMPLEX FUNCTION FA1A1P(XMSQ)
// C     ==================================================================
// C     complex form-factor for a1+a1prime.                       AJW 1/98
// C     ==================================================================

//       REAL XMSQ
//       REAL PKORB,WGA1
//       REAL XM1,XG1,XM2,XG2,XM1SQ,XM2SQ,GG1,GG2,GF,FG1,FG2
//       COMPLEX BET,F1,F2
//       INTEGER IFIRST/0/

//       IF (IFIRST.EQ.0) THEN
//         IFIRST = 1

// C The user may choose masses and widths that differ from nominal:
//         XM1 = PKORB(1,10)
//         XG1 = PKORB(2,10)
//         XM2 = PKORB(1,17)
//         XG2 = PKORB(2,17)
//         BET = CMPLX(PKORB(3,17),0.)
// C scale factors relative to nominal:
//         GG1 = XM1*XG1/(1.3281*0.806)
//         GG2 = XM2*XG2/(1.3281*0.806)

//         XM1SQ = XM1*XM1
//         XM2SQ = XM2*XM2
//       END IF

//       GF = WGA1(XMSQ)
//       FG1 = GG1*GF
//       FG2 = GG2*GF
//       F1 = CMPLX(-XM1SQ,0.0)/CMPLX(XMSQ-XM1SQ,FG1)
//       F2 = CMPLX(-XM2SQ,0.0)/CMPLX(XMSQ-XM2SQ,FG2)
//       FA1A1P = F1+BET*F2

//       RETURN
//       END
// C=======================================================================

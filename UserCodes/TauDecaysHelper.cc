#include "TauDecaysHelper.h"
#include <iostream>
TauDecaysHelper::TauDecaysHelper(){ }
TauDecaysHelper::TauDecaysHelper(vector<TLorentzVector> TauAndProd, string type){
  if(type=="muon" || type =="pion"){
    if(TauAndProd.size()!=2){
      std::cout<<" Warning!! Size of input vector  !=2  !! "<< " type:  "<<type<<std::endl;
    }
  }
  if(type=="rho"){
    if(TauAndProd.size()!=3){
      std::cout<<" Warning!! Size of input vector  !=3  !! "<< " type:  "<<type<<std::endl;
    }
  }

  Configure(TauAndProd,  type);
}
void 
TauDecaysHelper::Configure(vector<TLorentzVector> TauAndProd, string type){
  TauLV = TauAndProd.at(0);
  type_=type;
  if(type_ == "muon" || type_=="pion"){
    ProductLV= TauAndProd.at(1);
  }
  if(type_ == "rho" ){
    TauRhoPi  = TauAndProd.at(1);
    TauRhoPi0= TauAndProd.at(2);
    ProductLV = TauRhoPi+TauRhoPi0;
  }
  InvisibleLV = TauLV - ProductLV;
}
bool  
TauDecaysHelper::isConfigured(){
  if(TauLV.E()!=0 && ProductLV.E()!=0) return true; return false;
}
TauDecaysHelper::~TauDecaysHelper(){
}
TLorentzVector 
TauDecaysHelper::Boost(TLorentzVector pB, TLorentzVector frame){
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
TauDecaysHelper::getOmega(){
  double omega=-999;
  if(type_=="pion" || type_=="muon"){
    omega = 2*ProductLV.E()/TauLV.E() - 1;
  }
  if(type_=="rho"){
    omega = (TauRhoPi.E() -   TauRhoPi0.E())/(TauRhoPi.E() +  TauRhoPi0.E());
  }
  return omega;
}

TMatrixT<double> TauDecaysHelper::convertToMatrix(TVectorT<double> V){
  TMatrixT<double> M(V.GetNrows(),1);
  for(int i=0; i < M.GetNrows(); i++){
    M(i,0)=V(i);
  } return M;
}

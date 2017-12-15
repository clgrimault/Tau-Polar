#include "rhoHelper.h"
#include <iostream>

rhoHelper::rhoHelper(){
}

rhoHelper::rhoHelper(vector<TLorentzVector> TauRhoandProd, int helicity)
//  toplot(false)
{
  if(TauRhoandProd.size()!=3){
    std::cout<<" Warning!! Size of rho input vector != 4 !! "<<std::endl;
  }
  TLorentzVector fakeboost(0,0,0,0);
  Setup(TauRhoandProd,fakeboost,helicity);
}

rhoHelper::rhoHelper(vector<TLorentzVector> TauRhoandProd, TLorentzVector RefernceFrame, int helicity){
  if(TauRhoandProd.size()!=3){
    std::cout<<" Warning!! Size of Rho input vector != 3 !! "<<std::endl;
  }
  Setup(TauRhoandProd,RefernceFrame,helicity);
}

void 
rhoHelper::Setup(vector<TLorentzVector> TauRhoandProd, TLorentzVector ReferenceFrame, int helicity){
   mpi   = 0.13957018; // GeV 
   mpi0 = 0.1349766;   // GeV
   mtau = 1.77687; // GeV
   coscab = 0.975; 
   mrho = 0.773; // GeV
   debug=false;
   for(int i=0; i<TauRhoandProd.size(); i++){
     TauRhoandProd_RF.push_back(Boost(TauRhoandProd.at(i),ReferenceFrame));
   }



   TauHelicity = -helicity;







    // //   std::cout<<"tauy helicity  "<< helicity <<std::endl;
    //   TauLV         = TauRhoandProd.at(0);
    //   TauRhoPi      = TauRhoandProd.at(1);
    //   TauRhoPi0     = TauRhoandProd.at(2);
    //   ProductLV     = TauRhoPi+TauRhoPi0;
    //   // std::cout<<"Rho heloper  initial input: "<<std::endl;


    //   // std::cout<<"tau ";  TauLV.Print();
    //   // std::cout<<"pi ";  TauRhoPi.Print();
    //   // std::cout<<"pi0 ";  TauRhoPi0.Print();
    //   // std::cout<<"rho ";  ProductLV.Print();

    //   TVector3 Rot1 = TauLV.Vect();

    //   TLorentzVector TauLVR1    = TauLV;
    //   TLorentzVector TauRhoPiR1 = TauRhoPi;
    //   TLorentzVector TauRhoPi0R1= TauRhoPi0;
    //   TLorentzVector ProductLVR1= ProductLV;
 

      // std::cout<<"--------Rho heloper  first rotation tau along z: "<<std::endl;
      // TauLVR1.SetVect(Rotate(TauLVR1.Vect(),Rot1));
      // TauRhoPiR1.SetVect(Rotate(TauRhoPiR1.Vect(),Rot1));
      // TauRhoPi0R1.SetVect(Rotate(TauRhoPi0R1.Vect(),Rot1));
      // ProductLVR1.SetVect(Rotate(ProductLV.Vect(),Rot1));
      // TLorentzVector InvisibleR1= TauLVR1 - ProductLVR1;


      //      std::cout<<"helicity  "<< helicity << std::endl;
        // std::cout<<"tau ";  TauLVR1.Print();
	// std::cout<<"pi ";  TauRhoPiR1.Print();
	// std::cout<<"pi0 ";  TauRhoPi0R1.Print();
	// std::cout<<"rho ";  ProductLVR1.Print();
	// std::cout<<"neutrino ";  InvisibleR1.Print();

      // std::cout<<"========Rho heloper  boost to tauRF "<<std::endl;
      // TLorentzVector TauLVR1_TRF = Boost(TauLVR1,TauLVR1);
      // TLorentzVector TauRhoPiR1_TRF = Boost(TauRhoPiR1,TauLVR1);
      // TLorentzVector TauRhoPi0R1_TRF= Boost(TauRhoPi0R1,TauLVR1);
      // TLorentzVector ProductLVR1_TRF= Boost(ProductLVR1,TauLVR1);
      // TLorentzVector Invisible_TRF= Boost(InvisibleR1,TauLVR1);

     
      // //   std::cout<<"tau helicity  " << helicity << "   cos(theta)   " <<cos(Invisible_TRF.Theta())  << "  theta in grads  " <<  Invisible_TRF.Theta()*180./TMath::Pi()<<std::endl;
      // toplot=false;
      // if(Invisible_TRF.Theta()*180./TMath::Pi()< 90){toplot=true;
      // 	//	vec.RotateX(DPF_TauLV.Vect().Theta());
      // 	//	ProductLVR1_TRF.RotateX(TMath::Pi()-Invisible_TRF.Theta()*2);
      // 	//	Invisible_TRF.RotateX(TMath::Pi()-Invisible_TRF.Theta()*2);
      // }


      //  std::cout<<"tau ";  TauLVR1_TRF.Print();
// //  std::cout<<"pi ";  TauRhoPiR1_TRF.Print();
// //  std::cout<<"pi0 ";  TauRhoPi0R1_TRF.Print();
 // std::cout<<"rho ";  ProductLVR1_TRF.Print();
 // std::cout<<"neutrino ";  Invisible_TRF.Print();

   //   // TVector3 Rot2 = Invisible_TRF.Vect();

   //   // TLorentzVector TauLVR2 = TauLVR1_TRF;
   //   // TLorentzVector TauRhoPiR2 =TauRhoPiR1_TRF ;
   //   // TLorentzVector TauRhoPi0R2= TauRhoPi0R1_TRF;
   //   // TLorentzVector ProductLVR2= ProductLVR1_TRF;
   //   // TLorentzVector InvisibleR2=Invisible_TRF ;

   //   // DPF_TauRhoPi  =TauRhoPiR1_TRF;
   //   // DPF_TauRhoPi0 =TauRhoPi0R1_TRF;
   //   // InvisibleLV   =Invisible_TRF;
   //   // DPF_TauLV     =TauLVR1_TRF;
   //   // ProductLV     =ProductLVR1_TRF;



   // //  TauLVR2.SetVect(Rotate(TauLVR2.Vect(),Rot2));
   // //  TauRhoPiR2.SetVect(Rotate(TauRhoPiR2.Vect(),Rot2));
   // //  TauRhoPi0R2.SetVect(Rotate(TauRhoPi0R2.Vect(),Rot2));
   // //  ProductLVR2.SetVect(Rotate(ProductLVR2.Vect(),Rot2));
   // //  InvisibleR2.SetVect(Rotate(InvisibleR2.Vect(),Rot2));


   // //  // std::cout<<" >>>>>>>>>>> Rho heloper  second  rotation neutrino along z: "<<std::endl;
   // //  // std::cout<<"tau ";  TauLVR2.Print();
   // //  // std::cout<<"pi ";  TauRhoPiR2.Print();
   // //  // std::cout<<"pi0 ";  TauRhoPi0R2.Print();
   // //  // std::cout<<"rho ";  ProductLVR2.Print();
   // //  // std::cout<<"neutrino ";  InvisibleR2.Print();

   //   // TLorentzVector TauLVR1_TRF = Boost(TauLVR1,TauLVR1);
   //   // TLorentzVector TauRhoPiR1_TRF = Boost(TauRhoPiR1,TauLVR1);
   //   // TLorentzVector TauRhoPi0R1_TRF= Boost(TauRhoPi0R1,TauLVR1);
   //   // TLorentzVector ProductLVR1_TRF= Boost(ProductLVR1,TauLVR1);
   //   // TLorentzVector Invisible_TRF= Boost(InvisibleR1,TauLVR1);


      // TLorentzVector TauLVR2_HRF = Boost(TauLVR1_TRF,ProductLVR1_TRF);
      // TLorentzVector TauRhoPiR2_HRF = Boost(TauRhoPiR1_TRF,ProductLVR1_TRF);
      // TLorentzVector TauRhoPi0R2_HRF= Boost(TauRhoPi0R1_TRF,ProductLVR1_TRF);
      // TLorentzVector ProductLVR2_HRF= Boost(ProductLVR1_TRF,ProductLVR1_TRF);
      // TLorentzVector Invisible_HRF= Boost(Invisible_TRF,ProductLVR1_TRF);

   //    // std::cout<<"========Rho heloper  boost to HF "<<std::endl;
   //    // std::cout<<"tau ";  TauLVR2_HRF.Print();
   //    // std::cout<<"rho from prev  ";  ProductLVR2.Print();
   //    // std::cout<<"pi ";  TauRhoPiR2_HRF.Print();
   //    // std::cout<<"pi0 ";  TauRhoPi0R2_HRF.Print();
   //    // std::cout<<"rho ";  ProductLVR2_HRF.Print();
   //    // std::cout<<"neutrino ";  Invisible_HRF.Print();


       // DPF_TauRhoPi  =TauRhoPiR2_HRF;
       // DPF_TauRhoPi0 =TauRhoPi0R2_HRF;
       // InvisibleLV   =Invisible_HRF;
       // DPF_TauLV     =TauLVR2_HRF;
       // ProductLV     =ProductLVR1_TRF;

       // 
        // std::cout<<"rho in TRF  ";
        // ProductLV.Print();
	// Invisible_TRF.Print();

       // DPF_TauRhoPi.Print();
       // DPF_TauRhoPi0.Print();
   // //  DPF_TauRhoPi  =TauRhoPiR2;
   // //  DPF_TauRhoPi0 =TauRhoPi0R2;
   // //  InvisibleLV   =Invisible_HRF;
   // //  DPF_TauLV     =TauLVR1;
   // //  ProductLV     =ProductLVR1;





   // // // // std::cout<<"Angle beta 1:  "<< DPF_TauRhoPi.Angle(-ProductLV.Vect())<<std::endl;
   // // // // std::cout<<"Angle beta 2:  "<< DPF_TauRhoPi.Angle(DPF_TauLV.Vect())<<std::endl;






   // //--------------------------------------------------  test
   // // std::cout<<"helicity "<< helicity <<std::endl;
   //   TauLV         = TauRhoandProd.at(0);
   //   TauRhoPi      = TauRhoandProd.at(1);
   //   TauRhoPi0     = TauRhoandProd.at(2);
   //   ProductLV     = TauRhoPi+TauRhoPi0;
   //   TVector3 Rot1 = TauLV.Vect();

   //   TLorentzVector TauLVR1    = TauLV;
   //   TLorentzVector TauRhoPiR1 = TauRhoPi;
   //   TLorentzVector TauRhoPi0R1= TauRhoPi0;
   //   TLorentzVector ProductLVR1= ProductLV;
 

   //   // std::cout<<"--------Rho heloper  first rotation tau along z: "<<std::endl;
   //   TauLVR1.SetVect(Rotate(TauLVR1.Vect(),Rot1));
   //   TauRhoPiR1.SetVect(Rotate(TauRhoPiR1.Vect(),Rot1));
   //   TauRhoPi0R1.SetVect(Rotate(TauRhoPi0R1.Vect(),Rot1));
   //   ProductLVR1.SetVect(Rotate(ProductLV.Vect(),Rot1));
   //   TLorentzVector InvisibleR1= TauLVR1 - ProductLVR1;
   //   ProductLV=ProductLVR1;
   //   DPF_TauRhoPi  = Boost(TauRhoPiR1,ProductLVR1);
   //   DPF_TauRhoPi0 = Boost(TauRhoPi0R1,ProductLVR1);
   //   InvisibleLV = TauLV - ProductLVR1;
   //   DPF_TauLV=  Boost(TauLVR1,ProductLVR1);

    // std::cout<<" boosted  ";  DPF_TauLV.Print();




   //--------------------------------------------------  below is the previous configuration
   // std::cout<<"helicity "<< helicity <<std::endl;
    TauLV         = TauRhoandProd.at(0);
    TauRhoPi      = TauRhoandProd.at(1);
    TauRhoPi0     = TauRhoandProd.at(2);
    ProductLV     = TauRhoPi+TauRhoPi0;

    DPF_TauRhoPi  = Boost(TauRhoPi,ProductLV);
    DPF_TauRhoPi0 = Boost(TauRhoPi0,ProductLV);
    InvisibleLV = TauLV - ProductLV;
    DPF_TauLV=  Boost(TauLV,ProductLV);

    //    std::cout<<" boosted  ";  DPF_TauLV.Print();

    DPF_InvisibleLV=  Boost(InvisibleLV,ProductLV);

    TVector3 RotVector = DPF_TauLV.Vect();
    DPF_TauLV.SetVect(Rotate(DPF_TauLV.Vect(),RotVector));
    DPF_TauRhoPi.SetVect(Rotate(DPF_TauRhoPi.Vect(),RotVector));
    DPF_TauRhoPi0.SetVect(Rotate(DPF_TauRhoPi0.Vect(),RotVector));
    ProductLV.SetVect(Rotate(ProductLV.Vect(),RotVector)); //  rotate nL

    // std::cout<<"helicity  "<< helicity << std::endl;
    // std::cout<<"rho in TRF  "<< helicity << std::endl;

    // DPF_TauRhoPi.Print();
    // DPF_TauRhoPi0.Print();


}

void 
rhoHelper::Configure(vector<TLorentzVector> TauRhoandProd, int helicity){
  if(TauRhoandProd.size()!=3){
    std::cout<<" Warning!! Size of input vector != 4 !! "<<std::endl;
  }
  TLorentzVector fakeboost(0,0,0,0);
  Setup(TauRhoandProd,fakeboost,helicity);
}

void 
rhoHelper::Configure(vector<TLorentzVector> TauRhoandProd, TLorentzVector RefernceFrame, int helicity){
  if(TauRhoandProd.size()!=3){
    std::cout<<" Rho helper:  Warning!! Size of input vector != 3!   Size = "<< TauRhoandProd.size()<<std::endl;
  }
  Setup(TauRhoandProd,RefernceFrame,helicity);
}

bool
rhoHelper::isConfigured(){
  if(TauRhoandProd_RF.size()!=3){ std::cout<<"Error:   rhoHelper is not Configured! Check  the size of input vector!  Size =  "<< TauRhoandProd_RF.size() <<std::endl; return false;} return true;
}

rhoHelper::~rhoHelper(){
}

TLorentzVector 
rhoHelper::Boost(TLorentzVector pB, TLorentzVector frame){
   TMatrixT<double> transform(4,4);
   TMatrixT<double> result(4,1);
   TVectorT<double> vec(4); 
   TVector3 b;
   if(frame.Vect().Mag()==0){ std::cout<<"RH Boost is not set, perfrom calculation in the Lab Frame   "<<std::endl; return pB;}
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
TVector3
rhoHelper::Rotate(TVector3 LVec, TVector3 Rot){
  TVector3 vec = LVec;
  vec.RotateZ(0.5*TMath::Pi() - Rot.Phi());  // not 0.5, to avoid warnings about 0 pT
  vec.RotateX(Rot.Theta());
  return vec;
}

  double
  rhoHelper::getCosbetaRho(){
    double cb=-999;
    cb = (mrho/sqrt(mrho*mrho - 4*mpi*mpi))  *   (TauRhoPi.E() -   TauRhoPi0.E())/(TauRhoPi.E() +  TauRhoPi0.E()) ;
    if(fabs(cb) - 1 >0  && fabs(cb) < 1.01 &&  cb> 0  )  cb= 1;
    if(fabs(cb) - 1 >0  && fabs(cb) < 1.01 &&  cb< 0  )  cb=-1;
    if(fabs(cb) > 1 ){if(debug){std::cout<<"Warning! Cos beta > 1:  "<<cb <<std::endl;  }}
    return cb;
  }
 // double
 // rhoHelper::getCosbetaRho(){
 //   double cb=-999;
 //   cb = (mrho/sqrt(mrho*mrho - 4*mpi*mpi))  *   (DPF_TauRhoPi.E() -   DPF_TauRhoPi0.E())/(DPF_TauRhoPi.E() +  DPF_TauRhoPi0.E()) ;
 //   if(fabs(cb) - 1 >0  && fabs(cb) < 1.01 &&  cb> 0  )  cb= 1;
 //   if(fabs(cb) - 1 >0  && fabs(cb) < 1.01 &&  cb< 0  )  cb=-1;
 //   if(fabs(cb) > 1 ){if(debug){std::cout<<"Warning! Cos beta > 1:  "<<cb <<std::endl;  }}
 //   // std::cout<<"getCoseBeta "<< cb <<std::endl;
 //   return cb;
 // }



double
rhoHelper::getSinbetaRho(){
  double sb=-999;
  if(fabs(getCosbetaRho()) > 1 ){if(debug){std::cout<<"Warning! Cos beta > 1:  "<< getCosbetaRho()<<std::endl;  }return 0;}
  sb = sqrt(1- getCosbetaRho()*getCosbetaRho());
  return sb;
}

double
rhoHelper::getCosthetaRho(){
  double ct=-999;
  double QQ = ProductLV.M2();
  double x = ProductLV.E()/TauLV.E();
  double s = 4*TauLV.E()*TauLV.E();
  if(fabs(ct) - 1 >0  && fabs(ct) < 1.01 &&  ct> 0  )  ct = 1;
  if(fabs(ct) - 1 >0  && fabs(ct) < 1.01 &&  ct< 0  )  ct=-1;
  if( 1 - 4*mtau*mtau/s  <= 0 ){std::cout<<"Warning! In costheta root square <=0! return -999"<<std::endl; return ct;}
  ct= (2*x*mtau*mtau - mtau*mtau - QQ)/( (mtau*mtau - QQ)*sqrt(1 - 4*mtau*mtau/s) );
  if(fabs(ct) > 1 ){if(debug){std::cout<<"Warning! Cos theta > 1:  "<<ct<<std::endl; }}
  return ct;
}

double
rhoHelper::getSinthetaRho(){
  double st=-999;
   st = sqrt(1- getCosthetaRho()*getCosthetaRho());
   if(fabs(getCosthetaRho()) > 1 ){if(debug){std::cout<<"Warning! Cos theta > 1"<<std::endl; }}
  return st;
}

double 
rhoHelper::getUltrarel_cospsiLF(){
  double cos=-999;
  double QQ = ProductLV.M2();
  cos = (getCosthetaRho()*(mtau*mtau  + QQ)   + (mtau*mtau  - QQ))/(getCosthetaRho()*(mtau*mtau  - QQ)   + (mtau*mtau  + QQ));
  if(fabs(cos) - 1 >0  && fabs(cos) < 1.01 && cos > 0  ) cos = 1;
  if(fabs(cos) - 1 >0  && fabs(cos) < 1.01 && cos < 0  ) cos =-1;
  if(fabs(cos) > 1 )if(debug){std::cout<<"Warning! Cos psi > 1:  "<<cos<<std::endl; }
  return cos;
}

double 
rhoHelper::getSinpsiLF(){
  double sin = -999;
  sin = sqrt(1 - getUltrarel_cospsiLF()*getUltrarel_cospsiLF());
  if(getUltrarel_cospsiLF()*getUltrarel_cospsiLF() > 1  ){if(debug){std::cout<<"Warning! getultrarel_cospsiLF > 1"<<std::endl;}}
  return    sin;
}

TMatrixT<double> rhoHelper::convertToMatrix(TVectorT<double> V){
  TMatrixT<double> M(V.GetNrows(),1);
  for(int i=0; i < M.GetNrows(); i++){
    M(i,0)=V(i);
  } return M;
}
 
TLorentzVector
rhoHelper::sLV(){
  double QQ = ProductLV.M2();
  double l0 = 0.5*(mtau*mtau + QQ)/sqrt(QQ);
  double l   = 0.5*(mtau*mtau  - QQ)/sqrt(QQ);
  // std::cout<<"l "<< l<<" l0  "<< l0 << " costheta   "<<getCosthetaRho() <<std::endl;
  //  return TLorentzVector(getSinpsiLF(),0,-l0*getCosthetaRho()/mtau,-l*getCosthetaRho()/mtau);
  TLorentzVector tspin(getSinthetaRho(),0,-l0*TauHelicity*getCosthetaRho()/mtau,-l*getCosthetaRho()/mtau);

  // tspin.SetVect(Rotate(tspin.Vect(),DPF_TauLV.Vect()));
  return tspin;

   //  return TLorentzVector(getSinthetaRho(),0,l0*getCosthetaRho()/mtau,l*getCosthetaRho()/mtau);
}

TVector3
rhoHelper::nPerp(){
  return   DPF_TauRhoPi.Vect()*(1/DPF_TauRhoPi.Vect().Mag());
}


// TVector3
// rhoHelper::ns(){
//   return   sLV().Vect()*(1/sLV().Vect().Mag());
// }



TVector3
rhoHelper::ns(){  
  double QQ = ProductLV.M2();
  double l0 = 0.5*(mtau*mtau + QQ)/sqrt(QQ);
  double l   = 0.5*(mtau*mtau  - QQ)/sqrt(QQ);
  return   sLV().Vect()*(1/sLV().Vect().Mag());
}


TVector3
rhoHelper::nL(){
  return   -ProductLV.Vect()*(1/ProductLV.Vect().Mag());
}
 
TVector3
rhoHelper::nT(){
  return   DPF_TauLV.Vect()*(1/DPF_TauLV.Vect().Mag());
}
double rhoHelper::getCosBetaTest(){
  // std::cout<< toplot<<std::endl;
  //std::cout<<" nPerp()*nL()  "<< nPerp()*nL() <<std::endl;
  if(toplot) return nPerp()*nL();
  return -100.;
}
double 
rhoHelper::getSinBetaTest(){
  if(fabs(getCosBetaTest())> 1){std::cout<<" getCosBetaTest > 1"<< std::endl; return 0;}
  return sqrt(1-getCosBetaTest()*getCosBetaTest());
}



TVector3
rhoHelper::nTRotated(){
  TVector3 vec= DPF_TauLV.Vect();
  vec.RotateZ(0.5*TMath::Pi() - vec.Phi());
  vec.RotateX(vec.Theta());
  return   vec;
}
TVector3
rhoHelper::nPerpRotated(){
  TVector3 vec= DPF_TauRhoPi.Vect();
  vec.RotateZ(0.5*TMath::Pi() - DPF_TauLV.Vect().Phi());
  vec.RotateX(DPF_TauLV.Vect().Theta());
  return   vec;
}

double 
rhoHelper::TFK_cosbeta(){
  //std::cout<<" nT()*nPerp()  "<< nT()*nPerp() <<std::endl;

    return nT()*nPerp();
  //  return nT()*nPerp();
  //  return nL()*nPerp();
}
double 
rhoHelper::TFK_sinbeta(){
  if(fabs(TFK_cosbeta())> 1){std::cout<<" TFK_cosbeta > 1"<< std::endl; return 0;}
  return sqrt(1-TFK_cosbeta()*TFK_cosbeta());
}


double rhoHelper::TFK_costheta(){
  return -ns()*nT();
}
double rhoHelper::TFK_sintheta(){
  if(fabs(TFK_costheta())> 1){std::cout<<" TFK_costheta > 1"<< std::endl; return 0;}
  return sqrt(1-TFK_costheta()*TFK_costheta());
}

 
double  
rhoHelper::DPF_cosalpha(){
     TVector3 nTCrossns  = nT().Cross(ns());
     TVector3 nTCrossnPerp  = nT().Cross(nPerp());
     //  std::cout<<"-------- "<<std::endl;  
      // std::cout<<"nT "<<std::endl;  nT().Print();
      // std::cout<<"ns "<<std::endl;  ns().Print();
      // std::cout<<"nPerp "<<std::endl;  nPerp().Print();
      // (DPF_TauRhoPi0 + DPF_TauRhoPi).Print();
      // std::cout<<"rho Mass "<< (DPF_TauRhoPi0 + DPF_TauRhoPi).M() <<std::endl;

      // std::cout<<" Pi "<<std::endl;  DPF_TauRhoPi0.Print();
      // std::cout<<" Pi0 "<<std::endl;  DPF_TauRhoPi.Print();

       //  std::cout<<"beta "<<     TFK_cosbeta()<< "  " <<   getCosbetaRho() <<std::endl;





     //    std::cout<<"nT cross ns "<<std::endl;  nTCrossns.Print();
     //  std::cout<<"nT cross nPepr "<<std::endl;  nTCrossnPerp.Print();
     // std::cout<<""<<std::endl; std::cout<<std::endl; std::cout<<std::endl; 
     // std::cout<<"  TFK theta   "<< TFK_costheta()<< " costehta   "<< getCosthetaRho()<<std::endl; 
       // std::cout<<"  thete   "<< getCosthetaRho()<<std::endl; 
     // std::cout<<"     "<< sLV()*DPF_TauLV<<std::endl; 

     //  TFK_costheta() 

     //     sLV().Print();
     // DPF_TauLV.Print();

    if(nTCrossns.Mag() ==0 || nTCrossnPerp.Mag() ==0){if(debug){std::cout<<" Can not compute cos alpha, one denominator is 0, return DPF cos alpha =0  "<< std::endl; }return 0;}
   return nTCrossns.Dot(nTCrossnPerp)/nTCrossns.Mag()/nTCrossnPerp.Mag();
}


double  
rhoHelper::DPF_cosalphaTest(){
     TVector3 nTCrossns  = nL().Cross(nT());
     TVector3 nTCrossnPerp  = nL().Cross(nPerp());
     //  nL().Print();
    if(nTCrossns.Mag() ==0 || nTCrossnPerp.Mag() ==0){if(debug){std::cout<<" Can not compute cos alpha, one denominator is 0, return DPF cos alpha =0  "<< std::endl; }return 0;}
   return nTCrossns.Dot(nTCrossnPerp)/nTCrossns.Mag()/nTCrossnPerp.Mag();
}


 // double  
 // rhoHelper::DPF_cosalpha(){
 //      TVector3 nTCrossns  = nL().Cross(nT());
 //      TVector3 nTCrossnPerp  = nL().Cross(nPerp());
 //      // std::cout<<"-------- "<<std::endl;  
 //      // std::cout<<"nT "<<std::endl;  nT().Print();
 //      // std::cout<<"ns "<<std::endl;  ns().Print();
 //      // std::cout<<"nT cross ns "<<std::endl;  nTCrossns.Print();
 //      // std::cout<<""<<std::endl; std::cout<<std::endl; std::cout<<std::endl; 
 //        // std::cout<<"  TFK theta   "<< TFK_costheta()<<std::endl; 
 //        // std::cout<<"  thete   "<< getCosthetaRho()<<std::endl; 
 //      //   std::cout<<"     "<< sLV()*DPF_TauLV<<std::endl; 
 //      //     sLV().Print();
 //      // DPF_TauLV.Print();

 //     if(nTCrossns.Mag() ==0 || nTCrossnPerp.Mag() ==0){if(debug){std::cout<<" Can not compute cos alpha, one denominator is 0, return DPF cos alpha =0  "<< std::endl; }return 0;}
 //    return nTCrossns.Dot(nTCrossnPerp)/nTCrossns.Mag()/nTCrossnPerp.Mag();
 // }





double  
rhoHelper::DPF_sinalphaTest(){
    TVector3 nTCrossns  = nT().Cross(ns());
    TVector3 nTCrossnPerp  = nL().Cross(nPerp());

    if(nTCrossns.Mag() ==0 || nTCrossnPerp.Mag() ==0){if(debug){std::cout<<" Can not compute sin alpha, one denominator is 0, return DPF sin alpha =0  "<< std::endl; }return 0;}
    return -nT().Dot(nTCrossnPerp)/nTCrossns.Mag()/nTCrossnPerp.Mag();
 }





double  
rhoHelper::DPF_sinalpha(){
    TVector3 nTCrossns  = nT().Cross(ns());
    TVector3 nTCrossnPerp  = nT().Cross(nPerp());

    if(nTCrossns.Mag() ==0 || nTCrossnPerp.Mag() ==0){if(debug){std::cout<<" Can not compute sin alpha, one denominator is 0, return DPF sin alpha =0  "<< std::endl; }return 0;}
    return -ns().Dot(nTCrossnPerp)/nTCrossns.Mag()/nTCrossnPerp.Mag();
 }


double 
rhoHelper::getOmegaRho(){
  double omega=-999;
   omega = TFK_cosbeta();
  //  omega = getCosbetaRho();
  // omega = DPF_cosalpha();

  if(  isinf(fabs(omega)) ||  isnan(fabs(omega))) omega  = -999.;
  return omega;
}


// //----------------------  angles beta  + theta
//   double 
//   rhoHelper::getOmegaRhoBar(){
//     double omega=-999;
//     double QQ =  ProductLV.M2();
//     double Be = 0.5*(3*getCosbetaRho()*getCosbetaRho() -1);
//     double Ps = 0.5*(3*getUltrarel_cospsiLF() -1);
//     double RR = mtau*mtau/QQ;
//     double R = sqrt(RR);
//     //    omega = ((-2 + RR + 2*(1+RR)*Ps*Be)*getCosthetaRho() + 3*sqrt(RR)*Be*getSinthetaRho()*2*getUltrarel_cospsiLF()*getSinpsiLF())  /  ( 2 +RR - 2*(1-RR)*Ps*Be);
//     omega = ((-2 + RR + 2*(1+RR)*Ps*Be)*getCosthetaRho() + 3*sqrt(RR)*Be*getSinthetaRho()*2*getUltrarel_cospsiLF()*getSinpsiLF() -R*getSinthetaRho()*2*getCosbetaRho()*getSinbetaRho()*DPF_cosalphaTest() )  /  ( 2 +RR - 2*(1-RR)*Ps*Be);

//     //  std::cout<<"  "<< R*getSinthetaRho()*2*TFK_cosbeta()*TFK_sinbeta()*DPF_cosalphaTest() <<std::endl;


//     if(  isinf(fabs(omega)) ||  isnan(fabs(omega))) omega  = -999.;
//     return omega;
//   }




  



        double 
        rhoHelper::getOmegaRhoBar(){
          double omega=-999;
          double QQ =  ProductLV.M2();
          double RR = mtau*mtau/QQ;
   
        	 omega = (RR*TFK_costheta() - sqrt(RR)*TFK_sintheta()*2* TFK_sinbeta()*TFK_cosbeta()*DPF_cosalpha() -   TFK_costheta()* TFK_sinbeta()*TFK_sinbeta()*(1+RR) )  /  ( RR + (1-RR)*TFK_sinbeta()*TFK_sinbeta());


        		 // omega = ( - sqrt(RR)*TFK_sintheta()*2* TFK_sinbeta()*TFK_cosbeta()*DPF_cosalpha() -   TFK_costheta()* TFK_sinbeta()*TFK_sinbeta()*(1+RR) )  /  ( RR + (1-RR)*TFK_sinbeta()*TFK_sinbeta());
       
        	 // omega = (RR*getCosthetaRho() - sqrt(RR)*getSinthetaRho()*2* TFK_sinbeta()*TFK_cosbeta()*DPF_cosalpha() -   getCosthetaRho()* TFK_sinbeta()*TFK_sinbeta()*(1+RR) )  /  ( RR + (1-RR)*TFK_sinbeta()*TFK_sinbeta());
       

        	 // std::cout<<"  "<< RR*TFK_costheta()/( RR + (1-RR)*TFK_sinbeta()*TFK_sinbeta()) <<"  "<< (sqrt(RR)*TFK_sintheta()*2* TFK_sinbeta()*TFK_cosbeta()*DPF_cosalpha() )/  ( RR + (1-RR)*TFK_sinbeta()*TFK_sinbeta()) <<std::endl;
         if(  isinf(fabs(omega)) ||  isnan(fabs(omega))) omega  = -999. ;
          return omega;
        }



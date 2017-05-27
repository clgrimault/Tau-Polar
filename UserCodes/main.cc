#include "MultiplyNumbers.h"
#include "AddNumbers.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "TLorentzVector.h"
#include "a1Helper.h"
#include <vector>

 
void PrintLV(TLorentzVector a){
  std::cout<<" px:    "<< a.Px() << "  py:   "<< a.Py() <<"  pz:   "<< a.Pz() <<"  E:   "<< a.E() << "  mass:   "<< a.M()<< std::endl;
}

int main(int argc, const char* argv[]) {

  int a, b;
  a = atoi(argv[1]);
  b = atoi(argv[2]);
  AddNumbers ab;
  ab.setA(a);
  ab.setB(b);
  printf("%d + %d = %d\n", ab.getA(), ab.getB(), ab.getSum());


  MultiplyNumbers atimesb;
  atimesb.setA(a);
  atimesb.setB(b);
  printf("%d * %d = %d\n", atimesb.getA(), atimesb.getB(), atimesb.getProduct());

  TLorentzVector vec(2,0,0,5);
  std::cout<<" vec.M()    "<< vec.M() << "  vec.Px()   "<< vec.Px() <<std::endl;

  TLorentzVector os(-3.664335,0.233661,-5.183266,6.353556);
  TLorentzVector ss1(-7.876100,0.720873,-12.047975,14.412696);
  TLorentzVector ss2(-3.472285,0.243649,-4.660688,5.818730);
  TLorentzVector ta1(-26.788939,2.687771,-37.247365,45.993420);
  TLorentzVector a1(-15.012722,1.198184,-21.891933,26.584986);
  TLorentzVector rf(4,1,3,175);
  TLorentzVector ta2(2,-4,6,55);
  TLorentzVector zer(0,0,0,0);


  TLorentzVector rho(-13.491274,-4.967516,1.272212,14.454947);
  TLorentzVector pi0(-8.682081,-2.818832,0.750293,9.160063) ;
  TLorentzVector pi(-4.809193,-2.148683,0.521918,5.294883);


// ---------
// (x,y,z,t)=(-26.788939,2.687771,-37.247365,45.993420) (P,eta,phi,E)=(45.959086,-1.128328,3.041596,45.993420)
// (x,y,z,t)=(-15.012722,1.198184,-21.891933,26.584986) (P,eta,phi,E)=(26.572057,-1.168748,3.061950,26.584986)
// (x,y,z,t)=(-3.664335,0.233661,-5.183266,6.353556) (P,eta,phi,E)=(6.352023,-1.144735,3.077913,6.353556)
// (x,y,z,t)=(-7.876100,0.720873,-12.047975,14.412696) (P,eta,phi,E)=(14.412020,-1.207630,3.050320,14.412696)
// (x,y,z,t)=(-3.472285,0.243649,-4.660688,5.818730) (P,eta,phi,E)=(5.817056,-1.101985,3.071538,5.818730)
// ---------
// (x,y,z,t)=(-32.822804,29.511723,-12.752927,45.979044) (P,eta,phi,E)=(45.944700,-0.285049,2.409263,45.979044)
// (x,y,z,t)=(-16.949156,15.002583,-7.017585,23.721868) (P,eta,phi,E)=(23.698056,-0.305267,2.417042,23.721868)
// (x,y,z,t)=(-5.037089,4.346876,-2.412471,7.078636) (P,eta,phi,E)=(7.077260,-0.355084,2.429615,7.078636)
// (x,y,z,t)=(-0.923205,0.855447,-0.441767,1.341169) (P,eta,phi,E)=(1.333887,-0.344161,2.394271,1.341169)
// (x,y,z,t)=(-10.988859,9.800258,-4.163346,15.302059) (P,eta,phi,E)=(15.301423,-0.279118,2.413307,15.302059)
// ---------
// (x,y,z,t)=(-38.702364,-21.588178,10.946964,45.682771) (P,eta,phi,E)=(45.648203,0.244574,-2.632781,45.682771)
// (x,y,z,t)=(-23.317578,-13.107163,7.037263,27.685411) (P,eta,phi,E)=(27.659180,0.260141,-2.629496,27.685411)
// (x,y,z,t)=(-11.382100,-6.650994,3.389942,13.612462) (P,eta,phi,E)=(13.611746,0.254395,-2.612769,13.612462)
// (x,y,z,t)=(-3.030154,-1.805861,1.235359,3.740128) (P,eta,phi,E)=(3.737523,0.343422,-2.604147,3.740128)
// (x,y,z,t)=(-8.905323,-4.650308,2.411961,10.332820) (P,eta,phi,E)=(10.331878,0.237834,-2.660348,10.332820)
// ---------
// (x,y,z,t)=(16.862003,42.748406,0.932397,45.997607) (P,eta,phi,E)=(45.963276,0.020288,1.195086,45.997607)
// (x,y,z,t)=(13.909979,36.136450,1.181131,38.752321) (P,eta,phi,E)=(38.739200,0.030499,1.203349,38.752321)
// (x,y,z,t)=(3.471275,9.245016,0.396963,9.884186) (P,eta,phi,E)=(9.883201,0.040187,1.211609,9.884186)
// (x,y,z,t)=(5.788296,15.810430,0.630423,16.849064) (P,eta,phi,E)=(16.848486,0.037435,1.219846,16.849064)
// (x,y,z,t)=(4.650407,11.081000,0.153744,12.019067) (P,eta,phi,E)=(12.018256,0.012793,1.173446,12.019067)
// ---------
// (x,y,z,t)=(-18.658958,41.971850,1.745631,45.999970) (P,eta,phi,E)=(45.965641,0.037995,1.989116,45.999970)
// (x,y,z,t)=(-18.431158,41.073397,1.637837,45.066871) (P,eta,phi,E)=(45.049018,0.036373,1.992600,45.066871)
// (x,y,z,t)=(-8.572447,19.482827,0.453507,21.290668) (P,eta,phi,E)=(21.290210,0.021304,1.985303,21.290668)
// (x,y,z,t)=(-2.981254,7.248164,0.512469,7.855308) (P,eta,phi,E)=(7.854068,0.065342,1.961016,7.855308)
// (x,y,z,t)=(-6.877456,14.342405,0.671861,15.920894) (P,eta,phi,E)=(15.920282,0.042227,2.017925,15.920894)
// ---------
// (x,y,z,t)=(-21.363342,-31.640470,-25.580263,45.989332) (P,eta,phi,E)=(45.954995,-0.627947,-2.164677,45.989332)
// (x,y,z,t)=(-8.665741,-13.056724,-10.431898,18.858951) (P,eta,phi,E)=(18.825451,-0.624333,-2.156742,18.858951)
// (x,y,z,t)=(-0.701739,-1.383549,-0.893003,1.795432) (P,eta,phi,E)=(1.789999,-0.547820,-2.040189,1.795432)
// (x,y,z,t)=(-2.848837,-4.717945,-3.658006,6.616295) (P,eta,phi,E)=(6.614823,-0.622695,-2.114027,6.616295)
// (x,y,z,t)=(-5.115166,-6.955231,-5.880890,10.447225) (P,eta,phi,E)=(10.446292,-0.637162,-2.204915,10.447225)
// ---------




  std::vector<TLorentzVector> particles;
  particles.push_back(ta1);
  particles.push_back(os);
  particles.push_back(ss1);
  particles.push_back(ss2);

  // particles.at(0).Print();
  // particles.at(1).Print();
  // particles.at(2).Print();
  // particles.at(3).Print();

  a1Helper Helper(particles,a1);
   Helper.Configure(os,ss1,ss2,ta1,ta2);
   //   Helper.costheta();
   //   PrintLV(ss2);
   //PrintLV(Helper.Boost(ss2,ss2));

// Helper.getBoosted().at(0).Print();
// Helper.getBoosted().at(1).Print();
// Helper.getBoosted().at(2).Print();
// Helper.getBoosted().at(3).Print();


  // cout<<"rho Mass  "<<rho.M()<<endl;
  // cout<<"pi0 Mass  "<<pi.M()<<endl;
  // cout<<"pi Mass  "<<pi0.M()<<endl;
  // rho.Print();
  // (pi+pi0).Print();
  // Helper.nL().Print();
   std::cout<< "Helper.cosbetaLF()  "<< Helper.cosbetaLF() <<std::endl;
  std::cout<< " cosbeta =  "<<Helper.cosbeta()<<std::endl;
    std::cout<< " cos gamma  "<<Helper.cosgammaLF() <<std::endl;
    std::cout<< " sin gamma "<<Helper.singammaLF() <<std::endl;
  // std::cout<< " 1 =   "<<sqrt( Helper.cosgammaLF()*Helper.cosgammaLF()+ Helper.singammaLF()*Helper.singammaLF()) <<std::endl;




    std::cout<< "RF cos gamma  "<<Helper.cosgamma() <<std::endl;
    std::cout<< "RF sin gamma "<<Helper.singamma() <<std::endl;
    cout<<"getf  "<<Helper.getf()<<endl;
    cout<<"getg  "<<Helper.getg()<<endl;


    // cout<<"ultrarel_cospsiLF()  "<<Helper.ultrarel_cospsiLF()<<std::endl;
    // cout<<"cospsiLF()  "<<Helper.cospsiLF()<<std::endl;

    // cout<<"F1  "<<Helper.F1()<<endl;
    // cout<<"F2  "<<Helper.F2()<<endl;
    // cout<<"F4  "<<Helper.F4()<<endl;

//  cout<<"rho Mass  "<<Helper.Boost(rho,rho).M()<<endl;
//   cout<<"rho Mass  "<<Helper.Boost(pi,rho).M()<<endl;
// cout<<"rho Mass  "<<Helper.Boost(pi0,rho).M()<<endl;
// cout<<" ---------  "<<endl;

 // Helper.Boost(rho,pi).Print();
 // Helper.Boost(pi + pi0,pi).Print();
 // Helper.Boost(pi,pi).Print();
 // Helper.Boost(pi0,pi).Print();


   // std::cout<< "scalar  "<< Helper.Scalar(ss1,ss2) <<std::endl;
   // std::cout<< "mass   "<< Helper.Mass() <<std::endl;
   // std::cout<< "mass   "<< Helper.Mass("a1") <<std::endl;

   // std::cout<< "width    "<< Helper.Widths(4,"rho") <<std::endl;
   // std::cout<<  " width    "<< Helper.Widths(4) <<std::endl;

 


   //  ss2.Print();// Helper.Boost(ss2,ss2).Print();
   // particles.at(1).Print();

   cout<<"WA()   "<<      Helper.WA()<<endl;	 
   cout<<"WC();  "<<	 Helper.WC()<<endl;	 
   cout<<"WD();  "<<	 Helper.WD()<<endl;	 
   cout<<"WE();  "<<	 Helper.WE()<<endl;	 
   cout<<"WSA();  "<<	 Helper.WSA()<<endl; 
   cout<<"WSB();  "<<	 Helper.WSB()<<endl; 
   cout<<"WSD();  "<<	 Helper.WSD()<<endl; 
   cout<<"WSC();  "<<	 Helper.WSC()<<endl; 
   cout<<"WSE();   "<<	 Helper.WSE()<<endl; 


   cout<<     Helper.WA()*Helper.WA()<< "  = "<<Helper.WC()*Helper.WC()+  Helper.WD()*Helper.WD()+ Helper.WE()*Helper.WE()<<endl;	 


}

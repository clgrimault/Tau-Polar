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

  TLorentzVector os(1,1,1,sqrt(0.139*0.139 + 1*1 + 1*1+ 1*1));
  TLorentzVector ss1(2,2,2,sqrt(0.139*0.139 + 2*2 + 2*2+ 2*2));
  TLorentzVector ss2(3,3,3,sqrt(0.139*0.139 + 3*3 + 3*3+ 3*3));
  TLorentzVector ta1(6,6,6,sqrt(1.777*1.777 + 6*6 + 6*6+ 6*6));
  TLorentzVector rf(4,1,3,175);
  TLorentzVector ta2(2,-4,6,55);
  TLorentzVector zer(0,0,0,0);


  TLorentzVector rho(-13.491274,-4.967516,1.272212,14.454947);
  TLorentzVector pi0(-8.682081,-2.818832,0.750293,9.160063) ;
  TLorentzVector pi(-4.809193,-2.148683,0.521918,5.294883);






  std::vector<TLorentzVector> particles;
  particles.push_back(ta1);
  particles.push_back(os);
  particles.push_back(ss1);
  particles.push_back(ss2);

  // particles.at(0).Print();
  // particles.at(1).Print();
  // particles.at(2).Print();
  // particles.at(3).Print();

  a1Helper Helper(particles,ta1);
   Helper.Configure(os,ss1,ss2,ta1,ta2);
 
   //   PrintLV(ss2);
   //PrintLV(Helper.Boost(ss2,ss2));

// Helper.getBoosted().at(0).Print();
// Helper.getBoosted().at(1).Print();
// Helper.getBoosted().at(2).Print();
// Helper.getBoosted().at(3).Print();


  cout<<"rho Mass  "<<rho.M()<<endl;
  cout<<"pi0 Mass  "<<pi.M()<<endl;
  cout<<"pi Mass  "<<pi0.M()<<endl;
  rho.Print();
  (pi+pi0).Print();
cout<<" ---------  "<<endl;
 cout<<"rho Mass  "<<Helper.Boost(rho,rho).M()<<endl;
  cout<<"rho Mass  "<<Helper.Boost(pi,rho).M()<<endl;
cout<<"rho Mass  "<<Helper.Boost(pi0,rho).M()<<endl;
cout<<" ---------  "<<endl;

 Helper.Boost(rho,pi).Print();
 Helper.Boost(pi + pi0,pi).Print();
 Helper.Boost(pi,pi).Print();
 Helper.Boost(pi0,pi).Print();


   // std::cout<< "scalar  "<< Helper.Scalar(ss1,ss2) <<std::endl;
   // std::cout<< "mass   "<< Helper.Mass() <<std::endl;
   // std::cout<< "mass   "<< Helper.Mass("a1") <<std::endl;

   // std::cout<< "width    "<< Helper.Widths(4,"rho") <<std::endl;
   // std::cout<<  " width    "<< Helper.Widths(4) <<std::endl;

 


   //  ss2.Print();// Helper.Boost(ss2,ss2).Print();
   // particles.at(1).Print();


}

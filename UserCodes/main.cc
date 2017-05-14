#include "MultiplyNumbers.h"
#include "AddNumbers.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "TLorentzVector.h"
#include "a1Helper.h"



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

  TLorentzVector os(2,1,3,35);
  TLorentzVector ss1(1,4,6,25);
  TLorentzVector ss2(4,2.4,1.5,15);
  TLorentzVector ta1(4,1,3,75);
  TLorentzVector ta2(2,-4,6,55);
  TLorentzVector zer(0,0,0,0);

  a1Helper Helper;
  Helper.Configure(os,ss1,ss2,ta1,ta2);


  PrintLV(ss2);
  PrintLV(Helper.Boost(ss2,ss2));

  std::cout<< "scalar  "<< Helper.Scalar(ss1,ss2) <<std::endl;
  std::cout<< "mass   "<< Helper.Mass() <<std::endl;
  std::cout<< "mass   "<< Helper.Mass("a1") <<std::endl;

  std::cout<< "width    "<< Helper.Widths(4,"rho") <<std::endl;
  std::cout<<  " width    "<< Helper.Widths(4) <<std::endl;



  // ss2.Print(); Helper.Boost(ss2,ss2).Print();
 


}

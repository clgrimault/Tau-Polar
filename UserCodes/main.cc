#include "MultiplyNumbers.h"
#include "AddNumbers.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "TLorentzVector.h"



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

}

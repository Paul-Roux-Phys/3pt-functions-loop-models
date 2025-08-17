#include <print>
#include <iostream>
#include <complex>

#include "OnLoops.hpp"

using namespace OnLoops;

unsigned int L;
using Weight = std::complex<double>;

int main() {
  L = 4;
  Parameters<Weight> params(2, 2, 0.5);
  
  LinkPattern<Weight> l({EMPTY, EMPTY, OPEN, CLOSE});

  OnState<Weight> s(10);
  s.add(l, 1.0);
  std::cout << s << std::endl;
  return 0;
}

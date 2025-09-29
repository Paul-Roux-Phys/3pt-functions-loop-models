#include "../include/SparseVector.hpp"
#include <cassert>

int main() {
  Vector<double> v(10);
  v.add(0, 10.);

  assert(v[0] == 10.);
  assert(v[1] == 0  );

  v.add(0, 2.);
  assert(v[0] == 12);

  v.add(3, 2.0);

  double a = 0;
  for (auto it = v.begin(); it != v.end(); it++) {
    a += it->second;
  }

  assert(a == 14);

  double b = v.inner_product(v);
  assert(b == 12*12 + 2*2);

  Vector<double> v2(12);
  v2.add(0, 2.0);
  v2.add(1, 2.0);
  v2.add(3, 3.0);

  assert(v.inner_product(v2) == 30.);

  assert(v.size() == 2);
  assert(v.bucket_count() == 16);

  Vector<double> v3(100);
  assert(v3.bucket_count() == 128);

  v3.clear();
  for (auto it = v.begin(); it != v.end(); ) {
    v3.add(it->first, it->second);
    it = v.erase(it);
  }
  assert(v3[0] == 12);
  assert(v3[3] ==  2);
  assert(v.empty());

  v3 *= 1/12.0;

  assert(v3[0] == 1 && v3[3] == 1.0/6.0);

  return 0;
}

#include "../include/SparseVector.hpp"
#include "../include/Matrix.hpp"
#include "../examples/OnLoops/OnLoops.hpp"
#include <complex>

unsigned int L = 4;
size_t hashsize = 1e4;

using namespace std;
using namespace OnLoops;

using Weight = double;

int min_defects = 2;
Parameters<Weight> params(min_defects, 2, 0.5);

class OnTransferMatrix : public Matrix<Weight> {
public:
  using Matrix<Weight>::Matrix;
  void sparseMultMv(Vector& v1, Vector& v2) override {
    auto *o1 = static_cast<OnState<Weight> *>(&v1);
    auto *o2 = static_cast<OnState<Weight> *>(&v2);
    o1->transfer(*o2, params);
  }
};

int main () {
  OnState<Weight> v(hashsize), v2(hashsize);

  LinkPattern<Weight> l({EMPTY, EMPTY, OPEN, CLOSE, EMPTY, EMPTY});
  v.add(l);
  v.transfer(v2, params);
  cout << v << endl;
  v.transfer(v2, params);
  cout << v << endl;
  v.transfer(v2, params);
  cout << v << endl;
  // v.insert_aux_space(v2);

  // cout << v << endl;

  // v.r_matrix(v2, 0, params);

  // cout << v << endl;
  // v.r_matrix(v2, 1, params);
  // v.r_matrix(v2, 2, params);
  // v.r_matrix(v2, 3, params);

  // cout << v << endl;

  // v.contract_aux_space(v2, params);

  // cout << v << endl;

  // OnTransferMatrix M(v, v2);

  // M.generateSpace();
  
  return 0;
}

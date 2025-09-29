#include "../examples/OnLoops/OnLoops.hpp"
#include "../include/BigFloats.hpp"
#include "../include/Matrix.hpp"
// #include "../include/EigenSolver.hpp"
#include <iomanip>

unsigned int L = 4;
size_t hashsize = 1.5e4;

using namespace std;
using namespace OnLoops;
using namespace mp;

using Weight = mpfloat<1024U>;

Parameters<Weight> params(2, 2, 2, 0.5);

class OnTransferMatrix : public Matrix<LinkPattern<Weight>, Weight> {
public:
  using Matrix<LinkPattern<Weight>, Weight>::Matrix;
  void sparseMultMv(Vector<Weight> &v1, Vector<Weight> &v2) const override {
    auto *o1 = static_cast<OnState<Weight> *>(&v1);
    auto *o2 = static_cast<OnState<Weight> *>(&v2);
    o1->transfer(*o2, params);
  }
};

int main(void) {
  OnState<Weight> v(hashsize), v2(hashsize);

  LinkPattern<Weight> l({BOT, BOT, EMPTY, EMPTY});
  v.add(l, 2.34);

  Weight w = 2.34;

  int nb_multiplications = 10;
  auto start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < nb_multiplications; i++)
    v.transfer(v2, params);
  auto end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> duration = end - start;

  std::cout << "Average time over " << nb_multiplications
            << " multiplications of " << v.size() << "x" << v.size()
            << " matrix with 1024bits-float coefficients: "
            << duration.count() / nb_multiplications << std::endl;

  OnTransferMatrix M(v, v2);
  M.generateSpace();
  M.generateStateTable();

  Weight *vec = new Weight[10];
  Weight *vec2 = new Weight[10];

  vec[0] = 1;

  M.MultMv(vec, vec2);

  return 0;
}

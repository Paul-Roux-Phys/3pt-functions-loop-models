#include "SparseVector.hpp"
#include <assert.h>

template <typename Value> class Matrix {
public:
  using Vector = Vector<Value>;
private:
  using Key = uint64_t;
  std::vector<Key> state_table;
  Vector &v1, &v2;
  size_t dim;

public:
  /* computes v2 = M v1 then puts the result back in v1 */
  virtual void sparseMultMv(Vector &v1, Vector &v2) = 0;

  Matrix(Vector &_v1, Vector &_v2) : v1(_v1), v2(_v2) {}

  void generateSpace() {
    assert(v1.size() != 0 &&
           "v1 is zero, cannot generate the space of states from it.");
    size_t old_dim = 0, new_dim = 1;
    while (new_dim > old_dim) {
      sparseMultMv(v1, v2);
      old_dim = new_dim;
      new_dim = v1.size();
    }
    dim = new_dim;
  }

  void generateStateTable() {
    int i = 0;
    generateSpace();
    state_table.resize(dim);
    for (auto it : v1) {
      state_table[i] = it.first;
    }
  }
};

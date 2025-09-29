#include "SparseVector.hpp"
#include <assert.h>
#include <sstream>
#include <iomanip>

template <class Key, typename Value> class Matrix {
private:
  Vector<Value> &v1, &v2;

public:
  size_t dim;
  std::vector<Key> state_table;

  /* computes v2 = M v1 then puts the result back in v1 */
  virtual void sparseMultMv(Vector<Value> &v1, Vector<Value> &v2) const = 0;

  void MultMv(Value *x, Value *y) const {
    v1.clear();
    // place x in the hash table
    for (unsigned int i = 0; i < dim; i++) {
      v1.add(state_table[i], x[i]);
    }

    // multiply by matrix
    sparseMultMv(v1, v2);

    // put the result back into y
    for (unsigned int i = 0; i < dim; i++) {
      y[i] = v1[state_table[i]];
    }

    v1.clear();
  }

  Matrix(Vector<Value> &_v1, Vector<Value> &_v2) : v1(_v1), v2(_v2) {}

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
      state_table[i++] = it.first;
    }
  }
};

template <typename Key>
std::ostream &operator<<(std::ostream &os, const std::vector<Key> &v) {
  int i = 0;
  for (auto &x : v) {
    std::cout << i++ << ": " << x << std::endl;
  }
  return os;
}

template <class Key, typename Value>
std::ostream &operator<<(std::ostream &os, const Matrix<Key, Value> &M) {
  Value **coeffs = new Value *[M.dim];
  for (size_t i = 0; i < M.dim; i++) {
    coeffs[i] = new Value[M.dim](); // calls Value() constructor
  }

  Value *vec = new Value[M.dim]();
  vec[0] = 1;
  for (size_t i = 0; i < M.dim; i++) {
    M.MultMv(vec, coeffs[i]);
    vec[i] = 0;
    if (i < M.dim - 1) {
      vec[i + 1] = 1;
    }
  }
  // --- find max width ---
  size_t max_width = 0;
  for (size_t r = 0; r < M.dim; r++) {
    for (size_t c = 0; c < M.dim; c++) {
      std::ostringstream ss;
      ss << std::fixed << std::setprecision(8) << coeffs[c][r]; // use Value's operator<<
      max_width = std::max(max_width, ss.str().size());
    }
  }

  // --- print aligned ---
  for (size_t r = 0; r < M.dim; r++) {
    for (size_t c = 0; c < M.dim; c++) {
      os << std::fixed << std::setprecision(8) << coeffs[c][r] << "  ";
    }
    os << "\n";
  }

  for (size_t i = 0; i < M.dim; i++) {
    delete[] coeffs[i];
  }
  delete[] coeffs;
  return os;
}

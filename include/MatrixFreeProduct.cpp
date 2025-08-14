#include "SparseVector.hpp"

template <typename Value>
class MatrixFreeProduct {
private:
  using Vector = Vector<Value>;
protected:
  virtual void multMv(std::vector<Value> w, std::vector<Value> v);
};

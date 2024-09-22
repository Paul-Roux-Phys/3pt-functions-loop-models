template <typename T, typename... Params>
using RMatrix = std::function<void(Vector<T>&, BasisVector<T>&, Params...)>;

template <typename T, typename... Params>
void multiply(Vector<T>& v1, RMatrix<T, Params...> r, Vector<T>& v2, Params... params) {
    for (auto it = v2.draining_begin(); it != v2.draining_end(); ++it) {
        BasisVector<T> bv(*it);
        r(v1, bv, params...);
    }
}
using LargeFloat = std::pair<Real, int>;

void standardise_largefloat(LargeFloat& f) {
    int res_exp;

    f.first = std::frexp(f.first, &res_exp);
    f.second += res_exp;
}

LargeFloat sqrt(LargeFloat p) {
    LargeFloat res;
    res = std::make_pair(std::sqrt(p.first), p.second/2);
    if (p.second%2 == 1) {
        res.first *= std::sqrt(2);
    }
    standardise_largefloat(res);
    return res;
}

std::ostream& operator<<(std::ostream& os, LargeFloat p) {
    os << std::setprecision(std::numeric_limits<Real>::digits10);
    os << p.first << " x 2^" << p.second << endl;
    return os;
}

template<typename T>
class Vector : public std::unordered_map<Key, T, KeyHash> {
private:
    LargeFloat _norm_sq; // norm_squared in the form (mantissa, exponent)

public:
    // Using base class constructors
    using std::unordered_map<Key, T, KeyHash>::unordered_map;
    using iterator = typename std::unordered_map<Key, T, KeyHash>::iterator;
    
    Vector(std::size_t size=HASH_TABLE_SIZE) : 
    std::unordered_map<Key, T, KeyHash>(size),
    _norm_sq(std::make_pair(1.0, 0)) 
    {}

    Vector(Vector<T>& v, std::size_t size) :
    std::unordered_map<Key, T, KeyHash>(size),
    _norm_sq(std::make_pair(1.0, 0))
    {
        for (auto it = v.begin(); it != v.end(); it++) {
            BasisVector<T> b(*it);
            (*this)[b.key] = b.weight;
        }
    }

    Vector(VectorPair<T, BasisVector<T>, Vector<T>>, std::size_t size=HASH_TABLE_SIZE);

    Vector& operator+=(const std::pair<Key, T>& v) {
        // Add a basis vector to the Vector
        (*this)[v.first] += v.second;
        return *this;
    }

    Vector& operator+=(const BasisVector<T>& v) {
        (*this)[v.key] += v.weight;
        return *this;
    }

    Vector& operator+=(const Vector<T>& v) {
        for (auto it = v.begin(); it != v.end(); it++) {
            BasisVector<T> b(*it);
            (*this)[b.key] += b.weight;
        }
        return *this;
    }

    Vector& operator/=(const T w) {
        for (auto it = this->begin(); it != this->end(); it++) {
            BasisVector<T> b(*it);
            (*this)[b.key] /= w;
        }
        return *this;
    }

    void add_if(const BasisVector<T>& b, bool cond) {
        if (cond) {
            *this += b;
        }
    }

    void print() const {
        for (const auto& pair : *this) {
            cout << "Key: " << pair.first 
                 << std::setw(6) << " "
                 << "Value:  " << pair.second << endl;
        }
    }

    DrainingIterator<T> draining_begin() {
        return DrainingIterator<T>(this->begin(), this);
    }

    DrainingIterator<T> draining_end() {
        return DrainingIterator<T>(this->end(), this);
    }

    T inner_product(const Vector<T>& v) const {
        T res = 0;
        for (auto it = v.begin(); it != v.end(); it++) {
            BasisVector<T> b(*it);
            auto it2 = this->find(b.key);
            if (it2 != this->end()) {
                BasisVector<T> b2(*it2);
                res += b.weight * b2.weight;
            }
        }
        return res;
    }

    T inner_product(const VectorPair<T, BasisVector<T>, Vector<T>>& v);

    LargeFloat get_norm_sq() {
        _norm_sq.first *= inner_product(*this);
        standardise_largefloat(_norm_sq);
        return _norm_sq;
    }
    
    LargeFloat get_norm() {
        LargeFloat res = sqrt(get_norm_sq());
        standardise_largefloat(res);
        return res;
    }

    void factorise_norm() {
        Real tmp = 0.0;
        int tmp_exp;
        for (auto it = this->begin(); it != this->end(); it++) {
            BasisVector<T> b(*it);
            tmp += b.weight*b.weight;
        }
        *this /= std::sqrt(tmp);
        tmp = std::frexp(tmp, &tmp_exp);
        _norm_sq.first *= tmp;
        _norm_sq.second += tmp_exp;
    }
};
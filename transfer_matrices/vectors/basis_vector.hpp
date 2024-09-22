template <typename T>
class BasisVector {
public:
    Key key;
    T weight;

    BasisVector(Key k, T v) : key(k), weight(v) {}
    BasisVector(const BasisVector& other) : key(other.key), weight(other.weight) {}
    BasisVector(std::pair<Key, T> v) : key(v.first), weight(v.second) {}

    void print() const {
        cout << "Key: " << key 
             << std::setw(6) << " "
             << "Value:  " << weight << endl;
    }
};

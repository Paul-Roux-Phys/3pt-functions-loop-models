// Custom hash function using DJB2 algorithm
struct KeyHash {
    // Overload the function call operator to define the custom hash function
    std::size_t operator()(const Key& key) const {
        std::size_t hash = 5381;
        for (char c : key) {
            hash = ((hash << 5) + hash) + static_cast<unsigned char>(c); // hash * 33 + c
        }
        return hash;
    }
};

// Print keys
std::ostream& operator<<(std::ostream& os, const Key& k) {
    os << "(";
    for (const auto& e : k) {
        os << std::setw(3) << static_cast<int>(e);
    }
    os << "  )";
    return os;
}
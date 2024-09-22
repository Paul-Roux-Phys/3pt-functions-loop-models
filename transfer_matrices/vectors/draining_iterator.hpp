template <typename T>
class DrainingIterator {
public:
    using VectorType = Vector<T>;
    using IteratorType = typename Vector<T>::iterator;
    using ValueType = std::pair<const Key, T>;

    // Constructor
    DrainingIterator(IteratorType it, VectorType* vec) : current(it), vec(vec) {}

    // Dereference operator
    ValueType& operator*() {
        return *current;
    }

    // Arrow operator
    ValueType* operator->() {
        return &(*current); // Return pointer to the underlying pair
    }

    // Pre-increment
    DrainingIterator& operator++() {
        if (current != vec->end()) {
            auto next = std::next(current);
            vec->erase(current);
            current = next;
        }
        return *this;
    }

    DrainingIterator operator++(int) {
        DrainingIterator temp = *this;
        ++(*this);
        return temp;
    }

    // Comparison operators
    bool operator==(const DrainingIterator& other) const {
        return current == other.current;
    }

    bool operator!=(const DrainingIterator& other) const {
        return !(*this == other);
    }

private:
    IteratorType current;
    VectorType* vec;
};
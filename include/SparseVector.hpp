#pragma once

#include "unordered_dense.hpp"
#include <iostream>

struct custom_hash_uint64_t {
  using is_avalanching = void; // keep the optimization hint for ankerl map
  [[nodiscard]] auto operator()(uint64_t const &f) const noexcept -> uint64_t {
    uint64_t x = f;
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    return x;
  }
};
template <typename Value> class Vector {
public:
  using Key = uint64_t;

private:
  using HashMap = ankerl::unordered_dense::map<Key, Value, custom_hash_uint64_t,
                                               std::equal_to<>>;
  HashMap hashmap;

public:
  Vector(size_t _size) : hashmap(_size) {}
  ~Vector() { hashmap.clear(); }

  // Copy constructor
  Vector(const Vector& other) : hashmap(other.hashmap) {}

  // Copy assignment operator
  Vector& operator=(const Vector& other) {
    if (this != &other) {
      hashmap = other.hashmap;
    }
    return *this;
  }

  /* optimized addition:
     lookup key. if key not found, insert new (key, value) pair.
     If key present, directly add value in-place
   */
  void add(const Key key, const Value v = 1) {
    auto [it, inserted] = hashmap.try_emplace(key, v);
    if (!inserted) {
      // key already exists, add to current vue
      it->second += v;
    }
  }

  Value operator[](const Key k) const {
    auto it = hashmap.find(k);
    if (it != hashmap.end()) {
      return it->second;
    }
    return 0;
  }

  void swap(Vector &other) noexcept { hashmap.swap(other.hashmap); }

  friend void swap(Vector &a, Vector &b) noexcept { a.swap(b); }

  HashMap::iterator begin() { return hashmap.begin(); }
  HashMap::iterator end() { return hashmap.end(); }
  HashMap::const_iterator begin() const { return hashmap.begin(); }
  HashMap::const_iterator end() const { return hashmap.end(); }
  HashMap::iterator erase(HashMap::iterator it) { return hashmap.erase(it); }

  void clear() { hashmap.clear(); }
  bool empty() { return hashmap.empty(); }

  Value inner_product(Vector &other) {
    Value res = 0;
    if (this == &other) {
      // special case: inner product with itself
      for (const auto &[key, value] : hashmap) {
        res += value * value;
      }
      return res;
    }
    // General case
    Value b;
    for (const auto &[key, value] : hashmap) {
      b = other[key];
      if (b != 0) {
        res += value * b;
      }
    }
    return res;
  }

  auto bucket_count() { return hashmap.bucket_count(); }
  auto size() { return hashmap.size(); }

  // Fast draining: clears after iteration, frees memory along the way
  // auto drain() {
  //   struct DrainingIterator {
  //     using UnderlyingIter = typename HashMap::iterator;
  //     UnderlyingIter it;

  //     // Constructor
  //     DrainingIterator(UnderlyingIter iter) : it(iter) {}

  //     // Dereference returns reference to key-value pair
  //     auto &operator*() const { return *it; }
  //     auto *operator->() const { return &(*it); }

  //     // Pre-increment: call clear on value, then advance iterator
  //     DrainingIterator &operator++() {
  //       ++it;
  //       return *this;
  //     }

  //     // Equality/Inequality
  //     bool operator==(const DrainingIterator &other) const {
  //       return it == other.it;
  //     }
  //     bool operator!=(const DrainingIterator &other) const {
  //       return it != other.it;
  //     }
  //   };

  //   struct DrainingRange {
  //     HashMap &map;

  //     DrainingIterator begin() { return DrainingIterator{map.begin()}; }
  //     DrainingIterator end() { return DrainingIterator{map.end()}; }

  //     ~DrainingRange() {
  //       map.clear(); // clear whole map at the end as backup
  //     }
  //   };

  //   return DrainingRange{hashmap};
  // }

  friend std::ostream &operator<<(std::ostream &os, const Vector &v) {
    for (auto &it : v) {
      os << it.first << "  =>  " << it.second << std::endl;
    }
    return os;
  }

  void operator*=(Value v) {
    for (auto &it : *this) {
      it.second *= v;
    }
  }
};

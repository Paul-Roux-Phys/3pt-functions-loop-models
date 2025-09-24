#pragma once

#include "../../include/SparseVector.hpp"
#include <cstdint>
#include <initializer_list>
#include <iostream>
#include <math.h>

extern unsigned int L;

namespace OnLoops {

template <typename Weight> struct Parameters {
  int min_defects;
  int k1, k2, k3;
  Weight lambda;
  Weight n_loop;
  Weight n_ncloop;
  Weight w_empty;
  Weight w_turn;
  Weight w_straight;
  Weight w_full;

  Parameters(int _k1 = 0, int _k2 = 0, int _k3 = 0, Weight lambda = 0.5)
      : min_defects(_k3), k1(_k1), k2(_k2), k3(_k3) {
    n_loop = -2. * cos(4. * lambda);   // weight of loops
    n_ncloop = -2. * cos(4. * lambda); // weight of non-contractible loops
    w_empty = 1. + sin(lambda) + sin(3. * lambda) - sin(5. * lambda);
    w_turn = 2. * sin(2. * lambda) * sin((6. * lambda + M_PI) / 4.);
    w_straight = 1. + sin(3. * lambda);
    w_full = sin(lambda) + cos(2. * lambda);
  }
};

template <typename Weight> class OnState;

#define BITS_PER_SITE 3
uint64_t constexpr MASK = (1ULL << BITS_PER_SITE) - 1;
enum { EMPTY, OPEN, CLOSE, REF, BOT, MID, ERROR };

template <typename Weight> class LinkPattern {
private:

public:
  uint64_t bitfield;

  LinkPattern(uint64_t bits = 0) : bitfield(bits) {}

  LinkPattern(std::initializer_list<int> a) : bitfield(0) {
    int i = 0;
    for (int el : a)
      set_site(i++, el);
  }

  void clear_site(int i);
  void set_site(int i, int value);
  size_t operator[](int i) const;
  operator uint64_t() const;

  void shift_right();
  bool is_occupied(int i);
  bool is_defect(int i);
  bool is_ref(int i);
  bool is_bottom(int i);
  bool is_middle(int i);
  int get_nb_defects();
  int arch_end(int i);
  void contract_arch_defect(int i);
  void contract_arches(int i);
  bool can_contract_defects(int i, Parameters<Weight> &p);
  bool can_contract_sites(int i, Parameters<Weight> &p);
  Weight contract(int i, Weight w, Parameters<Weight> &p);
  void move_strand(int i);
  void put_arch(int i);
  void r_matrix_empty_empty(OnState<Weight> &v, int pos, Weight w, Parameters<Weight> &p);
  void r_matrix_empty_occ(OnState<Weight> &v, int pos, Weight w, Parameters<Weight> &p);
  void r_matrix_occ_occ(OnState<Weight> &v, int pos, Weight w, Parameters<Weight> &p);
  void r_matrix(OnState<Weight> &v, int pos, Weight w, Parameters<Weight> &p);
  void insert_aux_space(OnState<Weight> &v, Weight w);
  void contract_aux_space(OnState<Weight> &v, Weight w, Parameters<Weight> &p);
  void middle_op_ith_site(int max_down_def, int max_up_def, int pos);
  void insert_mid_op(OnState<Weight> &v, Weight w, Parameters<Weight> &p);
  void uncolor_defects(OnState<Weight> &v, Weight w);
};

template <typename Weight> class OnState : public Vector<Weight> {
public:
  using Vector<Weight>::Vector;
  // long double get_norm_sq();
  // long double get_norm();
  // void factorize_norm();
  void r_matrix(OnState &v, int i, Parameters<Weight> &p);
  void insert_aux_space(OnState &v);
  void contract_aux_space(OnState &v, Parameters<Weight> &p);
  void insert_mid_op(OnState &v);
  void transfer(OnState &v, Parameters<Weight> &p);
  void uncolor_defects(OnState &v);
};

template <typename Weight> void LinkPattern<Weight>::clear_site(int i) {
  bitfield &= ~(MASK << BITS_PER_SITE * i);
}

template <typename Weight>
void LinkPattern<Weight>::set_site(int i, int value) {
  clear_site(i);
  bitfield |= (static_cast<uint64_t>(value) << BITS_PER_SITE * i);
}

template <typename Weight> size_t LinkPattern<Weight>::operator[](int i) const {
  return (bitfield >> BITS_PER_SITE * i) & MASK;
}

template <typename Weight> LinkPattern<Weight>::operator uint64_t() const {
  return bitfield;
}

template <typename Weight> bool LinkPattern<Weight>::is_occupied(int i) {
  return (*this)[i] > EMPTY;
}
template <typename Weight> bool LinkPattern<Weight>::is_defect(int i) {
  return REF <= (*this)[i] && (*this)[i] <= MID;
}
template <typename Weight> bool LinkPattern<Weight>::is_ref(int i) {
  return (*this)[i] == REF;
}
template <typename Weight> bool LinkPattern<Weight>::is_bottom(int i) {
  return (*this)[i] == REF || (*this)[i] == BOT;
}
template <typename Weight> bool LinkPattern<Weight>::is_middle(int i) {
  return (*this)[i] == MID;
}
template <typename Weight> int LinkPattern<Weight>::get_nb_defects() {
  int res = 0;
  for (size_t i = 0; i < L; i++)
    if (is_defect(i))
      res += 1;
  return res;
}

template <typename Weight> void LinkPattern<Weight>::shift_right() {
  int last = (*this)[L + 1];
  for (int i = L; i >= 0; i--) {
    set_site(i + 1, (*this)[i]);
  }
  set_site(0, last);
}

template <typename Weight> int LinkPattern<Weight>::arch_end(int i) {
  unsigned int site = (*this)[i];
  if (site == EMPTY || is_defect(i)) {
    return i;
  }
  int direction = (site == OPEN) ? 1 : -1;
  int inner_arches = 1;
  int j = i;
  int limit = 0;
  do {
    if (limit++ > 1000) {
      std::cout << "Couldn't find end of arch at position " << i << " for key  "
                << *this << "  of size " << L + 2 << std::endl
                << "Check that the key is a valid link pattern" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    j = (j + direction + L + 2) % (L + 2);
    if ((*this)[j] == site) { // we are meeting a new arch of the same type
      inner_arches += 1;
    } else if ((*this)[j] ==
               site + direction) { // an arch of the same type is ending
      inner_arches -= 1;
    }
  } while (inner_arches > 0);
  return j;
}

template <typename Weight>
void LinkPattern<Weight>::contract_arch_defect(int i) {
  int ip1 = (i + 1) % (L + 2);
  int arch = is_defect(i) ? ip1 : i;
  int def = is_defect(i) ? i : ip1;
  set_site(arch_end(arch), (*this)[def]);
}

template <typename Weight> void LinkPattern<Weight>::contract_arches(int i) {
  int ip1 = (i + 1) % (L + 2);
  int site1 = (*this)[i], site2 = (*this)[ip1];
  if ((site1 != OPEN && site1 != CLOSE) || (site2 != OPEN && site2 != CLOSE)) {
    return;
  }
  // Contract arches at positions i and i+1
  int end1 = arch_end(i), end2 = arch_end(ip1);
  if (site1 == OPEN && site2 == OPEN) {
    set_site(end2, OPEN);
  } else if (site1 == CLOSE && site2 == CLOSE) {
    set_site(end1, CLOSE);
  }
}

template <typename Weight>
bool LinkPattern<Weight>::can_contract_defects(int i, Parameters<Weight> &p) {
  int ip1 = (i + 1) % (L + 2);
  return (
      !is_ref(i) && !is_ref(ip1) &&
      ((is_bottom(i) && is_middle(ip1)) || (is_bottom(ip1) && is_middle(i))) &&
      get_nb_defects() - 2 >= p.min_defects);
}

template <typename Weight>
bool LinkPattern<Weight>::can_contract_sites(int i, Parameters<Weight> &p) {
  int ip1 = (i + 1) % (L + 2);
  return (can_contract_defects(i, p) || !is_defect(i) || !is_defect(ip1));
}

template <typename Weight>
Weight LinkPattern<Weight>::contract(int i, Weight w, Parameters<Weight> &p) {
  int ip1 = (i + 1) % (L + 2);
  if ((*this)[i] == OPEN &&
      (*this)[ip1] == CLOSE) { // closing a contractible loop
    w *= p.n_loop;
  } else if ((*this)[i] == CLOSE && (*this)[ip1] == OPEN &&
             arch_end(i) == ip1) {
    // closing a loop winding around the cylinder
    w *= p.n_ncloop;
  }
  if (is_defect(i) != is_defect((i + 1) % (L + 2))) {
    // XOR: exactly one of the two sites is a defect
    contract_arch_defect(i);
  } else if (can_contract_defects(i, p)) {
    // nothing to do here
  } else {
    contract_arches(i);
  }
  return w;
}

template <typename Weight> void LinkPattern<Weight>::move_strand(int i) {
  int empty = ((*this)[i] == EMPTY) ? i : (i + 1) % (L + 2);
  int def = ((*this)[i] == EMPTY) ? (i + 1) % (L + 2) : i;
  set_site(empty, (*this)[def]);
  set_site(def, EMPTY);
}

template <typename Weight> void LinkPattern<Weight>::put_arch(int i) {
  // insert an arch between sites i and i+1;
  set_site(i, OPEN);
  set_site((i + 1) % (L + 2), CLOSE);
}

template <typename Weight>
void LinkPattern<Weight>::r_matrix_empty_empty(OnState<Weight> &v, int i, Weight w,
                                               Parameters<Weight> &p) {
  LinkPattern copy = *this;

  v.add(*this, w * p.w_empty);

  copy.put_arch(i);
  v.add(copy, w * p.w_turn);
}

template <typename Weight>
void LinkPattern<Weight>::r_matrix_empty_occ(OnState<Weight> &v, int i, Weight w,
                                             Parameters<Weight> &p) {
  LinkPattern copy = *this;

  v.add(*this, w * p.w_turn);

  copy.move_strand(i);
  v.add(copy, w * p.w_straight);
}

template <typename Weight>
void LinkPattern<Weight>::r_matrix_occ_occ(OnState<Weight> &v, int i, Weight w,
                                           Parameters<Weight> &p) {
  LinkPattern copy = *this;

  v.add(*this, w * p.w_full);

  if (can_contract_sites(i, p)) {
    w = copy.contract(i, w, p);
    LinkPattern copy2 = copy;
    copy.set_site(i, EMPTY);
    copy.set_site((i + 1) % (L + 2), EMPTY);
    v.add(copy, w * p.w_turn);

    copy2.put_arch(i);
    v.add(copy2, w * p.w_full);
  }
}

template <typename Weight>
void LinkPattern<Weight>::r_matrix(OnState<Weight> &v, int i, Weight w, Parameters<Weight> &p) {
  int ip1 = (i + 1) % (L + 2);
  if ((*this)[i] == EMPTY && (*this)[ip1] == EMPTY) {
    r_matrix_empty_empty(v, i, w, p);
  } else if ((*this)[i] == EMPTY || (*this)[ip1] == EMPTY) {
    r_matrix_empty_occ(v, i, w, p);
  } else {
    r_matrix_occ_occ(v, i, w, p);
  }
}

template <typename Weight>
void LinkPattern<Weight>::insert_aux_space(OnState<Weight> &v, Weight w) {
  shift_right();

  LinkPattern copy = *this;

  set_site(L + 1, EMPTY);
  set_site(0, EMPTY);
  v.add(*this, w);

  copy.set_site(L + 1, OPEN);
  copy.set_site(0, CLOSE);
  v.add(copy, w);
}

template <typename Weight>
void LinkPattern<Weight>::contract_aux_space(OnState<Weight> &v, Weight w,
                                             Parameters<Weight> &p) {
  int i1 = L, i2 = L + 1;
  if ((*this)[i1] == EMPTY && (*this)[i2] == EMPTY) {
    v.add(*this, w);
  } else if (can_contract_sites(i1, p) && is_occupied(i1) && is_occupied(i2)) {
    // both sites are occupied, at least one of them isn't a defect
    w = contract(i1, w, p);
    set_site(i1, EMPTY);
    set_site(i2, EMPTY);
    v.add(*this, w);
  }
}

// pos must be < max(max_down_def, max_up_def)
// e.g. if there are 1 up and 1 down defects then pos = 0
// this function handles the insertion of the middle operator at position pos
template <typename Weight>
void LinkPattern<Weight>::middle_op_ith_site(int nb_down_def, int nb_up_def,
                                             int pos) {
  if (is_bottom(pos)) {
    if (pos >= nb_down_def) {
      set_site(0, ERROR);
      return;
    }
    if (is_ref(pos)) {
      set_site(0, ERROR);
      return;
    }
    if (nb_up_def <= pos && pos < nb_down_def)
      set_site(pos, EMPTY);
    else
      set_site(pos, MID);
  } else if ((*this)[pos] == EMPTY) {
    if (nb_down_def <= pos && pos < nb_up_def)
      set_site(pos, MID);
    else {
      set_site(0, ERROR); // error
      return;
    }
  } else if (!is_defect(pos)) {
    int end_pos = arch_end(pos);
    if (pos >= nb_down_def) {
      set_site(0, ERROR); // error
      return;
    }
    if (end_pos < nb_down_def || end_pos < nb_up_def) {
      set_site(0, ERROR); // error
      return;
    }
    if (pos >= nb_up_def)
      set_site(pos, EMPTY);
    else
      set_site(pos, MID);
    set_site(end_pos, MID);
  }
}

template <typename Weight>
void LinkPattern<Weight>::insert_mid_op(OnState<Weight> &v, Weight w, Parameters<Weight> &p) {
  if (p.k2 % 2 == 0) {
    for (int i = 0; i < p.k2 / 2; i++)
      middle_op_ith_site(p.k2 / 2, p.k2 / 2, i);
    if (!((*this)[0] == ERROR))
      v.add(*this, w);
  } else {
    LinkPattern copy = *this;
    for (int i = 0; i < p.k2 / 2 + 1; i++) {
      // symmetrise up and down
      this->middle_op_ith_site(p.k2 / 2 + 1, p.k2 / 2, i);
      copy.middle_op_ith_site(p.k2 / 2, p.k2 / 2 + 1, i);
    }
    if (!((*this)[0] == ERROR))
      v.add(*this, w);
  }
}

template <typename Weight>
void LinkPattern<Weight>::uncolor_defects(OnState<Weight> &v, Weight w) {
  for (unsigned int i = 0; i < L; i++)
    if ((*this)[i] == MID)
      set_site(i, BOT);
  v.add(*this, w);
}

template <typename Weight>
std::ostream &operator<<(std::ostream &os, const LinkPattern<Weight> &l) {
  for (unsigned int i = 0; i < L + 2; i++) {
    os << l[i];
  }
  return os;
}

template <typename Weight>
std::ostream &operator<<(std::ostream &os, const OnState<Weight> &v) {
  for (auto &it : v) {
    os << LinkPattern<Weight>(it.first) << "  =>  " << it.second << std::endl;
  }
  return os;
};

// template <typename Weight> long double OnState<Weight>::get_norm_sq() {
//   long double tmp = 0.0;
//   for (auto &it : *this) {
//     tmp += it.second * it.second;
//   }
//   return norm_sq * tmp;
// }

// template <typename Weight> long double OnState<Weight>::get_norm() {
//   return sqrt(get_norm_sq());
// }

// template <typename Weight> void OnState<Weight>::factorize_norm() {
//   long double tmp = 0.0;
//   for (auto &it : *this) {
//     tmp += it.second * it.second;
//   }
//   long double invsqrt_tmp = 1/(sqrt(tmp));
//   for (auto &it : *this) {
//     it.second *= invsqrt_tmp;
//   }
//   norm_sq *= tmp;
// }

template <typename Weight>
void OnState<Weight>::r_matrix(OnState<Weight> &v, int i, Parameters<Weight> &p) {
  for (auto it = this->begin(); it != this->end(); ) {
    LinkPattern<Weight> lp(it->first);
    lp.r_matrix(v, i, it->second, p);
    it = this->erase(it);
  }
  this->swap(v);
}
template <typename Weight>
void OnState<Weight>::insert_aux_space(OnState<Weight> &v) {
  for (auto it = this->begin(); it != this->end(); ) {
    LinkPattern<Weight> lp(it->first);
    lp.insert_aux_space(v, it->second);
    it = this->erase(it);
  }
  this->swap(v);
}
template <typename Weight>
void OnState<Weight>::contract_aux_space(OnState<Weight> &v, Parameters<Weight> &p) {
  for (auto it = this->begin(); it != this->end(); ) {
    LinkPattern<Weight> lp(it->first);
    lp.contract_aux_space(v, it->second, p);
    it = this->erase(it);
  }
  this->swap(v);
}
template <typename Weight>
void OnState<Weight>::insert_mid_op(OnState<Weight> &v) {
  for (auto it = this->begin(); it != this->end(); ) {
    LinkPattern<Weight> lp(it->first);
    lp.insert_mid_op(v, it->second);
    it = this->erase(it);
  }
  this->swap(v);
}

template <typename Weight>
void OnState<Weight>::transfer(OnState<Weight> &v, Parameters<Weight> &p) {
  insert_aux_space(v);
  for (unsigned int i = 0; i < L; i++) {
    r_matrix(v, i, p);
  }
  contract_aux_space(v, p);
}
template <typename Weight> void uncolor_defects(OnState<Weight> &v);

template <typename Weight>
void transfer(OnState<Weight> *v1, OnState<Weight> *v2, Parameters<Weight> &p) {
  v1->transfer(*v2, p);
}

} // namespace OnLoops

#include "../include/SparseVector.hpp"
#include <initializer_list>
#include <iostream>
#include <print>
#include <complex>

unsigned int L;
#define SIZE (L + 2)
int min_defects;
int k2;

using Weight = std::complex<double>;

Weight lambda;
Weight n_loop, n_ncloop, w_empty, w_turn, w_straight, w_full;

void set_weights(Weight lambda) {
  n_loop = -2. * cos(4. * lambda);   // weight of loops
  n_ncloop = -2. * cos(4. * lambda); // weight of non-contractible loops
  w_empty = 1. + sin(lambda) + sin(3. * lambda) - sin(5. * lambda);
  w_turn = 2. * sin(2. * lambda) * sin((6. * lambda + M_PI) / 4.);
  w_straight = 1. + sin(3. * lambda);
  w_full = sin(lambda) + cos(2. * lambda);
}

class OnState;

#define BITS_PER_SITE 3
uint64_t constexpr MASK = (1ULL << BITS_PER_SITE) - 1;
enum { EMPTY, OPEN, CLOSE, REF, BOT, MID, ERROR };

class LinkPattern {
public:
  uint64_t bitfield;

  LinkPattern(uint64_t bits) : bitfield(bits) {}
  LinkPattern(std::initializer_list<int> a) : bitfield(0) {
    int i = 0;
    for (int el : a)
      set_site(i++, el);
  }

  void clear_site(int i) {
    bitfield &= ~(MASK << BITS_PER_SITE * i);
  }

  void set_site(int i, int value) {
    clear_site(i);
    bitfield |= (static_cast<uint64_t>(value) << BITS_PER_SITE * i);
  }

  size_t operator[](int i) const {
    return (bitfield >> BITS_PER_SITE * i) & MASK;
  }

  friend std::ostream &operator<<(std::ostream &os, const LinkPattern &l) {
    for (unsigned int i = 0; i < SIZE; i++) {
      os << l[i];
    }
    return os;
  }

  operator uint64_t() const { return bitfield; }

  void shift_right();
  bool is_occupied(int i) { return (*this)[i] > EMPTY; }
  bool is_defect(int i) { return REF <= (*this)[i] && (*this)[i] <= MID; }
  bool is_ref(int i) { return (*this)[i] == REF; }
  bool is_bottom(int i) { return (*this)[i] == REF || (*this)[i] == BOT; }
  bool is_middle(int i) { return (*this)[i] == MID; }
  int get_nb_defects() {
    int res = 0;
    for (size_t i = 0; i < SIZE; i++)
      if (is_defect(i))
        res += 1;
    return res;
  }
  int arch_end(int i);
  void contract_arch_defect(int i);
  void contract_arches(int i);
  bool can_contract_defects(int i);
  bool can_contract_sites(int i);
  Weight contract(int i, Weight w);
  void move_strand(int i);
  void put_arch(int i);
  void r_matrix_empty_empty(OnState &v, int pos, Weight w);
  void r_matrix_empty_occ(OnState &v, int pos, Weight w);
  void r_matrix_occ_occ(OnState &v, int pos, Weight w);
  void r_matrix(OnState &v, int pos, Weight w);
  void insert_aux_space(OnState &v, Weight w);
  void contract_aux_space(OnState &v, Weight w);
  void middle_op_ith_site(int max_down_def, int max_up_def, int pos);
  void insert_mid_op(OnState &v, Weight w);
  void uncolor_defects(OnState &v, Weight w);
};

class OnState : public Vector<Weight> {
public:
  using Vector<Weight>::Vector;
  friend std::ostream &operator<<(std::ostream &os, const OnState &v) {
    for (auto& it : v) {
      os << LinkPattern(it.first) << "  =>  " << it.second << std::endl;
    }
    return os;
  }
};

int main() {
  L = 2;
  LinkPattern l({EMPTY, EMPTY, OPEN, CLOSE});
  std::print("{:b}\n", l.bitfield);
  std::print("{}, {}, {}, {}\n", l[0], l[1], l[2], l[3]);
  std::cout << l << std::endl;

  OnState s(10);
  s.add(l, 10.1);
  std::cout << s << std::endl;
  return 0;
}

void LinkPattern::shift_right() {
  int last = (*this)[L + 1];
  for (int i = L; i >= 0; i--) {
    set_site(i + 1, (*this)[i]);
  }
  set_site(0, last);
}

int LinkPattern::arch_end(int i) {
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
                << *this << "  of size " << SIZE << std::endl
                << "Check that the key is a valid link pattern" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    j = (j + direction + SIZE) % (SIZE);
    if ((*this)[j] == site) { // we are meeting a new arch of the same type
      inner_arches += 1;
    } else if ((*this)[j] ==
               site + direction) { // an arch of the same type is ending
      inner_arches -= 1;
    }
  } while (inner_arches > 0);
  return j;
}

void LinkPattern::contract_arch_defect(int i) {
  int ip1 = (i + 1) % SIZE;
  int arch = is_defect(i) ? ip1 : i;
  int def = is_defect(i) ? i : ip1;
  set_site(arch_end(arch), (*this)[def]);
}

void LinkPattern::contract_arches(int i) {
  int ip1 = (i + 1) % SIZE;
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

bool LinkPattern::can_contract_defects(int i) {
  int ip1 = (i + 1) % SIZE;
  return (
      !is_ref(i) && !is_ref(ip1) &&
      ((is_bottom(i) && is_middle(ip1)) || (is_bottom(ip1) && is_middle(i))) &&
      get_nb_defects() - 2 >= min_defects);
}

bool LinkPattern::can_contract_sites(int i) {
  int ip1 = (i + 1) % SIZE;
  return (can_contract_defects(i) || !is_defect(i) || !is_defect(ip1));
}

Weight LinkPattern::contract(int i, Weight w) {
  int ip1 = (i + 1) % SIZE;
  if ((*this)[i] == OPEN &&
      (*this)[ip1] == CLOSE) { // closing a contractible loop
    w *= n_loop;
  } else if ((*this)[i] == CLOSE && (*this)[ip1] == OPEN &&
             arch_end(i) == ip1) {
    // closing a loop winding around the cylinder
    w *= n_ncloop;
  }
  if (is_defect(i) != is_defect((i + 1) % SIZE)) {
    // XOR: exactly one of the two sites is a defect
    contract_arch_defect(i);
  } else if (can_contract_defects(i)) {
    // nothing to do here
  } else {
    contract_arches(i);
  }
  return w;
}

void LinkPattern::move_strand(int i) {
  int empty = ((*this)[i] == EMPTY) ? i : (i + 1) % SIZE;
  int def = ((*this)[i] == EMPTY) ? (i + 1) % SIZE : i;
  set_site(empty, (*this)[def]);
  set_site(def, EMPTY);
}

void LinkPattern::put_arch(int i) {
  // insert an arch between sites i and i+1;
  set_site(i, OPEN);
  set_site((i + 1) % SIZE, CLOSE);
}

void LinkPattern::r_matrix_empty_empty(OnState &v, int i, Weight w) {
  LinkPattern copy = *this;

  v.add(*this, w * w_empty);

  copy.put_arch(i);
  v.add(copy, w * w_turn);
}

void LinkPattern::r_matrix_empty_occ(OnState &v, int i, Weight w) {
  LinkPattern copy = *this;

  v.add(*this, w * w_turn);

  copy.move_strand(i);
  v.add(copy, w * w_straight);
}

void LinkPattern::r_matrix_occ_occ(OnState &v, int i, Weight w) {
  LinkPattern copy = *this;

  v.add(*this, w * w_full);

  if (can_contract_sites(i)) {
    w = copy.contract(i, w);
    LinkPattern copy2 = copy;
    copy.set_site(i, EMPTY);
    copy.set_site((i + 1) % SIZE, EMPTY);
    v.add(copy, w * w_turn);

    copy2.put_arch(i);
    v.add(copy2, w * w_full);
  }
}

void LinkPattern::r_matrix(OnState &v, int i, Weight w) {
  int ip1 = (i + 1) % SIZE;
  if ((*this)[i] == EMPTY && (*this)[ip1] == EMPTY) {
    r_matrix_empty_empty(v, i, w);
  } else if ((*this)[i] == EMPTY || (*this)[ip1] == EMPTY) {
    r_matrix_empty_occ(v, i, w);
  } else {
    r_matrix_occ_occ(v, i, w);
  }
}

void LinkPattern::insert_aux_space(OnState &v, Weight w) {
  shift_right();

  LinkPattern copy = *this;

  set_site(SIZE - 1, EMPTY);
  set_site(0, EMPTY);
  v.add(*this, w);

  copy.set_site(SIZE - 1, OPEN);
  copy.set_site(0, CLOSE);
  v.add(copy, w);
}

void LinkPattern::contract_aux_space(OnState &v, Weight w) {
  int i1 = L, i2 = L + 1;
  if ((*this)[i1] == EMPTY && (*this)[i2] == EMPTY) {
    v.add(*this, w);
  } else if (can_contract_sites(i1) && is_occupied(i1) && is_occupied(i2)) {
    // both sites are occupied, at least one of them isn't a defect
    w = contract(i1, w);
    set_site(i1, EMPTY);
    set_site(i2, EMPTY);
    v.add(*this, w);
  }
}

// pos must be < max(max_down_def, max_up_def)
// e.g. if there are 1 up and 1 down defects then pos = 0
// this function handles the insertion of the middle operator at position pos
void LinkPattern::middle_op_ith_site(int nb_down_def, int nb_up_def, int pos) {
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

void LinkPattern::insert_mid_op(OnState &v, Weight w) {
  if (k2 % 2 == 0) {
    for (int i = 0; i < k2 / 2; i++)
      middle_op_ith_site(k2 / 2, k2 / 2, i);
    if (!((*this)[0] == ERROR))
      v.add(*this, w);
  } else {
    LinkPattern copy = *this;
    for (int i = 0; i < k2 / 2 + 1; i++) {
      // symmetrise up and down
      this->middle_op_ith_site(k2 / 2 + 1, k2 / 2, i);
      copy.middle_op_ith_site(k2 / 2, k2 / 2 + 1, i);
    }
    if (!((*this)[0] == ERROR))
      v.add(*this, w);
  }
}

void LinkPattern::uncolor_defects(OnState &v, Weight w) {
  for (unsigned int i = 0; i < L; i++)
    if ((*this)[i] == MID)
      set_site(i, BOT);
  v.add(*this, w);
}

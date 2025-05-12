#ifndef TRANSFER_HPP
#define TRANSFER_HPP

#include <algorithm>
#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <iostream>

#include "TransferMatrices.hpp"

// Lattice size
size_t L;
size_t M;
#define SIZE (L + 2)

// Hash table size
size_t HASH_SIZE = (2 << L);  // 2^{L+1}

// Number of iterations of transfer matrix
int iter = 0;

// Store keys on 64 bit integers. Store each site value on 3 bits.
#define bits_per_site 3
class Key : public key_64_bit_t<bits_per_site> {
   public:
    using key_64_bit_t<bits_per_site>::key_64_bit_t;
    void print () const { key_64_bit_t<bits_per_site>::print (SIZE); }
    void shift_right () {
        int last = operator[] (L + 1);
        for (int i = L; i >= 0; i--) {
            set (i + 1, operator[] (i));
        }
        set (0, last);
    }
};

using Hash = KeyHash<Key>;
// Complex numbers with 100 digits
using Weight = boost::multiprecision::cpp_complex_100;
class LinkPattern;
class OnState;

// Weights
Weight lambda;
Weight n_loop, n_ncloop, w_empty, w_turn, w_straight, w_full;
#define PI                                                                     \
    Weight (                                                                   \
        "3."                                                                   \
        "14159265358979323846264338327950288419716939937510582097494459230781" \
        "640628"                                                               \
        "620899862803482534211706798214808651328230")
Weight omega1, omega3;

// Fields parameters
int k1, k2, k3;     // 2r1, 2r2, 2r3: numbers of legs
int min_defects;    // number of defects that must remain at all times
int rs1, rs2, rs3;  // Spin of fields
Weight w1, w2, w3;  // weights of diagonal fields

void set_weights (Weight lambda) {
    n_loop = -2 * cos (4 * lambda);    // weight of loops
    n_ncloop = n_loop;  // weight of non-contractible loops
    w_empty = 1 + sin (lambda) + sin (3 * lambda) - sin (5 * lambda);
    w_turn = 2 * sin (2 * lambda) * sin ((6 * lambda + M_PI) / 4);
    w_straight = 1 + sin (3 * lambda);
    w_full = sin (lambda) + cos (2 * lambda);
}

/**
 * `REF` is the reference leg.
 * It has to go to the top of the cylinder
 */
enum { EMPTY, OPEN, CLOSE, REF, BOT, MID, ERROR };
// Weighted link patterns
class LinkPattern : public BasisVector<Key, Weight> {
   public:
    using BasisVector<Key, Weight>::BasisVector;  // use the base constructor
    LinkPattern (Key k, Weight v) : BasisVector<Key, Weight> (k, v) {}
    LinkPattern (Key k) : BasisVector<Key, Weight> (k, Weight (1)) {}

    void print () const { BasisVector<Key, Weight>::print (); }

    inline bool is_occupied (int i) { return key[i] > EMPTY; }
    inline bool is_defect (int i) { return REF <= key[i] && key[i] <= MID; }
    inline bool is_ref (int i) { return key[i] == REF; }
    inline bool is_bottom (int i) { return key[i] == REF || key[i] == BOT; }
    inline bool is_middle (int i) { return key[i] == MID; }

    int get_nb_defects () {
        int res = 0;
        for (size_t i = 0; i < SIZE; i++)
            if (is_defect (i)) res += 1;
        return res;
    }

    int arch_end (int i);
  bool arch_crosses_border (int i);
    void contract_arch_defect (int i);
    void contract_arches (int i);
    bool can_contract_defects (int i);
    bool can_contract_sites (int i);
    void contract (int i);
    void move_strand (int i);
    void put_arch (int i);
    void r_matrix_empty_empty (OnState& v, int pos);
    void r_matrix_empty_occ (OnState& v, int pos);
    void r_matrix_occ_occ (OnState& v, int pos);
    void r_matrix (OnState& v, int pos);
    void insert_aux_space (OnState& v);
    void contract_aux_space (OnState& v);
    void middle_op_ith_site (int max_down_def, int max_up_def, int pos);
    void insert_mid_op (OnState& v);
    void uncolor_defects (OnState& v);
};

class OnState : public Vector<LinkPattern, Hash> {
   public:
    OnState () : Vector<LinkPattern, Hash> (HASH_SIZE) {};
    OnState (size_t size) : Vector<LinkPattern, Hash> (size) {}
    OnState (size_t size, LinkPattern l) : Vector<LinkPattern, Hash> (size) {
        *this += l;
    }
    using Vector<LinkPattern, Hash>::operator+=;

    // OnState& operator+=(OnState& v) {
    //     Vector<LinkPattern, Hash>::operator+=(v);
    // }

    void print () const { Vector<LinkPattern, Hash>::print (); }

#define APPLY_TO_EACH(METHOD)                                      \
    for (auto it = draining_begin (); it != draining_end (); it++) \
    LinkPattern (*it).METHOD

    void r_matrix (OnState& v, int i) {
        APPLY_TO_EACH (r_matrix) (v, i);
        swap (v);  // pull the result to this instance (shallow copy)
    }

    void insert_aux_space (OnState& v) {
        APPLY_TO_EACH (insert_aux_space) (v);
        swap (v);
    }

    void contract_aux_space (OnState& v) {
        APPLY_TO_EACH (contract_aux_space) (v);
        swap (v);
    }

    void insert_mid_op (OnState& v) {
        APPLY_TO_EACH (insert_mid_op) (v);
        swap (v);
    }

    void transfer (OnState& v) {
        insert_aux_space (v);
        for (int i = 0; i < L; i++) {
            r_matrix (v, i);
        }
        contract_aux_space (v);
	iter += 1;
    }

    void uncolor_defects (OnState& v) {
        APPLY_TO_EACH (uncolor_defects) (v);
        swap (v);
    }

  void matrix_mul (OnState& v);
};

#endif  // TRANSFER_HPP

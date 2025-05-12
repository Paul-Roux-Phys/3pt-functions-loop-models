#ifndef TRANSFER_HPP
#define TRANSFER_HPP

#include <algorithm>
#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <iostream>

#include "TransferMatrices/Vectors.hpp"
#include "TransferMatrices/keys.hpp"

// Lattice size
size_t L;
#define SIZE (L + 2)

// Hash table size
size_t HASH_SIZE = (2 << L);  // 2^{L+1}

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
class UnState;

// Weights
Weight lambda;
Weight n_loop, n_ncloop;
#define PI                                                                     \
    Weight (                                                                   \
        "3."                                                                   \
        "14159265358979323846264338327950288419716939937510582097494459230781" \
        "640628"                                                               \
        "620899862803482534211706798214808651328230")

// Fields parameters
int k1, k2, k3;     // 2r1, 2r2, 2r3: numbers of legs
int min_defects;    // number of defects that must remain at all times
int rs1, rs2, rs3;  // Spin of fields

void set_weights (Weight lambda) {
    n_loop = -2 * cos (4 * lambda);    // weight of loops
    n_ncloop = n_loop;  // weight of non-contractible loops
}

/**
 * `REF` is the reference leg.
 * It has to go to the top of the cylinder
 */
enum { OPEN, CLOSE, REF, BOT, MID, ERROR };
// Weighted link patterns
class LinkPattern : public BasisVector<Key, Weight> {
   public:
    using BasisVector<Key, Weight>::BasisVector;  // use the base constructor
    LinkPattern (Key k, Weight v) : BasisVector<Key, Weight> (k, v) {}
    LinkPattern (Key k) : BasisVector<Key, Weight> (k, Weight (1)) {}

    void print () const { BasisVector<Key, Weight>::print (); }

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
    void contract_arch_defect (int i);
    void contract_arches (int i);
    bool can_contract_defects (int i);
    bool can_contract_sites (int i);
    void contract (int i);
    void move_strand (int i);
    void put_arch (int i);
    void r_matrix (UnState& v, int pos);
    void insert_aux_space (UnState& v);
    void contract_aux_space (UnState& v);
    void middle_op_ith_site (int max_down_def, int max_up_def, int pos);
    void insert_mid_op (UnState& v);
    void uncolor_defects (UnState& v);
};

class UnState : public Vector<LinkPattern, Hash> {
   public:
    UnState () : Vector<LinkPattern, Hash> (HASH_SIZE) {};
    UnState (size_t size) : Vector<LinkPattern, Hash> (size) {}
    UnState (size_t size, LinkPattern l) : Vector<LinkPattern, Hash> (size) {
        *this += l;
    }
    using Vector<LinkPattern, Hash>::operator+=;

    // UnState& operator+=(UnState& v) {
    //     Vector<LinkPattern, Hash>::operator+=(v);
    // }

    void print () const { Vector<LinkPattern, Hash>::print (); }

#define APPLY_TO_EACH(METHOD)                                      \
    for (auto it = draining_begin (); it != draining_end (); it++) \
    LinkPattern (*it).METHOD

    void r_matrix (UnState& v, int i) {
        APPLY_TO_EACH (r_matrix) (v, i);
        swap (v);  // pull the result to this instance (shallow copy)
    }

    void insert_aux_space (UnState& v) {
        APPLY_TO_EACH (insert_aux_space) (v);
        swap (v);
    }

    void contract_aux_space (UnState& v) {
        APPLY_TO_EACH (contract_aux_space) (v);
        swap (v);
    }

    void insert_mid_op (UnState& v) {
        APPLY_TO_EACH (insert_mid_op) (v);
        swap (v);
    }

    void transfer (UnState& v) {
        insert_aux_space (v);
        for (int i = 0; i < L; i++) {
            r_matrix (v, i);
        }
        contract_aux_space (v);
    }

    void contract (UnState& v, int i);

    void uncolor_defects (UnState& v) {
        APPLY_TO_EACH (uncolor_defects) (v);
        swap (v);
    }
};

#endif  // TRANSFER_HPP

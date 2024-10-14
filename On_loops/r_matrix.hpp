/* This file implements the transfer matrix that generates configurations
of the integrable O(n) model defined in
[[https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.69.710]
[Warnaar Nienhuis Seaton new construction of solvable 
lattice models including an Ising model in a field]]                       */

#pragma once

#include "TransferMatrices.hpp"
#include "link_patterns.hpp"
#include <stdexcept>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
using namespace transfer_matrices;
using boost::multiprecision::cpp_dec_float_100;
using boost::multiprecision::cpp_complex_100;

#ifndef LATTICE_SIZE
#define LATTICE_SIZE 4
#endif

constexpr std::size_t lattice_size = LATTICE_SIZE;
constexpr std::size_t size         = lattice_size + 2;
using K = OnKey<size>;
using V = double;
using BV = BasisVector<K, V>;
using Vec = Vector<BV, key_64_bit_hash_t<size>>;
using R = RMatrix<Vec, int>;

Vec v[2];
Vec initial(2);
VectorPair vp(&v[0], &v[1]);

#pragma region Weights
V lambda;
V n_loop, n_ncloop, w_empty, w_turn, w_straight, w_full;
void set_weights(V lambda) {
    n_loop = -2 * cos(4 * lambda);   // weight of loops
    n_ncloop = -2 * cos(4 * lambda); // weight of non-contractible loops
    w_empty = 1 + sin(lambda) + sin(3 * lambda) - sin(5 * lambda);
    w_turn = 2 * sin(2 * lambda) * sin((6 * lambda + M_PI) / 4);
    w_straight = 1 + sin(3 * lambda);
    w_full = sin(lambda) + cos(2 * lambda);
}

void print_weights() {
    cout << "lambda = " << lambda << endl;
    cout << "n_loop = " << n_loop << endl;
    // cout << "w_empty = " << w_empty << endl;
    // cout << "w_turn = " << w_turn << endl;
    // cout << "w_straight = " << w_straight << endl;
    // cout << "w_full = " << w_full << endl;
}
#pragma endregion

#pragma region R-matrices declarations
// Declare the R-matrix, auxiliary spaces;
void r_matrix(Vec* v, BV b, int i, int min_defects);
void insert_aux_space(Vec* v, BV b);
void contract_aux_space(Vec* v, BV b, int min_defects);
RMatrix<Vec, int, int> r_i = r_matrix;
RMatrix<Vec> ins_aux = insert_aux_space;
RMatrix<Vec, int> contr_aux = contract_aux_space;
#pragma endregion

#pragma region Implement the R matrix
void translate_lattice(K& k) {
    int tmp = k[lattice_size-1];
    for (int i = lattice_size-2; i >= 0; i--)
    {
        k.set(i+1, k[i]);
    }
    k.set(0, tmp);
}

void half_lattice_translation(K& k) {
    for (int i = 0; i < lattice_size/2; i++)
    {
        translate_lattice(k);
    }
}

void contract(BV& b, int i, int min_defects) {
    int ip1 = (i+1)%size;
    if (b.key[i] == OPENING && b.key[ip1] == CLOSING) 
    { // closing a contractible loop
        b.value *= n_loop;
    }
    else if (b.key[i] == CLOSING && b.key[ip1] == OPENING
             && b.key.arch_end(i) == ip1)
    { // closing a loop winding around the cylinder
        b.value *= n_ncloop;
    }
    b.key.contract(i, min_defects);
}

void r_matrix_empty_empty(Vec* v, BV b, int i) {
    BV b2 = b;

    b.value *= w_empty;
    *v += b;

    b2.key.put_arch(i);
    b2.value *= w_turn;
    *v += b2;
}

void r_matrix_empty_occ(Vec* v, BV b, int i) {
    BV b2 = b;

    b.value *= w_turn;
    *v += b;

    b2.key.move_strand(i);
    b2.value *= w_straight;
    *v += b2;
}

void r_matrix_occ_occ(Vec* v, BV b, int i, int min_defects) {
    BV b2 = b;
    
    b.value *= w_full;
    *v += b;

    if (b.key.can_contract_sites(i, min_defects))
    {
        contract(b2, i, min_defects);
        BV b3 = b2;
        b2.key.set(i, EMPTY);
        b2.key.set((i+1)%size, EMPTY);
        b2.value *= w_turn;
        *v += b2;

        b3.key.put_arch(i);
        b3.value *= w_full;
        *v += b3;
    }
}

void r_matrix(Vec* v, BV b, int i, int min_defects) {
    int ip1 = (i+1)%size;
    if (b.key[i] == EMPTY && b.key[ip1] == EMPTY)
    {
        r_matrix_empty_empty(v, b, i);
    }
    else if (b.key[i] == EMPTY || b.key[ip1] == EMPTY)
    {
        r_matrix_empty_occ(v, b, i);
    }
    else
    {
        r_matrix_occ_occ(v, b, i, min_defects);
    }
}
#pragma endregion

#pragma region Auxiliary spaces
void insert_aux_space(Vec* v, BV b) {
    b.key.shift_right();

    BV b2 = b;

    b.key.set(size-1, EMPTY);
    b.key.set(0, EMPTY);
    *v += b;

    b2.key.set(size-1, OPENING);
    b2.key.set(0, CLOSING);
    *v += b2;
}

void contract_aux_space(Vec* v, BV b, int min_defects) {
    int i1 = size-2, i2 = size-1;
    if (b.key[i1] == EMPTY && b.key[i2] == EMPTY)
    {
        *v += b;
    }
    else if (b.key.can_contract_sites(i1, min_defects)
             && b.key[i1] > EMPTY && b.key[i2] > EMPTY)
    { // both sites are occupied, at least one of them isn't a defect
        contract(b, i1, min_defects);
        b.key.set(i1, EMPTY);
        b.key.set(i2, EMPTY);
        *v += b;
    }
}


#pragma endregion

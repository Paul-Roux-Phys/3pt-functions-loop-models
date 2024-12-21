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
using boost::multiprecision::cpp_complex_double;

#ifndef LATTICE_SIZE
#define LATTICE_SIZE 4
#endif

constexpr std::size_t lattice_size = LATTICE_SIZE;
constexpr std::size_t size         = lattice_size + 2;
using K = OnKey<size>;
using V = cpp_complex_100;
using BV = BasisVector<K, V>;
using Vec = Vector<BV, key_64_bit_hash_t<size>>;
using R = RMatrix<Vec, int>;

Vec v[2];
Vec initial(2);
VectorPair vp(&v[0], &v[1]);

V lambda;
V n_loop, n_ncloop, w_empty, w_turn, w_straight, w_full;
V pi("3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230");
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


// Declare the R-matrix, auxiliary spaces;
void r_matrix(Vec* v, BV b, int i, int min_defects);
void insert_aux_space(Vec* v, BV b);
void contract_aux_space(Vec* v, BV b, int min_defects);
RMatrix<Vec, int, int> r_i = r_matrix;
RMatrix<Vec> ins_aux = insert_aux_space;
RMatrix<Vec, int> contr_aux = contract_aux_space;


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

void uncolor_bottom_legs(Vec* v, BV b)
{
    for (int i = 0; i < lattice_size; i++)
    {
        if (b.key.is_bottom(i))
            b.key.set(i, BOT_DEF);
    }

    *v += b;
}

RMatrix<Vec> uncolor_bottom = uncolor_bottom_legs;

void transfer(int k3)
{
    // insert aux space
    vp.mul<>(ins_aux);
    // multiply by r-matrices
    for (int i = 0; i < lattice_size; i++)
    {
        vp.mul<int, int>(r_i, i, k3);
    }
    // contract aux space
    vp.mul<int>(contr_aux, k3);
    vp.factorise_norm();
}

// find the position of the defect just before position pos
int find_previous_defect(BV b, int pos)
{
    int res = pos-1;
    while (b.key[res] < DEFECT && res > -2) res--;
    return res;
}

// cyclically permute the defects of b
BV pseudo_translation(BV b, int k, int rs, bool with_phase=true)
{
    if (k <= 1)
        return b;
    int current_defect_pos = find_previous_defect(b, lattice_size);
    int last_def = b.key[current_defect_pos];
    int previous_defect_pos;

    for (int i = 0; i < k-1; i++)
    {
        previous_defect_pos = find_previous_defect(b, current_defect_pos);
        b.key.set(current_defect_pos, b.key[previous_defect_pos]);
        current_defect_pos = previous_defect_pos;
    }
    previous_defect_pos = find_previous_defect(b, current_defect_pos);
    b.key.set(current_defect_pos, last_def);

    V imag_i(0, 1);
    V pi("3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230");
    V omega = exp(imag_i * 2 * pi * rs / k); // = e^{i\pi s}

    if (with_phase)
        b.value *= omega;
    return b;
}

BV permute_bottom_defect_labels(BV b, int k, int rs, bool with_phase=true)
{
    if (k <= 1 || rs == 0)
        return b;
    for (int j = 0; j < lattice_size; j++)
        if (b.key.is_bottom(j))
            b.key.set(j, (b.key[j] - BOT_DEF + 1) % k + BOT_DEF);

    if (with_phase)
    {
        V imag_i(0, 1);
        V omega = exp(imag_i * 2 * pi * rs / k); // = e^{i\pi s}
        b.value *= omega;
    }

    return b;
}

// add \sum_k e^{i\pi k s} u^k*b to v
void project_Zk_rep_bottom(Vec* v, BV b, int k, int rs, bool with_phase=true)
{
    if (k <= 1 || rs == 0)
        (*v)[b.key] = b.value;
    else for (int i = 0; i < k; i++)
    {
        (*v)[b.key] = b.value;
        b = permute_bottom_defect_labels(b, k, rs, with_phase);
    }
}

RMatrix<Vec, int, int, bool> proj_Zk_rep_bot = project_Zk_rep_bottom;

void project_Zk_rep_middle(Vec* v, BV b, int k, int rs, bool with_phase=true)
{
    // cout << "projecting state" << endl;
    // b.print();
    // cout << endl;

    if (k <= 1 || b.key.nb_defects() != k)
    {
        (*v)[b.key] = b.value;
        return;
    }

    for (int i = 0; i < k; i++)
    {
        (*v)[b.key] += b.value;
        b = pseudo_translation(b, k, -rs, with_phase);
    }

    // cout << "projected: ";
    // v->print();
    // cout << endl;
}

RMatrix<Vec, int, int, bool> proj_Zk_rep_middle = project_Zk_rep_middle;

void project_Zk_rep_top(Vec* v, BV b, int k, int rs, bool with_phase=true)
{
    if (k <= 1 || rs == 0 || b.key.nb_defects() != k)
        (*v)[b.key] = b.value;
    else for (int i = 0; i < k; i++)
    {
        (*v)[b.key] = b.value;
        b = pseudo_translation(b, k, rs, with_phase);
    }
}

RMatrix<Vec, int, int, bool> proj_Zk_rep_top = project_Zk_rep_top;

BV parity_reverse(BV b)
{
    int tmp;
    for (int i = 0; i < lattice_size/2; i++)
    {
        tmp = b.key[i];
        b.key.set(i, b.key[lattice_size - 1 - i]);
        b.key.set(lattice_size - 1 - i, tmp);
    }
    for (int i = 0; i < lattice_size; i++)
    { // restore the arches in the right orientation
        if (b.key[i] == OPENING)
            b.key.set(i, CLOSING);
        else if (b.key[i] == CLOSING)
            b.key.set(i, OPENING);
    }
    return b;
}

void project_parity(Vec* v, BV b, int parity)
{
    (*v)[b.key] += b.value;
    b = parity_reverse(b);
    (*v)[b.key] += parity*b.value;
}

RMatrix<Vec, int> proj_parity = project_parity;
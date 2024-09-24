/* This file implements the transfer matrix that generates configurations
of the integrable O(n) model defined in
[[https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.69.710]
[Warnaar Nienhuis Seaton new construction of solvable 
lattice models including an Ising model in a field]                       */

#include "TransferMatrices.hpp"
#include "link_patterns.hpp"
#include <stdexcept>
using namespace transfer_matrices;

#ifndef LATTICE_SIZE
#define LATTICE_SIZE 4
#endif

#pragma region Types
constexpr std::size_t lattice_size = LATTICE_SIZE;
constexpr std::size_t size         = lattice_size + 2;
using K = FKKey<size>;
using V = double;
using BV = BasisVector<K, V>;
using Vec = Vector<BV, key_64_bit_hash_t<size>>;
using R = RMatrix<Vec, int>;
#pragma endregion

Vec v[2];
Vec initial(2);
VectorPair vp(&v[0], &v[1]);

#pragma region Weights
V lambda     = 0.5;
V n_loop     = -2*cos(4*lambda); // weight of loops
V n_ncloop   = -2*cos(4*lambda); // weight of non-contractible loops
V w_empty    = 1 + sin(lambda) + sin(3*lambda) - sin(5*lambda);
V w_turn     = 2*sin(2*lambda)*sin((6*lambda+M_PI)/4);
V w_straight = 1 + sin(3*lambda);
V w_full     = sin(lambda) + cos(2*lambda);

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
void r_matrix(Vec* v, BV b, int i);
void insert_aux_space(Vec* v, BV b);
void contract_aux_space(Vec* v, BV b);
RMatrix<Vec, int> r_i = r_matrix;
RMatrix<Vec> ins_aux = insert_aux_space;
RMatrix<Vec> contr_aux = contract_aux_space;
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

void contract(BV& b, int i) {
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
    b.key.contract(i);
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

void r_matrix_occ_occ(Vec* v, BV b, int i) {
    BV b2 = b;
    
    b.value *= w_full;
    *v += b;

    if (!(b.key[i] == DEFECT && b.key[(i+1)%size] == DEFECT)) {
        contract(b2, i);
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

void r_matrix(Vec* v, BV b, int i) {
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
        r_matrix_occ_occ(v, b, i);
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

void contract_aux_space(Vec* v, BV b) {
    int i1 = size-2, i2 = size-1;
    if (b.key[i1] == EMPTY && b.key[i2] == EMPTY)
    {
        *v += b;
        return;
    }
    if (b.key[i1] != EMPTY && b.key[i2] != EMPTY \
             && !(b.key[i1] == DEFECT && b.key[i2] == DEFECT))
    { // both sites are occupied, at least one of them isn't a defect
        contract(b, i1);
        b.key.set(i1, EMPTY);
        b.key.set(i2, EMPTY);
        *v += b;
    }
}
#pragma endregion

#pragma region Transfer matrix
// filter states with too low weights
void filter(Vec* v, BV b) {
    if (abs(b.value) > 1e-14)
    {
        *v += b;
    }
}

void project(Vec* v, BV b) {
    (*v)[b.key] = b.value;
    half_lattice_translation(b.key);
    (*v)[b.key] = -b.value;
}

RMatrix<Vec> filt = filter;
RMatrix<Vec> proj = project;

// multiplies v by the transfer matrix.
template <typename Vec>
void VectorPair<Vec>::transfer() {
    // insert aux space
    mul<>(ins_aux);
    // multiply by r-matrices
    for (int i = 0; i < lattice_size; i++) {
        mul<int>(r_i, i);
    }
    // contract aux space
    mul<>(contr_aux);
    // mul<>(proj);
    factorise_norm();
    // mul<>(filt);
}
#pragma endregion

K initial_key(int nb_defects) {
    K res;
    for (int i = 0; i < size; i++)
    {
        res.set(i, EMPTY);
    }
    for (int i = 0 ; i < nb_defects; i++)
    {
        res.set(i, DEFECT);
    }
    return res;
}

void initialise_vector() {
    K k1 = initial_key(2);
    K k2 = k1;

    half_lattice_translation(k2);

    BV b1(k1, 1.0);
    BV b2(k2, -1.0);

    initial += b1;
    initial += b2;
    initial /= std::sqrt(2);
}

#pragma region main function
#define O OPENING
#define C CLOSING
#define E EMPTY
#define D DEFECT
int main() {
    initialise_vector();

    vp += initial;

    for (int i = 0; i < 200; i++)
    {
        vp.transfer();
    }

    LargeFloat res = vp.inner_product(initial);
    cout << res << endl;
    return 0;
}
#pragma endregion
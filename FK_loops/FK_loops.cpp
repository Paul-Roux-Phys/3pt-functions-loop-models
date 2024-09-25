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

void print_weights() {
    cout << "lambda = " << lambda << endl;
    cout << "n_loop = " << n_loop << endl;
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

void tl_generator(BV& b, int i) {
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
    b.key.tl_generator(i);
}

void r_matrix(Vec* v, BV b, int i) {
    BV b2 = b;
    
    *v += b;

    if (!(b2.key[i] == DEFECT && b2.key[(i+1)%size] == DEFECT)) {
        tl_generator(b2, i);
        *v += b2;
    }
}
#pragma endregion

#pragma region Auxiliary spaces
void insert_aux_space(Vec* v, BV b) {
    b.key.shift_right();
    b.key.set(size-1, OPENING);
    b.key.set(0, CLOSING);
    *v += b;
}

void contract_aux_space(Vec* v, BV b) {
    int i1 = size-2, i2 = size-1;

    if (!(b.key[i1] == DEFECT && b.key[i2] == DEFECT))
    { // both sites are occupied, at least one of them isn't a defect
        tl_generator(b, i1);
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
        if (i%2 == 0)
            res.set(i, OPENING);
        else
            res.set(i, CLOSING);
    }
    for (int i = 0 ; i < nb_defects; i++)
    {
        res.set(i, DEFECT);
    }
    return res;
}

void initialise_vector(int nb_defects) {
    K k1 = initial_key(nb_defects);
    K k2 = k1;

    half_lattice_translation(k2);

    BV b1(k1, 1.0);
    // BV b2(k2, -1.0);

    initial += b1;
    // initial += b2;
    // initial /= std::sqrt(2);
}

#pragma region main function
#define C CLOSING
#define E EMPTY
#define D DEFECT
int main() {
    print_weights();
    initialise_vector(2);

    vp += initial;
    
    vp.print();
    for (int i = 0; i < 2; i++)
    {
        vp.transfer();
        if (i > 20) cout << vp.inner_product(initial) << endl;
    }
    vp.print();
    return 0;

    // Matrix<VectorPair<Vec>> M(&vp, 10);
    // M.find_eigenvalues(true);
    // cout << M.sprint_eigenvalues();

}
#pragma endregion
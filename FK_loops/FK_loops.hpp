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

// #define AUX_SPACE

constexpr std::size_t lattice_size = LATTICE_SIZE;
constexpr std::size_t size         = lattice_size + 2;
using K = FKKey<size>;
using V = double;
using BV = BasisVector<K, V>;
using Vec = Vector<BV, key_64_bit_hash_t<size>>;
using R = RMatrix<Vec, int>;

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

void tl_generator(BV& b, int i, int min_defects) {
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
    b.key.tl_generator(i, min_defects);
}


void r_matrix(Vec* v, BV b, int i, int min_defects) {
    BV b2 = b;
    
    *v += b;

    if (b.key.can_contract_sites(i, min_defects))
    {
        tl_generator(b2, i, min_defects);
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

void contract_aux_space(Vec* v, BV b, int min_defects) {
    int i1 = size-2;
    if (b.key.can_contract_sites(i1, min_defects))
    {
        tl_generator(b, i1, min_defects);
        *v += b;
    }
}
#pragma endregion

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

class Amplitude : public VectorPair<Vec> {
private:
    int _bm, _bt, _mb, _mt;
public:
    int k1, k2, k3;

    // use base class constructors
    using VectorPair<Vec>::VectorPair;

    Amplitude(Vec* v1, Vec* v2, int k1, int k2, int k3)
    : VectorPair<Vec>(v1, v2), k1(k1), k2(k2), k3(k3)
    {
        _bm = (k1 + k2 - k3)/2;
        _bt = (k1 - k2 + k3)/2;
        _mb = _bm;
        _mt = (-k1 + k2 + k3)/2;
    }

    void insert_bottom_operator(int current=0);
    void transfer();
    LargeFloat sandwich_with_top_operator(int current=0);
    LargeFloat compute_three_point_fct(int half_size = 20*lattice_size);
};

void Amplitude::insert_bottom_operator(int current) {
    K key1;
    K key2;
    for (int i = 0; i < k1; i++)
        key1.set(i, BOTTOM+i);
    for (int i = 0; i < k1; i++)
        key2.set(i, BOTTOM+k1-1-i);

    BV b1(key1, 1.0);
    BV b2(key2, 1.0);
    *this += b1;
    *this += b2;
    // this->normalise();
}

void insert_middle_operator(Vec* v, BV b, int k2, int k3) {
    // assume k2=2 for now
    if (k2 == 0)
    {
        *v += b;
        return;
    }
    // glue to one up-facing and one down-facing defect in the middle
    if (b.key.is_bottom(0))
    {
        b.key.set(0, MIDDLE);
        *v += b;
    }
    if (b.key[0] < DEFECT)
    {
        int end_pos = b.key.arch_end(0);

        b.key.set(0, MIDDLE);
        b.key.set(end_pos, MIDDLE+1);
        *v += b;
    }
}

RMatrix<Vec, int, int> insert_mid = insert_middle_operator;

void tl_Rmat(Vec* v, BV b, int i, int min_defects) {
    if (b.key.can_contract_sites(i, min_defects))
    {
        tl_generator(b, i, min_defects);
        *v += b;
    }
}

RMatrix<Vec, int, int> tl_gen = tl_Rmat;

LargeFloat Amplitude::sandwich_with_top_operator(int current) {
    Vec top(2);
    K key1;
    K key2;
    for (int i = 0; i < k1; i++)
        key1.set(i, BOTTOM+i);
    for (int i = k1; i < k1 + k2; i++)
        key1.set(i, MIDDLE+i-k1);

    for (int i = 0; i < k2; i++)
        key2.set(i, MIDDLE+i);
    for (int i = k2; i < k1 + k2; i++)
        key2.set(i, BOTTOM+i-k2);

    BV b1(key1, 1.0);
    BV b2(key2, -1.0);
    top += b1;

    if (current)
        top += b2;

    return this->inner_product(top);
}

// multiplies v by the transfer matrix.
void Amplitude::transfer() {
    // for (int i = 0; i < lattice_size; i += 2)
    // {
    //     mul<int>(r_i, i);
    // }
    // for (int i = 1; i <= lattice_size; i += 2)
    // {
    //     mul<int>(r_i, i);
    // }
    // insert aux space
    mul<>(ins_aux);
    // multiply by r-matrices
    for (int i = 0; i < lattice_size; i++) {
        mul<int, int>(r_i, i, k3);
    }
    // contract aux space
    mul<int>(contr_aux, k3);
    // mul<>(proj);
    // factorise_norm();
    // mul<>(filt);
}

LargeFloat Amplitude::compute_three_point_fct(int half_size) {
    insert_bottom_operator();

    // lower half
    for (int i = 0; i < half_size; i++)
    {
        transfer();
    }

    // middle insertion
    mul<int, int>(insert_mid, _mb, _mt);

    // upper half
    for (int i = 0; i < half_size; i++)
    {
        transfer();
    }

    return sandwich_with_top_operator();    
}

LargeFloat compute_structure_constant() {
    Amplitude a1(&v[0], &v[1], 2, 2, 2);
    LargeFloat c222 = a1.compute_three_point_fct();
    vp.clear();

    Amplitude a2(&v[0], &v[1], 2, 2, 0);
    LargeFloat c220 = a1.compute_three_point_fct();
    vp.clear();

    Amplitude a3(&v[0], &v[1], 2, 0, 2);
    LargeFloat c202 = a1.compute_three_point_fct();
    vp.clear();

    Amplitude a4(&v[0], &v[1], 0, 0, 0);
    LargeFloat c000 = a4.compute_three_point_fct();
    vp.clear();

    return c222/c220*sqrt(c000/c202);
}
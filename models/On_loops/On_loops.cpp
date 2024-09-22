/* This file implements the transfer matrix that generates configurations
of the integrable O(n) model defined in
[[https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.69.710]
[Warnaar Nienhuis Seaton new construction of solvable 
lattice models including an Ising model in a field]                       */

#include "TransferMatrices.hpp"
#include "link_patterns.hpp"
using namespace transfer_matrices;

#pragma region Types
constexpr std::size_t lattice_size = 4;
constexpr std::size_t size         = lattice_size + 2;
using K = FKKey<size>;
using V = double;
using BV = BasisVector<K, V>;
using Vec = Vector<BV, key_64_bit_hash_t<size>>;
using R = RMatrix<Vec, int>;
#pragma endregion

Vec v[2];
int data_position = 0;

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
    cout << "w_empty = " << w_empty << endl;
    cout << "w_turn = " << w_turn << endl;
    cout << "w_straight = " << w_straight << endl;
    cout << "w_full = " << w_full << endl;
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
    if (b.key[i1] != EMPTY && b.key[i2] != EMPTY \
             && !(b.key[i1] == DEFECT && b.key[i2] == DEFECT))
    { // both sites are occupied, at least one of them isn't a defect
        contract(b, i1);
        b.key.set(i1, EMPTY);
        b.key.set(i2, EMPTY);
    }
    *v += b;
}
#pragma endregion

#pragma region Transfer matrix
// multiplies v by the transfer matrix.
// returns a pointer to the vector in which the result is stored.
void transfer() {
    // insert aux space
    multiply<Vec>(&v[1-data_position], ins_aux, &v[data_position]);
    data_position = 1 - data_position;
    // multiply by r-matrices
    for (int i = 0; i < lattice_size; i++) {
        multiply<Vec, int>(&v[1-data_position], r_i, &v[data_position], i);
        data_position = 1 - data_position;
    }
    // contract aux space
    multiply<Vec>(&v[1-data_position], contr_aux, &v[data_position]);
    data_position = 1 - data_position;
}
#pragma endregion

#pragma region main function
#define O OPENING
#define C CLOSING
#define E EMPTY
#define D DEFECT
int main() {
    print_weights();
    BV init({D, O, C, D, E, E}, 1.0);
    v[0] += init;

    transfer();

    v[data_position].print();
    // cout << "------------------" << endl << endl;
    // MatrixMul(&v1, &v2, );
    
    return 0;
}
#pragma endregion
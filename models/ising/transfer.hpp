#pragma once

#include "TransferMatrices.hpp"
using namespace transfer_matrices;

constexpr std::size_t size = 4;

using K = int8_key_t<size>;
using V = double;
using BV = BasisVector<K, V>;
using Vec = Vector<BV, KeyHash<K>>;
using R = RMatrix<Vec, int>;

#include "symmetrise.hpp"

Vec v1, v2;
Vec *p1 = &v1, *p2 = &v2; // pointers to the two vectors, for easy swap
void r_matrix(Vec*, BV, int); // forward declaration of r_matrix, implemented below
void symmetrise(Vec*, BV);
void antisymmetrise(Vec*, BV);
R r = r_matrix;
RMatrix<Vec> as = antisymmetrise;
RMatrix<Vec> s = symmetrise;

V J(log(1+sqrt(2))/2);
V expJ  = exp(J);
V expmJ = exp(-J);

void add_vertical_edge(Vec* v, BV b, int pos)
{
    BV b2 = b;

    b.value *= expJ;

    b2.key[pos] *= -1;
    b2.value *= expmJ;
    // project_sym(b2);
    // project_antisym(b2);

    *v += b;
    *v += b2;
}

void add_horizontal_edge(Vec* v, BV b, int pos1, int pos2)
{
    if (b.key[pos1] == b.key[pos2])
        b.value *= expJ;
    else 
        b.value *= expmJ;

    *v += b;
}

void r_matrix(Vec* v, BV b, int pos)
{
    if (pos%2 == 0)
        add_vertical_edge(v, b, pos/2);
    else
        add_horizontal_edge(v, b, pos/2, (pos/2+1)%size);
}

void swap () {
    Vec* tmp = p2;
    p2 = p1;
    p1 = tmp;
}

void transfer() {
    for (int pos = 0; pos < size; pos++) {
        multiply<Vec, int>(p1, r, p2, 2*pos);
        swap();
    }
    for (int pos = 0; pos < size; pos++) {
        multiply<Vec, int>(p1, r, p2, 2*pos+1);
        swap();
    }
    swap();
}
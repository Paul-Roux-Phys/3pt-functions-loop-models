#pragma once
#include "r_matrix.hpp"
#include "three_point_fcts.hpp"

void project_on_antisymmetric_part(Vec* v);

Vec two_defects_antisymmetrized_state()
{
    Vec res(2);
    K k;
    k.set(0, BOT_DEF);
    k.set(1, BOT_DEF+1);

    res += BV(k, 1);

    k.set(1, BOT_DEF);
    k.set(0, BOT_DEF+1);

    res += BV(k, -1);

    return res;
}

Vec two_defects_symmetrized_state()
{
    Vec res(2);
    K k;
    k.set(0, BOT_DEF);
    k.set(1, BOT_DEF+1);

    res += BV(k, 1);

    k.set(1, BOT_DEF);
    k.set(0, BOT_DEF+1);

    res += BV(k, 1);
    
    return res;
}

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

LargeFloat compute_two_point_function_current(int half_size)
{
    Vec init_end_state(2);
    init_end_state = two_defects_antisymmetrized_state();
    vp += init_end_state;

    vp /= 2;

    for (int i = 0; i < 2*half_size; i++)
    {
        transfer(2);
        project_on_antisymmetric_part(vp.get());
    }

    return vp.inner_product(init_end_state);
}

LargeFloat compute_two_point_function_2leg(int half_size)
{
    Vec init_end_state(2);
    init_end_state = two_defects_symmetrized_state();
    vp += init_end_state;

    vp /= 2;

    for (int i = 0; i < 2*half_size; i++)
    {
        transfer(2);
    }

    return vp.inner_product(init_end_state);
}
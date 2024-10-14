#pragma once
#include "r_matrix.hpp"
#include "two_point_fcts.hpp"

void insert_middle_operator(Vec *v, BV b, int k2, int k3);
RMatrix<Vec, int, int> insert_mid = insert_middle_operator;

void project_on_antisymmetric_part(Vec *v)
{
    for (auto it = v->begin(); it != v->end(); it++)
    {
        int pos_def1 = -1, pos_def2 = -1; // -1 means defect not found
        BV b(*it);
        for (int i = 0; i < lattice_size; i++)
        {
            if (b.key[i] == BOT_DEF)
                pos_def1 = i;
            if (b.key[i] == BOT_DEF + 1)
                pos_def2 = i;
        }

        K k_op = b.key;
        if (pos_def1 != -1)
            k_op.set(pos_def1, BOT_DEF + 1);
        if (pos_def2 != -1)
            k_op.set(pos_def2, BOT_DEF);

        (*v)[k_op] = -(*v)[b.key];
    }
}

void insert_middle_operator(Vec *v, BV b, int k2, int k3)
{
    // assume k2=0 or 2 for now
    // if (k2 == 0)
    // {
    //     *v += b;
    //     return;
    // }
    // glue to one up-facing and one down-facing defect in the middle
    if (b.key.is_bottom(0))
    {
        b.key.set(0, MID_DEF); 
        *v += b;
    }
    else if (b.key[0] == EMPTY)
    {
        return;
    }
    else if (b.key[0] < DEFECT)
    {
        int end_pos = b.key.arch_end(0);

        b.key.set(0, MID_DEF);
        b.key.set(end_pos, MID_DEF);
        *v += b;
    }
}

LargeFloat compute_three_point_function_current_2leg_2leg(int half_size)
{
    Vec init_state(2);
    init_state = two_defects_antisymmetrized_state();
    vp += init_state;
    
    vp /= 2;

    for (int i = 0; i < half_size; i++)
    {
        transfer(2);
        project_on_antisymmetric_part(vp.get());
    }

    // insert middle operator
    vp.mul<int, int>(insert_mid, 2, 2);

    for (int i = 0; i < half_size; i++)
    {
        transfer(2);
        project_on_antisymmetric_part(vp.get());
    }

    Vec end_state(4);
    K k_end;
    k_end.set(0, BOT_DEF);
    k_end.set(1, MID_DEF);
    end_state += BV(k_end, 1);
    k_end.set(0, MID_DEF);
    k_end.set(1, BOT_DEF);
    end_state += BV(k_end, 1);
    project_on_antisymmetric_part(&end_state);

    return vp.inner_product(end_state);
}

LargeFloat compute_three_point_function_2leg_2leg_id(int half_size)
{
    Vec init_state(2);
    init_state = two_defects_symmetrized_state();
    vp += init_state;

    vp /= 2;

    for (int i = 0; i < half_size; i++)
    {
        transfer(0);
    }

    // insert middle operator
    vp.mul<int, int>(insert_mid, 2, 0);

    for (int i = 0; i < half_size; i++)
    {
        transfer(0);
    }

    Vec top_state(1);
    K k;
    top_state += BV(k, 1);
    return vp.inner_product(top_state);
}

LargeFloat compute_three_point_function_2leg_2leg_2leg(int half_size)
{
    Vec init_state(2);
    init_state = two_defects_symmetrized_state();
    vp += init_state;

    for (int i = 0; i < half_size; i++)
    {
        transfer(2);
    }

    // insert middle operator
    vp.mul<int, int>(insert_mid, 2, 2);

    for (int i = 0; i < half_size; i++)
    {
        transfer(2);
    }

    Vec end_state(2);
    K k;
    k.set(0, BOT_DEF);
    k.set(1, MID_DEF);
    end_state += BV(k, 1);
    k.set(0, MID_DEF);
    k.set(1, BOT_DEF);
    end_state += BV(k, 1);

    return vp.inner_product(end_state);
}

LargeFloat compute_partition_function(int half_size)
{
    Vec init_state(1);
    K k;
    init_state += BV(k, 1);
    vp += init_state;

    for (int i = 0; i < 2*half_size; i++)
    {
        transfer(0);
    }

    return vp.inner_product(init_state);
}
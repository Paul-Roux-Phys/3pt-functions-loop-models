#include "r_matrix.hpp"

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

void project_on_antisymmetric_part()
{
    int pos_def1, pos_def2;
    for (auto it = vp.begin(); it != vp.end(); it++)
    {
        BV b(*it);
        for (int i = 0; i < lattice_size; i++)
        {
            if (b.key[i] == BOT_DEF)
                pos_def1 = i;
            if (b.key[i] == BOT_DEF+1)
                pos_def2 = i;
        }

        K k_op = b.key;
        k_op.set(pos_def1, BOT_DEF+1);
        k_op.set(pos_def2, BOT_DEF);

        (*v)[k_op] = -(*v)[b.key];
    }
}

void transfer()
{
    // insert aux space
    vp.mul<>(ins_aux);
    // multiply by r-matrices
    for (int i = 0; i < lattice_size; i++)
    {
        vp.mul<int, int>(r_i, i, 2);
    }
    // contract aux space
    vp.mul<int>(contr_aux, 2);
    vp.factorise_norm();
}

LargeFloat compute_current_two_point_function(int cylinder_size)
{
    Vec init_end_state(2);
    init_end_state = two_defects_antisymmetrized_state();
    vp += init_end_state;

    for (int i = 0; i < cylinder_size; i++)
    {
        transfer();
        project_on_antisymmetric_part();
        if (i >= cylinder_size / 2)
            cout << i + cylinder_size << "\t"
                 << std::setprecision(15) << vp.inner_product(init_end_state)
                 << endl;
    }

    return vp.inner_product(init_end_state);
}
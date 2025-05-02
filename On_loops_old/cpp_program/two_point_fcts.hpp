#pragma once
#include "r_matrix.hpp"
#include "three_point_fcts.hpp"

// insert k legs carrying the representation e^{2i\pi s/L} of Z_k. The second argument
// is r*s where r = k/2
Vec bottom_state(int k, int rs, int starting_position, V P1=0)
{
    if (P1 != 0)
        n_ncloop = 2*cos(4*sqrt(pi*lambda)*P1);
    Vec res(k);
    if (k == 0)
    {
        K key;
        res += BV(key, 1);
    }
    else
    {
        BV b;
        b.value = 1;
        for (int i = 0; i < k; i++)
            b.key.set((i+starting_position)%lattice_size, BOT_DEF + (rs != 0) * i);
        project_Zk_rep_bottom(&res, b, k, rs);
        if (k > 1 && rs != 0)
            res /= k;
    }

    return res;
}

LargeFloat compute_2pt_function(int half_size, int k, int rs, int starting_position)
{
    vp += bottom_state(k, rs, starting_position);
    // vp.mul<int>(proj_parity, 1);
    Vec end_state_no_parity = bottom_state(k, -rs, 0);
    Vec end_state;
    end_state = end_state_no_parity;
    // multiply(&end_state, proj_parity, &end_state_no_parity, 1);

    for (int i = 0; i < 2*half_size; i++)
    {
        transfer(k);
        // vp.mul<int>(proj_parity, 1);
        // if (rs != 0)
        //     vp.mul<int, int>(proj_Zk_rep_bot, k, rs, true); // project on Z_k irrep
    }

    LargeFloat res = vp.inner_product(end_state);
    vp.clear();
    return res;
}
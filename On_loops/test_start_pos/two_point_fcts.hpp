#pragma once
#include "r_matrix.hpp"
#include "three_point_fcts.hpp"

// insert k legs carrying the representation e^{2i\pi s/L} of Z_k. The second argument
// is r*s where r = k/2
Vec bottom_state(int k1, int k2, int k3, int pos1, int pos2)
{
    int mid = (k1 + k2 - k3) / 2, top = (k1 - k2 + k3) / 2;
    Vec res(k1);
    K key;
    if (mid == 0 && top == 0)
    {
        // do nothing
    }
    else if (mid == 0 && top == 2)
    {
        key.set(pos1, BOT_DEF + 1);
        key.set(pos2, BOT_DEF + 1);
    }
    else if (mid == 1 && top == 1)
    {
        key.set(pos1, BOT_DEF);
        key.set(pos2, BOT_DEF + 1);
        res += BV(key, 1);
        // res /= k;
    }
    else if (mid == 2 && top == 0)
    {
        key.set(pos1, BOT_DEF);
        key.set(pos2, BOT_DEF);
    }
    else
    {
        fprintf(stderr, "Error: the current bottom_state() function only works for 0 or 2 defects");
    }

    res += BV(key, 1);

    return res;
}

LargeFloat compute_2pt_function(int half_size, int k, int pos1, int pos2, int pos3, int pos4)
{
    vp += bottom_state(k, 0, k, pos1, pos2);
    Vec end_state = bottom_state(k, 0, k, pos3, pos4);

    // cout << "initial state" << endl;
    // vp.print();
    // cout << "end state: " << endl;
    // end_state.print();

    for (int i = 0; i < 2 * half_size; i++)
    {
        transfer(k);
    }

    LargeFloat res = vp.inner_product(end_state);
    vp.clear();
    return res;
}
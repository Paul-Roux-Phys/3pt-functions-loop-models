#pragma once
#include "r_matrix.hpp"
#include "two_point_fcts.hpp"

// pos must be < max(max_down_def, max_up_def)
// this function handles the insertion of the middle operator at position pos
// void insert_middle_defect(BV& b, int max_down_def, int max_up_def, int pos)
// {
//     // cout << "inserting "; b.print(); 
//     if (b.key.is_bottom(pos))
//     {
//         if (pos > max_down_def)
//         {
//             b.key.set(0, ERROR); // error
//             return;
//         }
//         if (max_up_def < pos && pos < max_down_def)
//             b.key.set(pos, EMPTY);
//         else
//             b.key.set(pos, MID_DEF);
//     }
//     else if (b.key[pos] == EMPTY)
//     {
//         if (max_down_def < pos && pos <= max_up_def)
//             b.key.set(pos, MID_DEF);
//         else
//         {
//             b.key.set(0, ERROR); // error
//             return;
//         }
//     }
//     else if (b.key[pos] < DEFECT)
//     {
//         int end_pos = b.key.arch_end(pos);
//         if (pos > max_down_def)
//         {
//             b.key.set(0, ERROR); // error
//             return;
//         }
//         if (end_pos < max_down_def || end_pos < max_up_def)
//         {
//             b.key.set(0, ERROR); // error
//             return;
//         }
//         if (pos > max_up_def)
//             b.key.set(pos, EMPTY);
//         else
//             b.key.set(pos, MID_DEF);
//         b.key.set(end_pos, MID_DEF);
//     }
//     // cout << "inserted "; b.print(); cout << endl;
// }

// void insert_middle_operator(Vec *v, BV b, int k2, int k3, V P3=0)
// {
//     if (P3 != 0)
//         n_ncloop = 2*cos(4*sqrt(pi*lambda)*P3);
//     if (k2 % 2 == 0)
//     {
//         for (int i = 0; i < k2 / 2; i++)
//             insert_middle_defect(b, k2 / 2 - 1, k2 / 2 - 1, i);
//         if (!(b.key[0] == ERROR))
//             *v += b;
//     }
//     else
//     {
//         for (int i = 0; i < k2 / 2 + 1; i++)
//             insert_middle_defect(b, k2 / 2 - 1, k2/2 , i);
//         if (!(b.key[0] == ERROR))
//             *v += b;
//     }
// }

void insert_middle_operator(Vec *v, BV b, int k2, int k3)
{
    if (b.key[0] == BOT_DEF)
    {
        b.key.set(0, MID_DEF);
    }
    else if (b.key[0] == BOT_DEF+1)
    {
        b.key.set(0, ERROR);
    }
    else if (b.key[0] == EMPTY)
    {
        b.key.set(0, ERROR);
    }
    else if (b.key[0] == OPENING || b.key[0] == CLOSING)
    {
        b.key.set(b.key.arch_end(0), MID_DEF);
        b.key.set(0, MID_DEF);
    }
    if (b.key[0] != ERROR)
        *v += b;
}

RMatrix<Vec, int, int> insert_mid = insert_middle_operator;

// TODO: allow for non-zero s in the middle
// BV permute_top_defect_labels(BV b, int k1, int rs1, int k3, int rs3)
// {
//     if (k3 <= 1)
//         return b;
//     for (int j = 0; j < lattice_size; j++)
//         if (rs1 == 0)
//         {
//             if (b.key.is_bottom(j))
//                 b.key.set(j, MID_DEF);
//             else if (b.key.is_middle(j))
//                 b.key.set(j, BOT_DEF);
//         }
//         else
//         {
//             if (b.key.is_bottom(j) && b.key[j] < BOT_DEF + k1 - 1)
//                 b.key.set(j, b.key[j] + 1);
//             else if (b.key.is_bottom(j) && b.key[j] == BOT_DEF + k1 - 1)
//                 b.key.set(j, MID_DEF);
//             else if (b.key.is_middle(j))
//                 b.key.set(j, BOT_DEF);
//         }

//     V imag_i(0, 1);
//     V pi("3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230");
//     V omega = exp(imag_i * 2 * pi * rs3 / k3); // = e^{i\pi s}

//     b.value *= omega;

//     return b;
// }

// Vec end_state(int k1, int k2, int k3, int rs1, int rs2, int rs3)
// {
//     Vec res(10);
//     // put the right number of legs coming from bot or mid, with correct labels
//     int bot = (k1 - k2 + k3) / 2;
//     int mid = (-k1 + k2 + k3) / 2;
//     BV b;
//     b.value = 1;
//     for (int i = 0; i < bot; i++)
//         b.key.set(i, BOT_DEF + (rs1 != 0) * i);
//     for (int i = bot; i < bot + mid; i++)
//         b.key.set(i, MID_DEF + (rs2 != 0) * (i-bot));

//     project_Zk_rep_bottom(&res, b, k1, rs1, true);

//     for (int i = 0; i < k3; i++)
//     {
//         project_Zk_rep_bottom(&res, b, k1, rs1, true);
//         b = pseudo_translation(b, k3, rs3);
//     }
//     // if (k3 > 1 && rs3 != 0)
//     //     res /= k3;

//     return res;
// }

Vec end_state(int k1, int k2, int k3, int pos1, int pos2)
{
    Vec res(2);

    int bot = (k1 - k2 + k3) / 2;
    int mid = (-k1 + k2 + k3) / 2;

    int label1 = EMPTY, label2 = EMPTY;

    if (bot == 2 && mid == 0)
    {
        label1 = BOT_DEF;
        label2 = BOT_DEF+1;
    } else if (bot == 1 && mid == 1)
    {
        label1 = BOT_DEF+1;
        label2 = MID_DEF;
    }
    else if (bot == 0 && mid == 2)
    {
        label1 = MID_DEF;
        label2 = MID_DEF+1;
    }
    K key;
    key.set(pos1, label1);
    key.set(pos2, label2);

    res += BV(key, 1);
    return res;
}

LargeFloat compute_3pt_function(int half_size, int k1, int k2, int k3, int pos1, int pos2, int pos3, int pos4)
{
    if (k2 == 0)
    {
        return compute_2pt_function(half_size, k1, pos1, pos2, pos3, pos4);
    }

    // else

    // initial state
    vp += bottom_state(k1, k2, k3, pos1, pos2);

    // cout << "initial state" << endl;
    // vp.print();
    // end state
    Vec end = end_state(k1, k2, k3, pos3, pos4);
    // cout << "end state" << endl;
    // end.print();

    LargeFloat res;
    // vp.print();

    for (int i = 0; i < half_size; i++)
    {
        transfer(k3);
    }

    
    // cout << "before gluing" << endl;
    // vp.print();
    vp.mul<int, int>(insert_mid, k2, k3);
    // vp.mul<>(uncolor_bottom);
    // cout << "after gluing" << endl;
    // vp.print();

    for (int i = 0; i < half_size; i++)
    {
        transfer(k3);
    }
    // vp.print();

    res = vp.inner_product(end);
    vp.clear();
    return res;
}
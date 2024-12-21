#include "r_matrix.hpp"
#include <complex>
#include <boost/multiprecision/cpp_dec_float.hpp>
using boost::multiprecision::cpp_dec_float_100;
cpp_dec_float_100 pi = boost::math::constants::pi<cpp_dec_float_100>();

class Amplitude : public VectorPair<Vec>
{
private:
    int _bm, _bt, _mb, _mt;

public:
    int k1, k2, k3;

    using VectorPair<Vec>::VectorPair; // use base class constructors

    Amplitude(Vec *v1, Vec *v2, int k1, int k2, int k3)
        : VectorPair<Vec>(v1, v2), k1(k1), k2(k2), k3(k3)
    {
        _bm = (k1 + k2 - k3) / 2;
        _bt = (k1 - k2 + k3) / 2;
        _mb = _bm;
        _mt = (-k1 + k2 + k3) / 2;
    }

    void insert_bottom_operator();
    void transfer();
    LargeFloat sandwich_with_top_operator(int current = 0);
    LargeFloat compute_two_point_amplitude_current_with_phases();
    LargeFloat compute_three_point_fct(int half_size = 20 * lattice_size);
};

// multiplies v by the transfer matrix.
void Amplitude::transfer()
{
    // insert aux space
    mul<>(ins_aux);
    // multiply by r-matrices
    for (int i = 0; i < lattice_size; i++)
    {
        mul<int, int>(r_i, i, k3);
    }
    // contract aux space
    mul<int>(contr_aux, k3);
    factorise_norm();
}

void Amplitude::insert_bottom_operator()
{
    K key1;
    K key2;
    for (int i = 0; i < k1; i++)
        key1.set(i, BOT_DEF + i);
    for (int i = 0; i < k1; i++)
        key2.set(i, BOT_DEF + k1 - 1 - i);

    BV b1(key1, 1.0);
    BV b2(key2, -1.0);
    *this += b1;
    *this += b2;
    this->normalise();
}

void insert_middle_operator(Vec *v, BV b, int k2, int k3)
{
    // assume k2=0 or 2 for now
    if (k2 == 0)
    {
        *v += b;
        return;
    }
    // glue to one up-facing and one down-facing defect in the middle
    if (b.key.is_bottom(0))
    {
        b.key.set(0, MID_DEF);
        *v += b;
    }
    if (b.key[0] < DEFECT)
    {
        int end_pos = b.key.arch_end(0);

        b.key.set(0, MID_DEF);
        b.key.set(end_pos, MID_DEF + 1);
        *v += b;
    }
}

RMatrix<Vec, int, int> insert_mid = insert_middle_operator;

LargeFloat Amplitude::sandwich_with_top_operator(int current)
{
    Vec top(2);
    K key1;
    K key2;
    for (int i = 0; i < k1; i++)
        key1.set(i, BOT_DEF + i);
    for (int i = k1; i < k1 + k2; i++)
        key1.set(i, MID_DEF + i - k1);

    key2 = key1;
    key2.permute_defects();

    BV b1(key1, 1.0);
    BV b2(key2, -1.0);
    top += b1;

    if (current)
        top += b2;

    // top.normalise();

    return this->inner_product(top);
}

Vec two_defect_state_with_phases(int s)
{
    Vec top(100);

    for (int i = 0; i < lattice_size - 1; i++)
    {
        for (int j = i + 1; j < lattice_size; j++)
        {
            K k, k2;
            k.set(i, BOT_DEF);
            k.set(j, BOT_DEF + 1);
            // the -1 ensures that when i=0, j=1 the angle is 0 (state inserted at the bottom)
            double total_angle = M_PI / lattice_size * s * (i + j - 1);
            V phase(0, total_angle);

            top += BV(k, exp(phase));

            k2.set(i, BOT_DEF + 1);
            k2.set(j, BOT_DEF);
            total_angle += M_PI;
            phase = V(0, total_angle);

            top += BV(k2, exp(phase));
        }
    }

    // top.normalise();

    return top;
}

LargeFloat Amplitude::compute_three_point_fct(int half_size)
{
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
        cout << half_size + i << "\t" << sandwich_with_top_operator(1) << endl;
    }

    return sandwich_with_top_operator();
}

LargeFloat Amplitude::compute_two_point_amplitude_current_with_phases()
{
    Vec bottom = two_defect_state_with_phases(1);
    Vec top = two_defect_state_with_phases(-1);
    bottom.normalise();
    top.normalise();

    *this += bottom;

    for (int i = 0; i < 20*lattice_size; i++)
    {
        transfer();
        if (i >= 10*lattice_size)
            cout << 10*lattice_size+i << "\t" << inner_product(top) << endl;
    }

    return this->inner_product(top);
}
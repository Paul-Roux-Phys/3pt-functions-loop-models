#include "FK_loops.hpp"

#pragma region main function
#define C CLOSING
#define E EMPTY
#define D DEFECT
int main() {
    // cout << compute_structure_constant() << endl;

    Amplitude a1(&v[0], &v[1], 2, 0, 2);
    a1.insert_bottom_operator();
    a1.transfer();
    a1.transfer();
    a1.transfer();
    a1.transfer();
    a1.transfer();
    a1.mul<int, int>(insert_mid, a1.k2, a1.k3);
    a1.transfer();
    a1.transfer();
    a1.transfer();
    a1.transfer();
    a1.transfer();
    a1.print();

    a1.clear();
    cout << a1.compute_three_point_fct(5);
    a1.print();

    // a1.mul<int, int>(insert_mid, 2, 0);
    // a1.print();
    // a1.transfer();
    // a1.transfer();
    // a1.print();
    // cout << a1.compute_three_point_fct(2) << endl;
    // a1.clear();

    // Amplitude a2(&v[0], &v[1], 0, 2, 2);
    // cout << a2.compute_three_point_fct(2) << endl;
    // a2.clear();

    return 0;
}
#pragma endregion
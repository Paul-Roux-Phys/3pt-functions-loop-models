#include "three_point_fcts.hpp"
// #include "amplitudes.hpp"

#define O OPENING
#define C CLOSING
#define E EMPTY
int main()
{
    cout << "lambda" << "\t" << std::setw(20) << "n" << "\t" << "omega" << endl;
    for (lambda = 0.43; lambda < std::sqrt(1.5); lambda += 0.04)
    {
        set_weights(lambda);
        int half_size = 20 * lattice_size;
        LargeFloat c123 = compute_three_point_function_current_2leg_2leg(half_size);
        vp.clear();
        LargeFloat c101 = compute_two_point_function_current(half_size);
        vp.clear();
        LargeFloat c220 = compute_three_point_function_2leg_2leg_id(half_size);
        vp.clear();
        LargeFloat c000 = compute_partition_function(half_size);
        vp.clear();

        // cout << "c123 = " << c123 << endl;
        // cout << "c101 = " << c101 << endl;
        // cout << "c220 = " << c220 << endl;
        // cout << "c000 = " << c000 << endl;

        cout << lambda << "\t" << std::setw(20) << std::setprecision(15) << n_loop << "\t"
             << std::setprecision(15) << c123 / c220 * sqrt(abs(c000 / c101)) << endl;
    }
    cout << endl;
    return 0;
}
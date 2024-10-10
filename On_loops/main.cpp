#include "two_point_fcts.hpp"
#include "amplitudes.hpp"

#define O OPENING
#define C CLOSING
#define E EMPTY
int main()
{
    Amplitude a(&v[0], &v[1], 2, 0, 2);
    a.compute_two_point_amplitude_current_with_phases();
    // compute_current_two_point_function(20*lattice_size);
    return 0;
}
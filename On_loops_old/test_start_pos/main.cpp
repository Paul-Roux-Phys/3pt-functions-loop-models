#include "three_point_fcts.hpp"
#include <unistd.h>
#include <stdio.h>

#define O OPENING
#define C CLOSING
#define E EMPTY

int main(int argc, char *argv[])
{
    int k1, rs1, k2, rs2, k3, rs3;
    int pos1, pos2, pos3, pos4;

    lambda = std::stod(argv[1]);
    k1 = std::stoi(argv[2]);
    k2 = std::stoi(argv[3]);
    k3 = std::stoi(argv[4]);
    pos1 = std::stoi(argv[5]);
    pos2 = std::stoi(argv[6]);
    pos3 = std::stoi(argv[7]);
    pos4 = std::stoi(argv[8]);
    set_weights(lambda);

    int half_size = 10 * lattice_size;
    cout << std::setprecision(15) << compute_3pt_function(half_size, k1, k2, k3, pos1, pos2, pos3, pos4)
         << endl;

    return 0;
}
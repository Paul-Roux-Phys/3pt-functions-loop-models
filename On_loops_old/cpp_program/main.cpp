#include "three_point_fcts.hpp"
#include <unistd.h>
#include <stdio.h>

#define O OPENING
#define C CLOSING
#define E EMPTY
int main(int argc, char* argv[])
{
    int flags, opt;
    int diagonals = 0;

    flags = 0;

    int k1, k2, k3, rs1, rs2, rs3, starting_position;
    V P1, P3;

    // while ((opt = getopt(argc, argv, "d:")) != -1)
    // {
    //     switch (opt)
    //     {
    //     case 'd':
    //         diagonals = atoi(optarg);
    //         break;
    //     default: /* '?' */
    //         lambda = std::stod(argv[1]);
    //         k1 = std::stoi(argv[2]);
    //         rs1 = std::stoi(argv[3]);
    //         k2 = std::stoi(argv[4]);
    //         rs2 = std::stoi(argv[5]);
    //         k3 = std::stoi(argv[6]);
    //         rs3 = std::stoi(argv[7]);
    //         set_weights(lambda);
    //     }
    // }

    lambda = std::stod(argv[1]);
    k1 = std::stoi(argv[2]);
    rs1 = std::stoi(argv[3]);
    k2 = std::stoi(argv[4]);
    rs2 = std::stoi(argv[5]);
    k3 = std::stoi(argv[6]);
    rs3 = std::stoi(argv[7]);
    set_weights(lambda);

    // switch (diagonals)
    // {
    // case 0:
    //     cout << "no diagonal fields" << endl;
    //     break;
    // case 1:
    //     cout << "diagonal at position 1" << endl;
    //     break;
    // case 3:
    //     cout << "diagonal at position 3" << endl;
    //     break;
    // case 13:
    //     cout << "diagonal at positions 1 and 3" << endl;
    //     break;
    // default:
    //     cout << "error in -d options" << endl;
    //     return 1;
    // }

    // cout << "argv" << argv << endl;
    // cout << "optind = " << optind << endl;
    // cout << "argc = " << argc << endl;

    int half_size = 10 * lattice_size;
    cout << std::setprecision(15) << compute_3pt_function(half_size, k1, k2, k3, rs1, rs2, rs3, 0).real()
         << endl;

    return 0;
}
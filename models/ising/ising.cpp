#include "transfer.hpp"

K init_key() {
    K k;
    std::fill(k.begin(), k.end(), 1);
    return k;
}

// std::string print_eigvals(VP* v, int nev, const char* which ="LM", bool no_mult=false, int nev_to_print=0) {
//     Matrix<Weight, VP> M(v, nev, which);

//     M.find_eigenvalues(no_mult);
//     return M.sprint_eigenvalues(nev_to_print);
// }

LargeFloat two_point_amplitude(int M, Vec& init, Vec& end) {
    V res=0;

    p1->clear();
    *p1 += init;

    cout << "init: "; init.print();

    for (int i = 0; i < M; i++) {
        transfer();
        p1->factorise_norm();
        // v.mul<>(S); // make sure the state remains exactly symmetric even with numerical errors
    }

    res = p1->inner_product(end);

    return std::make_pair(res*p1->get_norm().first, p1->get_norm().second);
}

int main(int argc, char* argv[]) {
    
    Vec init(2);
    BV b1(init_key(), 1);
    BV b2(init_key(), -1);
    key_flip(b2.key);
    init += b1;
    init += b2;

    LargeFloat two_point = two_point_amplitude(1000, init, init);
    cout << two_point;


    // if (argc < 2) {
    //     cout << "Expected argument: number of eigenvalues to compute";
    //     return 1;
    // }
    // int nev_to_print = atoi(argv[1]);
    
    // v += BV(init_key(), 1);

    // int nev = 20;
    // const char* which = "LM";
    // bool no_mult=true;

    // cout << print_eigvals(&v, nev, which, no_mult, nev_to_print);


    return 0;
}
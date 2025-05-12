#include "transfer.hpp"

int LinkPattern::arch_end (int i) {
    int site = key[i];
    if (site == EMPTY || is_defect (i)) {
        return i;
    }
    int direction = (site == OPEN) ? 1 : -1;
    int inner_arches = 1;
    int j = i;
    int limit = 0;
    do {
        if (limit++ > 1000) {
            std::cout << "Couldn't find end of arch at position " << i
                      << " for key  ";
            key.print ();
            std::cout << "  of size " << SIZE << std::endl;
            std::cout << "Check that the key is a valid link pattern"
                      << std::endl;
            std::exit (EXIT_FAILURE);
        }
        j = (j + direction + SIZE) % (SIZE);
        if (key[j] == site) {  // we are meeting a new arch of the same type
            inner_arches += 1;
        } else if (key[j] ==
                   site + direction) {  // an arch of the same type is ending
            inner_arches -= 1;
        }
    } while (inner_arches > 0);
    return j;
}

void LinkPattern::contract_arch_defect (int i) {
    int ip1 = (i + 1) % SIZE;
    int arch = is_defect (i) ? ip1 : i;
    int def = is_defect (i) ? i : ip1;
    key.set (arch_end (arch), key[def]);
}

void LinkPattern::contract_arches (int i) {
    int ip1 = (i + 1) % SIZE;
    int site1 = key[i], site2 = key[ip1];
    if ((site1 != OPEN && site1 != CLOSE) ||
        (site2 != OPEN && site2 != CLOSE)) {
        return;
    }
    // Contract arches at positions i and i+1
    int end1 = arch_end (i), end2 = arch_end (ip1);
    if (site1 == OPEN && site2 == OPEN) {
        key.set (end2, OPEN);
    } else if (site1 == CLOSE && site2 == CLOSE) {
        key.set (end1, CLOSE);
    }
}

bool LinkPattern::can_contract_defects (int i) {
    int ip1 = (i + 1) % SIZE;
    return (!is_ref (i) && !is_ref (ip1) &&
            ((is_bottom (i) && is_middle (ip1)) ||
             (is_bottom (ip1) && is_middle (i))) &&
            get_nb_defects () - 2 >= min_defects);
}

bool LinkPattern::can_contract_sites (int i) {
    int ip1 = (i + 1) % SIZE;
    return (can_contract_defects (i) || !is_defect (i) || !is_defect (ip1));
}

void LinkPattern::contract (int i) {
    int ip1 = (i + 1) % SIZE;
    if (key[i] == OPEN && key[ip1] == CLOSE) {  // closing a contractible loop
        value *= n_loop;
    } else if (key[i] == CLOSE && key[ip1] == OPEN && arch_end (i) == ip1) {
        // closing a loop winding around the cylinder
        value *= n_ncloop;
    }
    if (is_defect (i) != is_defect ((i + 1) % SIZE)) {
        // XOR: exactly one of the two sites is a defect
        contract_arch_defect (i);
    } else if (can_contract_defects (i)) {
        // nothing to do here
    } else {
        contract_arches (i);
    }
}

void LinkPattern::move_strand (int i) {
    int empty = (key[i] == EMPTY) ? i : (i + 1) % SIZE;
    int def = (key[i] == EMPTY) ? (i + 1) % SIZE : i;
    key.set (empty, key[def]);
    key.set (def, EMPTY);
}

void LinkPattern::put_arch (int i) {
    // insert an arch between sites i and i+1;
    key.set (i, OPEN);
    key.set ((i + 1) % SIZE, CLOSE);
}

void LinkPattern::r_matrix_empty_empty (OnState& v, int i) {
    LinkPattern copy = *this;

    value *= w_empty;
    v += *this;

    copy.put_arch (i);
    copy.value *= w_turn;
    v += copy;
}

void LinkPattern::r_matrix_empty_occ (OnState& v, int i) {
    LinkPattern copy = *this;

    value *= w_turn;
    v += *this;

    copy.move_strand (i);
    copy.value *= w_straight;
    v += copy;
}

void LinkPattern::r_matrix_occ_occ (OnState& v, int i) {
    LinkPattern copy = *this;

    value *= w_full;
    v += *this;

    if (can_contract_sites (i)) {
        copy.contract (i);
        LinkPattern copy2 = copy;
        copy.key.set (i, EMPTY);
        copy.key.set ((i + 1) % SIZE, EMPTY);
        copy.value *= w_turn;
        v += copy;

        copy2.put_arch (i);
        copy2.value *= w_full;
        v += copy2;
    }
}

void LinkPattern::r_matrix (OnState& v, int i) {
    int ip1 = (i + 1) % SIZE;
    if (key[i] == EMPTY && key[ip1] == EMPTY) {
        r_matrix_empty_empty (v, i);
    } else if (key[i] == EMPTY || key[ip1] == EMPTY) {
        r_matrix_empty_occ (v, i);
    } else {
        r_matrix_occ_occ (v, i);
    }
}

void LinkPattern::insert_aux_space (OnState& v) {
    key.shift_right ();

    LinkPattern copy = *this;

    key.set (SIZE - 1, EMPTY);
    key.set (0, EMPTY);
    v += *this;

    copy.key.set (SIZE - 1, OPEN);
    copy.key.set (0, CLOSE);
    v += copy;
}

void LinkPattern::contract_aux_space (OnState& v) {
    int i1 = L, i2 = L + 1;
    if (key[i1] == EMPTY && key[i2] == EMPTY) {
        v += *this;
    } else if (can_contract_sites (i1) && is_occupied (i1) &&
               is_occupied (i2)) {
        // both sites are occupied, at least one of them isn't a defect
        contract (i1);
        key.set (i1, EMPTY);
        key.set (i2, EMPTY);
        v += *this;
    }
}

// pos must be < max(max_down_def, max_up_def)
// e.g. if there are 1 up and 1 down defects then pos = 0
// this function handles the insertion of the middle operator at position pos
void LinkPattern::middle_op_ith_site (int nb_down_def, int nb_up_def, int pos) {
    if (is_bottom (pos)) {
        if (pos >= nb_down_def) {
            key.set (0, ERROR);
            return;
        }
        if (is_ref (pos)) {
            key.set (0, ERROR);
            return;
        }
        if (nb_up_def <= pos && pos < nb_down_def)
            key.set (pos, EMPTY);
        else
            key.set (pos, MID);
    } else if (key[pos] == EMPTY) {
        if (nb_down_def <= pos && pos < nb_up_def)
            key.set (pos, MID);
        else {
            key.set (0, ERROR);  // error
            return;
        }
    } else if (!is_defect (pos)) {
        int end_pos = arch_end (pos);
        if (pos >= nb_down_def) {
            key.set (0, ERROR);  // error
            return;
        }
        if (end_pos < nb_down_def || end_pos < nb_up_def) {
            key.set (0, ERROR);  // error
            return;
        }
        if (pos >= nb_up_def)
            key.set (pos, EMPTY);
        else
            key.set (pos, MID);
        key.set (end_pos, MID);
    }
}

void LinkPattern::insert_mid_op (OnState& v) {
    if (k2 % 2 == 0) {
        for (int i = 0; i < k2 / 2; i++) middle_op_ith_site (k2 / 2, k2 / 2, i);
        if (!(key[0] == ERROR)) v += *this;
    } else {
        LinkPattern copy = *this;
        for (int i = 0; i < k2 / 2 + 1; i++) {
            // symmetrise up and down
            this->middle_op_ith_site (k2 / 2 + 1, k2 / 2, i);
            copy.middle_op_ith_site (k2 / 2, k2 / 2 + 1, i);
        }
        if (!(key[0] == ERROR)) v += *this;
    }
}

void LinkPattern::uncolor_defects (OnState& v) {
    for (int i = 0; i < L; i++)
        if (key[i] == MID) key.set (i, BOT);
    v += *this;
}

LinkPattern defects (std::vector<int> positions, Weight phase,
                     std::vector<int> labels) {
    LinkPattern p;
    for (int i = 0; i < positions.size (); i++) {
        p.key.set (positions[i] % L, labels[i]);
    }
    p.value = phase;
    return p;
}

OnState psi (int k) {
    OnState res;
    LinkPattern p (0, 1);
    for (int i = 0; i < k; i++) {
        p.key.set (i, BOT);
    }
    res += p;
    return res;
}

OnState bottom_state () { return psi (k1); }

OnState top_state () { return psi (k3); }

Weight compute_3pt (size_t half_size) {
    OnState s, s2;  // s2 is used as an intermediary variable
    s += bottom_state ();

    for (size_t i = 0; i < half_size; i++) {
        s.transfer (s2);
    }

    s.insert_mid_op (s2);

    for (size_t i = 0; i < half_size; i++) {
        s.transfer (s2);
    }

    s.uncolor_defects (s2);
    return s.inner_product (top_state ());
}

int main (int argc, char* argv[]) {
    const char usage[] = "Usage: %s L lambda k1 k2 k3\n";
    int opt;
    // parse command-line options
    while ((opt = getopt (argc, argv, "h")) != -1) {
        switch (opt) {
            case 'h':
                printf (usage, argv[0]);
                return 0;
            case '?':
                printf (usage, argv[0]);
                break;
        }
    }

    // parse other arguments
    if (argc - optind != 5) {
        printf (usage, argv[0]);
        return 1;
    }
    L = std::atoi (argv[optind++]);
    lambda = Weight (std::atof (argv[optind++]));
    k1 = std::atoi (argv[optind++]);
    rs1 = 0;
    // rs1 = std::atoi (argv[optind++]);
    k2 = std::atoi (argv[optind++]);
    rs2 = 0;
    // rs2 = std::atoi (argv[optind++]);
    k3 = std::atoi (argv[optind++]);
    rs3 = 0;
    // rs3 = std::atoi (argv[optind++]);

    set_weights (lambda);

    std::cout << compute_3pt (10 * L).real () << std::endl;

    return 0;
}

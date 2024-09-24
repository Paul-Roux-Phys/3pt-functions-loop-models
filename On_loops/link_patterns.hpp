/* We encode link states as sequences of 0 for arches
openings 1 for closings, and numbers >= 2 for defects.
Such a sequence can be encoded on 64 bits for lattice
sizes of up to 16 if there are at most 14 defects   */

#pragma once

#include <cstdint>
#include <iostream>
#include <iomanip>
#include <initializer_list>

enum SITE {
    EMPTY,
    OPENING,
    CLOSING,
    DEFECT
};

template<std::size_t size=4>
class FKKey : public key_64_bit_t<size> {
public:
    using key_64_bit_t<size>::key_64_bit_t; // use base class constructors
    using key_64_bit_t<size>::set;

    int arch_end(int i) {
        int site = (*this)[i];
        if (site == EMPTY || site == DEFECT) {
            return i;
        }
        int direction = (site == OPENING) ? 1 : -1;
        int inner_arches = 1;
        int j = i;
        int limit = 0;
        do 
        {
            if (limit++ > 1000)
            {
                std::cerr << "Couldn't find end of arch at position " << i << " for key  " << *this
                          << "  of size " << size << endl;
                std::cerr << "Check that the key is a valid link pattern" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            j = (j+direction+size)%(size);
            if ((*this)[j] == site)
            { // we are meeting a new arch of the same type
                inner_arches += 1;
            }
            else if ((*this)[j] == site+direction)
            { // an arch of the same type is ending
                inner_arches -= 1;
            }
        }
        while (inner_arches > 0) ;
        return j;
    }

    int arch_opening(int i) {
        int end = arch_end(i);
        return (end == OPENING) ? end : i;
    }

    int arch_closing(int i) {
        int end = arch_end(i);
        return (end == CLOSING) ? end : i;
    }

    void contract_arches(int i) {
        int ip1 = (i+1)%size;
        int site1 = (*this)[i], site2 = (*this)[ip1];
        if ((site1 != OPENING && site1 != CLOSING) \
            || (site2 != OPENING && site2 != CLOSING))
        {
            return;
        }
        // Contract arches at positions i and i+1
        int end1 = arch_end(i), end2 = arch_end(ip1);
        if (site1 == OPENING && site2 == OPENING)
        {
            set(end2, OPENING);
        }
        else if (site1 == CLOSING && site2 == CLOSING)
        {
            set(end1, CLOSING);
        }
    }

    void contract_arch_defect(int i) {
        // Contract an arch and a defect between positions i and i+1
        int ip1 = (i+1)%size;
        int arch = ((*this)[i] == DEFECT) ? ip1 : i;
        set(arch_end(arch), DEFECT);
    }

    void contract(int i) {
        // Contract occupied sites i and i+1
        if (((*this)[i] == DEFECT) != ((*this)[(i+1)%size] == DEFECT)) // XOR: one of the two sites is a defect
        {
            contract_arch_defect(i);
        } else
        {
            contract_arches(i);
        }
    }

    void move_strand(int i) {
        int empty = ((*this)[i] == EMPTY) ? i : (i+1)%size;
        int def   = ((*this)[i] == EMPTY) ? (i+1)%size : i;
        set(empty, (*this)[def]);
        set(def, EMPTY);
    }

    void put_arch(int i) {
        // insert an arch between sites i and i+1;
        set(i, OPENING);
        set((i+1)%size, CLOSING);
    }
};
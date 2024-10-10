/* We encode link states as sequences of 1 for arches
openings, 2 for closings, and numbers >= 3 for defects.
If this can be encoded on 3 bits (i.e. 6 types of
defects at most), such a sequence can be encoded on 64
bits for lattice sizes of up to 21 */

#pragma once

#include <cstdint>
#include <iostream>
#include <iomanip>
#include <initializer_list>
#define XOR !=
enum SITE {
    OPENING,
    CLOSING,
    BOTTOM,
    MIDDLE=5
};
#define DEFECT BOTTOM

template<std::size_t size=4>
class FKKey : public key_64_bit_t<size> {
public:
    using key_64_bit_t<size>::key_64_bit_t; // use base class constructors
    using key_64_bit_t<size>::set;

    FKKey() {
        for (int i = 0; i < size; i++)
        {
            if (i%2 == 0)
                set(i, OPENING);
            else
                set(i, CLOSING);
        }
    }

    int arch_end(int i) {
        int site = (*this)[i];
        if (site >= DEFECT) {
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
        int arch = ((*this)[i] >= DEFECT) ? ip1 : i;
        int defect = ((*this)[i] >= DEFECT) ? i : ip1;
        set(arch_end(arch), (*this)[defect]);
    }

    int nb_defects() {
        int res = 0;
        for (int i = 0; i < size; i++)
        {
            if ((*this)[i] >= DEFECT)
                res += 1;
        }
        return res;
    }

    bool is_bottom(int i) {
        return (BOTTOM <= (*this)[i] && (*this)[i] < MIDDLE);
    }

    bool is_middle(int i) {
        return MIDDLE <= (*this)[i];
    }

    bool can_contract_defects(int i, int min_defects) {
        int ip1 = (i+1)%size;
        return (((is_bottom(i) && is_middle(ip1)) \
            || (is_bottom(ip1) && is_middle(i))) \
            && nb_defects()-2 >= min_defects);
    }

    bool can_contract_sites(int i, int min_defects) {
        // sites can be contracted if one of them isn't a defect or
        // if the two defects are allowed to join
        return (can_contract_defects(i, min_defects) || \
                !((*this)[i] >= DEFECT && (*this)[(i+1)%size] >= DEFECT));
    }

    void tl_generator(int i, int min_defects) {
        // Contract sites i and i+1
        if (can_contract_defects(i, min_defects))
        { 
            // nothing to do here
        }
        else if (((*this)[i] >= DEFECT) != ((*this)[(i+1)%size] >= DEFECT)) // XOR
        {
            contract_arch_defect(i);
        }
        else
        {
            contract_arches(i);
        }

        // insert an arch between sites i and i+1;
        set(i, OPENING);
        set((i+1)%size, CLOSING);
    }
};
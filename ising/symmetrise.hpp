#pragma once

#include "transfer.hpp"

void key_flip(K& k) {
    for (int i=0; i < size; i++) {
        k[i] *= -1;
    }
}

// we can represent a pair ((state, weight), (flipped_state, -weight))
// by only keeping the state with a +1 in first position.
void project_antisym(BV& b) {
    if (b.key[0] == -1) {
        key_flip(b.key);
        b.value *= -1;
    }
}

void project_sym(BV& b) {
    if (b.key[0] == -1) {
        key_flip(b.key);
    }
}

void symmetrise(Vec* v, BV b) {
    if (b.key[0] == 1) {
        key_flip(b.key);
        (*v)[b.key] = b.value;
    }
}

void antisymmetrise(Vec* v, BV b) {
    if (b.key[0] == 1) {
        key_flip(b.key);
        (*v)[b.key] = -b.value;
    }
}
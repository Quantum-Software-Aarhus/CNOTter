#ifndef TREE_H
#define TREE_H

#include <cstdint>
#include "hashset.h"
#include "matrixN.h"

using rootset = HashSet<uint64_t, Linear, MurmurHash>;
using pairset = HashSet<uint32_t, Linear, MurmurHash>;
using leafset = HashSet<uint32_t, Linear, MurmurHash>;

pairset intermediate;
leafset leaves;

// We need to fold NR leaves (uint64_t) into one uint64_t root
// We build a tree as follows:
// - The leaves are stored in idxs [0,NR)
// - The intermediate nodes are stored in idxs [NR,2*NR-2)
// - The root is not stored

// Test if vetctor x is in the table.
// Early termination if a subvector of x is unknown. 
bool CONTAINS(Matrix x, rootset &table) {
    uint32_t idxs[2*NR-2];

    // Check if leaf indices are known and store them
    for (uint8_t i=0; i<NR; i++) {
        idxs[i] = leaves.contains(x._bits[i]);
        if (!(idxs[i])) return false;
    }

    // Check if intermediate indices are known and store them
    // node NR+i has children 2*i and 2*i+1
    for (uint8_t i=0; i<NR-2; i++) {
        uint64_t pair = (uint64_t)idxs[i*2]<<32 | idxs[i*2+1];
        idxs[NR+i] = intermediate.contains(pair);
        if (!(idxs[NR+i])) return false;
    }
    uint64_t root = (uint64_t)idxs[2*NR-4]<<32 | idxs[2*NR-3];
    return table.contains(root);
}

bool INSERT(Matrix x, rootset &table) {
    uint32_t idxs[2*NR-2];

    // insert the leaves
    for (uint8_t i=0; i<NR; i++)
        idxs[i] = leaves.findOrPut(x._bits[i]);

    // insert intermediate nodes
    for (uint8_t i=NR; i<2*NR-2; i++) {
        uint64_t pair = (uint64_t)idxs[i*2]<<32 | idxs[i*2+1];
        idxs[NR+i] = intermediate.findOrPut(pair);
    }

    mat_idx root = (uint64_t)idxs[2*NR-4]<<32 | idxs[2*NR-3];
    
    // insert root
    return table.insert(root);
}

Matrix GET(uint64_t root) {
    uint32_t idxs[2*NR-2];
    idxs[2*NR-3] = root;
    idxs[2*NR-4] = root >> 32;
    for (uint8_t i=NR-3; i>=0; i--) {
        uint64_t next = intermediate.get(idxs[NR+i]);
        idxs[2*i+1] = intermediate.get(next);
        idxs[2*i] = intermediate.get(next >> 32);
    }
    Matrix result;
    for (uint8_t i=0; i<NR; i++) {
        result._bits[i] = leaves.get(idxs[i]);
    }
    return result;
}

#endif
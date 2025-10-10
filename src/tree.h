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
constexpr int PairSize=32;

// We need to fold NR leaves (uint64_t) into one uint64_t root
// We build a tree as follows:
// - The leaves are stored in idxs [0,NR)
// - The intermediate nodes are stored in idxs [NR,2*NR-2)
// - The root is not stored

// Test if vetctor x is in the table.
// Early termination if a subvector of x is unknown. 
bool CONTAINS(Matrix x, rootset &table) {
    uint64_t root;
    if (NR==1) { // no need for compression
        root = x._bits[0];
    }
    else {
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
        root = (uint64_t)idxs[2*NR-4]<<32 | idxs[2*NR-3];
    }
    return table.contains(root);
}

Matrix GET(mat_idx root) {
    Matrix result;
    if (NR==1) { // no need for compression
        result._bits[0] = root;
        return result;
    }
    else {
        uint32_t idxs[2*NR-2];
        idxs[2*NR-3] = root;
        idxs[2*NR-4] = root >> 32;
        for (uint8_t i=NR-3; i < NR; i--) { // terminate when "negative"
            uint64_t next = intermediate.get(idxs[NR+i]);
            idxs[2*i+1] = next;
            idxs[2*i] = next >> 32;
        }
        for (uint8_t i=0; i<NR; i++) {
            result._bits[i] = leaves.get(idxs[i]);
        }
    }
    return result;
}

bool INSERT(Matrix x, rootset &table) {
    mat_idx root;

    if (NR==1) { // no need for compression
        root = x._bits[0];
    }
    else {
        uint32_t idxs[2*NR-2];
        // insert the leaves
        for (uint8_t i=0; i<NR; i++)
            idxs[i] = leaves.findOrPut(x._bits[i]);
        // insert intermediate nodes
        for (uint8_t i=0; i<NR-2; i++) {
            uint64_t pair = (uint64_t)idxs[i*2]<<32 | idxs[i*2+1];
            idxs[NR+i] = intermediate.findOrPut(pair);
        }
        root = (uint64_t)idxs[2*NR-4]<<32 | idxs[2*NR-3];
    }
    return table.insert(root);
}

#endif
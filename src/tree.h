#include <cstdint>
#include "hashset.h"
#include "matrixN.h"

using rootset = HashSet<uint64_t, Linear, MurmurHash>;
using pairset = HashSet<uint32_t, Linear, MurmurHash>;
using leafset = HashSet<uint32_t, Linear, MurmurHash>;

pairset intermediate;
leafset leaves;

// We need to fold NR leaves into one root
bool CONTAINS(Matrix x, rootset &table) {
    uint32_t leave_idx[2*NR-1];

    for (uint8_t i=0; i<NR; i++)
        if (!(leave_idx[i] = leaves.contains(x._bits[i]))) 
            return false;
    for (uint8_t i=NR, j=0; i<2*NR-2; i++, j+=2)
        if (!(leave_idx[i] = intermediate.contains((uint64_t)x._bits[j]<<32 | x._bits[j+1]))) 
            return false;

    uint64_t root = (uint64_t)leave_idx[2*NR-2]<<32 | leave_idx[2*NR-1];
    return table.contains(root);
}

bool INSERT(Matrix x, rootset &table) {
    uint32_t leave_idx[2*NR-1];

    for (uint8_t i=0; i<NR; i++)
        leave_idx[i] = leaves.findOrPut(x._bits[i]);

    for (uint8_t i=NR, j=0; i<2*NR-2; i++, j+=2)
        leave_idx[i] = intermediate.findOrPut((uint64_t)x._bits[j]<<32 | x._bits[j+1]);

    return table.insert((uint64_t)leave_idx[2*NR-2]<<32 | leave_idx[2*NR-1]);
}

Matrix GET(uint64_t root) {
    uint32_t leave_idx[2*NR-1];
    leave_idx[2*NR-1] = root >> 32;
    leave_idx[2*NR-2] = root;
    for (uint8_t i=2*NR-1, j=2*NR-3; i>=NR; i--, j-=2) {
        uint64_t next = intermediate.get(leave_idx[i]);
        leave_idx[j] = intermediate.get(next >> 32);
        leave_idx[j-1] = intermediate.get(next);
    }
    Matrix result;
    for (uint8_t i=0; i<NR; i++) {
        result._bits[i] = leaves.get(leave_idx[i]);
    }
    return result;
}

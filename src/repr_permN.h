#ifndef REPR_H
#define REPR_H

#include <algorithm>
#include "matrixN.h"

typedef std::array<byte,4> finger_t;
    // Fingerprint of an index i:
    // 0: NEGATION of diagonal [i][i]
    // 1: sum of row [i]
    // 2: sum of column [i]
    // 3: store diagonal i  // this is only used locally to sort the matrix

inline void fingerprint(const Matrix &x, finger_t finger[N]) {
    for (byte i=0; i<N; i++) {
        finger[i][0] = 1;
        finger[i][1] = 0;
        finger[i][2] = 0;
        finger[i][3] = i;
    }
    for (byte i=0; i<N; i++)
        for (byte j=0; j<N; j++) // traverse all bits in x
            if (x.get(i,j)) {
                if (i == j)     
                    finger[i][0] = 0; // 1 on diagonal
                else {
                    finger[i][1]++; // count 1 on row i
                    finger[j][2]++; // count 1 on col j
}           }   }

// compute and sort the finger-print, return the normalized matrix and permutation
inline Matrix normalize(const Matrix &x, finger_t finger[N], perm pi) {
    fingerprint(x, finger);
    std::sort(finger, finger+N);
    for (byte i=0; i<N; i++) pi[i] = finger[i][3];
    return x.permute(pi);
}

// generate next permutation of list[] as specified by cycles[]
// cycles[] is 0-terminated, and sum(cycles) = len(list)
inline bool next_cycle_perm(const byte cycles[], byte list[]) {
    if (cycles[0]==0) return false; 
    return std::next_permutation(list, list + cycles[0])
            || next_cycle_perm(cycles+1, list + cycles[0]);
}

// This function normalizes y, initializes cycles, and returns the first essential index
inline byte compute_cycles(byte cycles[], Matrix &y) {
    finger_t finger[N]; // finger print
    perm pi;
    y = normalize(y, finger, pi); // this also initializes fingerprint and pi

    // skip and count inessential indices
    byte i = 0;
    while (i<N && !(finger[i][0] || finger[i][1] || finger[i][2])) 
        i++;
    byte essential = i;

    // identify the cycles that must be permuted
    byte c = 0; // index in cycles
    while (i<N) {
        // count constant fragment with the same fingerprint
        byte j = i+1;
        while (j<N  && finger[i][0] == finger[j][0]
                    && finger[i][1] == finger[j][1]
                    && finger[i][2] == finger[j][2])
            j++;
        cycles[c++] = j-i;
        i = j;
    }
    cycles[c] = 0; //terminator
    return essential;
}

// Assume y is normalized
// Update y to the smallest representative
// Return the number of "essential" stabilizers
inline counter explore_orbit(Matrix &y, byte cycles[], byte essential) {
    perm pi; id_perm(pi);

    // traverse all nested permutations (orbit) and count stabilizers
    counter stabilizers = 1; // the id is surely a stabilizer
    Matrix smallest = y;      // maintain the smallest matrix in the orbit
    while (next_cycle_perm(cycles, pi + essential)) { // we counted id already
        Matrix z = y.permute(pi);
        if (z == y)
            stabilizers++;
        else if (z < smallest)
            smallest = z;
    }

    y = smallest; // set y to the smallest representative
    return stabilizers;
}

inline counter representative(Matrix &y) {
    byte cycles[N+1];  // cycles for permutation
    byte essential = compute_cycles(cycles, y); // now y is normalized
    counter stab = explore_orbit(y, cycles, essential);
    return stab * fac[essential];
}

// return the permutation from x to its representative
void representativePerm(const Matrix &x, perm pi) {
    finger_t finger[N]; // finger print
    perm pi1;                     // permutation from x to y
    Matrix y = normalize(x, finger, pi1);

    // Similar to representative
    perm pi2;                     // permutation from y to smallest
    for (byte i=0; i<N; i++) 
        pi[i] = pi2[i] = i;       // initialize "running" permutation pi
    byte cycles[N+1];             // cycles for permutation
    byte essential = compute_cycles(cycles, y);
    Matrix smallest = y;          // maintain the smallest matrix in the orbit
    while (next_cycle_perm(cycles, pi + essential)) { // we counted id already
        Matrix z = y.permute(pi);
        if (z < smallest) {
            smallest = z;
            for (byte i=0; i<N; i++) pi2[i]=pi[i];
        }
    }
    compose_perm(pi2, pi1, pi);
    assert(x.permute(pi) == smallest);
}

// Assume that x is normalized
// Compute the first essential index
inline byte countEss(const Matrix &x) {
    for (byte ess=0; ess<N; ess++)
        for (byte i=0; i<N; i++)
            if (i==ess) {
                if (!(x.get(i,i))) return ess;
            }
            else {
                if (x.get(ess,i) || x.get(i,ess)) return ess;
            }
    return N;
}

void pretty_finger(finger_t finger[N]) {
    for (byte i=0; i<N; i++)
        printf("(%u,%u,%u) ",finger[i][0], finger[i][1], finger[i][2]);
    printf("\n");
}

void pretty_cycles(byte essential, byte cycles[N+1]) {
    printf("Non-essential indices: %u. Cycles: ",essential);
    byte c=0;
    byte i=essential;
    while (cycles[c]) {
        printf("( ");
        for (byte k=0;k<cycles[c];k++)
            printf("%u ",i+k);
        i+=cycles[c];
        printf(") ");
        c++;
    }
    printf("\n");
}

void investigate(const Matrix &x) {
    printf("Original matrix:\n");
    x.print();
    finger_t finger[N];
    perm pi;
    byte cycles[N+1];
    fingerprint(x,finger);
    printf("Original Fingerprint: ");
    pretty_finger(finger);
    Matrix y = normalize(x,finger,pi);
    printf("Sorted Fingerprint  : ");
    pretty_finger(finger);
    printf("Permutation:\n");
    pretty_perm(pi);
    printf("Normalized matrix:\n");
    y.print();
    byte essential = compute_cycles(cycles, y);
    pretty_cycles(essential, cycles);
    counter stabilizer = explore_orbit(y, cycles, essential);
    counter orbits = fac[N]/(stabilizer * fac[essential]);
    printf("Minimized matrix:\n");
    y.print();
    representativePerm(x,pi);
    printf("Permutation:\n");
    pretty_perm(pi);
    assert(x.permute(pi) == y);
    printf("Represents %lu matrices.\n\n", orbits);
}

#endif
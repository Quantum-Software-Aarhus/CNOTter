#ifndef REPR_H
#define REPR_H

#include "options.h"
#include "nauty.h"
#include "nautinv.h"
#include "matrixN.h"

const byte m=1; // nauty
const byte n=N; // 
static DEFAULTOPTIONS_DIGRAPH(options);

void matrix2nauty(const Matrix &y, graph g[m*n]) { // now always the same g
    EMPTYGRAPH(g,m,n);
    for (byte i=0; i<N; i++) {
        div_t wordi = div(i,NC);
        uint64_t row = (y._bits[wordi.quot] >> (N * wordi.rem)) & ((1ULL << N) - 1);
        // Nauty uses the most-significant bits: 
        // - shift all bits to the left
        // - our row i is Nauty's row N-i-1
        g[N-i-1] = row << (WORDSIZE-N);
        // printf("row %u = %lu; g[%u]: %lu\n",i,row,i,g[i]);
    }
}

Matrix nauty2matrix(const graph* g) {
    Matrix y(0); // initialize
    for (byte i=0; i<N; i++) {
        // Nauty uses most-significant bits:
        // - shift bits to the right again
        // - Nauty's row N-i-1 is our row i
        div_t wordi = div(i,NC);
        uint64_t row = g[N-i-1] >> (WORDSIZE-N); // convert to 64-bit before shift-left
        y._bits[wordi.quot] |= (row << (N*wordi.rem));
    }
    return y;
}

inline uint64_t representative(Matrix &y) {
    graph g[m*n];
    graph h[m*n];
    int lab[n], ptn[n], orbits[n];
    statsblk stats;
    matrix2nauty(y,g);
    densenauty(g,lab,ptn,orbits,&options,&stats,m,n,h);
    uint64_t mysize = stats.grpsize;
    // if (stats.grpsize2) {printf("Large group!\n"); exit(-1); }
    y=nauty2matrix(h); // this value is returned
    return fac[N] / mysize;
}

// return the permutation from x to its representative
void representativePerm(const Matrix &x, perm pi) {
    graph g[m*n];
    graph h[m*n];
    int lab[n], ptn[n], orbits[n];
    statsblk stats;
    matrix2nauty(x,g);
    densenauty(g,lab,ptn,orbits,&options,&stats,m,n,h);
    for (byte i=0; i<N; i++)
        pi[N-1-i] = N-1-lab[i];       // revert N-1..0 to 0..N-100
    assert(x.permute(pi) == nauty2matrix(h));
}

void investigate(const Matrix &x) {
    Matrix y = x;
    printf("Original matrix:\n");
    y.print();
    perm pi;
    representativePerm(y,pi);
    printf("Permutation:\n");
    pretty_perm(pi);
    uint64_t orbit = representative(y);
    printf("Canonical matrix:\n");
    y.print();
    printf("Represents %lu matrices.\n\n",orbit);
}

#endif
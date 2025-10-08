#ifndef REPR_H
#define REPR_H

#include "options.h"
#include "nauty.h"
#include "nautinv.h"
#include "matrix.h"

const byte m=1; // nauty
const byte n=N; // 
static DEFAULTOPTIONS_DIGRAPH(options);

void matrix2nauty(const matrix &y, graph g[m*n]) { // now always the same g
    EMPTYGRAPH(g,m,n);
    for (byte i=0; i<N; i++) {
        matrix row = ((y >> N*i) & ((1<<N) - 1));
        // Nauty uses the most-significant bits: 
        // - shift all bits to the left
        // - our row i is Nauty's row N-i-1
        g[N-i-1] = row << (WORDSIZE-N);
        // printf("row %u = %lu; g[%u]: %lu\n",i,row,i,g[i]);
    }
}

matrix nauty2matrix(const graph* g) {
    matrix y=0LL;
    for (byte i=0; i<N; i++) {
        // Nauty uses most-significant bits:
        // - shift bits to the right again
        // - Nauty's row N-i-1 is our row i
        matrix row = g[N-i-1] >> (WORDSIZE-N); 
        // converted to 64-bit before shift-left
        y |= row << (N*i);
    }
    return y;
}

inline uint64_t representative(matrix &y) {
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
void representativePerm(matrix x, perm pi) {
    graph g[m*n];
    graph h[m*n];
    int lab[n], ptn[n], orbits[n];
    statsblk stats;
    matrix2nauty(x,g);
    densenauty(g,lab,ptn,orbits,&options,&stats,m,n,h);
    for (byte i=0; i<N; i++)
        pi[N-1-i] = N-1-lab[i];       // revert N-1..0 to 0..N-100
    assert(permute(x, pi) == nauty2matrix(h));
}

void investigate(matrix x) {
    printf("Original matrix:\n");
    pretty_matrix(x);
    perm pi;
    representativePerm(x,pi);
    printf("Permutation:\n");
    pretty_perm(pi);
    uint64_t orbit = representative(x);
    printf("Canonical matrix:\n");
    pretty_matrix(x);
    printf("Represents %lu matrices.\n\n",orbit);
}

#endif
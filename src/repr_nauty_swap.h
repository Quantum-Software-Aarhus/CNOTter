#ifndef REPR_H
#define REPR_H

#include "nauty.h"
#include "nautinv.h"
#include "matrix.h"

const byte m=1; // nauty
const byte n=2*N; // 
static DEFAULTOPTIONS_DIGRAPH(options);

void matrix2nauty(const matrix &y, graph g[m*n]) { // now always the same g
    EMPTYGRAPH(g,m,n);

    /* 
    We will construct the following bipartite graph:

    ( 0 0 )
    ( M 0 )

    Note: In Nauty the bits start on the left, so we shift WORDSIZE.
    */
    for (byte i=0; i<N; i++) {
        g[N-i-1] = 0;
        matrix M_row = ((y >> N*i) & ((1<<N) - 1));
        g[2*N-i-1] = M_row << (WORDSIZE-N);
    }
}

matrix nauty2matrix(const graph* g) {
    matrix y=0LL;

    /* Assumption: nauty will give back a matrix of the form
        (0  0)
        (M' 0)
        The assumption seems to hold for N=1..8.
        But it might fail, so we test that rows N..2N-1 are non-zero
    */

    for (byte i=N; i<2*N; i++) {
        if (g[i] == 0) {
            printf("\nProblem: Assumption on Nauty failed (maybe input was not full-rank?)\n");
            exit(-1);
        }
        matrix row = g[2*N-(i-N)-1] >> (WORDSIZE-N); // convert to 64-bit before shift-left
        y |= row << (N*(i-N));
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
    y=nauty2matrix(h); // this value is returned
    return fac[N] * fac[N] / mysize;
}

// return the permutation from x to its representative
void representativePerm2(matrix x, perm pi1, perm pi2) {
    graph g[m*n];
    graph h[m*n];
    int lab[n], ptn[n], orbits[n];
    statsblk stats;
    matrix2nauty(x,g);
    densenauty(g,lab,ptn,orbits,&options,&stats,m,n,h);
    matrix y=nauty2matrix(h);
    for (byte i=0; i<N; i++) {
        pi2[N-1-i] = N-1-lab[i];       // revert N-1..0 to 0..N-100
        pi1[N-1-i] = 2*N-1-lab[N+i]; // revert N-1..0 to 0..N-100
    }
    assert(permute2(x, pi1, pi2) == y);
}

// assuming m1 and m2 are equivalent, find sig, tau such that (sig,tau) . m1 = m2
void equiv_perm2(matrix m1, matrix m2, perm sig, perm tau) {
    perm sig1, tau1, sig2, tau2;
    representativePerm2(m1, sig1, tau1);   // repr = (sig1,tau1) . m1
    representativePerm2(m2, sig2, tau2);   // repr = (sig2,tau2) . m2
    compose_inv_perm(sig2, sig1, sig);
    compose_inv_perm(tau2, tau1, tau);
    assert(permute2(m1, sig1, tau1) == permute2(m2, sig2, tau2)); // both are repr
    assert(permute2(m1, sig, tau) == m2); 
}

void investigate(matrix x) {
    printf("Original matrix:\n");
    pretty_matrix(x);
    perm pi1, pi2;
    representativePerm2(x, pi1, pi2);
    printf("rows:\n"); pretty_perm(pi1);
    printf("cols:\n"); pretty_perm(pi2);
    uint64_t stabilizers = representative(x);
    printf("Canonical matrix:\n");
    pretty_matrix(x);
    printf("Represents %lu matrices.\n\n",stabilizers);
}

#endif
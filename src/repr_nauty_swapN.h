#ifndef REPR_H
#define REPR_H

#include "options.h"
#include "nauty.h"
#include "nautinv.h"
#include "matrixN.h"

const byte m=1; // nauty
const byte n=2*N; // 
static DEFAULTOPTIONS_DIGRAPH(options);

// NOTE: The current encoding only supports N up to 32, since we encode graphs in 64-bit words.
// For larger N, we need to increase m.

void matrix2nauty(const Matrix &y, graph g[m*n]) { // now always the same g
    EMPTYGRAPH(g,m,n);

    /* 
    We will construct the following bipartite graph:

    ( 0 0 )
    ( M 0 )

    Note: In Nauty the bits start on the left, so we shift WORDSIZE.
    */
    for (byte i=0; i<N; i++) {
        g[N-i-1] = 0;
        div_t wordi = div(i,NC);
        uint64_t row = ((y._bits[wordi.quot] >> N*wordi.rem) & ((1ULL << N) - 1));
        g[2*N-i-1] = row << (WORDSIZE-N);
    }
}

Matrix nauty2matrix(const graph* g) {
    Matrix y(0);

    /* Assumption: nauty will give back a matrix of the form
        (0  0)
        (M' 0)
        The assumption seems to hold for N=1..8.
        But it might fail, so we test that rows N..2N-1 are non-zero
    */

    for (byte i=0; i<N; i++) {
        if (g[i+N] == 0) {
            printf("\nProblem: Assumption on Nauty failed (maybe input was not full-rank?)\n");
            exit(-1);
        }
        div_t wordi = div(i,NC);
        uint64_t row = g[2*N-i-1] >> (WORDSIZE-N); // convert to 64-bit before shift-left
        y._bits[wordi.quot] |= row << (N*wordi.rem);
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
    y=nauty2matrix(h); // this value is returned
    return fac[N] * fac[N] / mysize;
}

// return the permutation from x to its representative
void representativePerm2(const Matrix &x, perm pi1, perm pi2) {
    graph g[m*n];
    graph h[m*n];
    int lab[n], ptn[n], orbits[n];
    statsblk stats;
    matrix2nauty(x,g);
    densenauty(g,lab,ptn,orbits,&options,&stats,m,n,h);
    for (byte i=0; i<N; i++) {
        pi2[N-1-i] = N-1-lab[i];       // revert N-1..0 to 0..N-100
        pi1[N-1-i] = 2*N-1-lab[N+i]; // revert N-1..0 to 0..N-100
    }
    Matrix y=nauty2matrix(h);
    assert(x.permute2(pi1, pi2) == y);
}

// assuming m1 and m2 are equivalent, find sig, tau such that (sig,tau) . m1 = m2
void equiv_perm2(const Matrix &m1, const Matrix &m2, perm sig, perm tau) {
    perm sig1, tau1, sig2, tau2;
    representativePerm2(m1, sig1, tau1);   // repr = (sig1,tau1) . m1
    representativePerm2(m2, sig2, tau2);   // repr = (sig2,tau2) . m2
    compose_inv_perm(sig2, sig1, sig);
    compose_inv_perm(tau2, tau1, tau);
    assert(m1.permute2(sig1, tau1) == m2.permute2(sig2, tau2)); // both are repr
    assert(m1.permute2(sig, tau) == m2); 
}

void investigate(Matrix &x) {
    printf("Original matrix:\n");
    x.print();
    perm pi1, pi2;
    representativePerm2(x, pi1, pi2);
    printf("rows:\n"); pretty_perm(pi1);
    printf("cols:\n"); pretty_perm(pi2);
    uint64_t stabilizers = representative(x);
    printf("Canonical matrix:\n");
    x.print();
    printf("Represents %lu matrices.\n\n",stabilizers);
}

#endif
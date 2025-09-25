#include "nauty.h"
#include "nautinv.h"

const byte m=1; // nauty
const byte n=2*N; // 
static DEFAULTOPTIONS_DIGRAPH(options);

void matrix2nauty(const matrix &y, graph g[m*n]) { // now always the same g
    EMPTYGRAPH(g,m,n);

    /* 
    We will construct the following bipartite graph:

    ( 0 M )
    ( 0 0 )

    Note: In Nauty the bits start on the left, so we shift WORDSIZE.
    */
    for (byte i=0; i<N; i++) {
        matrix M_row = ((y >> N*i) & ((1<<N) - 1));
        g[N-i-1] = M_row << (WORDSIZE-2*N);
    }
    for (byte i=0; i<N; i++) {
        g[2*N-i-1] = 0;
    }

}

matrix nauty2matrix(const graph* g, int lab[2*N]) {
    matrix y=0LL;
    byte rows[N];
    byte cols[N];
    byte r=0,c=0;

    /* Assumption: nauty will give back a matrix of the form
        (0  0)
        (M' 0)
        The assumption seems to hold for N=1..8.
        But it might fail, so we test that rows N..2N-1 are non-zero
    */

    for (byte i=N; i<2*N; i++) {
        if (g[i] == 0) {
            printf("assumption on Nauty failed\n");
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
    y=nauty2matrix(h,lab); // this value is returned
    return fac[N] * fac[N] / mysize;
}

// return the permutation from x to its representative
void representativePerm(matrix x, byte pi[N]) {
    printf("Permutation: Not yet adapted to swap-free\n");
    exit(-1);
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
    byte pi[N];
    // representativePerm(x,pi);
    printf("Permutation: skipped\n");
    // pretty_perm(pi);
    uint64_t stabilizers = representative(x);
    printf("Canonical matrix:\n");
    pretty_matrix(x);
    printf("Represents %lu matrices.\n\n",stabilizers);
}

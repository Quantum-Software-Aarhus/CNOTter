#ifndef OPTIONS_H
#define OPTIONS_H

#include <assert.h>
#include <iostream>
#include <fstream>
#include <cstdint>

using counter = uint64_t;

// Global definition of factorials up to 20!
const counter fac[] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 
    39916800, 479001600, 6227020800, 87178291200, 1307674368000, 20922789888000, 
    355687428096000, 6402373705728000, 121645100408832000, 2432902008176640000};

#ifndef N 
#define N 6 // Number of qubits, can be at most 8, set with -DN=7
#endif

#ifndef E 
#define E 1 // extra bits added to log of hash table size, set with -DE=2
#endif

#ifndef MAX
#define MAX 34 // maximum allocated table, set with -DMAX=36
#endif

#ifndef SWAP
#define SWAP 0 // SWAPS for free, enable with -DSWAP=1
#endif

#if SWAP==1
#define NAUTY 1  // SWAP requires Nauty
#endif

#ifndef NAUTY
#define NAUTY 0 // using Nauty instead of permutations, disable with -DNAUTY=0
#endif

#ifndef POLY
#define POLY 0 // compute polynomial coefficients, enable with -DPOLY=1
#endif

#ifndef BEAT
#define BEAT 0 // frequency of lifebeat in seconds (0 if no lifebeat), set with -DBEAT=60
#endif

#endif

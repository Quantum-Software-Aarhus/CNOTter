#ifndef OPTIONS_H
#define OPTIONS_H

#include <assert.h>
#include <iostream>
#include <fstream>

// Global definition of factorials up to 8!
const uint64_t fac[] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320};

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
#define POLY 0  // SWAP is incompatible with POLY
#endif

#ifndef NAUTY
#define NAUTY 1 // using Nauty instead of permutations, disable with -DNAUTY=0
#endif

#ifndef POLY
#define POLY 0 // compute polynomial coefficients, enable with -DPOLY=1
#endif

#ifndef BEAT
#define BEAT 0 // frequency of lifebeat in seconds (0 if no lifebeat), set with -DBEAT=60
#endif

#endif
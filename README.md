# Optimal CNOT synthesis

Computes the minimal CNOT circuit to implement a given parity matrix.

## Features:
- Based on BFS search, can enumerate all parity matrices on 1-8 qubits
- Efficient due to parallel hash-tables, symmetries, use of graph-isomorphism
- Can search for a goal-matrix with (bounded) bi-directional search
- Can handle the case where Swaps-are-for-free
- Can compute polynomial coefficients to predict group sizes, see [1]

## Usage and Options:
(see also first-time build below)

Usage: `sh matrix_cnot.sh [options] [goal]`: optimal CNOT synthesis  

Compile-time Options:
```
  -Q qubits  : number of Qubits (default 6)
  -S swap    : swaps-are-for-free, requires nauty (0 no, 1 yes) (default 0)
  -P poly    : compute polynomial coefficients (0 no, 1 yes) (default 0)
  -T threads : number of OpenMP threads to use (default: number of cores)
  -B beat    : heart-beat every BEAT seconds (default 0: no beat)
  -E extra   : size of hash-tables for levels + EXTRA bits (default 1)
  -M max     : max table size 2^MAX (default 34)
  -N nauty   : using Nauty (0 no, 1 yes) (default 1)
  -h         : this help
```
Run-time Options:
```
  -D dist    : maximum distance to try (default -1: unlimited)
  goal       : filename for Goal matrix (see 'Inputs/' for examples)
```

Note: all options are compile-time options, except for Dist and Goal.
Compiles a binary `./matrix_cnot<Q>.exe` (where Q=QUBITS) and runs it.
To rerun with the same Q on a different goal/distance, use 

    ./matrix_cnot<Q>.exe -<dist> <goal>

## Examples

Enumerate all matrices on 5 Qubits:

```
    sh matrix_cnot.sh -Q5
```
Total size: 9999360 (85411 orbits), completed at depth 12

---

Enumerate all matrices on 5 Qubits, where swaps are free:
```
    sh matrix_cnot.sh -Q5 -S1
```
Total size: 9999360 (885 orbits), completed at depth 8

---

Compute polynomial for distance 4 (needs 8 qubits):
```
    sh matrix_cnot.sh -Q8 -D4 -P1
```
Polynomial coefficients (8/4): [0, 0, 0, 60, 1818, 9990, 13200, 7560, 1680]  
Interpretation: The number of *different* CNOT circuits of depth 4 on $n$ qubits is precisely:

>$60\binom{n}{3} + 1818\binom{n}{4} + 9990\binom{n}{5} + 13200\binom{n}{6} + 7560\binom{n}{7} + 1680\binom{n}{8}$

---

Synthesize circuit for a long cyclic permutation on 6 Qubits with Bi-directional search:
```
    sh matrix_cnot.sh -Q6 Inputs/cycle6.txt
```
Found at distance 15 (9 + 6)

---

Synthesize circuit for the inverse matrix on (again) 6 Qubits:
```
    ./matrix_cnot6.exe Inputs/inverse6.txt
```
Found at distance 13 (7 + 6)  
**NOTE:** The binary could only be reused because other options didn't change.


## First-time build

This code depends on Nauty (see [README_NAUTY_MODIFIED.md](nauty/README_NAUTY_MODIFIED.md)). On the first-time build, type:

```
    sh install-nauty.sh
```

## Literature

The code is based on the following paper:

>Jens Emil Christensen, Søren F. Jørgensen, Andreas Pavlogiannis, Jaco van de Pol,  
**On Exact Sizes of Minimal CNOT Circuits.**  
In: *Reversible Computing*, RC’25, Odense, DK, 2025  
[[DOI]](https://doi.org/10.1007/978-3-031-97063-4_6) [[arXiv]](https://arxiv.org/abs/2503.01467)

# Nauty - modified subset

This is only a part of the original Nauty + Traces software version 2.9.1,
downloaded from https://pallini.di.uniroma1.it/nauty2_9_1.tar.gz.
See `LICENSE-2.0.txt` for the original license of Nauty + Traces.
Nauty + Traces is described in:

>McKay, B.D. and Piperno, A.,  
Practical Graph Isomorphism, II,  
Journal of Symbolic Computation, 60 (2014), pp. 94-112,  
<https://doi.org/10.1016/j.jsc.2013.09.003>


## Modifications made to Nauty

The following patch has already been applied to the source from `nauty2_9_1` (see `grpsize.patch`)
This modification is made by Jaco van de Pol, Aarhus University, Denmark, September 2025.

- The structure `double *grpsize1, int *grpsize2` is replaced everywhere by `uint64_t grpsize`
Reason: I needed exact counting for fairly large numbers.
**Warning:** This may not suit your purposes, since overflow can happen.

- Affected files:
    - gtnauty.c
    - nauty-h.in
    - nauty.c
    - schreier.c
    - schreier.h
    - traces.c
    - traces.h

##  Installation:

```
    ./configure --enable-tls
    make nautyW1.a
```
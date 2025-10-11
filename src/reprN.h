#ifndef REPRwrap_H
#define REPRwrap_H
#include "matrixN.h"

#if SWAP==1
#include "repr_nauty_swapN.h"
#elif NAUTY==1
#include "repr_nautyN.h"
#else
#include "repr_permN.h"
#endif


//TODO: move elsewhere (repr.h)

#if SWAP==0
// assuming m1 and m2 are equivalent, find pi such that pi . m1 = m2
void equiv_perm(const Matrix &m1, const Matrix &m2, perm pi) {
    perm pi1, pi2;
    representativePerm(m1, pi1);    // repr = pi1 . m1
    representativePerm(m2, pi2);    // repr = pi2 . m2
    compose_inv_perm(pi2, pi1, pi);
    assert(m1.permute(pi1) == m2.permute(pi2)); // both are repr
    assert(m1.permute(pi) == m2);  // pi . m1 = m2
}
#endif

#endif
#if SWAP==1
#include "repr_nauty_swap.h"
#elif NAUTY==1
#include "repr_nauty.h"
#else
#include "repr_perm.h"
#endif

#if SWAP==0
// assuming m1 and m2 are equivalent, find pi such that pi . m1 = m2
void equiv_perm(matrix m1, matrix m2, perm pi) {
    perm pi1, pi2;
    representativePerm(m1, pi1);    // repr = pi1 . m1
    representativePerm(m2, pi2);    // repr = pi2 . m2
    compose_inv_perm(pi2, pi1, pi);
    assert(permute(m1, pi1) == permute(m2, pi2)); // both are repr
    assert(permute(m1, pi) == m2);  // pi . m1 = m2
}
#else

#endif
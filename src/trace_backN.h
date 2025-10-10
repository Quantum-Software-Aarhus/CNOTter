#ifndef TRACE_BACK_H
#define TRACE_BACK_H

#include "matrixN.h"
#include "tree.h"
#include "repr_nautyN.h" // TODO: other reprs
#include <vector>

using trace = std::vector<std::pair<byte,byte>>;
using hashset = HashSet<uint64_t, Linear, MurmurHash>;

inline bool find_level(const Matrix &m, hashset &level) {
    Matrix n=m;
    representative(n);
    return CONTAINS(n, level);
}

Matrix step_back(const Matrix &x, hashset &level, trace &tr) {
    for (byte i=0; i<N; i++) {
        for (byte j=0; j<N; j++) // add to row j
            if (i != j) {
                Matrix prev = x.addrow(i,j);
                if (find_level(prev, level)) {
                    tr.push_back(std::pair<byte,byte>(i,j));
                    return prev;
    }   }       }
    assert(false && "Predecessor not found");
    exit(-1);
}

Matrix trace_back(const Matrix &goal, hashset levels[], int depth, trace &tr) {
    Matrix next = goal;
    for (int d=depth - 1; d>0; d--)
        next = step_back(next, levels[d], tr);
    return next;
}

// apply pi to all elements in trace tr
trace permute_trace(const perm pi, const trace &tr) {

    trace result;
    for (auto pair = tr.begin(); pair < tr.end(); pair++) 
        result.push_back(std::pair<byte,byte>(pi[pair->first], pi[pair->second]));
    
    return result;
}

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

trace trace_back_middle(Matrix &id, Matrix &middle, Matrix &goal, hashset bfs_fwd[], hashset bfs_bwd[], int fdepth, int bdepth, perm pi) {
    trace fwd_trace, bwd_trace, result;
    Matrix id_found, goal_found;
    id_found = trace_back(middle, bfs_fwd, fdepth, fwd_trace);
    goal_found = trace_back(middle, bfs_bwd, bdepth, bwd_trace);

    // REVERSE the forward trace and concatenate backward trace
    std::reverse(fwd_trace.begin(),fwd_trace.end());
    fwd_trace.insert(fwd_trace.end(), bwd_trace.begin(), bwd_trace.end()); 
    // NOTE: now trace runs from id_found to goal_found

#if SWAP==0
    assert(id==id_found);
    // Reconstruct the proper permutation (we want "goal" instead of "found")
    equiv_perm(goal, goal_found, pi);
    assert(goal.permute(pi) == goal_found);  // goal_found = pi . goal 
    // PERMUTE the whole trace to the true goal
    result = permute_trace(pi, fwd_trace);
    id_perm(pi); // return the identity permutation
#else
    // Reconstruct the proper permutations for the identity and goal
    perm pi1, pi2, pi3, pi4;
    equiv_perm2(id, id_found, pi1, pi2);
    equiv_perm2(goal_found, goal, pi3, pi4);
    assert(permute2(id, pi1, pi2) == id_found);      // id_found = (pi1,pi2) . id 
    assert(permute2(goal_found, pi3, pi4) == goal);  // goal = (pi3,pi4) . goal_found

    // construct q1 := pi1 ; pi2^-1 ; pi4^-1 and apply to the trace
    perm q0, q1, p4inv;
    inv_perm(pi4, p4inv);
    compose_inv_perm(pi2, p4inv, q0);
    compose_perm(pi1, q0, q1);
    result = permute_trace(q1, fwd_trace);

    // construct the final permutation pi := pi3 ; q1
    compose_perm(pi3, q1, pi);
#endif

    return result;
}

void print_trace(Matrix m, Matrix goal, const trace &tr, perm pi) {
    printf("\nOPENQASM 2.0;\n");
    printf("include \"qelib1.inc\";\n");
    printf("qreg q[%u];\n\n",N);
    for (std::pair<byte,byte> pair : tr) {
        byte i = pair.first, j = pair.second;
        printf("cx q[%u],q[%u];\n",i,j);
        m = m.addrow(i,j);
    }
    printf("\nResult of the circuit:\n");
    m.print();
#if SWAP==1
    printf("\nRow permutation:\n");
    pretty_perm(pi);
    printf("\nPermuted Result:\n");
    byte id[N];
    id_perm(id);
    m=permute2(m,pi,id);
    pretty_matrix(m);
#endif
    if (m==goal)
        printf("The result is correct!\n");
    else
        printf("Error: result is incorrect!\n");
}

#endif
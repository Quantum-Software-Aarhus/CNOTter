using trace = std::vector<std::pair<byte,byte>>;

matrix step_back(matrix x, hashset &level, trace &tr) {
    matrix mask = (1 << N) - 1; // to select row i
    for (byte i=0; i<N; i++) {
        matrix row = (x & mask) >> N*i;
        mask <<= N;
        for (byte j=0; j<N; j++, row <<= N) // add to row j
            if (i != j) {
                matrix prev = x ^ row; 
                if (find_level(prev, level)) {
                    tr.push_back(std::pair<byte,byte>(i,j));
                    return prev;
    }   }       }
    assert(false && "Predecessor not found");
    exit(-1);
}

matrix trace_back(matrix goal, hashset levels[], int depth, trace &tr, byte pi[N]) {
    for (int d=depth - 1; d>0; d--)
        goal = step_back(goal, levels[d], tr);
    id_perm(pi); // return identity permutation
    return goal;
}

// apply pi to all elements in trace tr
trace permute_trace(byte pi[N], trace &tr) {

    trace result;
    for (auto pair = tr.begin(); pair < tr.end(); pair++) 
        result.push_back(std::pair<byte,byte>(pi[pair->first], pi[pair->second]));
    
    return result;
}

trace trace_back_middle(matrix id, matrix middle, matrix goal, int fdepth, int bdepth, byte pi[N]) {
    trace fwd_trace, bwd_trace, result;
    matrix id_found, goal_found;
    id_found = trace_back(middle, bfs_fwd, fdepth, fwd_trace, pi); // dummy pi
    goal_found = trace_back(middle, bfs_bwd, bdepth, bwd_trace, pi); // dummy pi

    // REVERSE the forward trace and concatenate backward trace
    std::reverse(fwd_trace.begin(),fwd_trace.end());
    fwd_trace.insert(fwd_trace.end(), bwd_trace.begin(), bwd_trace.end()); 

#if SWAP==1
    // Reconstruct the proper permutations for the identity and goal
    byte pi1[N], pi2[N], pi3[N], pi4[N];
    equiv_perm(id, id_found, pi1, pi2);
    equiv_perm(goal_found, goal, pi3, pi4);
    assert(permute2(id, pi1, pi2) == id_found);      // id_found = (pi1,pi2) . id 
    assert(permute2(goal_found, pi3, pi4) == goal);  // goal = (pi3,pi4) . goal_found

    // construct q1 := pi1 ; pi2^-1 ; pi4^-1 and apply to the trace
    byte q0[N], q1[N],p4i[N];
    inv_perm(pi4, p4i);
    compose_inv_perm(pi2, p4i, q0);
    compose_perm(pi1, q0, q1);
    result = permute_trace(q1, fwd_trace);

    // construct the final permutation pi := pi3 ; q1
    compose_perm(pi3, q1, pi);

#else
    assert(id==id_found);
    // Reconstruct the proper permutation (we want "goal" instead of "found")
    byte pi_goal[N], pi_goal_found[N];
    representativePerm(goal_found, pi_goal_found); // goal_found = pi_found . found
    representativePerm(goal, pi_goal);   // goal_found = pi_goal . goal
    for (byte i=0; i<N; i++) pi[pi_goal_found[i]] = pi_goal[i];
    assert(permute(goal, pi) == goal_found);  // found = pi . goal 
    // PERMUTE the whole trace to the true goal
    result = permute_trace(pi, fwd_trace);
    id_perm(pi); // return the identity permutation
#endif

    return result;
}

void print_trace(matrix m, matrix goal, const trace &tr, byte pi[N]) {
    printf("\nOPENQASM 2.0;\n");
    printf("include \"qelib1.inc\";\n");
    printf("qreg q[%u];\n\n",N);
    matrix mask = (1 << N) - 1;
    for (std::pair<byte,byte> pair : tr) {
        byte i = pair.first, j = pair.second;
        printf("cx q[%u],q[%u];\n",i,j);
        matrix row = (m & (mask << N*i)) >> N*i;
        m ^= row << N*j;
    }
    printf("\nResult of the circuit:\n");
    pretty_matrix(m);
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

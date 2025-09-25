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

matrix trace_back(matrix goal, hashset levels[], int depth, trace &tr) {
    for (int d=depth - 1; d>0; d--)
        goal = step_back(goal, levels[d], tr);
    return goal;
}

// permute trace, with permutation that brings goal to other
trace permute_trace(matrix goal, matrix other, trace &tr) {

    // Reconstruct the proper permutation (we want "goal" instead of "other")
    byte pi_goal[N], pi_other[N], pi[N];
    representativePerm(other, pi_other); // goal_repr = pi_other . other
    representativePerm(goal, pi_goal);   // goal_repr = pi_goal . goal
    for (byte i=0; i<N; i++) pi[pi_other[i]] = pi_goal[i];
    assert(permute(goal, pi) == other);  // other = pi . goal 

    trace result;
    // REVERSED and PERMUTED forward trace (first half)
    for (auto pair = tr.begin(); pair < tr.end(); pair++) 
        result.push_back(std::pair<byte,byte>(pi[pair->first], pi[pair->second]));
    
    return result;
}

trace trace_back_middle(matrix id, matrix middle, matrix goal, int fdepth, int bdepth) {
    trace fwd_trace, bwd_trace, result;
    matrix other;
    other = trace_back(middle, bfs_fwd, fdepth, fwd_trace);
    assert(other==id);
    other = trace_back(middle, bfs_bwd, bdepth, bwd_trace);

    // REVERSE the forward trace and concatenate backward trace
    std::reverse(fwd_trace.begin(),fwd_trace.end());
    fwd_trace.insert(fwd_trace.end(), bwd_trace.begin(), bwd_trace.end()); 
    
    // PERMUTE the whole trace to the true goal
    result = permute_trace(goal, other, fwd_trace);
    return result;
}

void print_trace(matrix m, matrix goal, const trace &tr) {
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
    printf("\nResult:\n");
    pretty_matrix(m);
    assert(m==goal && "Result is incorrect");
}

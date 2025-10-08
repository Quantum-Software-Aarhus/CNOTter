// Jaco van de Pol, Aarhus University, 
// September 2024, September 2025
// main idea by Jens Emil Christensen

// Compile: 
// Assume nauty is built with "./configure --enable-tls; make" (flag is needed for thread safety)
// g++ -o matrix_cnot matrix_cnot.cpp -fopenmp -DN=6 -DE=1 -DNAUTY=1 -DSWAP=0 -O3 -DNDEBUG -march=native -Inauty/ nauty/nautyW1.a -DWORDSIZE=32 -DMAXN=WORDSIZE

#include <array>
#include <vector>
#include <omp.h>
#include "hashset.h" // thread-safe hash set from dtree project
#include "options.h" // defines N,E,MAX,SWAP,NAUTY,POLY,BEAT, see also matrix_cnot.sh
#include "timing.h"
#include "matrix.h"
#include "repr.h"
#include "trace_back.h"

// precalculated 2-log of the orbit level sizes (0-terminated)
// NOTE: the size depends on if SWAPs are free or not

#if SWAP==0
const std::array<std::vector<byte>,9> levelSizes = {{
    // first number is for level depth=2 (externally: Depth=1)
    {}, //0
    {0}, // 1
    {0,0,0,0}, // 2
    {0,3,4,4,3,0,0}, // 3
    {0,3,5,7,8,9,8,5,0,0}, // 4
    {0,3,5,8,11,13,14,15,15,13,8,0,0}, // 5
    {0,3,6,8,11,14,17,19,22,23,24,23,20,11,0,0}, //6
    {0,3,6,8,11,15,18,21,24,27,30,32,33,34,33,29,17,0,0}, //7
    {0,3,6,8,11,15,18,22,25,29,32,35, /* guess from here on */ 37,38,40,41,40,38,36,34,0,0} //8
}};
#else
const std::array<std::vector<byte>,9> levelSizes = {{
    // first number is for level depth=2 (externally: Depth=1)
    {}, //0
    {0}, // 1
    {0,0,0,0}, // 2
    {0,3,4,4,3,0,0}, // 3
    {0,3,5,5,3,0,0,0,0,0}, // 4
    {0,3,5,7,9,9,7,3,0,0,0,0,0}, // 5
    {0,3,5,8,10,13,14,15,13,10,3,0,0,0,0,0}, //6
    {0,3,5,8,11,14,16,19,21,22,22,20,13,2,0,0,0,0,0}, //7
    {0,3,5,8,11,14,17,20,23,26,28,30,31,30,28,21,3,0,0,0,0,0} //8
}};
#endif

#if POLY==1
std::array<std::atomic<uint64_t>,N+1>poly; // coefficients of the polynomial at distance N/2
#endif

hashset bfs_levels[3*N];    // for one-directional BFS
hashset bfs_fwd[(3*N+2)/2]; // for bi-directional BFS
hashset bfs_bwd[(3*N+1)/2];

void Add(const matrix &x, byte i, byte j, 
                hashset *prev, hashset *current, hashset *next, int depth,
                uint64_t &level, uint64_t &count) {
    uint64_t mask = (1UL<<N*(i+1)) - (1UL<<N*i);
    uint64_t row = (x & mask) >> i*N;
    matrix y = x ^ (row << j*N);
    uint64_t Orbit = representative(y);
    if (!prev->contains(y) && !current->contains(y) && next->insert(y)) {
        // only insert and count if new; 
        level += Orbit;
        count++;
#if POLY==1
        if (2*(depth-1)==N) {
            byte ess = countEssential(y);
            poly[N-ess] += Orbit * (fac[ess] * fac[N-ess]) / fac[N];
        }
#endif
    }
}

uint64_t init_level(hashset levels[], matrix start) {
    levels[0] = hashset(); // level 0 (prev)
    levels[0].init(3);
    levels[1] = hashset(); // level 1 (current)
    levels[1].init(3);
    uint64_t Orbit = representative(start); // modifies start
    levels[1].insert(start);
    return Orbit;
}

// explore and count all successors of the current level
uint64_t next_level(uint64_t &size, hashset levels[], uint32_t depth) { 
    std::atomic<uint64_t> level(0);
    std::atomic<uint64_t> count(0);

    // current and prev are accessed read-only
    // next is modified (extended) concurrently

    auto prev = &levels[depth-2];
    auto current = &levels[depth-1];
    auto next = &levels[depth];

    current->parallelForAll(
        [&](matrix x){
            uint64_t loc_level=0, loc_count=0;
            for (byte i=0; i<N; i++)
                for (byte j=0; j<N; j++) // add to row j
                    if (i != j) Add(x, i, j, prev, current, next, depth, loc_level, loc_count);
        if (loc_level > 0) {
            level += loc_level;
            count += loc_count;
        }
#if BEAT>0
        size_t worker = omp_get_thread_num();
        if (passedTime(lifeTime[worker]) >= BEAT) { // every minute
            # pragma omp critical
            {
                lifeBeat(worker, loc_level, loc_count);
            }
            lifeTime[worker] = system_clock::now();
        }
#endif
    });
    size = count;
    return level;
}

int generate_bfs(matrix start, matrix goal, byte limit, hashset bfs_levels[]) {

    // initialize Breadth-First Search
    byte depth = 1, tableSize = 3;
    uint64_t level, levels, orbit, orbits;
    orbit = orbits = 1;

    printf("Depth 0 (2^3): "); fflush(stdout);
    levels = level = init_level(bfs_levels, start);

    while (orbit) {
        report(level, orbit);
        if (goal)
            { if (find_level(goal, bfs_levels[depth])) return -depth; }
        else 
            { if (depth > 1) bfs_levels[depth-2].deinit(); }
        if (depth-1 == limit) return depth;
        depth++;
        tableSize = std::min(std::max(levelSizes[N][depth-2] + E, 3), MAX);
        bfs_levels[depth] = hashset();
        bfs_levels[depth].init(tableSize);
        printf("Depth %u (2^%u): ", depth-1, tableSize); fflush(stdout);
        levels += level = next_level(orbit, bfs_levels, depth);
        orbits += orbit;
    }
    depth--;
    printf("--\n");
/*  // print matrices at the last level (out of curiosity)
    bfs_levels[depth].forAll(
        [](matrix c) { pretty_matrix(c); }
    );
*/
    printf("Total size: %lu (%lu orbits), completed at depth %u\n", levels, orbits, depth-1);
    return depth;
}

matrix intersect(hashset &L1, hashset &L2) {
    std::atomic<matrix> joint(0);
    L1.parallelForAll([&](matrix x){
        if (L2.contains(x)) joint=x; // How to terminate when found?
    });
    return joint;
}

// Bidirectional search yields a matrix in the intersection of Fwd(start) and Bwd(goal)
// We also return the depths of the fwd and bwd search (fdepth,bdepth)
// We return (0,fdepth,bdepth) if start and goal are not connected

using triple = std::pair<matrix,std::pair<byte,byte>>;
inline triple Triple(matrix m, byte d1, byte d2) {
    return std::pair<matrix,std::pair<byte,byte>>(m, std::pair<byte,byte>(d1, d2));
}

triple bidirectional(matrix start, matrix goal, byte limit, hashset bfs_fwd[], hashset bfs_bwd[]) {

    // initialize Bidirectional fwd/bwd Search
    byte fdepth = 1, bdepth=1, tableSize;
    uint64_t level, forbit, borbit, levels, orbits;
    forbit = borbit = 1; orbits = 2;
    levels = level = init_level(bfs_fwd, start);
    printf("Fwd Depth 0 (2^3): "); report(level, forbit);
    levels += level = init_level(bfs_bwd, goal);
    printf("Bwd Depth 0 (2^3): "); report(level, borbit);
    matrix m = intersect(bfs_fwd[fdepth], bfs_bwd[bdepth]);
    if (m) return Triple(m, fdepth, bdepth);

    while (fdepth + bdepth - 2 < 3*(N-1)) { // expand the smallest level
        if (fdepth+bdepth-2 == limit) return Triple(m, fdepth, bdepth);
        if (forbit <= borbit) {
            fdepth++; 
            tableSize = std::min(std::max(levelSizes[N][fdepth-2] + E, 3), MAX);
            printf("Fwd Depth %u (2^%u): ", fdepth-1, tableSize); fflush(stdout);
            bfs_fwd[fdepth] = hashset();
            bfs_fwd[fdepth].init(tableSize);
            levels += level = next_level(forbit, bfs_fwd, fdepth);
            orbits += forbit;
            report(level, forbit);
        }
        else {
            bdepth++;
            // Note: this Bwd level is smaller than next Fwd one
            // Problem: Bwd's successor can still be larger than Fwd's successor (hence 10)
            tableSize = std::min(std::max(levelSizes[N][fdepth-1] + E, 10), MAX); 
            printf("Bwd Depth %u (2^%u): ", bdepth-1, tableSize); fflush(stdout);
            bfs_bwd[bdepth] = hashset();
            bfs_bwd[bdepth].init(tableSize);
            levels += level = next_level(borbit, bfs_bwd, bdepth);
            orbits += borbit;
            report(level, borbit);
        }
        m = intersect(bfs_fwd[fdepth], bfs_bwd[bdepth]);
        if (m) return Triple(m, fdepth, bdepth);
    }
    printf("Not found at distance %u+%u (%lu, %lu)\n", fdepth-1, bdepth-1, levels, orbits);
    return Triple(0, fdepth, bdepth);
}

/*
 * Main: if first arg is -K, set K as limit. If last argument is not -*, set as goal
 */

int main(int argc, char const *argv[]) {
#if NAUTY==1
    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
    options.getcanon=true;   // we want the canonical graph
    options.defaultptn=true; // default coloring
#endif
    if (N<1 || N>8) {
        printf("N={%u} not supported, only N=1..8\n", N);
        exit(-1);
    }
    if (POLY==1 && SWAP==1) {
        printf("Polynomial coefficients are not supported with SWAP\n");
        exit(-1);
    }
    printf("Handling matrices of size N = %u\n", N);
    printf("Using DTree + %u extra bits, max-size %u\n", E, MAX);
    printf("Use Nauty: %u. Swaps-for-free: %u. Polynomial: %u\n", NAUTY, SWAP, POLY);
    #if defined(_OPENMP)
        printf("Running with %d OpenMP threads\n",omp_get_max_threads());
    #endif

    matrix id=1; // compute identity matrix
    for (byte i=1; i<N; i++) id = (id << (N+1)) | 1;
    matrix goal=0; // search for goal: set with last argument "filename"
    byte  limit=-1; // search limit when >=0: set with argument "-<limit>"

    if (argc>1 && argv[1][0]=='-') {
        limit = atoi(argv[1]+1); // skip the leading '-'
        if (limit!=(byte)-1)     // unsigned, so this is 255
            printf("Cutting off at maximum distance: %d\n", limit);
    }
    if (argc>1 && argv[argc-1][0]!='-') {
        goal = read_matrix(argv[argc-1]);
        //investigate(goal);
        assert(goal!=0 && "0-matrix cannot be generated");
    }
    if (goal) {
        triple m = bidirectional(id, goal, limit, bfs_fwd, bfs_bwd);
        matrix middle = m.first;
        int fdepth = m.second.first;
        int bdepth = m.second.second;
        if (m.first) {
            printf("Found at distance %u (%u + %u)\n", fdepth + bdepth - 2, fdepth-1, bdepth-1);
            perm pi;
            trace concat = trace_back_middle(id, middle, goal, bfs_fwd, bfs_bwd, fdepth, bdepth, pi);
            print_trace(id, goal, concat, pi);
        } else {
            printf("Goal not found after %d steps: \n", fdepth+bdepth-2);
            pretty_matrix(goal);
        }
    } else {
        int depth = generate_bfs(id, goal, limit, bfs_levels); 
        if (goal) { // currently unreachable, since bidirectional is preferred
            if (depth < 0) { // negative means goal is found 
                depth = -depth;
                printf("Goal found at level %d\n", depth-1);
                trace bfs_trace;
                matrix other = trace_back(goal, bfs_levels, depth, bfs_trace);
                assert(other==id);
                std::reverse(bfs_trace.begin(), bfs_trace.end());
                perm pi; id_perm(pi);
                print_trace(other, goal, bfs_trace, pi);
            }
            else { // currently unreachable
                printf("Goal not found after %d steps: \n", depth-1);
                pretty_matrix(goal);
            }
        }
#if POLY==1
        if (!(N&1) && (depth-1)*2 >= N) {
            printf("Polynomial coefficients (%u/%u): [", N, N/2);
            for (byte i=0; i<=N; i++)
                printf("%lu%c ", poly[N-i].load(std::memory_order_relaxed), (i<N ? ',' : ']'));
            printf("\n");
        }
#endif
    }
    std::cout << std::setprecision(std::numeric_limits<double>::digits10)
              << "Total time: " << currentTime() << "s" << std::endl;
}

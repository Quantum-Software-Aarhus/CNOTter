#include <fstream>

typedef uint8_t byte;
typedef uint64_t matrix; // store at most 8x8 Booleans

void pretty_perm(byte pi[N]) {
    for (byte i=0; i<N; i++)
        printf("%3u", i);
    printf("\n");
    for (byte i=0; i<N; i++)
        printf("%3u", pi[i]);
    printf("\n");
}

void pretty_matrix(matrix x) {
    std::string delimiter(N*2-1,'=');
    std::cout << delimiter << std::endl;
    for (byte i=0; i<N; i++) {
        for (byte j=0; j<N; j++, x >>= 1)
            printf("%lu ", x & 1);
        printf("\n");
    }
    std::cout << delimiter << std::endl;
}

matrix read_matrix(std::string filename) {
    std::ifstream input(filename, std::ios_base::in);
    if (!input.is_open()) { 
        std::cerr << "Could not open input file: " << filename << "\n";
        exit(-1); 
    }
    matrix result=0;
    byte idx=0;
    for (byte i=0; i<N; i++)
        for (byte j=0; j<N; j++, idx++) {
            char c=0;
            do {
                input.get(c);
            } while (!input.eof() && (c==' ' || c=='\n' || c=='\t' || c=='\r'));
            if (c=='1') result ^= 1UL << idx;
            else assert(c=='0' && "Expected input 0 or 1");
        }
    return result;
}

inline matrix permute(matrix x, const byte pi[N]) { // y[i][j] := x[pi[i]][pi[j]]
    matrix y = 0;
    for (byte i=N-1; i<N; i--)
        for (byte j=N-1; j<N; j--) {
            y <<= 1;
            y |= (x >> (pi[i]*N + pi[j])) & 1;
        }
    return y;
}

inline bool testEssential(matrix x, byte i) {
    if (!(x & 1UL<<(N+1)*i))
        return true;
    for (byte j=0; j<N; j++)
        if (j!=i && (x & 1UL<<(N*j+i) || x & 1UL<<(N*i+j)))
            return true;
    return false;
}

// Count inessential indices (=isolated vertices in nauty)
inline byte countEssential(matrix x) {
    byte ess=0;
    for (byte i=0; i<N; i++)
        if (testEssential(x,i)) ess++;
    return ess;
}
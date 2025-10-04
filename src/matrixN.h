#include "options.h"

#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <array>


#define N 18
#define NC (64 / N)        // Columns: number of rows in a single word.
#define NR ((N-1) / NC +1) // Rows: number of words needed (last word may not be filled)


class Matrix {

public:
    std::array<uint64_t,NR> _bits = {{ }}; // initalizes to 0

public:

    Matrix(bool diag=1) { // create identity matrix
        for (uint8_t i=0; i<N; i++)
            set(i,i,diag);
    }

    inline bool get(uint8_t i, uint8_t j) const {
        div_t wordi = div(i,NC);
        uint64_t res = _bits[wordi.quot] & (1UL << (N * wordi.rem + j));
        return (res ? true : false);
    }

    inline void set(uint8_t i, uint8_t j, bool val) {
        div_t wordi = div(i,NC);
        uint64_t mask = (1UL << (N * wordi.rem + j));
        if (val)
            _bits[wordi.quot] |= mask;
        else
            _bits[wordi.quot] &= ~ mask;
    }

    bool operator==(const Matrix &other) const { // assume unused bits are 0
        for (uint8_t i=NR-1; i<NR; i--) {
            if (_bits[i] != other._bits[i]) return false;
        }
        return true;
    }
    
    bool operator<(const Matrix &other) const { // assume unused bits are 0
        for (uint8_t i=NR-1; i<NR; i--) {
            if (_bits[i] < other._bits[i]) return true;
            if (_bits[i] > other._bits[i]) return false;
        }
        return false;
    }

    /* pretty print matrix */
    void print() const {
        std::string delimiter(N*2-1,'=');
        std::cout << delimiter << std::endl;
            for (uint8_t i=0; i<N; i++) {
                for (uint8_t j=0; j<N; j++)
                    printf("%c ", (get(i,j) ? '1' : '0'));
                printf("\n");   
            }
        std::cout << delimiter << std::endl;
    };

    static Matrix read(std::string filename) {
        std::ifstream input(filename, std::ios_base::in);
        if (!input.is_open()) { 
            std::cerr << "Could not open input file: " << filename << "\n";
            exit(-1); 
        }
        Matrix result;
        for (u_int8_t i=0; i<N; i++)
            for (u_int8_t j=0; j<N; j++) {
                char c=0;
                do {
                    input.get(c);
                } while (!input.eof() && (c==' ' || c=='\n' || c=='\t' || c=='\r'));
                if (c=='1') result.set(i,j,true);
                else {
                    assert(c=='0' && "Expected input 0 or 1");
                    result.set(i,j,false);
                }
        }
        return result;
    }

    /* Add row i to row j*/
    Matrix addrow(uint8_t i, uint8_t j) const {
        assert(i!=j && i<N && j<N);
        Matrix result=*this;
        //printf("Adding row %d to row %d\n",i,j);
        for (uint8_t k=0; k<N; k++)
            result.set(j,k,get(i,k) != get(j,k)); // xor
        return result;
    };

    Matrix permute(const uint8_t pi[N]) const { // other[i][j] := this[pi[i]][pi[j]]
        Matrix other;
        for (uint8_t i=0; i<N; i++) {
            for (uint8_t j=0; j<N; j++) {
                other.set(i,j, get(pi[i],pi[j]));
            }
        }
        return other;
    }

};

#include "options.h"

#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <array>


#define N 18
#define NC (64 / N)        // Columns: number of rows in a single word.
#define NR ((N-1) / NC +1) // Rows: number of words needed (last word may not be filled)

typedef uint8_t byte;

// TODO: move to permutation.h

typedef byte perm[N];       // permutation of N elements

void pretty_perm(const perm pi) {
    for (byte i=0; i<N; i++)
        printf("%3u", i);
    printf("\n");
    for (byte i=0; i<N; i++)
        printf("%3u", pi[i]);
    printf("\n");
}

// return the identity permutation in pi
inline void id_perm(perm pi) {
    for (byte i=0; i<N; i++)
        pi[i] = i;
}

// return the inverse permutation in pi_inv
inline void inv_perm(const perm pi, perm pi_inv) {
    for (byte i=0; i<N; i++)
        pi_inv[pi[i]] = i;
}

// return the composition pi = pi1 . pi2
inline void compose_perm(const perm pi1, const perm pi2, perm pi) {
    for (byte i=0; i<N; i++)
        pi[i] = pi2[pi1[i]]; // non-standard, since we permute indices
}

// return the composition pi = pi1^-1 . pi2
inline void compose_inv_perm(const perm pi1, const perm pi2, perm pi) {
    for (byte i=0; i<N; i++)
        pi[pi1[i]] = pi2[i]; // non-standard, since we permute indices
}

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

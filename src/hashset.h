/*
 * Jaco van de Pol, Aarhus University, September 2024
 * Copied and modified from Freark van der Berg, https://github.com/bergfi/dtree
 */

/*
 * Dtree - a concurrent compression tree for variable-length vectors
 * Copyright © 2018-2021 Freark van der Berg
 *
 * This file is part of Dtree.
 *
 * Dtree is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * Dtree is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dtree.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include <atomic>
#include <algorithm>
#include <cstdio>
#include <sys/mman.h>
#include <thread>
#include <assert.h>
#include <map>


template<typename TO_TYPE, typename TREE>
class Linear {
public:
    TREE& tree;
    TO_TYPE& e;

    Linear(TREE& tree, TO_TYPE& e): tree(tree), e(e) {}

    __attribute__((always_inline))
    void next() {
        e = (e+1) & tree._entriesMask;
    }
};

template<typename TO_TYPE, typename TREE>
class QuadLinear {
public:
    TREE& tree;
    TO_TYPE& e;
    TO_TYPE eBase;
    TO_TYPE eOrig;
    TO_TYPE inc;

    QuadLinear(TREE& tree, TO_TYPE& e): tree(tree), e(e), eBase(e & ~0x7ULL), eOrig(e & 0x7ULL), inc(1ULL) {}

    __attribute__((always_inline))
    void next() {
        e = (e + 1) & 0x7ULL;
        if(e == eOrig) {
            TO_TYPE diff = (inc*2);
            diff -= __builtin_popcountll(diff); // TODO: SHOULD DEPEND ON TO_TYPE
            eBase = (eBase + diff*8) & tree._entriesMask;

            inc++;
            if(inc==1000) {
                printf("Quad Linear 1000 inc++ exceeded\n");
            }
        }
        e += eBase;
    }
};

// Adapted from: https://github.com/aappleby/smhasher/blob/master/src/MurmurHash2.cpp
//-----------------------------------------------------------------------------
// MurmurHash2, 64-bit versions, by Austin Appleby
// 64-bit hash for 64-bit platforms
// adapted: just one 64-bit word

inline uint64_t MurmurHash64 ( uint64_t k ) {
    const uint64_t m = 0xc6a4a7935bd1e995;
    const int r = 47;
    if (k==0) return 0;
    uint64_t h = 8*m; // ^ seed

    k ^= k*m >> r; 
    h ^= k*m;
    h ^= h*m >> r;
    h ^= h*m >> r;

    return h;
}

template<typename T=uint64_t>
struct MurmurHash {
    __attribute__((always_inline))
    size_t hash( const T& k ) const {
        return MurmurHash64(k);
    }
};


template<typename T>
struct HashCompare {
    __attribute__((always_inline))
    size_t hash( const T& k ) const {
        return k;
    }
};

template< typename TO_TYPE
        , template<typename,typename> typename BUCKETFINDER = Linear
        , template<typename> typename HASH = HashCompare
        >
class HashSet {
public:
    using Bucketfinder = BUCKETFINDER<TO_TYPE, HashSet<TO_TYPE, BUCKETFINDER, HASH>>;
    friend Bucketfinder;
public:

    HashSet(): _scale(0), _buckets(0), _entriesMask(0), _map(nullptr) {
        assert (HASH<uint64_t>().hash(0) == 0 && "0 should be hashed to 0");
    }

    HashSet& init(size_t scale=28ULL) {
        assert(scale>2 && "scale should be at least 3");
        assert(scale <= sizeof(TO_TYPE)*8 && "scale to large for entry type");
        _scale = scale;
        _buckets = 1ULL << _scale;
        _entriesMask = (_buckets - 1);
        
        assert(!_map && "map already in use");
        _map = (decltype(_map))mmap(nullptr, _buckets * sizeof(uint64_t), PROT_READ|PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
        assert(_map && "failed to mmap data");
        //printf("......allocated table %ld: %ld bytes\n", _scale, _buckets * sizeof(uint64_t));
        return *this;
    }

    void deinit() {
        if(_map) {
            munmap(_map, _buckets * sizeof(uint64_t));
            //printf("......deallocated table %ld: %ld bytes\n", _scale, _buckets * sizeof(uint64_t));
            _map = nullptr;
        }
        //else printf("......huh???\n");
    }

    ~HashSet() { deinit(); }

    HashSet& operator=(const HashSet& other) = delete;

    HashSet& operator=(HashSet&& other) {
        deinit();
        _scale = other._scale;
        _buckets = other._buckets;
        _entriesMask = other._entriesMask;
        _map = other._map;
        other._map = nullptr;
        return *this;
    }

    TO_TYPE entry(uint64_t key) {
        TO_TYPE h = HASH<uint64_t>().hash(key);
        return h & _entriesMask;
    }

    __attribute__((always_inline))
    constexpr uint64_t newlyInserted(uint64_t v) const {
        return v | 0x8000000000000000ULL;
    }

    template<int INSERT>
    TO_TYPE insertOrContains(uint64_t key, bool &is_new) {
        assert(_map && "storage not initialized");
        assert(key && "cannot store 0");
        assert(_buckets == 1ULL << _scale);
        TO_TYPE e = entry(key);
        e += e == 0;
        Bucketfinder searcher(*this, e);
        std::atomic<uint64_t>* current = &_map[e];

        size_t probeCount = 1;

        while(probeCount < _buckets) {
            uint64_t k = current->load(std::memory_order_relaxed);
            if(k == 0ULL) {

                // Make sure we do not use the 0th index
                if(e == 0) goto findnext;
                if (!INSERT) {
                    return 0;
                }
                if(current->compare_exchange_strong(k, key, std::memory_order_release, std::memory_order_relaxed)) {
                    // printf("...insert %ld at %ld\n",key,e);
                    is_new = true;
                    return e;
                }
            }
            if(k == key) {
                is_new = false;
                return e;
            }
            findnext:
            searcher.next();
            current = &_map[e];
            probeCount++;
        }
        printf("Hash map full! (2^%ld=%ld buckets)\n",_scale,_buckets);
        exit(-1);
    }

    __attribute__((always_inline))
    TO_TYPE findOrPut(uint64_t key) {
        bool dummy;
        return insertOrContains<1>(key, dummy);
    }

    __attribute__((always_inline))
    bool insert(uint64_t key) {
        bool is_new;
        insertOrContains<1>(key, is_new);
        return is_new;
    }

    __attribute__((always_inline))
    TO_TYPE contains(uint64_t key) {
        bool dummy;
        return insertOrContains<0>(key, dummy);
    }

    __attribute__((always_inline))
    uint64_t get(TO_TYPE idx) {
        assert(0 <= idx);
        assert(idx < _buckets);
        uint64_t result = _map[idx].load(std::memory_order_relaxed);
        return result;
    }

    void stats() {
        printf("...table 2^%ld: ",_scale);
        std::atomic<size_t> count(0);
        #pragma omp parallel for
        for (uint64_t i=0; i<_buckets; i++)
            if (_map[i]) count++;
        printf("Size: %ld\n",count.load());
    }

    void statistics(bool verbose=false) {
        printf("...table 2^%ld: ",_scale);
        size_t count=0;
        size_t run=0;
        bool last=false;
        std::map<uint64_t,uint64_t> frequency;

        for (uint64_t i=0; i<_buckets; i++) {
            // printf("%d",get(i)>0);
            if (_map[i]) {
                count++;
                run++;
                last=true;
            }
            else {
                if (last && verbose) {
                    frequency[run]++;
                    last=false;
                    run=0;
                }
            }
        }
        printf("Size: %ld\n",count);
        for (auto f : frequency) {
            if (verbose) printf("      probes of length %ld : %ldx\n",f.first,f.second);
        }
    }

    template<typename FUNC>
    void forAll(FUNC&& func) {
        for(TO_TYPE idx = 0; idx < _buckets; ++idx) {
            uint64_t value = get(idx);
            if(value) {
                func(value);
            }
        }
    }

    template<typename FUNC>
    void parallelForAll(FUNC&& func) {
        #pragma omp parallel for
        for(TO_TYPE idx = 0; idx < _buckets; ++idx) {
            uint64_t value = get(idx);
            if(value) {
                func(value);
            }
        }
    }

public:
    size_t _scale;
    size_t _buckets;
    size_t _entriesMask;
    std::atomic<uint64_t>* _map;
};

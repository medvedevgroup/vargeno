#ifndef GENERATEDICT_H
#define GENERATEDICT_H
#include <iostream>
#include <string>
#include <sdsl/bit_vectors.hpp>

#include <vector>
#include <unordered_set>
#ifdef __cplusplus
extern "C" {
#endif

#include <assert.h>


#ifdef __cplusplus
}
#endif
#include "util.h"
#include "allsome_util.h"

using namespace std;

class GenomeSeq {
public:
	string name;
	string seq;
	int size;
	GenomeSeq() {
		name = "";
		seq = "";
		size = -1;
	}
};

// simple version of bloom filter
// modified based on https://github.com/medvedevgroup/bloomtree-allsome
class BloomFilter {
public:

    static const uint64_t BIT32 = 4294967296;
    static const uint64_t BIT40 = 1099511627776;
	
	BloomFilter(int _value_range) {
		bv = NULL;
		value_range = _value_range;
	}

	BloomFilter(uint64_t num_filter_bits, int _value_range) {
		if (num_filter_bits > 0) {
			bv = new sdsl::bit_vector(num_filter_bits);
			value_range = _value_range;
			bv_size = bv->size();
		}
		else {
			std::cerr << "[BloomFilter] Error: bloom filter size <= 0" << std::endl;
		}
	}
	
	~BloomFilter() {
		if (bv != NULL) {
			delete bv;
			bv = NULL;
		}
	}
	
	void load(string & filename) {
		assert(bv == NULL);
		//auto start_time = get_wall_time();
		bv = new sdsl::bit_vector();
        if (!sdsl::load_from_file(*bv, filename)) {
			std::cerr << "WARNING: Attempt to load \"" << filename << "\" failed!" << std::endl;
		}
		//auto elapsed_time = elapsed_wall_time(start_time);
		//std::cerr << "[BloomFilter load] " << elapsed_time << " secs " << filename << std::endl;
		bv_size = bv->size();
	}
	
	uint64_t num_filter_bits() const {
		return bv_size;
	}

	void save(string & filename) {
		//std::cerr << "Saving uncompressed BF to " << filename << std::endl;
		//auto start_time = get_wall_time();
		sdsl::store_to_file(*bv, filename);
		//auto elapsed_time = elapsed_wall_time(start_time);
		//std::cerr << "[BloomFilter save] " << elapsed_time << " secs " << filename << std::endl;
	}
	
	int check_bit(uint64_t pos) const {
		//std::cout << "UncompressedBF::[" << pos << "]" << std::endl;
		return (*bv)[pos];
	}

	void set_bit(uint64_t p) {
		(*bv)[p] = 1;
	}

	void set_value(uint64_t v) {
		if (value_range == 32) {
			set_bit(hash32(v) % bv_size);
		}
		else if (value_range == 40){
			set_bit(hash40(v) % bv_size);
		}
		else {
			std::cerr << "[BloomFilter set_value] Error: unrecognized value_range" << std::endl;
		}
	}

    int check_value(uint64_t v){
        if(value_range == 32){
            return check_bit(hash32(v) % bv_size);
        }else if (value_range == 40){
            return check_bit(hash40(v) % bv_size);
        }else{
            std::cerr << "[BloomFilter check_value] Error: unrecognized value_range" << std::endl;
            return 0;
        }
    }

	// this function is directly borrowed from :
	// http://stackoverflow.com/questions/664014/what-integer-hash-function-are-good-that-accepts-an-integer-hash-key/12996028#12996028
    static uint32_t hash32(uint32_t x) {
		x = ((x >> 16) ^ x) * 0x45d9f3b;
		x = ((x >> 16) ^ x) * 0x45d9f3b;
		x = (x >> 16) ^ x;
		return x;
	}

	static uint32_t hash24(uint32_t) {
		std::cerr << "hash24 not avaialble" << std::endl;
		return 0;
	}

	static uint64_t hash40(uint64_t x) {
		x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
		x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
		x = x ^ (x >> 31);
        return x;
	}

    /*
    static uint32_t hash32(uint32_t x) {
        return x;
        //return (x*2654435761) % BIT32;
    }
    
    static uint64_t hash40(uint64_t x) {
        return x;
        //return (x*2654435761) % BIT40;
    }
    */
    

	uint64_t count_ones() {
		uint64_t* data = bv->data();
		// RSH: using size() instead of capacity() misses the bits in the final
		// .. chunk if size%64 is non-zero; caveat: using capacity() risks counting
		// .. extra bits in the final chunk
		sdsl::bit_vector::size_type len = bv->capacity() >> 6;
		uint64_t count = 0;
		for (sdsl::bit_vector::size_type p = 0; p < len; ++p) {
			int pc = popcount(*data);
			assert(pc == __builtin_popcountl(*data));
			data++;
			count += pc;
		}

		sdsl::rank_support_v<> rbv(bv);
		assert(rbv.rank(bv->size()) == count);
		return count;
	}

private:

	int popcount(uint64_t b) {
		int count = 0;
		for (int i = 0; i < 64; i++) {
			if (b % 2) count++;
			b >>= 1;
		}
		return count;
	}

	sdsl::bit_vector* bv;
	int value_range;
	uint64_t bv_size;
};

class BFGenerator {
public:
	BFGenerator();
	~BFGenerator();
	int readFasta(string & fasta_filename);
	int constructBfFromGenomeseq(string bf_filename, bool is_canonical);
	int constructBfFromVcf(const string & vcf_filename, string bf_filename, bool is_canonical);
	int constructBfFromUcsc(const string & text_filename, string bf_filename, bool is_canonical);
    int constructBfFromEncode(const string & encode_filename, string bf_filename, bool is_canonical);
    static const int64_t REF_BF_BYTES = 1200000000; // in terms of bytes, 1 byte = 8 bit
    static const int64_t HUGE_LITE_BF_BYTES = 2300000000; // this is only for temp usage
    static const int64_t SNP_BF_BYTES = 140000000; // in terms of bytes
    static const int64_t REF_LITE_BF_BYTES = 280000000; // in terms of bytes
	// double the size of bloom filter
    //static const int64_t REF_BF_BYTES = 2400000000; // in terms of bytes, 1 byte = 8 bit
    //static const int64_t REF_LITE_BF_BYTES = 4600000000;
    //static const int64_t SNP_BF_BYTES = 280000000; // in terms of bytes
	static const int REF_BF_RANGE = 32;
	static const int SNP_BF_RANGE = 40;

protected:
	vector<GenomeSeq> genome_vector;
	static const int KMER_SIZE = 32; // only right half kmers
	static const int HASH_NUM = 1;

	char rev(const char c);

};

#endif /* GENERATEDICT_H */

#ifndef BF_H
#define BF_H

#include <string>
#include <sys/mman.h>
#include <sdsl/bit_vectors.hpp>
#include <jellyfish/mer_dna_bloom_counter.hpp>
#include "Kmers.h"
#include "sbt_util.h"
#include "roaring.hh"

using HashPair = jellyfish::hash_pair<jellyfish::mer_dna>;

using sdsl_rrr_vector = sdsl::rrr_vector<255>;

// a kmer bloom filter
class BF {
public:
    BF(const std::string & filename, HashPair hp, int nh);
    virtual ~BF();

    virtual void load(uint64_t partitions = 1);
    virtual void save();

    virtual int operator[](uint64_t pos) const;
    virtual void set_bit(uint64_t p);
    virtual uint64_t num_filter_bits() const;
    virtual uint64_t hash_function_range() const;
    void update_hash_function_range();
    void set_filter_partitions(uint64_t partitions);

    void add(const jellyfish::mer_dna & m);

    bool contains(const jellyfish::mer_dna & m) const;
    bool contains(const std::string & str) const;
    bool contains(const uint64_t pos) const;

    uint64_t kmer_position(const jellyfish::mer_dna & m) const;

    virtual uint64_t similarity(const BF* other, int type) const;
    virtual std::tuple<uint64_t, uint64_t> b_similarity(const BF* other) const;
    virtual BF* union_with(const std::string & new_name, const BF* f2) const;
    virtual BF* intersection_with(const std::string & new_name, const BF* f2) const;
    virtual void union_into(const BF* f2);
    virtual BF* append(const std::string & new_name, const BF* f2) const;
    virtual uint64_t count_ones() const;
    virtual void compress_rrr();
    virtual void compress_wraprrr(const BF* f2);
    virtual void compress_roar();
    virtual void uncompress();
    virtual void split(BF* child0_bf_all, const BF* child0_bf_some,
                       BF* child1_bf_all, const BF* child1_bf_some,
                       const std::string& f_all,  BF*& bf_all,
                       const std::string& f_some, BF*& bf_some);
    virtual std::string compression_name();
//protected: //needed to unprotect them for "uncompress()" (Rayan)
    std::string filename;
    sdsl_rrr_vector* bits;

    HashPair hashes;
    unsigned long num_hash;
    uint64_t hash_function_modulus;
    uint64_t filter_partitions;      // number of filters stored in the bits (usually 1)
};


class WrapRrrBF: public BF  {
public:
    WrapRrrBF(const std::string & filename, HashPair hp, int nh, uint64_t num_filter_bits = 0);
    virtual ~WrapRrrBF();

    virtual void load(uint64_t partitions = 1);
    virtual void save();
    virtual int operator[](uint64_t pos) const;

    virtual uint64_t num_filter_bits() const;

    virtual BF* append(const std::string & new_name, const BF* f2) const;
    virtual void compress_rrr();
    virtual void compress_wraprrr(const BF* f2);
    virtual void compress_roar();
    virtual void uncompress();

    virtual std::string compression_name();

	int bits_constant;       // 0 or 1, valid only if bits == nullptr
	int bits2_constant;      // 0 or 1, valid only if bits2 == nullptr
    sdsl_rrr_vector* bits2;
    uint64_t total_filter_bits;
};


class RoarBF: public BF  {
public:
    RoarBF(const std::string & filename, HashPair hp, int nh, uint64_t num_filter_bits = 0);
    virtual ~RoarBF();

    virtual void load(uint64_t partitions = 1);
    virtual void save();

    virtual int operator[](uint64_t pos) const;
    virtual void set_bit(uint64_t p);
    virtual uint64_t num_filter_bits() const;
    virtual uint64_t num_compressed_bytes() const;

    void add(const jellyfish::mer_dna & m);

    virtual uint64_t similarity(const BF* other, int type) const;
    virtual std::tuple<uint64_t, uint64_t> b_similarity(const BF* other) const;
    virtual BF* union_with(const std::string & new_name, const BF* f2) const;
    virtual BF* intersection_with(const std::string & new_name, const BF* f2) const;
    virtual void union_into(const BF* f2);
    virtual BF* append(const std::string & new_name, const BF* f2) const;
    virtual uint64_t count_ones() const;
    virtual void compress_rrr();
    virtual void compress_wraprrr(const BF* f2);
    virtual void compress_roar();
    virtual void uncompress();
    virtual void split(BF* child0_bf_all, const BF* child0_bf_some,
                       BF* child1_bf_all, const BF* child1_bf_some,
                       const std::string& f_all,  BF*& bf_all,
                       const std::string& f_some, BF*& bf_some);
    virtual std::string compression_name();
//protected: //needed to unprotect them for "uncompress()" (Rayan)
    roaring_bitmap_t *roar;
protected:
    uint64_t size_in_bits;
};


class UncompressedBF : public BF {
public:
    UncompressedBF(const std::string & filename, HashPair hp, int nh, uint64_t num_filter_bits = 0);

    virtual ~UncompressedBF();

    virtual void load(uint64_t partitions = 1);
    virtual void save();

    virtual int operator[](uint64_t pos) const;
    virtual void set_bit(uint64_t p);
    virtual uint64_t num_filter_bits() const;
    virtual uint64_t similarity(const BF* other, int type) const;
    virtual std::tuple<uint64_t, uint64_t> b_similarity(const BF* other) const;
    virtual BF* union_with(const std::string & new_name, const BF* f2) const;
    virtual BF* intersection_with(const std::string & new_name, const BF* f2) const;
    virtual void union_into(const BF* f2);
    virtual BF* append(const std::string & new_name, const BF* f2) const;
    virtual uint64_t count_ones() const;
    virtual void compress_rrr();
    virtual void compress_wraprrr(const BF* f2);
    virtual void compress_roar();
    virtual void uncompress();
    virtual void split(BF* child0_bf_all, const BF* child0_bf_some,
                       BF* child1_bf_all, const BF* child1_bf_some,
                       const std::string& f_all,  BF*& bf_all,
                       const std::string& f_some, BF*& bf_some);

    virtual std::string compression_name();
//protected: // same as BF above (Rayan)
    sdsl::bit_vector* bv;
};


sdsl::bit_vector* union_bv_fast(const sdsl::bit_vector & b1, const sdsl::bit_vector& b2);
int64_t bv_bit_constant(const sdsl::bit_vector* bv);
BF* load_bf_from_file(const std::string & fn, HashPair hp, int nh);

UncompressedBF* uncompress_bf(const BF *in);
roaring_bitmap_t *roar_from_file(std::string filename, uint64_t &size_in_bits);
void roar_to_file(roaring_bitmap_t *roar_bits, uint64_t num_filter_bits, std::string filename);

#endif

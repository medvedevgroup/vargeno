#ifndef BLOOMTREE_H
#define BLOOMTREE_H

#include <string>
#include <queue>
#include <mutex>
#include "Heap.h"
#include "BF.h"
#include "ThreadPool.h"

// this is the max number of BF allowed in memory at once.
extern int BF_INMEM_LIMIT;

class BloomTree {
public:
    BloomTree(const std::string & f, HashPair hp, int nh);
    virtual ~BloomTree();
    std::string name() const;
    virtual std::string name_all() const; // (used for split nodes)
    virtual bool is_split() const;
    virtual bool is_single_filter() const;

    BloomTree* child(int which) const;
    void set_child(int which, BloomTree* c);
    int num_children() const;
    void set_parent(const BloomTree* p);
    const BloomTree* get_parent() const;
    uint64_t similarity(BloomTree* other, int type) const;
    std::tuple<uint64_t, uint64_t> b_similarity(BloomTree* other) const;
    BF* bf() const;
    virtual BF* bf_all() const; // (used for split nodes)
    BF* bf_inter_leaves() const; // get the bloom filter corresponding to the intersection of all leaves below this node

    BloomTree* union_bloom_filters(const std::string & new_name, BloomTree* f2);
    void union_into(const BloomTree* other);

    int usage() const;
    void increment_usage() const;
    static void protected_cache(bool b);
    static bool is_cache_protected();
    static int cache_size();

    // additions w.r.t original SBT implementation (Rayan)
    static std::mutex btLock; // may or may not be used in the current implementation
    static bool caching;
    static ThreadPool pool;

    // RSH: subclass SplitBloomTree needs access to several of these
//private:
    std::string compression_name();
    bool load() const;
    void unload() const;

    static Heap<const BloomTree> bf_cache;
    static void drain_cache();

    std::string filename;
    std::string filename_all; // (used for split nodes)
    HashPair hashes;
    int num_hash;
    mutable BF* bloom_filter;
    mutable BF* bloom_filter_all; // (used for split nodes)
    mutable BF* bloom_filter_inter_leaves;
    mutable Heap<const BloomTree>::heap_reference* heap_ref;

    BloomTree* children[2];
    BloomTree* parent;
    mutable int usage_count;
    mutable bool dirty;
};


class SplitBloomTree : public BloomTree {
public:
    SplitBloomTree(const std::string & f_all, const std::string & f_some, HashPair hp, int nh);
    virtual ~SplitBloomTree();

    virtual std::string name_all() const;
    virtual bool is_split() const;
    virtual bool is_single_filter() const;
    virtual BF* bf_all() const;
};


class AllsomeBloomTree : public BloomTree {
public:
    AllsomeBloomTree(const std::string & f, HashPair hp, int nh);
    virtual ~AllsomeBloomTree();

    virtual bool is_split() const;
    virtual bool is_single_filter() const;
};


SplitBloomTree* split_bt(BloomTree* root);
void rebuild_bt(BloomTree* root);
HashPair* get_hash_function(const std::string & matrix_file, int & nh);
BloomTree* read_bloom_tree(const std::string & filename, bool read_hashes=true);
void write_bloom_tree(const std::string & outfile, BloomTree* root, const std::string & matrix_file);
void write_compressed_bloom_tree(const std::string & outfile, const std::string & compression_name, BloomTree* root, const std::string & matrix_file);
void write_split_bloom_tree(const std::string & outfile, const std::string & compression_name, BloomTree* root, const std::string & matrix_file);
void dump_bloom_tree_bit_vector_stats(BloomTree* root, int level = 0);
void dump_bloom_tree_bit_vectors(BloomTree* root, int spacing = 0);
#endif

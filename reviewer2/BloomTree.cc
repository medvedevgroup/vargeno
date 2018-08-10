#include "BloomTree.h"
#include "util.h"
#include "BF.h"
#include "gzstream.h"

#include <fstream>
#include <list>
#include <cassert>
#include <jellyfish/file_header.hpp>

std::mutex BloomTree::btLock;
Heap<const BloomTree> BloomTree::bf_cache;
int BF_INMEM_LIMIT = 1;
bool BloomTree::caching;
ThreadPool BloomTree::pool(8);

// construct a bloom filter with the given filter backing.
BloomTree::BloomTree(
    const std::string & f,
    HashPair hp,
    int nh
) :
    filename(f),
    filename_all(""),
    hashes(hp),
    num_hash(nh),
    bloom_filter(0),
    bloom_filter_all(0),
    heap_ref(nullptr),
    parent(0),
    usage_count(0),
    dirty(false)
{
    children[0] = nullptr;
    children[1] = nullptr;
    //std::cerr << "BTTRACK BloomTree(" << f << ")" << std::endl;
}

// free the memory for this node.
// RSH: but not for its children(!)
BloomTree::~BloomTree() {
    unload();
}

std::string BloomTree::name() const {
    return filename;
}

std::string BloomTree::name_all() const {
    return "";  // only split subsclasses have a non-empty second name
}

bool BloomTree::is_split() const {
    return false;  // because we have only one backing file
}

bool BloomTree::is_single_filter() const {
    return true;  // because we have only one bloom filter
}

// Return the node for the given child
BloomTree* BloomTree::child(int which) const {
    assert(which >= 0 || which < 2);
    return children[which];
}

// Set the given child
void BloomTree::set_child(int which, BloomTree* c) {
    assert(which >= 0 || which < 2);
    if (c != nullptr) c->parent = this;
    children[which] = c;
}

int BloomTree::num_children() const {
    return ((children[0]==nullptr)?0:1) + ((children[1]==nullptr)?0:1);
}

const BloomTree* BloomTree::get_parent() const {
    return this->parent;
}

void BloomTree::set_parent(const BloomTree* p) {
    parent = const_cast<BloomTree*>(p);
}

// return the bloom filter, loading first if necessary
BF* BloomTree::bf() const {
    load();
    return bloom_filter;
}

BF* BloomTree::bf_all() const {
    return nullptr;  // only split subsclasses have a second BF populated
}

BF* BloomTree::bf_inter_leaves() const {
    // $$$ this may be incomplete; as of Oct/2016 we probably don't need it anyway
    // stripped down copy of load()
    if (bloom_filter_inter_leaves)
    {
        // strip the .roar extension
        size_t lastindex = filename.find_last_of(".");
        std::string rawname = filename.substr(0, lastindex);

        // read the BF file and set bloom_filter
        bloom_filter_inter_leaves = load_bf_from_file(rawname + ".inter_leaves.roar", hashes, num_hash);
        bloom_filter_inter_leaves->load();
    }

    return bloom_filter_inter_leaves;
}


// return the number of times this bloom filter has been used.
int BloomTree::usage() const {
    return usage_count;
}

// increment the usage counter, and update the filter's location in
// the heap if needed.
void BloomTree::increment_usage() const {
    if (!caching) return;

    //std::lock_guard<std::mutex> lock(BloomTree::btLock); // actually put that same lock guard below in compress_roar
    usage_count++;
    // if we're in the cache, let the cache know we've been used.
    if (heap_ref != nullptr) {
        bf_cache.increase_key(heap_ref, usage_count);
    }
}

std::string BloomTree::compression_name() {
    size_t suffix_ix = filename.rfind(".bf.bv");
    if (suffix_ix == std::string::npos) return "none";
    suffix_ix += 6;
    if (suffix_ix == filename.size()) return "none";
    if (filename[suffix_ix] == '.') suffix_ix += 1;
    return filename.substr(suffix_ix);
    }

// Loads the bloom filtering into memory
// nota bene: this always returns true -- what was the return value supposed to mean?
bool BloomTree::load() const {
    bool bf_resident     = (bloom_filter != nullptr);
    bool bf_all_resident = (filename_all.empty()) || (bloom_filter_all != nullptr);
    if (bf_resident && bf_all_resident) {
        //std::cerr << "BTTRACK load BT " << filename << " already resident" << std::endl;
        return true;
    }

    //std::cerr << "BTTRACK load BT " << filename << std::endl;

    // figure out how many partitions are in the BF

    uint64_t filter_partitions = 1;
    if ((!is_single_filter()) && (!is_split())) filter_partitions = 2;

    // if the cache isn't protected from deleting elements, remove enough
    // elements so that there is 1 cache spot free (if the cache is
    // protected, we're allowed to go over the cache limit)
    if (caching && (!bf_cache.is_protected())) BloomTree::drain_cache();

    if (!bf_resident) {
        // read the BF file and set bloom_filter
        //std::cerr << "BTTRACK loading BF " << filename << std::endl;
        bloom_filter = load_bf_from_file(filename, hashes, num_hash);
        bloom_filter->load(filter_partitions);
    }

    if (!bf_all_resident) {
        // read the BF file and set bloom_filter
        //std::cerr << "BTTRACK loading BF " << filename_all << std::endl;
        bloom_filter_all = load_bf_from_file(filename_all, hashes, num_hash);
        bloom_filter_all->load(1);
    }

    if (caching) heap_ref = bf_cache.insert(this, usage());
    dirty = false;

    // since we had to load, we bump up the usage to make it less likely we
    // load again in the near future.
    increment_usage();

    return true;
}

// Frees the memory associated with the bloom filter
void BloomTree::unload() const {
    // you can't unload something until you remove it from the cache
    //std::cerr << "BTTRACK Unloading BT: " << filename << std::endl;

    // free the memory
    if (bloom_filter != nullptr) {
        //std::cerr << "BTTRACK Unloading BF: " << filename << std::endl;
        if (dirty) {
            bloom_filter->save();
        }
        delete bloom_filter;
        bloom_filter = nullptr;
    }

    if (bloom_filter_all != nullptr) {
        //std::cerr << "BTTRACK Unloading BF: " << filename_all << std::endl;
        if (dirty) {
            bloom_filter_all->save();
        }
        delete bloom_filter_all;
        bloom_filter_all = nullptr;
    }

    dirty = false;
}

void BloomTree::drain_cache() {
    if (!caching) return;

    // while the cache is too big, pop an item and unload it
    while (bf_cache.size() >= BF_INMEM_LIMIT && !bf_cache.is_protected()) {
        // toss the bloom filter with the lowest usage
        const BloomTree* loser = bf_cache.pop();
        //std::cerr << "BTTRACK loser BT: " << loser->filename
        //          << " cache size = " << bf_cache.size() << std::endl;
        loser->heap_ref = nullptr;
        loser->unload();
    }
}

void BloomTree::protected_cache(bool b) {
    if (!caching) return;
    bf_cache.set_protected(b);
    if (!b) {
        BloomTree::drain_cache();
    }
}

bool BloomTree::is_cache_protected() {
    if (!caching) {
        return true;
    } else {
        return bf_cache.is_protected();
    }
}

int BloomTree::cache_size() {
    if (!caching) {
        return 0;
    } else {
        return bf_cache.size();
    }
}

uint64_t BloomTree::similarity(BloomTree* other, int type) const {
    // $$$ this doesn't properly consider split nodes
    protected_cache(true);
    uint64_t sim = this->bf()->similarity(other->bf(), type);
    protected_cache(false);
    return sim;
}

std::tuple<uint64_t,uint64_t> BloomTree::b_similarity(BloomTree* other) const {
    // $$$ this doesn't properly consider split nodes
    protected_cache(true);
    std::cerr << "Before \n";
    std::tuple<uint64_t, uint64_t> sim = this->bf()->b_similarity(other->bf());
    protected_cache(false);
    std::cerr << "After \n";
    return sim;
}


// Create a new node that is the union of the bloom filters
// in two other nodes;
BloomTree* BloomTree::union_bloom_filters(const std::string & new_name, BloomTree* f2) {
    // $$$ this doesn't properly consider split nodes

    // move the union op into BloomTree?
    BloomTree* bt = new BloomTree(new_name, hashes, num_hash);

    protected_cache(true);
    bt->bloom_filter = bf()->union_with(new_name, f2->bf());

    bt->set_child(0, this);
    bt->set_child(1, f2);
    //bf_cache.insert(bt, bt->usage());
    bt->dirty = true;
    bt->unload();

    protected_cache(false);
    return bt;
}

void BloomTree::union_into(const BloomTree* other) {
    protected_cache(true);
    bf()->union_into(other->bf());
    dirty = true;
    protected_cache(false);
}

/*
BloomTree* create_union_node(BloomTree* T, BloomTree* N) {
    assert(T != nullptr && N != nullptr);
    assert(T->hashes == N->hashes);
    assert(T->num_hashes == N->num_hashes);

    std::string new_name = "union";

    BloomTree* unionNode = new BloomTree(new_name, T->hashes, T->num_hash);
    unionNode->bloom_filter = T->bloom_filter->union_with(new_name, N->bloom_filter);
    unionNode->set_child(0, T);
    unionNode->set_child(1, N);
    return unionNode;
}
*/

/*============================================*/

// construct a split bloom filter with the given filter-pair backing, both names
// provided
SplitBloomTree::SplitBloomTree(
    const std::string & f_all,
    const std::string & f_some,
    HashPair hp,
    int nh
) :
    BloomTree(f_some, hp, nh)
{
    filename_all = f_all;
    //std::cerr << "BTTRACK SplitBloomTree(" << f_all << ","  << f_some << ")" << std::endl;
}

// free the memory for this node.
SplitBloomTree::~SplitBloomTree() {
    // call to base destructor happens automatically
}

std::string SplitBloomTree::name_all() const {
    return filename_all;
}

bool SplitBloomTree::is_split() const {
    return true;  // because we have two backing files
}

bool SplitBloomTree::is_single_filter() const {
    return false;  // because we have two bloom filters
}

BF* SplitBloomTree::bf_all() const {
    load();
    return bloom_filter_all;
}

/*============================================*/

// construct a bloom filter with the given one-file filter-pair backing
AllsomeBloomTree::AllsomeBloomTree(
    const std::string & f,
    HashPair hp,
    int nh
) :
    BloomTree(f, hp, nh)
{
    //std::cerr << "BTTRACK AllsomeBloomTree(" << f << ")" << std::endl;
}

// free the memory for this node.
AllsomeBloomTree::~AllsomeBloomTree() {
    // call to base destructor happens automatically
}

bool AllsomeBloomTree::is_split() const {
    return false;  // because we have only one backing file
}
bool AllsomeBloomTree::is_single_filter() const {
    return false;  // because we have two bloom filters (but stored in one vector)
}

/*============================================*/

static unsigned split_bt_count;

SplitBloomTree* split_bt_worker (BloomTree* root) {
    // note that the incoming tree is destroyed
    if (root == nullptr) return nullptr;

    //DEBUG
    split_bt_count++;
    std::cerr << "splitting node #" << split_bt_count << ", " << root->name() << std::endl;

    // split the subtrees (if they exist)

    SplitBloomTree* child0 = nullptr;
    if (root->child(0) != nullptr) child0 = split_bt_worker(root->child(0));

    SplitBloomTree* child1 = nullptr;
    if (root->child(1) != nullptr) child1 = split_bt_worker(root->child(1));

    // create a split node that will take the place of the old root, but as yet
    // unpopulated

    std::string root_name = root->name();
    size_t bv_ix = root_name.rfind(".bv");
    if (bv_ix == std::string::npos) {
        bv_ix = root_name.length();
    }

    std::string filename_some = root_name.substr(0,bv_ix) + "-some" + root_name.substr(bv_ix);
    std::string filename_all  = root_name.substr(0,bv_ix) + "-all"  + root_name.substr(bv_ix);

    SplitBloomTree* new_root = new SplitBloomTree(filename_all, filename_some, root->hashes, root->num_hash);

    // make sure the children are resident; note that we have to do this
    // *after* both subtrees have been split, otherwise we'd risk the children
    // being unloaded from memory, which would invalidate the bf pointers

    BF* child0_bf_all  = nullptr;
    BF* child0_bf_some = nullptr;
    if (child0 != nullptr) {
        child0_bf_all  = child0->bf_all();
        child0_bf_some = child0->bf();
    }

    BF* child1_bf_all  = nullptr;
    BF* child1_bf_some = nullptr;
    if (child1 != nullptr) {
        child1_bf_all  = child1->bf_all();
        child1_bf_some = child1->bf();
    }

    bool was_protected = BloomTree::is_cache_protected();
    BloomTree::protected_cache(true);

    // create a new split filter from the children's filters (the children are
    // split nodes);  note that if there are no children a simple split filter
    // is created from the old root's filter's bits
    // $$$ when we call split() below, the call to root->bf()loads root->bv;
    //     we really only need root-bv if root is a leaf, so there's some waste
    //     in loading it for internal nodes

    BF* new_bf_all;
    BF* new_bf_some;
    root->bf()->split(child0_bf_all, child0_bf_some, child1_bf_all, child1_bf_some,
                      filename_all, new_bf_all, filename_some, new_bf_some);

    if (child0 != nullptr) {
        //std::cerr << "force unload " << child0->filename << std::endl;
        child0->unload();   // this node was never involved in the node cache
    }

    if (child1 != nullptr) {
        //std::cerr << "force unload " << child1->filename << std::endl;
        child1->unload();   // this node was never involved in the node cache
    }

    // hook up the children to the new root, and get rid of the old root

    new_root->bloom_filter_all = new_bf_all;   // $$$ change the way this is set
    new_root->bloom_filter     = new_bf_some;  // $$$ change the way this is set
    new_root->dirty = true;
    new_root->set_child(0,child0);
    new_root->set_child(1,child1);

    BloomTree::protected_cache(was_protected);

    //DEBUG
    //std::cerr << "new SplitBloomTree([" << new_root << "]) named \"" << new_root->name() << "\""
    //          << " replaces [" << root  << "] " << root->name() << std::endl;
    //std::cerr << std::endl;

// $$$ maybe it has entries in the cache?
//  delete root;

    return new_root;
}

SplitBloomTree* split_bt (BloomTree* root) {
    // note that the incoming tree is destroyed

    DIE_IF(!root->is_single_filter(), "splitting of trees with multi-filter nodes is not supported");

    BloomTree::caching = true; // make sure bloom filter caching is on
    if (BF_INMEM_LIMIT < 2)    // .. and that we will cache enough nodes
        BF_INMEM_LIMIT = 2;

    split_bt_count = 0;
    assert (!BloomTree::is_cache_protected());
    SplitBloomTree* new_root = split_bt_worker(root);

    if (new_root != nullptr) {
        //std::cerr << "force unload " << new_root->filename << std::endl;
        new_root->unload();   // this node was never involved in the node cache
    }

    return new_root;
}

/*============================================*/

static unsigned rebuild_bt_count;

void rebuild_bt_worker (BloomTree* root) {
    // don't bother to rebuild leaves

    if (root == nullptr) return;

    BloomTree* child0 = root->child(0);
    BloomTree* child1 = root->child(1);

    if ((child0 == nullptr) && (child0 == nullptr))
        return;

    // rebuild the subtrees (if they exist)

    //DEBUG
    rebuild_bt_count++;
    std::cerr << "rebuilding internal node #" << rebuild_bt_count << ", " << root->name() << std::endl;

    if (child0 != nullptr) rebuild_bt_worker(child0);
    if (child1 != nullptr) rebuild_bt_worker(child1);

    bool was_protected = BloomTree::is_cache_protected();
    BloomTree::protected_cache(true);

    // union the children to form this node

    assert (root->bloom_filter == nullptr);
    BF* union_bf = child0->bf()->union_with(root->name(), child1->bf());
    union_bf->save();
    delete union_bf;

    // we're done with the children, unload them

    if (child0 != nullptr) {
        //std::cerr << "force unload " << child0->filename << std::endl;
        child0->unload();
    }

    if (child1 != nullptr) {
        //std::cerr << "force unload " << child1->filename << std::endl;
        child1->unload();
    }

    BloomTree::protected_cache(was_protected);
}

void rebuild_bt (BloomTree* root) {
    DIE_IF(!root->is_single_filter(), "rebuilding of trees with multi-filter nodes is not supported");
    if (root->compression_name() == "rrr")
        std::cerr << "WARNING: rebuilding an rrr tree can be very slow" << std::endl;

    BloomTree::caching = true; // make sure bloom filter caching is on
    if (BF_INMEM_LIMIT < 3)    // .. and that we will cache enough nodes
        BF_INMEM_LIMIT = 3;

    assert (!BloomTree::is_cache_protected());
    rebuild_bt_count = 0;
    rebuild_bt_worker(root);
    root->unload();
}

/*============================================*/

HashPair* get_hash_function(const std::string & matrix_file, int & nh) {
    std::cerr << "Loading hashes from " << matrix_file << std::endl;
    igzstream in(matrix_file.c_str(), std::ios::in | std::ios::binary);
    jellyfish::file_header header(in);
    DIE_IF(!in.good(), "Couldn't parse bloom filter header!");
    HashPair * hp = new HashPair(header.matrix(1), header.matrix(2));
    in.close();

    nh = header.nb_hashes();
    std::cerr << "# Hash applications=" << nh << std::endl;

    jellyfish::mer_dna::k(header.key_len() / 2);
    std::cerr << "Read hashes for k=" << jellyfish::mer_dna::k() << std::endl;
    return hp;
}


/* Read a file that defines the bloom tree structure.

   There are (currently) three diffferent types of tree:
     BloomTree:        one filter and one file per node
     SplitBloomTree:   two filters and two files per node ("all" and "some") (it's our previous attempt at allsome)
     AllsomeBloomTree: two filters but one file per node

   The file has lines of the form:
     Root,HashFile
     *Child1
     ***Child3
     ****Child4
     *Child2
   where the '*' indicates the level and where "Root" and "Child1" etc is
   either a single filename or two filenames separated by a comma. The latter
   case indicates a SplitBloomTree tree with nodes split into "all" and "some"
   files (the "all" filename comes first). In the single file case, if the
   filename contains the substring "-allsome" the node is an AllsomeBloomTree;
   otherwise it is a BloomTree.  All nodes in the tree must be the same type.

   If a BF filename (or HashFile) doesn't contain a path, we copy the path
   from the incoming filename.  The idea is that if the structure file doesn't
   contain paths, we assume the bloom filter (and hash function) files are in
   the same directory as the structure file.  This lets the user keep the bloom
   tree self-contained in its own directory, with a structure file that
   contains no paths, while performing operations from another directory.

   This function will return a pointer to the root of the constructed bloom
   tree.
*/

BloomTree* read_bloom_tree(const std::string & filename, bool read_hashes) {
    std::ifstream in(filename.c_str());

    if (!in.good()) {
        std::cerr << ">>> Failed to open \"" << filename << "\" for reading <<<" << std::endl;
        return nullptr;
    }

    std::string file_path = "";
    size_t slashix = filename.find_last_of("/");
    if (slashix != std::string::npos) {
        file_path = filename.substr(0,slashix+1);
    }

    std::list<BloomTree*> path;
    BloomTree* tree_root = 0;
    int n = 0;
    // if read_hashes is false, you must promise never to access the bloom filters
    HashPair* hashes = new HashPair; // useless hashpair used if read_hashes is false
    int num_hashes = 0;
    unsigned int filters_per_node = 0;
    unsigned int files_per_node = 0;
    std::string tree_description = "";

    std::string node_info;
    while (getline(in, node_info)) {
        node_info = Trim(node_info);
        if (node_info.size() == 0) continue;
        size_t level = node_info.find_first_not_of("*");
        node_info.erase(0, level);

        // each node info is a comma separated list
        std::vector<std::string> fields;
        SplitString(node_info, ',', fields);
        if (files_per_node == 0) {
            DIE_IF(fields.size() < 2, "Must specify hash file for root.");
            DIE_IF(fields.size() > 3, "Too many fields in first line of bloom tree file.");
            files_per_node = fields.size() - 1;
        } else {
            DIE_IF(fields.size() != files_per_node, "Inconsistent number of fields in bloom tree file.");
        }

        if (filters_per_node == 0) {
            if ((files_per_node == 2) || (fields[0].find("-allsome") != std::string::npos))
                filters_per_node = 2;
            else
                filters_per_node = 1;
        }

        if (files_per_node == 2) {
            DIE_IF(fields[0].find("-allsome") != std::string::npos, "can't have -allsome nodes in a split bloom tree file.");
            DIE_IF(fields[1].find("-allsome") != std::string::npos, "can't have -allsome nodes in a split bloom tree file.");
        } else if (filters_per_node == 1) {
            DIE_IF(fields[0].find("-allsome") != std::string::npos, "cant' mix -allsome nodes with regular nodes in a bloom tree file.");
        }

        std::string bf_filename = fields[0];
        std::string bf_all_filename = "";
        if ((!file_path.empty()) && (bf_filename.find_last_of("/") == std::string::npos)) {
            bf_filename = file_path + bf_filename;
        }

        if (files_per_node == 2) {
            bf_all_filename = bf_filename;
            bf_filename = fields[1];
            if ((!file_path.empty()) && (bf_filename.find_last_of("/") == std::string::npos)) {
                bf_filename = file_path + bf_filename;
            }
        }

        n++;

        BloomTree* bn = nullptr;

        // if we're at the root
        if (path.size() == 0) {
            DIE_IF(level != 0, "Root must start in column 0");
            DIE_IF(tree_root != 0, "Can't set root twice!");

            // set the hash function up
            if (read_hashes) {
                DIE_IF(fields.size() < files_per_node+1, "Must specify hash file for root.");
                std::string hash_filename = fields[files_per_node];
                if ((!file_path.empty()) && (hash_filename.find_last_of("/") == std::string::npos)) {
                    hash_filename = file_path + hash_filename;
                }
                hashes = get_hash_function(hash_filename, num_hashes);
            }

            // create the root node

            if (files_per_node == 2) {
                tree_description = "split allsome Bloom Tree";
                bn = new SplitBloomTree(bf_all_filename, bf_filename, *hashes, num_hashes);
            } else if (filters_per_node == 2) {
                tree_description = "concatenated allsome Bloom Tree";
                bn = new AllsomeBloomTree(bf_filename, *hashes, num_hashes);
            } else {
                tree_description = "union Bloom Tree";
                bn = new BloomTree(bf_filename, *hashes, num_hashes);
            }
            tree_root = bn;

        // if we're adding a child
        } else {
            if (files_per_node == 2) {
                bn = new SplitBloomTree(bf_all_filename, bf_filename, *hashes, num_hashes);
            } else if (filters_per_node == 2) {
                bn = new AllsomeBloomTree(bf_filename, *hashes, num_hashes);
            } else {
                bn = new BloomTree(bf_filename, *hashes, num_hashes);
            }

            while (path.size() > level) {
                path.pop_back();
            }
            DIE_IF(level != path.size(), "Must increase level by <= 1");

            if (path.back()->child(0) == nullptr) {
                path.back()->set_child(0, bn);
            } else if (path.back()->child(1) == nullptr) {
                path.back()->set_child(1, bn);
            } else {
                DIE("Tried to add >= 2 children to a node.");
            }
        }
        path.push_back(bn);
    }
    delete hashes;

    std::cerr << "Read " << n << " nodes in " << tree_description << std::endl;

    return tree_root;
}

// write the bloom tree file format in a way that can be read by read_bloom_tree()
void write_bloom_tree_helper(std::ostream & out, BloomTree* root, int level=1) {
    std::string lstr(level, '*');

    for (int i = 0; i < 2; i++) {
        if (root->child(i) != nullptr) {
            out << lstr << root->child(i)->name() << std::endl;
            write_bloom_tree_helper(out, root->child(i), level+1);
        }
    }
}

void write_bloom_tree(
    const std::string & outfile,
    BloomTree* root,
    const std::string & matrix_file
) {
    std::cerr << "Writing to " << outfile << std::endl;
    std::ofstream out(outfile.c_str());
    out << root->name() << "," << matrix_file << std::endl;
    write_bloom_tree_helper(out, root);
    std::cerr << "Done writing " << outfile << "." << std::endl;
}

std::string strip_compression_from_name(std::string filename) {
    size_t suffix_ix = filename.rfind(".bv");
    if (suffix_ix == std::string::npos) return filename;
    suffix_ix += 3;
    return filename.substr(0,suffix_ix);
    }

// write the bloom tree file format in a way that can be read by read_bloom_tree()
void write_compressed_bloom_tree_helper(std::ostream & out, const std::string & compression_name, BloomTree* root, int level=1) {
    std::string lstr(level, '*');

    for (int i = 0; i < 2; i++) {
        if (root->child(i) != nullptr) {
            out << lstr << strip_compression_from_name(root->child(i)->name()) << "." << compression_name << std::endl;
            write_compressed_bloom_tree_helper(out, compression_name, root->child(i), level+1);
        }
    }
}

void write_compressed_bloom_tree(
    const std::string & outfile,
    const std::string & compression_name,
    BloomTree* root,
    const std::string & matrix_file
) {
    std::cerr << "Writing to " << outfile << std::endl;

    std::ofstream out(outfile.c_str());
    out << strip_compression_from_name(root->name()) << "." << compression_name << "," << matrix_file << std::endl;
    write_compressed_bloom_tree_helper(out, compression_name, root);
    std::cerr << "Done writing " << outfile << "." << std::endl;
}


// write the bloom tree file format in a way that can be read by read_bloom_tree()
void write_split_bloom_tree_helper(
    std::ostream & out,
    const std::string & compression_ext,
    BloomTree* root,
    int level=1
) {
    std::string lstr(level, '*');

    for (int i = 0; i < 2; i++) {
        if (root->child(i) != nullptr) {
            out << lstr << strip_compression_from_name(root->child(i)->name_all()) << compression_ext
                << "," << strip_compression_from_name(root->child(i)->name()) << compression_ext
                << std::endl;
            write_split_bloom_tree_helper(out, compression_ext, root->child(i), level+1);
        }
    }
}

void write_split_bloom_tree(
    const std::string & outfile,
    const std::string & compression_name,  // empty string implies uncompressed
    BloomTree* root,
    const std::string & matrix_file
) {
    std::cerr << "Writing to " << outfile << std::endl;

    std::string compression_ext;
    if ((compression_name.empty()) || (compression_name == "none")) {
        compression_ext = "";
    } else {
        compression_ext = "." + compression_name;
    }

    std::ofstream out(outfile.c_str());
    out << strip_compression_from_name(root->name_all()) << compression_ext
        << "," << strip_compression_from_name(root->name()) << compression_ext
        << "," << matrix_file << std::endl;
    write_split_bloom_tree_helper(out, compression_ext, root);
    std::cerr << "Done writing " << outfile << "." << std::endl;
}

/*============================================*/

static unsigned dump_bit_vector_stats_count;

void dump_bit_vector_stats(BF* bf, std::string name, bool is_leaf, int level) {
    std::ios::fmtflags save_flags = std::cout.flags();
    int save_precision = std::cout.precision();
    std::cout.setf(std::ios::fixed, std::ios::floatfield);

    unsigned long nb_ones = bf->count_ones();
    double occupancy = ((double) nb_ones) / bf->num_filter_bits();

    std::cout        << std::setw(4) << dump_bit_vector_stats_count;
    std::cout << " " << std::setw(10) << nb_ones;
    std::cout << " " << std::setprecision(6) << std::setw(9) << occupancy;
    if (is_leaf) std::cout << " @ ";
            else std::cout << " . ";
    std::string level_str(level, '*');
    std::cout << level_str <<name << std::endl;

    std::cout.flags(save_flags);
    std::cout.precision(save_precision);
}

void dump_bloom_tree_bit_vector_stats(BloomTree* root, int level) {
    bool is_leaf = (root->num_children() == 0);

    if (level == 0) {
        dump_bit_vector_stats_count = 1;
    } else {
        dump_bit_vector_stats_count++;
    }

    if (root->is_split()) {
        dump_bit_vector_stats (root->bf_all(), root->name_all(), is_leaf, level);
    }
    dump_bit_vector_stats (root->bf(), root->name(), is_leaf, level);

    for (int i = 0; i < 2; i++) {
        if (root->child(i) != nullptr) {
            dump_bloom_tree_bit_vector_stats(root->child(i), level+1);
        }
    }
}

/*============================================*/

static uint64_t max_bit_vector_dump = 512;

// write the bloom tree bit vectors to stdout, for debugging
void dump_bloom_tree_bit_vector_header(BF* bf, int spacing) {
    uint64_t num_bits = bf->num_filter_bits();
    if (num_bits > max_bit_vector_dump) num_bits = max_bit_vector_dump;

    for (size_t i=0; i<num_bits; i++) {
        if ((i != 0) && (i % spacing == 0)) std::cout << " ";
        if (i % 10 == 0) std::cout << ((i/10) % 10);
                    else std::cout << " ";
    }
    std::cout << std::endl;

    for (size_t i=0; i<num_bits; i++) {
        if ((i != 0) && (i % spacing == 0)) std::cout << " ";
        std::cout << (i % 10);
    }
    std::cout << std::endl;
}

void dump_bloom_tree_bit_vector(BF* bf, std::string name, bool is_leaf, int level, int spacing) {
    uint64_t num_bits = bf->num_filter_bits();
    if (num_bits > max_bit_vector_dump) num_bits = max_bit_vector_dump;
    for (size_t i=0; i<num_bits; i++) {
        if ((i != 0) && (i % spacing == 0)) std::cout << " ";
        if ((*bf)[i] == 1) std::cout << "+";
                      else std::cout << "-";
    }
    if (num_bits <  bf->num_filter_bits()) std::cout << " ...";

    if (is_leaf) std::cout << " @ ";
            else std::cout << " . ";
    std::string level_str(level, '*');
    std::cout << level_str <<name << std::endl;
}

void dump_bloom_tree_bit_vectors_helper(BloomTree* root, int level, int spacing) {
    bool is_leaf = (root->num_children() == 0);

    if (root->is_split()) {
        dump_bloom_tree_bit_vector (root->bf_all(), root->name_all(), is_leaf, level, spacing);
    }
    dump_bloom_tree_bit_vector (root->bf(), root->name(), is_leaf, level, spacing);

    for (int i = 0; i < 2; i++) {
        if (root->child(i) != nullptr) {
            dump_bloom_tree_bit_vectors_helper(root->child(i), level+1, spacing);
        }
    }
}

void dump_bloom_tree_bit_vectors(BloomTree* root, int spacing) {
    dump_bloom_tree_bit_vector_header(root->bf(),spacing);
    dump_bloom_tree_bit_vectors_helper(root,0,spacing);
}

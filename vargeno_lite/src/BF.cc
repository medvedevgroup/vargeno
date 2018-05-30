#include "BF.h"
#include "Kmers.h"
#include "sbt_util.h"

#include <iomanip> //for std::setprecision
#include <iostream> //for std::ios::fixed

#include <jellyfish/file_header.hpp>
#include <sdsl/sfstream.hpp>


//#define USE_SIMPLE_TIMER   // if this is defined, we use a second-resolution
//                           // timer;  otherwise, we use something better

#ifdef USE_SIMPLE_TIMER
#define get_wall_time()          time(NULL)
#define elapsed_wall_time(start) (get_wall_time()-start)
#else
#include <chrono>
#define get_wall_time()          std::chrono::system_clock::now()
#define elapsed_wall_time(start) ((std::chrono::duration_cast<std::chrono::nanoseconds>(get_wall_time()-start).count()) / 1000000000.0)
#endif

// File header for wrapped rrr files

struct WrapRrrHeader {
    uint32_t num_bit_vectors;  // (expected to be 2)
    uint32_t padding1;         // (expected to be 0)
    uint64_t num_bv1_bits;     // count of bits in uncompressed version of bit vector 1
    uint32_t bv1_compressor;   // 0=>all zeros; 1=>all ones; 2=>rrr
    uint32_t padding2;         // (expected to be 0)
    uint64_t num_bv2_bits;     // count of bits in uncompressed version of bit vector 2
    uint32_t bv2_compressor;   // 0=>all zeros; 1=>all ones; 2=>rrr
    uint32_t padding3;         // (expected to be 0)
};

/*============================================*/
// Rrr BF

BF::BF(const std::string & f, HashPair hp, int nh) :
    filename(f),
    bits(nullptr),
    hashes(hp),
    num_hash(nh),
    hash_function_modulus(0),
    filter_partitions(1)
{ 
}

BF::~BF() {
    if (bits != nullptr) {
        delete bits;
    }
}

void BF::add(const jellyfish::mer_dna & m) {
    jellyfish::mer_dna can(m);
    can.canonicalize();

    if (num_hash == 1) {
        uint64_t h0 = hashes.m1.times(can);
        const size_t pos = h0 % hash_function_range();
        this->set_bit(pos);
    } else {
        uint64_t h0 = hashes.m1.times(can);
        uint64_t h1 = hashes.m2.times(can);
        const size_t base = h0 % hash_function_range();
        const size_t inc = h1 % hash_function_range();

        for (unsigned long i = 0; i < num_hash; ++i) {
            const size_t pos = (base + i * inc) % hash_function_range();
            this->set_bit(pos);
        }
    }
}

// set the bit at position 1 (can't use operator[] b/c we'd need
// to return a proxy, which is more trouble than its worth)
void BF::set_bit(uint64_t p) {
    DIE("Compressed BF are not mutable!");
}

// returns true iff the bloom filter contains the given kmer
bool BF::contains(const jellyfish::mer_dna & m) const {
    //std::cout << "BF::contains( \"" << m.to_str() << "\")" << std::endl;
    std::string temp = m.to_str();
    jellyfish::mer_dna n = jellyfish::mer_dna(temp);
    n.canonicalize();

    if (num_hash == 1) {
        uint64_t h0 = hashes.m1.times(n);
        const size_t pos = h0 % hash_function_range();
        //DEBUG
        //if ((*this)[pos] != 0) { std::cout << "testing pos " << pos << std::endl; }
        //                  else { std::cout << "testing pos " << pos << ", absent" << std::endl; }
        if ((*this)[pos] == 0) return false;
    } else {
        uint64_t h0 = hashes.m1.times(n);
        uint64_t h1 = hashes.m2.times(n);

        const size_t base = h0 % hash_function_range();
        const size_t inc = h1 % hash_function_range();

        for (unsigned long i = 0; i < num_hash; ++i) {
            const size_t pos = (base + i * inc) % hash_function_range();
            //DEBUG: std::cerr << "pos=" << pos << std::endl;
            if ((*this)[pos] == 0) return false;
        }
    }
    return true;
}

// convenience function
bool BF::contains(const std::string & str) const {
    //jellyfish::mer_dna temp = jellyfish::mer_dna(str);
    //temp.canonicalize();
    //std::cout << "checking: " << str << "\n";
    //std::cout << "canonical: " << temp.to_str() << "\n";
    return contains(jellyfish::mer_dna(str));
}

// convenience function
bool BF::contains(const uint64_t pos) const {
    //std::cout << "BF::contains( " << pos << ")" << std::endl;
    return ((*this)[pos] != 0);
}

uint64_t BF::num_filter_bits() const {
    return bits->size();
}

uint64_t BF::hash_function_range() const {
    return hash_function_modulus;
}

void BF::update_hash_function_range() {
    // we assume (without checking) that num_filter_bits is a multiple of filter_partitions
    hash_function_modulus = num_filter_bits() / filter_partitions;
}

void BF::set_filter_partitions(uint64_t partitions) {
    filter_partitions = partitions;
}

int BF::operator[](uint64_t pos) const {
    //std::cout << "BF::[" << pos << "]" << std::endl;
    return (*bits)[pos];
}

// returns the position of a given kmer in the bloom filter
uint64_t BF::kmer_position(const jellyfish::mer_dna & m) const {
// $$$ this only supports the first hash function!
    std::string temp = m.to_str();
    jellyfish::mer_dna n = jellyfish::mer_dna(temp);
    n.canonicalize();

    assert(num_hash == 1);
    uint64_t h0 = hashes.m1.times(n);
    return h0 % hash_function_range();
}

// read the bit vector and the matrices for the hash functions.
void BF::load(uint64_t partitions) {
    //std::cerr << "Loading rrr BF from " << filename << std::endl;
    // read the actual bits
    bits = new sdsl_rrr_vector();
    auto start_time = get_wall_time();
    sdsl::load_from_file(*bits, filename);
    auto elapsed_time = elapsed_wall_time(start_time);
    uint64_t fsize= filesize(filename);
    std::cerr << std::setiosflags(std::ios::fixed) << std::setprecision(4)  << "[BF load] " << elapsed_time << " secs " << filename << " " << fsize/1024/1024 << " MB, " << fsize / elapsed_time /1024/1024 << " MB/sec" << std::endl;

    set_filter_partitions(partitions);
    update_hash_function_range();
}

void BF::save() {
    //std::cerr << "Saving compressed BF to " << filename << std::endl;
    auto start_time = get_wall_time();
    sdsl::store_to_file(*bits, filename);
    auto elapsed_time = elapsed_wall_time(start_time);
    std::cerr << "[BF save] " << elapsed_time << " secs " << filename << std::endl;
}


// create a new RRR bloom filter that is the union of this BF and the given BF.
// Will re-use the hashes from this and both BFs must use exactly the same hash
// (not checked).
BF* BF::union_with(const std::string & new_name, const BF* f2) const {
    assert(bits->size() == f2->num_filter_bits());

    // create an uncompressed version of this bf
    sdsl::bit_vector b(bits->size(), 0);

    std::cerr << "Performing OR... (num_filter_bits " << b.size() << ")" << std::endl;

    // union it with f2
    for (unsigned long i = 0; i < b.size(); i++) {
        b[i] = (*bits)[i] | (*f2->bits)[i];
        if (i % 1000000 == 0) std::cerr << "i=" << i << std::endl;
    }

    // create a new BF wrapper for the new BF
    std::cerr << "Building BF object..." << std::endl;
    BF* out = new BF(new_name, hashes, num_hash);

    // "load" the new BF by converting it to a RRR
    std::cerr << "Building RRR vector..." << std::endl;
    out->bits = new sdsl_rrr_vector(b);
    return out;
}

BF* BF::intersection_with(const std::string & new_name, const BF* f2) const {
    DIE("RRR intersection_with not yet implemented");
    return nullptr;
}

// modified by Rayan
// used to be not implemented, I've added a dynamic uncompressed step to support it
uint64_t BF::similarity(const BF* other, int type) const {
    UncompressedBF *u_this = uncompress_bf(this);
    UncompressedBF *u_other = uncompress_bf(other);
    return u_this->similarity(u_other,type);
}

std::tuple<uint64_t, uint64_t> BF::b_similarity(const BF* other) const {
    DIE("b_similarity not yet implemented");
    return std::make_tuple(0,0);
}

void BF::union_into(const BF* other) {
    DIE("union_into not yet implemented");
}

BF* BF::append(const std::string & new_name, const BF* f2) const {
    DIE("appending to RRR not yet implemented");
    return nullptr;
}

uint64_t BF::count_ones() const {
    DIE("count_ones not yet implemented for rrr");
    return 0;
}

void BF::compress_rrr() {
    DIE("Cant compress rrr further with existing code base");
}

void BF::compress_wraprrr(const BF* f2) {
    DIE("Cant compress rrr w/ wraprrr with existing code base");
}

void BF::compress_roar() {
    DIE("Cant compress rrr w/ roar further with existing code base");
}


void BF::uncompress() {
    auto t = sdsl::bit_vector(bits->size(), 0);
    std::cerr << "uncompressing RRR bitvector " << filename
              << " (" << bits->size() << " bits)" << std::endl;
    //unsigned long sz = bits->size();
    //for (unsigned long i = 0; i <  sz; i++) 
    //    t[i] = (*bits)[i];
    /*
    auto it = bits->begin(); 
    for (unsigned long i = 0; i <  sz; i++) {
        if (*it != (*bits)[i])
            std::cerr << "boo!" << std::endl;
        it++;
    }*/

    unsigned long i = 0;
    for (auto it = bits->begin(); it != bits->end(); it++)
    {
        t[i] = (*it);
        i++;
    }

    // strip the .rrr extension
    size_t lastindex = filename.find_last_of("."); 
    std::string rawname = filename.substr(0, lastindex); 
    std::cerr << "Compressed RRR vector would be " << sdsl::size_in_mega_bytes(*bits) << "MB, storing uncompressed to " << rawname << std::endl;

    sdsl::store_to_file(t,rawname);
    std::cerr << "done" << std::endl;
}

void BF::split(BF* child0_bf_all, const BF* child0_bf_some,
               BF* child1_bf_all, const BF* child1_bf_some,
               const std::string& f_all,  BF*& _bf_all,
               const std::string& f_some, BF*& _bf_some) {
    DIE("Can't split an rrr-compressed BF");
}

std::string BF::compression_name() {
    return "rrr";
}


/*============================================*/
// WrapRrr BF -- wrapper around RRR, storing multiple RRR vectors in one file

WrapRrrBF::WrapRrrBF(const std::string & f, HashPair hp, int nh, uint64_t num_filter_bits) :
    BF(f, hp, nh),
    bits_constant(0),
    bits2_constant(0),
    bits2(nullptr),
    total_filter_bits(0)
{
    //std::cerr << "Constructing wraprrr BF " << f << std::endl;
}

WrapRrrBF::~WrapRrrBF() {
    // call to base destructor happens automatically
    if (bits2 != nullptr) {
        delete bits2;
    }
}

void WrapRrrBF::load(uint64_t partitions) {
    // $$$ add load time report

    //std::cerr << "Loading wraprrr BF from " << filename << std::endl;
    DIE_IF(partitions!=2, "wraprrr is not implemented for anything other than partitions" );

    WrapRrrHeader header;

    sdsl::isfstream in(filename, std::ios::binary | std::ios::in);
    DIE_IF(!in, "failed to load " + filename);

    auto start_time = get_wall_time();
    in.read((char*)&header, sizeof(header));
    auto elapsed_time = elapsed_wall_time(start_time);

    //std::cerr << "  num_bit_vectors = " << header.num_bit_vectors << std::endl;
    //std::cerr << "  num_bv1_bits    = " << header.num_bv1_bits << std::endl;
    //std::cerr << "  num_bv2_bits    = " << header.num_bv2_bits << std::endl;
    //std::cerr << "  bv1_compressor  = " << header.bv1_compressor << std::endl;
    //std::cerr << "  bv2_compressor  = " << header.bv2_compressor << std::endl;
    DIE_IF(header.num_bv1_bits != header.num_bv2_bits, "mismatched filter sizes in " + filename);

    if (header.bv1_compressor < 2) {
        bits = nullptr;
        bits_constant = header.bv1_compressor;
    } else if (header.bv1_compressor == 2) {
        bits = new sdsl_rrr_vector();
        auto start_time2 = get_wall_time();
        sdsl::load(*bits,in);
        elapsed_time += elapsed_wall_time(start_time2);
        //std::cerr << "  loaded " << bits->size() << " bit filter for -all component" << std::endl;
    } else {
        DIE("bad compression identifier for filter 1 in " + filename);
    }

    if (header.bv2_compressor < 2) {
        bits2 = nullptr;
        bits2_constant = header.bv2_compressor;
    } else if (header.bv2_compressor == 2) {
        bits2 = new sdsl_rrr_vector();
        auto start_time2 = get_wall_time();
        sdsl::load(*bits2,in);
        elapsed_time += elapsed_wall_time(start_time2);
        //std::cerr << "  loaded " << bits2->size() << " bit filter for -some component" << std::endl;
    } else {
        DIE("bad compression identifier for filter 2 in " + filename);
    }

    in.close();

    std::cerr << std::setiosflags(std::ios::fixed) << std::setprecision(4)  << "[BF load] " << elapsed_time << " secs " << filename << " " << std::endl;

    total_filter_bits = header.num_bv1_bits + header.num_bv2_bits;
    set_filter_partitions(2);
    update_hash_function_range();
}

void WrapRrrBF::save() {
    DIE("saving WrapRRR not yet implemented");
}

uint64_t WrapRrrBF::num_filter_bits() const {
    return total_filter_bits;
}

int WrapRrrBF::operator[](uint64_t pos) const {
    //std::cout << "WrapRrrBF::[" << pos << "]" << std::endl;
    if (pos < hash_function_modulus) {
        if (bits != nullptr) return (*bits)[pos];
        return bits_constant;
    } else {
        if (bits2 != nullptr) return (*bits2)[pos-hash_function_modulus];
        return bits2_constant;
    }
}

BF* WrapRrrBF::append(const std::string & new_name, const BF* f2) const {
    DIE("appending to WrapRRR not yet implemented");
    return nullptr;
}

void WrapRrrBF::compress_rrr() {
    DIE("Cant compress wraprrr w/ rrr with existing code base");
}

void WrapRrrBF::compress_wraprrr(const BF* f2) {
    DIE("Cant compress wraprrr further with existing code base");
}

void WrapRrrBF::compress_roar() {
    DIE("Cant compress wraprrr w/ roar with existing code base");
}

void WrapRrrBF::uncompress() {
    DIE("Cant uncompress wraprrr with existing code base");
}

std::string WrapRrrBF::compression_name() {
    return "wraprrr";
}

/*============================================*/

UncompressedBF::UncompressedBF(const std::string & f, HashPair hp, int nh, uint64_t num_filter_bits) :
    BF(f, hp, nh),
    bv(nullptr)
{
    hash_function_modulus = num_filter_bits;
    if (num_filter_bits > 0) {
        bv = new sdsl::bit_vector(num_filter_bits);
    }
}


UncompressedBF::~UncompressedBF() {
    // call to base destructor happens automatically
    if (bv != nullptr) {
        delete bv;
    }
}

void UncompressedBF::load(uint64_t partitions) {
    // read the actual bits
    assert(bv == nullptr);
    // std::cerr << "Loading uncompressed BF from " << filename << std::endl;
    bv = new sdsl::bit_vector();
    auto start_time = get_wall_time();
    sdsl::load_from_file(*bv, filename);
    auto elapsed_time = elapsed_wall_time(start_time);
    std::cerr << "[BF load] " << elapsed_time << " secs " << filename << std::endl;

    set_filter_partitions(partitions);
    update_hash_function_range();
}

void UncompressedBF::save() {
    //std::cerr << "Saving uncompressed BF to " << filename << std::endl;
    auto start_time = get_wall_time();
    sdsl::store_to_file(*bv, filename);
    auto elapsed_time = elapsed_wall_time(start_time);
    std::cerr << "[BF save] " << elapsed_time << " secs " << filename << std::endl;
}

uint64_t UncompressedBF::num_filter_bits() const {
    return bv->size();
}

int UncompressedBF::operator[](uint64_t pos) const {
    //std::cout << "UncompressedBF::[" << pos << "]" << std::endl;
    return (*bv)[pos];
}

void UncompressedBF::set_bit(uint64_t p) {
    (*bv)[p] = 1;
}

BF* UncompressedBF::union_with(const std::string & new_name, const BF* f2) const {
    assert(num_filter_bits() == f2->num_filter_bits());
    const UncompressedBF* b = dynamic_cast<const UncompressedBF*>(f2);
    if (b == nullptr) {
        DIE("Can only union two uncompressed BF");
    }
    UncompressedBF* out = new UncompressedBF(new_name, hashes, num_hash);
    out->bv = union_bv_fast(*this->bv, *b->bv);
    return out;
}

BF* UncompressedBF::intersection_with(const std::string & new_name, const BF* f2) const {
    DIE("Uncompressed intersection_with not implemented");
    return nullptr;
}

void UncompressedBF::union_into(const BF* f2) {
    assert(num_filter_bits() == f2->num_filter_bits());

    const UncompressedBF* b = dynamic_cast<const UncompressedBF*>(f2);
    if (b == nullptr) {
        DIE("Can only union two uncompressed BF");
    }

    uint64_t* b1_data = this->bv->data();
    const uint64_t* b2_data = b->bv->data();

    // RSH: using num_filter_bits() instead of capacity() misses the bits in
    // .. the final chunk if size%64 is non-zero
    // .. replaced sdsl::bit_vector::size_type len = num_filter_bits()>>6;
    sdsl::bit_vector::size_type len = this->bv->capacity()>>6;
    for (sdsl::bit_vector::size_type p = 0; p < len; ++p) {
        *b1_data = (*b1_data) | (*b2_data++);
        b1_data++;
    }
}

BF* UncompressedBF::append(const std::string & new_name, const BF* f2) const {
    const UncompressedBF* b2 = dynamic_cast<const UncompressedBF*>(f2);
    if (b2 == nullptr) {
        DIE("Can only append two uncompressed BF");
    }

    uint64_t new_filter_bits = num_filter_bits() + b2->num_filter_bits();
    UncompressedBF* out = new UncompressedBF(new_name, hashes, num_hash, new_filter_bits);

    // copy the first bit vector

    const uint64_t* b1_data = bv->data();
    const uint64_t* b2_data = b2->bv->data();
    uint64_t* b3_data = out->bv->data();

    sdsl::bit_vector::size_type len1 = bv->capacity()>>6;
    sdsl::bit_vector::size_type len2 = b2->bv->capacity()>>6;
    sdsl::bit_vector::size_type len3 = out->bv->capacity()>>6;
    sdsl::bit_vector::size_type copyLen2 = len3 - len1;

    for (sdsl::bit_vector::size_type p = 0; p < len1; ++p) {
        *(b3_data++) = *(b1_data++);
    }

    // copy the second bit vector, shifting bits if necessary
    // RSH: this makes serious assumptions about how sdsl stores bit vectors,
    // specifically the order of bits within the vectors, that the vectors are
    // 64-bits chunks, and that the first bit into a 64-bit chunk is the
    // least-significant bit

    uint32_t fillShift = bv->capacity() - num_filter_bits();
    assert(fillShift < 64);
    if (fillShift == 0) {
        for (sdsl::bit_vector::size_type p = 0; p < len1; ++p) {
            *(b3_data++) = *(b2_data++);
        }
    } else {
        b3_data--;
        for (sdsl::bit_vector::size_type p = 0; p < copyLen2; ++p) {
            uint64_t chunk = *(b2_data++);
            *(b3_data++) |= (chunk << (64-fillShift));
            *b3_data     =  (chunk >> fillShift);
        }
        if (copyLen2 < len2) {
            uint64_t chunk = *b2_data;
            *b3_data |= (chunk << (64-fillShift));
        }
    }

    return out;
}

uint64_t UncompressedBF::similarity(const BF* other, int type) const {
    assert(other->num_filter_bits() == num_filter_bits());

    const uint64_t* b1_data = this->bv->data();
    const UncompressedBF* o = dynamic_cast<const UncompressedBF*>(other);
    if (o == nullptr) {
        DIE("Can only compute similarity on same type of BF.");
    }

    const uint64_t* b2_data = o->bv->data();
    if (type == 1) {
        uint64_t xor_count = 0;
        uint64_t or_count = 0;
    
        // RSH: using num_filter_bits() instead of capacity() misses the bits in
        // .. the final chunk if size%64 is non-zero; caveat: using capacity()
        // .. risks counting extra bits in the final chunk
        sdsl::bit_vector::size_type len = this->bv->capacity()>>6;
        for (sdsl::bit_vector::size_type p = 0; p < len; ++p) {
            xor_count += __builtin_popcountl((*b1_data) ^ (*b2_data));
            or_count += __builtin_popcountl((*b1_data++) | (*b2_data++));
            //count += __builtin_popcountl((*b1_data++) ^ (*b2_data++));
        }
   
        return uint64_t(float(or_count - xor_count) / float(or_count) * num_filter_bits() );
    }
    else if (type == 0) {
        uint64_t count = 0;
        // RSH: using num_filter_bits() instead of capacity() misses the bits in
        // .. the final chunk if size%64 is non-zero; caveat: using capacity()
        // .. risks counting extra bits in the final chunk
        sdsl::bit_vector::size_type len = this->bv->capacity()>>6;
        for (sdsl::bit_vector::size_type p = 0; p < len; ++p) {
            count += __builtin_popcountl((*b1_data++) ^ (*b2_data++));
        }
        return num_filter_bits() - count;
    }

    DIE("ERROR: ONLY TWO TYPES IMPLEMENTED");
    return 0;
}



std::tuple<uint64_t, uint64_t> UncompressedBF::b_similarity(const BF* other) const {
    assert(other->num_filter_bits() == num_filter_bits());

    const uint64_t* b1_data = this->bv->data();
    const UncompressedBF* o = dynamic_cast<const UncompressedBF*>(other);
    if (o == nullptr) {
        DIE("Can only compute similarity on same type of BF.");
    }

    const uint64_t* b2_data = o->bv->data();

    uint64_t and_count = 0;
    uint64_t or_count = 0;
    // RSH: using num_filter_bits() instead of capacity() misses the bits in
    // .. the final chunk if size%64 is non-zero; caveat: using capacity()
    // .. risks counting extra bits in the final chunk
    sdsl::bit_vector::size_type len = this->bv->capacity()>>6;
    for (sdsl::bit_vector::size_type p = 0; p < len; ++p) {
        and_count += __builtin_popcountl((*b1_data) & (*b2_data)); // Rayan: there was a bug here
        or_count += __builtin_popcountl((*b1_data++) | (*b2_data++));
        //count += __builtin_popcountl((*b1_data++) ^ (*b2_data++));
    }
    return std::make_tuple(and_count, or_count);
//num_filter_bits() - count;
}


// 00 = 0
// 01 = 1
// 10 = 1
// 11 = 0

int popcount(uint64_t b) {
    int count = 0;
    for (int i = 0; i < 64; i++) {
        if (b % 2) count++;
        b >>= 1;
    }
    return count;
}

// return the # of 1s in the bitvector
uint64_t UncompressedBF::count_ones() const {
    uint64_t* data = bv->data();
    // RSH: using size() instead of capacity() misses the bits in the final
    // .. chunk if size%64 is non-zero; caveat: using capacity() risks counting
    // .. extra bits in the final chunk
    sdsl::bit_vector::size_type len = bv->capacity()>>6;
    uint64_t count = 0;
    for (sdsl::bit_vector::size_type p = 0; p < len; ++p) {
        int pc = popcount(*data);
        DIE_IF(pc != __builtin_popcountl(*data), "popcountl and us disagree about popcount");
        data++;
        count += pc;
    }

    sdsl::rank_support_v<> rbv(bv);
    DIE_IF(rbv.rank(bv->size()) != count, "SDSL and us disagree about number of 1s");
    return count;
}

void UncompressedBF::compress_rrr() {
    sdsl_rrr_vector rrr(*bv);
    //std::cerr << "Compressed RRR vector is " << sdsl::size_in_mega_bytes(rrr)
    //          << "MB,  storing it to " << filename+".rrr" << std::endl;
    sdsl::store_to_file(rrr,filename+".rrr");
}

void UncompressedBF::compress_wraprrr(const BF* f2) {
    const UncompressedBF* b2 = dynamic_cast<const UncompressedBF*>(f2);
    if (b2 == nullptr) {
        DIE("Can only use wraprrr with two uncompressed BF");
    }

    WrapRrrHeader header;
    header.padding1        = 0;
    header.padding2        = 0;
    header.padding3        = 0;
    header.num_bit_vectors = 2;
    header.num_bv1_bits    = bv->size();
    header.num_bv2_bits    = b2->bv->size();
    header.bv1_compressor  = 2;  // assume compressor is rrr unless we
    header.bv2_compressor  = 2;  // .. determine otherwise

    int64_t bv1_constant_val = bv_bit_constant(bv);
    int64_t bv2_constant_val = bv_bit_constant(b2->bv);

    if (bv1_constant_val >= 0) {
        //std::cerr << "  " << filename << " is constant " << bv1_constant_val << std::endl;
        header.bv1_compressor = bv1_constant_val;
    }

    if (bv2_constant_val >= 0) {
        //std::cerr << "  " << f2->filename << " is constant " << bv2_constant_val << std::endl;
        header.bv2_compressor = bv2_constant_val;
    }

    std::string wrapRrrFilename = filename;
    size_t suffix_ix = wrapRrrFilename.rfind(".bf-all.bv");
    if (suffix_ix != std::string::npos) {
         wrapRrrFilename = wrapRrrFilename.substr(0,suffix_ix) + ".bf-allsome.bv" + wrapRrrFilename.substr(suffix_ix+10);
    }
    wrapRrrFilename = wrapRrrFilename + ".wraprrr";

    std::cerr << "Creating " << wrapRrrFilename << std::endl;
    //std::cerr << "  bv->size() = " << bv->size() << std::endl;
    //std::cerr << "  b2->bv->size() = " << b2->bv->size() << std::endl;

    sdsl::osfstream out(wrapRrrFilename, std::ios::binary | std::ios::trunc | std::ios::out);
    out.write((char*)&header, sizeof(header));
    if (header.bv1_compressor == 2) {
        sdsl_rrr_vector rrr1(*bv);
        //std::cerr << "  adding " << sdsl::size_in_bytes(rrr1) << " byte component for " << filename << std::endl;
        rrr1.serialize(out, nullptr, "");
    }
    if (header.bv2_compressor == 2) {
        sdsl_rrr_vector rrr2(*b2->bv);
        //std::cerr << "  adding " << sdsl::size_in_bytes(rrr2) << " byte component for " << f2->filename << std::endl;
        rrr2.serialize(out, nullptr, "");
    }
    out.close();
}

void UncompressedBF::uncompress() {
    DIE("Cant uncompress an already uncompressed BF");
}

void UncompressedBF::compress_roar() {
    roaring_bitmap_t *roar_bits = roaring_bitmap_create();
    unsigned long nb_ones = 0;
    for (unsigned long i = 0; i < num_filter_bits(); i++) {
        if ((*bv)[i]) {
            roaring_bitmap_add(roar_bits, i);
            nb_ones ++;
        }
    }
    std::cerr << "Creating a roaring bitvector for file " << filename << std::endl;
    std::cerr << "Saturation: " << nb_ones << "/" << num_filter_bits()
              << " = " << (double)nb_ones/(double)num_filter_bits()
              << std::endl;
    roar_to_file(roar_bits, num_filter_bits(), filename + ".roar");
    roaring_bitmap_free(roar_bits);

    // union timing code
    /*  
    clock_t begin = clock();
    roaring_bitmap_t *r1_2_3 = roaring_bitmap_or(r1, r1);
          clock_t end = clock();
      double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
      std::cout << " time to union " << elapsed_secs << std::endl;

      begin=clock();
    roaring_bitmap_or_inplace(r1_2_3, r1);
    end=clock();
      elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
      std::cout << " time to inplace-union " << elapsed_secs << std::endl;

            begin=clock();
    this->union_into(this);
    end=clock();
      elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
      std::cout << " time to uncompressed union " << elapsed_secs << std::endl;
    */
}

void UncompressedBF::split(BF* child0_bf_all, const BF* child0_bf_some,
                           BF* child1_bf_all, const BF* child1_bf_some,
                           const std::string& f_all,  BF*& _bf_all,
                           const std::string& f_some, BF*& _bf_some) {
    // note that the child bf_all bit vectors may be modified as a side effect

    // assure either both children are NULL, or neither is
    assert ((child0_bf_all == nullptr) == (child0_bf_some == nullptr));
    assert ((child1_bf_all == nullptr) == (child1_bf_some == nullptr));
    assert ((child0_bf_all == nullptr) == (child1_bf_all  == nullptr));
    bool is_leaf = (child0_bf_all == nullptr);

    //DEBUG
    //if (is_leaf) {
    //    std::cerr << "UncompressedBF[" << filename << "]::split(null,null)" << std::endl;
    //} else {
    //    std::cerr << "UncompressedBF[" << filename << "]::split("
    //              << child0_bf_all->filename
    //              << "," << child1_bf_all->filename << ")" << std::endl;
    //}

    UncompressedBF* bf_all  = new UncompressedBF(f_all,  hashes, num_hash, num_filter_bits());
    UncompressedBF* bf_some = new UncompressedBF(f_some, hashes, num_hash, num_filter_bits());
    _bf_all  = bf_all;
    _bf_some = bf_some;

    if (is_leaf) {
        // at leaf L,
        //   BFall(L)  = BF(L)
        //   BFsome(L) = all zero (empty set)
        *(bf_all->bv)  =  *bv;
        *(bf_some->bv) =  *bv;
        *(bf_some->bv) ^= *bv;
    } else {
        // at node N with children C0 and C1,
        UncompressedBF*       uchild0_bf_all  = (UncompressedBF*)       child0_bf_all;
        UncompressedBF*       uchild1_bf_all  = (UncompressedBF*)       child1_bf_all;
        const UncompressedBF* uchild0_bf_some = (const UncompressedBF*) child0_bf_some;
        const UncompressedBF* uchild1_bf_some = (const UncompressedBF*) child1_bf_some;
        //   bits that are in both children's ALL are percolated up to the
        //   parent, and are removed from the children
        //     BFall(N)  = BFall(C0) and BFall(C1)
        //     BFall(C0) = BFall(C0) and not BFall(N)
        //     BFall(C1) = BFall(C1) and not BFall(N)
        *(bf_all->bv)         =  *(uchild0_bf_all->bv);
        *(bf_all->bv)         &= *(uchild1_bf_all->bv);
        *(uchild0_bf_all->bv) ^= *(bf_all->bv);
        *(uchild1_bf_all->bv) ^= *(bf_all->bv);
        //   bits that are in either children's SOME are percolated up to the
        //   parent, as well as any bits in a child's ALL that did not make it
        //   to the parent ALL
        //     BFsome(N) = BFsome(C0) or BFsome(C1)
        //                 or [BFall(C0) and not BFall(N)]
        //                 or [BFall(C1) and not BFall(N)]
        *(bf_some->bv) =  *(uchild0_bf_all->bv);
        *(bf_some->bv) |= *(uchild0_bf_some->bv);
        *(bf_some->bv) |= *(uchild1_bf_all->bv);
        *(bf_some->bv) |= *(uchild1_bf_some->bv);
    }
}

std::string UncompressedBF::compression_name() {
    return "none";
}

// union using 64bit integers
sdsl::bit_vector* union_bv_fast(const sdsl::bit_vector & b1, const sdsl::bit_vector& b2) {
    assert(b1.size() == b2.size());

    sdsl::bit_vector* out = new sdsl::bit_vector(b1.size(), 0);
    uint64_t* out_data = out->data();

    const uint64_t* b1_data = b1.data();
    const uint64_t* b2_data = b2.data();
    // RSH: using size() instead of capacity() misses the bits in the final
    // .. chunk if size%64 is non-zero
    // .. replaced sdsl::bit_vector::size_type len = b1.size()>>6;
    sdsl::bit_vector::size_type len = b1.capacity()>>6;
    for (sdsl::bit_vector::size_type p = 0; p < len; ++p) {
        (*out_data++) = (*b1_data++) | (*b2_data++);
    }
    return out;
}

int64_t bv_bit_constant(const sdsl::bit_vector* bv) {
    // returns:
    //   0  if the segment contains all zeros
    //   1  if the segment contains all ones
    //   -1 if the segment contains a mix of zeros and ones
    // empty bit vectors are considered to be all zeros.

    sdsl::bit_vector::size_type bvLen = bv->size();
    if (bvLen == 0) return -1;
    sdsl::bit_vector::size_type whole_chunks = bvLen >> 6;
    sdsl::bit_vector::size_type extra_bits   = bvLen - (whole_chunks << 6);

    uint64_t all_ones = (uint64_t) -1;

    // get the constant from the first chunk, which may be incomplete

    uint64_t alien_bits = 0;
    if ((whole_chunks == 0) && (extra_bits != 0)) {
        // the vector exists entirely in a single chunk, and is missing
        // some upper bits in that chunk
        alien_bits = all_ones << extra_bits;  // set 64-N most-significant bits
        whole_chunks = 1;
        extra_bits   = 0;
    }

    const uint64_t* data = bv->data();
    uint64_t constant = *(data++);

    if      ((constant & ~alien_bits) == 0)        constant = 0;
    else if ((constant |  alien_bits) == all_ones) constant = all_ones;
    else return -1;

    //std::cerr << "    bv_bit_constant: initial constant = " << std::hex << std::uppercase << constant << std::endl;

    // check that all the remaining whole chunks match the constant

    for (sdsl::bit_vector::size_type p = 1; p < whole_chunks; ++p) {
        if (*(data++) != constant) {
            //std::cerr << "    bv_bit_constant: chunk " << p << " is " << std::hex << std::uppercase << data[-1] << std::endl;
            return -1;
        }
    }

    // check that the final partial chunk, if there is one, matches the constant

    if (extra_bits != 0) {
        // the vector extends into another chunk, and is missing some upper
        // bits in that chunk
        alien_bits = all_ones << extra_bits;  // set 64-N most-significant bits
        if (constant == 0) {
            if ((*data & ~alien_bits) != 0) return -1;
        } else /* if (constant == all_ones) */ {
            if ((*data |  alien_bits) != all_ones) return -1;
        }
    }

    // return the constant as 0 or 1

    //std::cerr << "    bv_bit_constant: final constant = " << std::hex << std::uppercase << constant << std::endl;

    return constant & 1;
}


// RSH: misleading name; bf is not actually loaded, just created and bound to a filename
BF* load_bf_from_file(const std::string & fn, HashPair hp, int nh) {
    if (ends_with(fn, ".rrr")) {
        return new BF(fn, hp, nh);
    } else if (ends_with(fn, ".roar")) {
        return new RoarBF(fn, hp, nh);
    }  else if (ends_with(fn, ".bv")) {
        return new UncompressedBF(fn, hp, nh);
    } else if (ends_with(fn, ".wraprrr")) {
        return new WrapRrrBF(fn, hp, nh);
    } else {
        DIE("unknown bloom filter filetype (make sure extension is .rrr, .roar, .bv or .wraprrr)");
        return nullptr;
    }
}


// added by Rayan
// create an in-memory uncompressed representation of a compressed bloom filter
// when calling similarity() on compressed bloom filters, for example
UncompressedBF* uncompress_bf(const BF *in) 
{
    UncompressedBF *out = new UncompressedBF("dummy_file", in->hashes, in->num_hash);
    // create an uncompressed version of this bf
    // it's slow.
    out->bv = new sdsl::bit_vector(in->bits->size(), 0);

    for (unsigned long i = 0; i < in->bits->size(); i++) 
        (*(out->bv))[i] = (*(in->bits))[i];
    
    return out;
}



/*============================================*/
// Roar BF

RoarBF::RoarBF(const std::string & f, HashPair hp, int nh, uint64_t num_filter_bits) :
    BF(f, hp, nh),
    roar(nullptr)
{
    //std::cerr << "Constructing roar BF " << f << std::endl;
    size_in_bits = num_filter_bits;
    //roar = roaring_bitmap_create();  this defeats the assert(roar == nullptr) in load()
}

RoarBF::~RoarBF() {
    // call to base destructor happens automatically
    if (roar != nullptr) {
        roaring_bitmap_free(roar);
    }
}

void RoarBF::load(uint64_t partitions) {
    // read the actual bits
    assert(roar == nullptr);
    size_in_bits = 0;

    // std::cerr << "Loading roar BF from " << filename << std::endl;
    auto start_time = get_wall_time();
    roar = roar_from_file(filename, size_in_bits);
    auto elapsed_time = elapsed_wall_time(start_time);
    uint64_t fsize= filesize(filename);
    std::cerr << std::setiosflags(std::ios::fixed) << std::setprecision(4)  << "[BF load] " << elapsed_time << " secs " << filename << " " << fsize/1024/1024 << " MB, " << fsize / elapsed_time /1024/1024 << " MB/sec" << std::endl;

    set_filter_partitions(partitions);
    update_hash_function_range();

    //fprintf (stderr, "loaded %s size=%d allocation_size=%d\n",
    //                 filename.c_str(),
    //                 roar->high_low_container->size,
    //                 roar->high_low_container->allocation_size);
    //roaring_array_t* ra = roar->high_low_container;
    //for (int i=0; i<ra->size; ++i) {
    //    fprintf (stderr, "  container [%d] key=%d type=%s card=%d\n",
    //             i,
    //             (int) ra->keys[i],
    //             get_full_container_name(ra->containers[i], ra->typecodes[i]),
    //             container_get_cardinality(ra->containers[i], ra->typecodes[i]));
    //    if (ra->typecodes[i] == BITSET_CONTAINER_TYPE_CODE) {
    //        bitset_container_s* bs = (bitset_container_s*) ra->containers[i];
    //        fprintf (stderr, "  bitset card=%d\n", bs->cardinality);
    //        fprintf (stderr, "  bitset array=");
    //        for (int j=0; j<16; ++j)
    //            fprintf (stderr, " %02X", ((unsigned char*) bs->array)[j]);
    //        fprintf (stderr, " ..\n");
    //    }
    //}
}

void RoarBF::save() {
    //std::cerr << "Saving Roar BF to " << filename << std::endl;
    assert(roar != nullptr);
    auto start_time = get_wall_time();
    roar_to_file(roar, num_filter_bits(), filename);
    auto elapsed_time = elapsed_wall_time(start_time);
    std::cerr << "[BF save] " << elapsed_time << " secs " << filename << std::endl;
}

uint64_t RoarBF::num_filter_bits() const {
    return size_in_bits;
}

uint64_t RoarBF::num_compressed_bytes() const {
    assert(roar != nullptr);
    return roaring_bitmap_portable_size_in_bytes(roar);
}

int RoarBF::operator[](uint64_t pos) const {
    //std::cout << "RoarBF::[" << pos << "]" << std::endl;
    assert(roar != nullptr);
    return roaring_bitmap_contains(roar, pos);
}

void RoarBF::set_bit(uint64_t p) {
    assert(roar != nullptr);
    roaring_bitmap_add(roar, p);
}

BF* RoarBF::union_with(const std::string & new_name, const BF* f2) const {
    assert(roar != nullptr);
    assert(num_filter_bits() == f2->num_filter_bits());
    const RoarBF* o = dynamic_cast<const RoarBF*>(f2);
    if (o == nullptr) {
        DIE("Can only Roar-union two Roar BF");
    }

    //std::cerr << "Performing RoarBF union_with of " << this->filename << " with " << f2->filename << std::endl;
    RoarBF* out = new RoarBF(new_name, hashes, num_hash, num_filter_bits());
    out->roar = roaring_bitmap_or(this->roar, o->roar);

    // let's compute some stats also!
    //std::cerr << "Union" << std::endl;
    //std::cerr << out->count_ones() << std::endl;

    //RoarBF* tmp = new RoarBF(new_name+".dummytmp", hashes, num_hash, num_filter_bits());
    //tmp->roar = roaring_bitmap_and(this->roar, o->roar);
    //std::cerr << "Intersection" << std::endl;
    //std::cerr << tmp->count_ones() << std::endl;
    //delete tmp;

    return out;
}

BF* RoarBF::intersection_with(const std::string & new_name, const BF* f2) const {
    assert(roar != nullptr);
    assert(num_filter_bits() == f2->num_filter_bits());
    const RoarBF* o = dynamic_cast<const RoarBF*>(f2);
    if (o == nullptr) {
        std::cerr << "Problem with " << f2->filename << std::endl;
        DIE("Can only Roar-intersect two Roar BF");
    }

    //std::cerr << "Performing RoarBF intersection_with..." << std::endl;
    RoarBF* out = new RoarBF(new_name, hashes, num_hash, num_filter_bits());
    out->roar = roaring_bitmap_and(this->roar, o->roar);

    //std::cerr << "Intersection of " << this->filename << " with " << f2->filename << std::endl;
    //std::cerr << out->count_ones() << std::endl;
    
    return out;
}

void RoarBF::union_into(const BF* f2) {
    assert(roar != nullptr);
    assert(num_filter_bits() == f2->num_filter_bits());
    const RoarBF* o = dynamic_cast<const RoarBF*>(f2);
    if (o == nullptr) {
        DIE("Can only Roar-union two Roar BF");
    }

    roaring_bitmap_or_inplace(this->roar, o->roar);
}

BF* RoarBF::append(const std::string & new_name, const BF* f2) const {
    DIE("appending to ROAR not yet implemented");
    return nullptr;
}

uint64_t RoarBF::similarity(const BF* other, int type) const {
    assert(roar != nullptr);
    assert(other->num_filter_bits() == num_filter_bits());

    const RoarBF* o = dynamic_cast<const RoarBF*>(other);
    if (o == nullptr) {
        DIE("Can only compute similarity on same type of BF.");
    }

    if (type == 1) {
        uint64_t xor_count = 0;
        uint64_t or_count = 0;

        roaring_bitmap_t *roar_or= roaring_bitmap_or(this->roar, o->roar);
        or_count = roaring_bitmap_get_cardinality(roar_or);
        roaring_bitmap_free(roar_or);

        roaring_bitmap_t *roar_xor= roaring_bitmap_xor(this->roar, o->roar);
        xor_count = roaring_bitmap_get_cardinality(roar_xor);
        roaring_bitmap_free(roar_xor);

        return uint64_t(float(or_count - xor_count) / float(or_count) * num_filter_bits() );
    }
    else if (type == 0) {
        uint64_t count = 0;
        roaring_bitmap_t *roar_xor= roaring_bitmap_xor(this->roar, o->roar);
        count = roaring_bitmap_get_cardinality(roar_xor);
        roaring_bitmap_free(roar_xor);
        return num_filter_bits() - count;
    }

    DIE("ERROR: ONLY TWO TYPES IMPLEMENTED");
    return 0;
}



std::tuple<uint64_t, uint64_t> RoarBF::b_similarity(const BF* other) const {
    assert(roar != nullptr);
    assert(other->num_filter_bits() == num_filter_bits());

    const RoarBF* o = dynamic_cast<const RoarBF*>(other);
    if (o == nullptr) {
        DIE("Can only compute similarity on same type of BF.");
    }

    uint64_t and_count = 0;
    uint64_t or_count = 0;

    roaring_bitmap_t *roar_or= roaring_bitmap_or(this->roar, o->roar);
    or_count = roaring_bitmap_get_cardinality(roar_or);
    roaring_bitmap_free(roar_or);

    roaring_bitmap_t *roar_and= roaring_bitmap_and(this->roar, o->roar);
    and_count = roaring_bitmap_get_cardinality(roar_and);
    roaring_bitmap_free(roar_and);

    return std::make_tuple(and_count, or_count);
}

// return the # of 1s in the bitvector
uint64_t RoarBF::count_ones() const {
    assert(roar != nullptr);
    unsigned long count = roaring_bitmap_get_cardinality(this->roar);
    return count;
}

void RoarBF::compress_rrr() {
    // decompress roar to a bit vector
    auto t = sdsl::bit_vector(num_filter_bits(), 0);
    std::cerr << "Converting ROAR bitvector " << filename
              << " to RRR (" << num_filter_bits() << " bits)" << std::endl;

    for (unsigned long pos = 0; pos < num_filter_bits(); pos++) {
        if (roaring_bitmap_contains(roar,pos)) {
            t[pos] = 1;
        }
    }

    // compress bit vector to rrr
    sdsl_rrr_vector rrr(t);

    // replace the .roar extension with .rrr
    size_t lastindex = filename.find_last_of("."); 
    std::string rawname = filename.substr(0, lastindex) + ".rrr"; 
    //std::cerr << "Compressed RRR vector is " << sdsl::size_in_mega_bytes(rrr)
    //          << "MB,  storing it to " << rawname << std::endl;
    sdsl::store_to_file(rrr,rawname);
}

void RoarBF::compress_wraprrr(const BF* f2) {
    DIE("Cant compress roar w/ wraprrr with existing code base");
}

void RoarBF::compress_roar() {
    DIE("Cant compress an already compressed RoarBF");
}

void RoarBF::uncompress() {
    auto t = sdsl::bit_vector(num_filter_bits(), 0);
    std::cerr << "Uncompressing ROAR bitvector " << filename
              << " (" << num_filter_bits() << " bits)" << std::endl;

    for (unsigned long pos = 0; pos < num_filter_bits(); pos++) {
        if (roaring_bitmap_contains(roar,pos)) {
            t[pos] = 1;
        }
    }

    // strip the .roar extension
    size_t lastindex = filename.find_last_of("."); 
    std::string rawname = filename.substr(0, lastindex); 
    double megaBytes = roaring_bitmap_portable_size_in_bytes(roar) / 1000000.0;
    std::cerr << "Compressed Roar vector would be " << megaBytes << "MB, storing uncompressed to " << rawname << std::endl;

    sdsl::store_to_file(t,rawname);
}

void RoarBF::split(BF* child0_bf_all, const BF* child0_bf_some,
                   BF* child1_bf_all, const BF* child1_bf_some,
                   const std::string& f_all,  BF*& _bf_all,
                   const std::string& f_some, BF*& _bf_some) {
    // note that the child bf_all bit vectors may be modified as a side effect

    // assure either both children are NULL, or neither is
    assert ((child0_bf_all == nullptr) == (child0_bf_some == nullptr));
    assert ((child1_bf_all == nullptr) == (child1_bf_some == nullptr));
    assert ((child0_bf_all == nullptr) == (child1_bf_all  == nullptr));
    bool is_leaf = (child0_bf_all == nullptr);

    //DEBUG
    //if (is_leaf) {
    //    std::cerr << "RoarBF[" << filename << "]::split(null,null)" << std::endl;
    //} else {
    //    std::cerr << "RoarBF[" << filename << "]::split("
    //              << child0_bf_all->filename
    //              << "," << child1_bf_all->filename << ")" << std::endl;
    //}

    RoarBF* bf_all  = new RoarBF(f_all,  hashes, num_hash, num_filter_bits());
    RoarBF* bf_some = new RoarBF(f_some, hashes, num_hash, num_filter_bits());
    _bf_all  = bf_all;
    _bf_some = bf_some;

    if (is_leaf) {
        // at leaf L,
        //   BFall(L)  = BF(L)
        //   BFsome(L) = all zero (empty set)
        bf_all->roar  = roaring_bitmap_copy(roar);
        bf_some->roar = roaring_bitmap_create_with_capacity(num_filter_bits());
    } else {
        // at node N with children C0 and C1,
        RoarBF*       uchild0_bf_all  = (RoarBF*)       child0_bf_all;
        RoarBF*       uchild1_bf_all  = (RoarBF*)       child1_bf_all;
        const RoarBF* uchild0_bf_some = (const RoarBF*) child0_bf_some;
        const RoarBF* uchild1_bf_some = (const RoarBF*) child1_bf_some;
        //   bits that are in both children's ALL are percolated up to the
        //   parent, and are removed from the children
        //     BFall(N)  = BFall(C0) and BFall(C1)
        //     BFall(C0) = BFall(C0) and not BFall(N)
        //     BFall(C1) = BFall(C1) and not BFall(N)
        bf_all->roar = roaring_bitmap_and(uchild0_bf_all->roar, uchild1_bf_all->roar);
        roaring_bitmap_xor_inplace(uchild0_bf_all->roar, bf_all->roar);
        roaring_bitmap_xor_inplace(uchild1_bf_all->roar, bf_all->roar);
        //   bits that are in either children's SOME are percolated up to the
        //   parent, as well as any bits in a child's ALL that did not make it
        //   to the parent ALL
        //     BFsome(N) = BFsome(C0) or BFsome(C1)
        //                 or [BFall(C0) and not BFall(N)]
        //                 or [BFall(C1) and not BFall(N)]
        bf_some->roar = roaring_bitmap_or(uchild0_bf_all->roar, uchild0_bf_some->roar);
        roaring_bitmap_or_inplace(bf_some->roar, uchild1_bf_all->roar);
        roaring_bitmap_or_inplace(bf_some->roar, uchild1_bf_some->roar);
    }
}

std::string RoarBF::compression_name() {
    return "roar";
}


/*============================================*/

roaring_bitmap_t *roar_from_file(std::string filename, uint64_t &size_in_bits) {
    // RSH: note that the roar file contains a small header before roaring's serialized data
    //from https://gist.github.com/sanmarcos/991042
    int fd = open(filename.c_str(), O_RDONLY, (mode_t)0660);
    struct stat fileInfo = {0};
    if (fstat(fd, &fileInfo) == -1)
    {
        fprintf(stderr, "problem with roar file %s\n", filename.c_str());
        perror("Can't get file size, probably a missing file, assumed to be an empty BF. Error code (should be bad file descriptor):");
        //exit(EXIT_FAILURE);
        return roaring_bitmap_create();
        
    }
    std::size_t headerbytes = sizeof(uint64_t);
    if (fileInfo.st_size == 0)
    {
        fprintf(stderr, "File %s is empty, nothing to load\n", filename.c_str());
        return roaring_bitmap_create();
    }
    //printf("File size is %ji\n", (intmax_t)fileInfo.st_size);
    if (fileInfo.st_size < (int64_t) headerbytes)
    {
        fprintf(stderr, "File %s has an incomplete header, nothing to load\n", filename.c_str());
        return roaring_bitmap_create();
    }
    if (fileInfo.st_size == (int64_t) headerbytes)
    {
        fprintf(stderr, "File %s is empty (except for the header), nothing to load\n", filename.c_str());
        return roaring_bitmap_create();
    }
    char *map = (char*) mmap(0, fileInfo.st_size, PROT_READ, MAP_SHARED, fd, 0);
    if (map == MAP_FAILED)
    {
        close(fd);
        perror("Error mmapping the file");
        exit(EXIT_FAILURE);
    }
    size_in_bits = *((uint64_t*) map);
    roaring_bitmap_t *roar_bits = roaring_bitmap_portable_deserialize(map + headerbytes);
    if (munmap(map, fileInfo.st_size) == -1)
    {
        close(fd);
        perror("Error un-mmapping the file");
        exit(EXIT_FAILURE);
    }
    close(fd);

    return roar_bits;
}


void roar_to_file(roaring_bitmap_t *roar_bits, uint64_t num_filter_bits, std::string filename) {
    uint32_t expected_size_basic = roaring_bitmap_portable_size_in_bytes(roar_bits);
    /*roaring_bitmap_run_optimize(roar_bits);
    uint32_t expected_size_run = roaring_bitmap_portable_size_in_bytes(roar_bits);
    printf("size before run optimize %d bytes, and after %d bytes\n",
                   expected_size_basic, expected_size_run);
    */ // optimization doesn't seem to change the bitmap size?!

    // RSH: we write a small header followed by roaring's serialized data; the
    //      header is simply the number of bits in the bloom filter, that info
    //      is not otherwise recoverable when the file is read back in

    // from https://gist.github.com/sanmarcos/991042
    int fd = open(filename.c_str(), O_RDWR | O_CREAT | O_TRUNC, (mode_t)0664);
    if (fd == -1) {
        perror("Error opening file for writing");
        exit(EXIT_FAILURE);
    }
    std::size_t headerbytes = sizeof(uint64_t);            // uint64 is overkill, but that's the type of num_filter_bits
    unsigned long sz = headerbytes + expected_size_basic;  // warning: possible overflow
    if (lseek(fd, sz-1, SEEK_SET) == -1) {
        close(fd);
        perror("Error calling lseek() to 'stretch' the file");
        exit(EXIT_FAILURE);
    }
    if (write(fd, "", 1) == -1) {
        close(fd);
        perror("Error writing last byte of the file");
        exit(EXIT_FAILURE);
    }
    char *map = (char*)mmap(0, sz, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0); 
    if (map == MAP_FAILED) {
        close(fd);
        perror("Error mmapping the file");
        exit(EXIT_FAILURE);
    }
    char *serializedbytes = (char *)malloc(sz);
    *((uint64_t*) serializedbytes) = num_filter_bits;
    roaring_bitmap_portable_serialize(roar_bits, serializedbytes + headerbytes);
    memcpy(map,serializedbytes,sz); // serialize failed on the map directly, so i'm trying this indirect route which should work
    if (msync(map, sz, MS_SYNC) == -1)
        perror("Could not sync the file to disk");
    if (munmap(map, sz) == -1) {
        close(fd);
        perror("Error un-mmapping the file");
        exit(EXIT_FAILURE);
    }
    close(fd);
   
    free(serializedbytes);
}



#include <iostream>
#include <string>
#include <jellyfish/mer_dna_bloom_counter.hpp>
#include <jellyfish/file_header.hpp>
#include <sdsl/bit_vectors.hpp>

#include <jellyfish/mer_dna.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/hash_counter.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/mer_overlap_sequence_parser.hpp>
#include <jellyfish/mer_iterator.hpp>

#include <vector>
#ifdef __cplusplus
extern "C" {
#endif

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <assert.h>
#include "fasta_parser.h"
#include "dictgen.h"
#include "dict_filt.h"
#include "util.h"
#include "vartype.h"
#include <stdint.h>

#ifdef __cplusplus
}
#endif
#include "BF.h"
#include "Count.h"
#include "Kmers.h"
#include "BloomTree.h"
#include "allsome_util.h"

using namespace std;

int main(int argc, char* argv[]){

    string hash_filename = argv[1];
    string bf_filename = argv[2];
	int number_of_hashes;
	jellyfish::mer_dna::k(32); // Set length of mers (k=32)
	HashPair* hash_pair = get_hash_function(hash_filename, number_of_hashes);
	/*read bloom filter file*/
	UncompressedBF* bloom_filter = new UncompressedBF(bf_filename, *hash_pair, number_of_hashes);
	/*please remember to delete hash_pair; delete bloom_filter;*/
	/*load bloom filter into memory*/
	bloom_filter->load();

    if(bloom_filter->contains("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")){
        cout << "yet" << endl;
    }else{
        cout << "no" << endl;
    }

    delete bloom_filter;
    delete hash_pair;

}

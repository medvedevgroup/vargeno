#ifndef KMERCOUNT_H
#define KMERCOUNT_H
#include <iostream>
#include <string>
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
#include "generate_bf.h"

using namespace std;

class KmerCounter {
public:
	KmerCounter();
	~KmerCounter();
	int readFasta(string & fasta_filename);
	int constructBfFromGenomeseq(string bf_filename, bool is_canonical);
	int constructBfFromVcf(const string & vcf_filename, string bf_filename, bool is_canonical);
	int constructBfFromUcsc(const string & text_filename, string bf_filename, bool is_canonical);
    int constructBfFromEncode(const string & encode_filename, string bf_filename, bool is_canonical);
	static const int REF_BF_RANGE = 32;
	static const int SNP_BF_RANGE = 40;

protected:
	vector<GenomeSeq> genome_vector;
	static const int KMER_SIZE = 32; // only right half kmers
	static const int HASH_NUM = 1;

	char rev(const char c);

};

#endif /* KMERCOUNT_H */

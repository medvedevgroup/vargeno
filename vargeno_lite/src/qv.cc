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
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <unordered_map>
#include <set>
#include <bitset>
#include "generate_bf.h"
using namespace std;

#if PCOMPACT
  #include "pileup.h"
#endif

/*global variable declearation*/
/*
uint64_t odd_mask_64 = 0x5555555555555555;
uint64_t even_mask_64 = 0xAAAAAAAAAAAAAAAA;
uint32_t odd_mask_32 = 0x55555555;
uint32_t even_mask_32 = 0xAAAAAAAA;
*/

// bloom filter for ref and snp kmer hi part
BloomFilter* ref_bf;
BloomFilter* snp_bf;

// bit mask declearation, initialization is in main function.
uint64_t odd_mask_64 = 0;
uint64_t even_mask_64 = 0;
uint32_t odd_mask_32 = 0;
uint32_t even_mask_32 = 0;

/* === calcualte diff base dictionary === */
// if two kmers have only one different base, this dictionary is to directly return the diff base position within the kmer
unordered_map<kmer_t, int> diff_base_dictionary;

/*
 * We use the following structures for quickly finding the index that
 * the most k-mers in a read (including all hamming neighbors) agree
 * on.
 */

#define INDEX_TABLE_SLOT_COUNT  1009
#define INDEX_TABLE_ENTRY_DEPTH  500  /* enough for 5 k-mers */

typedef struct {
	uint32_t index;
	uint8_t freq;
} IndexTableEntry;

typedef struct {
	unsigned short count;
	IndexTableEntry entries[INDEX_TABLE_ENTRY_DEPTH];
} IndexTableSlot;

typedef struct {
	IndexTableEntry *best;  // highest frequency
	bool ambiguous;         // `best` ambiguous?
	IndexTableSlot table[INDEX_TABLE_SLOT_COUNT];
} IndexTable;

void index_table_clear_index(IndexTable *index_table, uint32_t index);
void index_table_clear(IndexTable *index_table);
//void index_table_add(IndexTable *index_table, uint32_t index);

void index_table_clear_index(IndexTable *index_table, uint32_t index)
{
	index_table->table[index % INDEX_TABLE_SLOT_COUNT].count = 0;
}

void index_table_clear(IndexTable *index_table)
{
	index_table->best = NULL;
	index_table->ambiguous = false;

	for (size_t i = 0; i < INDEX_TABLE_SLOT_COUNT; i++) {
		index_table_clear_index(index_table, i);
	}
}

void index_table_add(IndexTable *index_table, uint32_t index, uint32_t kmer_pos, unordered_map<uint32_t, unordered_set<uint32_t>> & index_2_kmer_pos_set, bool is_neighbor=true)
{
    if(is_neighbor){
        if (index_2_kmer_pos_set.find(index) == index_2_kmer_pos_set.end()){
            // we need to make sure that at least one origin kmer support current position
            return;
        }
    }

	size_t slot_index = index % INDEX_TABLE_SLOT_COUNT;
	IndexTableSlot *slot = &index_table->table[slot_index];
	IndexTableEntry *target = NULL;

	for (unsigned short i = 0; i < slot->count; i++) {
		IndexTableEntry *e = &slot->entries[i];
		if (e->index == index) {
			++e->freq;
			target = e;
			break;
		}
	}

	if (target == NULL) {
#if DEBUG
		assert(slot->count < INDEX_TABLE_ENTRY_DEPTH);
#endif
		slot->entries[slot->count] = (IndexTableEntry){.index = index, .freq = 1};
		target = &slot->entries[slot->count];
		++slot->count;
	}

	index_2_kmer_pos_set[index].insert(kmer_pos);
	// if no two kmer pos, do not consider it as a match
	if (index_2_kmer_pos_set[index].size() <= 1) return;

	if (index_table->best == NULL) {
		index_table->best = target;
		index_table->ambiguous = false;
	} else if (target == index_table->best) {
		index_table->ambiguous = false;
	} else if (target->freq == index_table->best->freq) {
		index_table->ambiguous = true;
	} else if (target->freq > index_table->best->freq) {
		index_table->best = target;
		index_table->ambiguous = false;
	}
}

/* --- */

#define PILEUP_TABLE_INIT_SIZE (1 << 25)

static struct kmer_entry *query_ref_dict(kmer_t key,
                                         uint32_t *ref_jumpgate,
                                         struct kmer_entry *ref_dict,
                                         const size_t ref_dict_size);

static struct snp_kmer_entry *query_snp_dict(kmer_t key,
                                             uint32_t *snp_jumpgate,
                                             struct snp_kmer_entry *snp_dict,
                                             const size_t snp_dict_size);

static int ref_dict_entry_cmp(const void *key, const void *entry)
{
#if REF_LITE
	const uint64_t kmer_lo = *(uint64_t *)key;
	const uint64_t entry_lo = ((struct kmer_entry *)entry)->kmer_lo40;
#else
	const uint32_t kmer_lo = *(uint32_t *)key;
	const uint32_t entry_lo = ((struct kmer_entry *)entry)->kmer_lo;
#endif
	return (kmer_lo > entry_lo) - (kmer_lo < entry_lo);
}

static struct kmer_entry *query_ref_dict(kmer_t key,
                                         uint32_t *ref_jumpgate,
                                         struct kmer_entry *ref_dict,
                                         const size_t ref_dict_size)
{
#if REF_LITE
	const uint32_t kmer_hi = HI24(key);
	const uint64_t kmer_lo = LO40(key);
#else
	const uint32_t kmer_hi = HI(key);
	const uint32_t kmer_lo = LO(key);
#endif

	const uint32_t lo = ref_jumpgate[kmer_hi];

	if (lo == ref_dict_size) {
		return NULL;
	}

#if REF_LITE
	const uint32_t hi = (kmer_hi == 0xFFFFFF ? ref_dict_size : ref_jumpgate[kmer_hi + 1]);
#else
	const uint32_t hi = (kmer_hi == 0xFFFFFFFF ? ref_dict_size : ref_jumpgate[kmer_hi + 1]);
#endif

	if (lo == hi) {
		return NULL;
	}

#if DEBUG
	assert(hi > lo);
#endif
	struct kmer_entry *target = (struct kmer_entry *)bsearch(&kmer_lo, &ref_dict[lo], (hi - lo), sizeof(*ref_dict), ref_dict_entry_cmp);
	return target;
}

int check_block_size(kmer_t key, uint32_t *ref_jumpgate, const size_t ref_dict_size){

#if REF_LITE
		const uint32_t kmer_hi = HI24(key);
		const uint64_t kmer_lo = LO40(key);
#else
		const uint32_t kmer_hi = HI(key);
		const uint32_t kmer_lo = LO(key);
#endif

		const uint32_t lo = ref_jumpgate[kmer_hi];

		if (lo == ref_dict_size) {
			return 0;
		}

#if REF_LITE
		const uint32_t hi = (kmer_hi == 0xFFFFFF ? ref_dict_size : ref_jumpgate[kmer_hi + 1]);
#else
		const uint32_t hi = (kmer_hi == 0xFFFFFFFF ? ref_dict_size : ref_jumpgate[kmer_hi + 1]);
#endif
		return hi - lo;
}

// this might not be the best way, an easy way is to do binary search??
inline bool one_hamming_distance_32(uint32_t a, uint32_t b, int & diff_base_pos){
  uint32_t x = a ^ b;
  if(x == 0) return false; // k and l direct match

  if ((x & (x - 1)) == 0) {
	  diff_base_pos = diff_base_dictionary[x];
	  return true; // difference between k and l are only one bits
  }
  // next, we need to check difference between k and l are only two consecutive bits, higher bit should be in odd position, and lower bit should be in even position
  // do not need to check if y == 0 or z == 0, as the final step will check this.
  uint32_t y = x & odd_mask_32; // get odd position
  if((y & (y-1)) != 0) return false; // only one bit in odd position

  uint32_t z = x & even_mask_32; // get even position,
  if((z & (z-1)) != 0) return false; // only one bit in even position

  if (y == (z << 1)) {
	  diff_base_pos = diff_base_dictionary[x];
	  return true; // y and z should be consecutive
  }
  else return false;
}

// this might not be the best way, an easy way is to do binary search??
inline bool one_hamming_distance_64(uint64_t a, uint64_t b, int & diff_base_pos){
  uint64_t x = a ^ b;
  if(x == 0) return false; // k and l direct match

  if ((x & (x - 1)) == 0) {
	  diff_base_pos = diff_base_dictionary[x];
	  return true; // difference between k and l are only one bits
  }
  // next, we need to check difference between k and l are only two consecutive bits, higher bit should be in odd position, and lower bit should be in even position
  // do not need to check if y ==0 and z == 0, as the final step will check this.
  uint64_t y = x & odd_mask_64; // get odd position
  if((y & (y-1)) != 0) return false; // only one bit in odd position

  uint64_t z = x & even_mask_64; // get even position,
  if((z & (z-1)) != 0) return false; // only one bit in even position

  if (y == (z << 1)) {
	  diff_base_pos = diff_base_dictionary[x];
	  return true; // y and z should be consecutive
  }
  else return false;
}



void iterate_ref_dict(kmer_t key,
				uint32_t *ref_jumpgate,
				struct kmer_entry *ref_dict,
				const size_t ref_dict_size,
				struct kmer_entry* ref_hit_set[],
				int & ref_hit_set_size,
				uint64_t neighbors[],
				int ref_diff_base_pos_list[])
{
#if REF_LITE
	const uint32_t kmer_hi = HI24(key);
	const uint64_t kmer_lo = LO40(key);
#else
	const uint32_t kmer_hi = HI(key);
	const uint32_t kmer_lo = LO(key);
#endif

	// lo and hi are low and high position of the searching block
	ref_hit_set_size = 0;
	const uint32_t lo = ref_jumpgate[kmer_hi];

	if (lo == ref_dict_size) {
		return;
	}

#if REF_LITE
	const uint32_t hi = (kmer_hi == 0xFFFFFF ? ref_dict_size : ref_jumpgate[kmer_hi + 1]);
#else
	const uint32_t hi = (kmer_hi == 0xFFFFFFFF ? ref_dict_size : ref_jumpgate[kmer_hi + 1]);
#endif

	if (lo == hi) {
		return;
	}

#if DEBUG
	assert(hi > lo);
#endif
	
	const void * pointer_base = &ref_dict[lo];
	char* pointer_base_helper = (char*) pointer_base;
	size_t element_size = sizeof(*ref_dict);

#if REF_LITE
	for (uint32_t i = lo; i < hi; i++) {
		char* pointer_location_helper = pointer_base_helper + (i - lo) * element_size;
		//const void * pointer_location = pointer_base + (i - lo) * element_size;
		const void * pointer_location = (void *) pointer_location_helper;
		const uint64_t entry_lo = ((struct kmer_entry *)pointer_location)->kmer_lo40;
		//[TODO] here we need to make sure pointer_location->kmer_lo contains 40 bits	
		int diff_base_pos = -1;
		if (one_hamming_distance_64(kmer_lo, entry_lo, diff_base_pos)) {
			ref_hit_set[ref_hit_set_size] = &(ref_dict[i]);
			// construct the whole kmer
			neighbors[ref_hit_set_size] = kmer_hi;
			neighbors[ref_hit_set_size] <<= 40;
			neighbors[ref_hit_set_size] |= entry_lo;
			// for debug
			//assert(diff_base_pos >= 0);
			ref_diff_base_pos_list[ref_hit_set_size] = diff_base_pos;
			ref_hit_set_size++;
		}
	}
#else
	for (uint32_t i = lo; i < hi; i++) {
		char* pointer_location_helper = pointer_base_helper + (i - lo) * element_size;
		//const void * pointer_location = pointer_base + (i - lo) * element_size;
		const void * pointer_location = (void *) pointer_location_helper;
		const uint32_t entry_lo = ((struct kmer_entry *)pointer_location)->kmer_lo;		
		int diff_base_pos = -1;
		if (one_hamming_distance_32(kmer_lo, entry_lo, diff_base_pos)) {
			ref_hit_set[ref_hit_set_size] = &(ref_dict[i]);
			neighbors[ref_hit_set_size] = kmer_hi;
			neighbors[ref_hit_set_size] <<= 32;
			neighbors[ref_hit_set_size] |= entry_lo;
			// for debug
			//assert(diff_base_pos >= 0);
			ref_diff_base_pos_list[ref_hit_set_size] = diff_base_pos;
			ref_hit_set_size++;
		}
	}
#endif

	return;
}

static int snp_dict_entry_cmp(const void *key, const void *entry)
{
	const uint64_t kmer_lo = *(uint64_t *)key;
	const uint64_t entry_lo = ((struct snp_kmer_entry *)entry)->kmer_lo40;
	return (kmer_lo > entry_lo) - (kmer_lo < entry_lo);
}

static struct snp_kmer_entry *query_snp_dict(kmer_t key,
                                             uint32_t *snp_jumpgate,
                                             struct snp_kmer_entry *snp_dict,
                                             const size_t snp_dict_size)
{
  // so HI24 is specially designed for small dictionary
	const uint32_t kmer_hi = HI24(key);
	const uint64_t kmer_lo = LO40(key);

	const uint32_t lo = snp_jumpgate[kmer_hi];

	if (lo == snp_dict_size) {
		return NULL;
	}

	const uint32_t hi = (kmer_hi == 0xFFFFFF ? snp_dict_size : snp_jumpgate[kmer_hi + 1]);

	if (lo == hi) {
		return NULL;
	}

#if DEBUG
	assert(hi > lo);
#endif
	struct snp_kmer_entry *target = (struct snp_kmer_entry *)bsearch(&kmer_lo, &snp_dict[lo], (hi - lo), sizeof(*snp_dict), snp_dict_entry_cmp);
	return target;
}

void iterate_snp_dict(kmer_t key,
                      uint32_t *snp_jumpgate,
                      struct snp_kmer_entry *snp_dict,
                      const size_t snp_dict_size,
                      struct snp_kmer_entry* snp_hit_set[],
                      int & snp_hit_set_size,
                      uint64_t neighbors[],
					  int snp_diff_base_pos_list[])
{
  // so HI24 is specially designed for small dictionary
	const uint32_t kmer_hi = HI24(key);
	const uint64_t kmer_lo = LO40(key);

	snp_hit_set_size = 0;

	const uint32_t lo = snp_jumpgate[kmer_hi];

	if (lo == snp_dict_size) {
		return;
	}

	const uint32_t hi = (kmer_hi == 0xFFFFFF ? snp_dict_size : snp_jumpgate[kmer_hi + 1]);

	if (lo == hi) {
		return;
	}

#if DEBUG
	assert(hi > lo);
#endif

	const void * pointer_base = &snp_dict[lo];
	char* pointer_base_helper = (char*) pointer_base;
	size_t element_size = sizeof(*snp_dict);

  for (int i = lo; i < hi; i++) {
  	char* pointer_location_helper = pointer_base_helper + (i - lo) * element_size;
	//const void * pointer_location = pointer_base + (i - lo) * element_size;
	const void * pointer_location = (void *) pointer_location_helper;
	const uint64_t entry_lo = ((struct snp_kmer_entry *)pointer_location)->kmer_lo40;
	int diff_base_pos = -1;
    if (one_hamming_distance_64(kmer_lo, entry_lo, diff_base_pos)) {
      snp_hit_set[snp_hit_set_size] = &(snp_dict[i]);
      neighbors[snp_hit_set_size] = kmer_hi;
      neighbors[snp_hit_set_size] <<= 40;
      neighbors[snp_hit_set_size] |= entry_lo;
#if DEBUG
	  assert(diff_base_pos >= 0);
#endif
	  snp_diff_base_pos_list[snp_hit_set_size] = diff_base_pos;
      snp_hit_set_size++;
    }
  }
  return;
}

/* --- */

struct call { int genotype; double confidence; };
#define CALL(g, c) ((struct call){.genotype = (g), .confidence = (c)})
static inline struct call choose_best_genotype(const int ref_cnt,
                                               const int alt_cnt,
                                               const uint8_t ref_freq_enc,
                                               const uint8_t alt_freq_enc);

static void genotype(FILE *refdict_file, FILE *snpdict_file, FILE *fastq_file, FILE *chrlens_file, FILE *out)
{
	clock_t begin, end;
	double time_spent;
	begin = clock();

	/* Load chrlens file */
	struct { char name[32]; size_t len; } chrlens[128];
	size_t num_chrs = 0;

	char chrlen_buf[256];
	while (fgets(chrlen_buf, sizeof(chrlen_buf), chrlens_file)) {
		size_t i = 0;
		while (!isspace(chrlen_buf[i]) && i < sizeof(chrlens[0].name)) {
			chrlens[num_chrs].name[i] = chrlen_buf[i];
			++i;
		}
		chrlens[num_chrs].name[i] = '\0';

		while (isspace(chrlen_buf[i])) ++i;

		chrlens[num_chrs].len = atol(&chrlen_buf[i]);

		++num_chrs;
	}

	uint32_t *ref_jumpgate;
	struct kmer_entry *ref_dict;
	struct aux_table *ref_aux_table;
	uint32_t *snp_jumpgate;
	struct snp_kmer_entry *snp_dict;
	struct snp_aux_table *snp_aux_table;

#if PCOMPACT
	PileupTable ptable;
#else
	struct packed_pileup_entry *pileup_table;
#endif

	uint32_t last_hi;
	uint32_t max_pos = 0;

	fprintf(stderr, "Initializing...\n");

	/* === Reference Dictionary Construction === */
	const size_t ref_dict_size = read_uint64(refdict_file);
	const size_t ref_aux_table_size = read_uint64(refdict_file);

	if (ref_dict_size > POW_2_32) {
		fprintf(stderr, "Reference dictionary is too large (limit: %lu 32-mers)\n", POW_2_32);
		exit(EXIT_FAILURE);
	}

#if REF_LITE
	ref_jumpgate = (uint32_t*)malloc(POW_2_24 * sizeof(*ref_jumpgate));
#else
	ref_jumpgate = (uint32_t*)malloc(POW_2_32 * sizeof(*ref_jumpgate));
#endif
	assert(ref_jumpgate);
	ref_dict = (struct kmer_entry*)malloc(ref_dict_size * sizeof(*ref_dict));
	assert(ref_dict);
	ref_aux_table = (struct aux_table*)malloc(ref_aux_table_size * sizeof(*ref_aux_table));
	assert(ref_aux_table);

	ref_jumpgate[0] = 0;
	last_hi = 0;
	for (size_t i = 0; i < ref_dict_size; i++) {
		const kmer_t kmer = read_uint64(refdict_file);
		const uint32_t pos = read_uint32(refdict_file);
		const uint8_t ambig_flag = read_uint8(refdict_file);

#if REF_LITE
		ref_dict[i].kmer_lo40 = LO40(kmer);
#else
		ref_dict[i].kmer_lo = LO(kmer);
#endif
		ref_dict[i].pos = pos;
		ref_dict[i].ambig_flag = ambig_flag;

		if (pos > max_pos)
			max_pos = pos;

#if REF_LITE
		const uint32_t hi = HI24(kmer);
#else
		const uint32_t hi = HI(kmer);
#endif

		if (hi != last_hi) {
#if DEBUG
			assert(hi > last_hi);
#endif
			for (size_t j = (last_hi + 1); j <= hi; j++)
				ref_jumpgate[j] = i;

			last_hi = hi;
		}
	}

	// test printing all kmers_t and corresponding hi and lo
/*	ofstream tempfile;
	tempfile.open ("kmer_info.txt");
	for (size_t i = 0; i < ref_dict_size; i++) {
		
		const kmer_t kmer = read_uint64(refdict_file);
		const uint32_t pos = read_uint32(refdict_file);
		const uint8_t ambig_flag = read_uint8(refdict_file);

#if REF_LITE
		ref_dict[i].kmer_lo40 = LO40(kmer);
#else
		ref_dict[i].kmer_lo = LO(kmer);
#endif
		ref_dict[i].pos = pos;
		ref_dict[i].ambig_flag = ambig_flag;

		if (pos > max_pos)
			max_pos = pos;

#if REF_LITE
		const uint32_t hi = HI24(kmer);
#else
		const uint32_t hi = HI(kmer);
#endif
		// print here
		tempfile << kmer << "," << hi << "," << ref_dict[i].kmer_lo << endl;
	}
	tempfile.close();
*/
	// end test printing

#if REF_LITE
	if (last_hi != 0xFFFFFF) {
		for (size_t j = (last_hi + 1); j < POW_2_24; j++)
			ref_jumpgate[j] = ref_dict_size;
	}
#else
	if (last_hi != 0xFFFFFFFF) {
		for (size_t j = (last_hi + 1); j < POW_2_32; j++)
			ref_jumpgate[j] = ref_dict_size;
	}
#endif

	for (size_t i = 0; i < ref_aux_table_size; i++) {
		for (size_t j = 0; j < AUX_TABLE_COLS; j++) {
			ref_aux_table[i].pos_list[j] = read_uint32(refdict_file);
		}
	}

	/* === Pileup Table Initialization === */
#if PCOMPACT
	ptable_init(&ptable, PILEUP_TABLE_INIT_SIZE);
#else
	/*
	 * We assume that the maximum position encountered in the
	 * reference dictionary will not be smaller than the
	 * maximum position that will be encountered in the SNP
	 * dictionary. If not, we will reallocate our pileup table.
	 */
	size_t pileup_size = (size_t)max_pos + 32 + 1;
	pileup_table = (struct packed_pileup_entry*)calloc(pileup_size, sizeof(*pileup_table));
#endif

	/* === SNP Dictionary Construction === */
	const size_t snp_dict_size = read_uint64(snpdict_file);
	const size_t snp_aux_table_size = read_uint64(snpdict_file);

	if (snp_dict_size > POW_2_32) {
		fprintf(stderr, "SNP dictionary is too large (limit: %lu 32-mers)\n", POW_2_32);
		exit(EXIT_FAILURE);
	}

	snp_jumpgate = (uint32_t*)malloc(POW_2_24 * sizeof(*snp_jumpgate));
	assert(snp_jumpgate);
	snp_dict = (struct snp_kmer_entry*)malloc(snp_dict_size * sizeof(*snp_dict));
	assert(snp_dict);
	snp_aux_table = (struct snp_aux_table*)malloc(snp_aux_table_size * sizeof(*snp_aux_table));
	assert(snp_aux_table);

	snp_jumpgate[0] = 0;
	last_hi = 0;
	for (size_t i = 0; i < snp_dict_size; i++) {
		const kmer_t kmer = read_uint64(snpdict_file);
		const uint32_t pos = read_uint32(snpdict_file);
		const snp_info snp = read_uint8(snpdict_file);
		const uint8_t ambig_flag = read_uint8(snpdict_file);
		const uint8_t ref_freq = read_uint8(snpdict_file);
		const uint8_t alt_freq = read_uint8(snpdict_file);

		snp_dict[i].kmer_lo40 = LO40(kmer);
		snp_dict[i].pos = pos;
		snp_dict[i].snp = snp;
		snp_dict[i].ambig_flag = ambig_flag;

		const unsigned snp_info_ref = SNP_INFO_REF(snp);

		if ((snp_info_ref & BASE_N) == 0 &&  // i.e. if the reference base is A, C, G or T
		     pos != POS_AMBIGUOUS &&
		     ambig_flag == FLAG_UNAMBIGUOUS) {

			const unsigned snp_info_pos = SNP_INFO_POS(snp);  // relative to k-mer
			const uint32_t snp_pos = pos + snp_info_pos;      // relative to reference

#if PCOMPACT
			ptable_add(&ptable, snp_pos, snp_info_ref, kmer_get_base(kmer, snp_info_pos), ref_freq, alt_freq);
#else
			if (snp_pos >= pileup_size) {
				pileup_size = (snp_pos + 1) * sizeof(*pileup_table);
				printf("Re-allocing pileup table to %lu entries...\n", pileup_size);
				pileup_table = (struct packed_pileup_entry*)realloc(pileup_table, pileup_size);
				assert(pileup_table);
			}
			pileup_table[snp_pos].ref = snp_info_ref;
			pileup_table[snp_pos].alt = kmer_get_base(kmer, snp_info_pos);
			pileup_table[snp_pos].ref_freq = ref_freq;
			pileup_table[snp_pos].alt_freq = alt_freq;
#endif
		}

		const uint32_t hi = HI24(kmer);

		if (hi != last_hi) {
#if DEBUG
			assert(hi > last_hi);
#endif
			for (size_t j = (last_hi + 1); j <= hi; j++)
				snp_jumpgate[j] = i;

			last_hi = hi;
		}
	}

	if (last_hi != 0xFFFFFF) {
		for (size_t j = (last_hi + 1); j < POW_2_24; j++)
			snp_jumpgate[j] = snp_dict_size;
	}

	for (size_t i = 0; i < snp_aux_table_size; i++) {
		const kmer_t kmer = read_uint64(snpdict_file);
		UNUSED(kmer);

		for (size_t j = 0; j < AUX_TABLE_COLS; j++) {
			const uint32_t pos = read_uint32(snpdict_file);
			const snp_info snp = read_uint8(snpdict_file);
			const uint8_t ref_freq = read_uint8(snpdict_file);
			const uint8_t alt_freq = read_uint8(snpdict_file);
			UNUSED(ref_freq);
			UNUSED(alt_freq);

			snp_aux_table[i].pos_list[j] = pos;
			snp_aux_table[i].snp_list[j] = snp;

			/*
			if (pos != 0) {
				const unsigned snp_info_ref = SNP_INFO_REF(snp);
				const unsigned snp_info_pos = SNP_INFO_POS(snp);  // relative to k-mer
				const size_t snp_pos = pos + snp_info_pos;        // relative to reference

				if (snp_pos >= pileup_size) {
					pileup_size = (snp_pos + 1) * sizeof(*pileup_table);
					printf("Re-allocing pileup table to %lu entries...\n", pileup_size);
					pileup_table = realloc(pileup_table, pileup_size);
					assert(pileup_table);
				}

				pileup_table[snp_pos].ref = snp_info_ref;
				pileup_table[snp_pos].alt = kmer_get_base(kmer, snp_info_pos);
				pileup_table[snp_pos].ref_freq = ref_freq;
				pileup_table[snp_pos].alt_freq = alt_freq;
			}
			*/
		}
	}

bool is_ref_lite = false;
#if REF_LITE
	is_ref_lite = true;
#endif

	/* === Walk FASTQ File === */
#define BUF_SIZE 1024
	char id[BUF_SIZE];
	char read[BUF_SIZE];
	char read_revcompl[BUF_SIZE];
	char sep[BUF_SIZE];
	char qual[BUF_SIZE];

	kmer_t kmers[BUF_SIZE];

#define MAX_HITS 2000
#define NO_MODIFICATION 10086
	//struct kmer_entry *ref_hits[MAX_HITS];
	//struct snp_kmer_entry *snp_hits[MAX_HITS];

	size_t n_ref_hits;
	size_t n_snp_hits;

	/* convenient way to store k-mer information */
	typedef struct {
		kmer_t kmer;
		uint32_t position;  // 1-based position of read based on kmer hit
		uint32_t kmer_pos;  // 1-based position of k-mer
		uint32_t modified_pos;
#if DEBUG
		bool is_neighbor;
#endif
	} kmer_context;

	kmer_context ref_hit_contexts[MAX_HITS];  // same count as `ref_hits` (i.e. `n_ref_hits`)
	kmer_context snp_hit_contexts[MAX_HITS];  // same count as `snp_hits` (i.e. `n_snp_hits`)
#undef MAX_HITS
#undef BUF_SIZE

	IndexTable index_table;
	index_table_clear(&index_table);

#if DEBUG
	size_t total_count = 0;
	size_t match_count = 0;
	size_t multi_count = 0;
	size_t nohit_count = 0;

	size_t good_reads = 0;  // give us SNP information
	size_t bad_reads = 0;   // don't give us anything

	size_t ambig_hits = 0;
	size_t unambig_hits = 0;

	size_t ref_covs = 0;
	size_t alt_covs = 0;
	size_t non_ref_or_alt_covs = 0;
#endif

	fprintf(stderr, "Processing...\n");

#if DEBUG
	FILE *read_data = fopen("read_data.txt", "w");
	assert(read_data);
#endif

	while (fgets(id, sizeof(id), fastq_file)) {
		char *unused = fgets(read, sizeof(read), fastq_file);
		unused =       fgets(sep,  sizeof(sep),  fastq_file);
		unused =       fgets(qual, sizeof(qual), fastq_file);
		UNUSED(unused);
#if DEBUG
		assert(!ferror(fastq_file));
#endif

		bool revcompl = false;
		unordered_map<uint32_t, unordered_set<uint32_t>> index_2_kmer_pos_set;

		/*
		 * We process reads in 32-base chunks, so we trim off
		 * any remainder if the read length is not a multiple
		 * of 32. Also, we assume this is a valid FASTQ file,
		 * so the read length is equal to the quality string
		 * length.
		 */
		const size_t read_len_true = strlen(read) - 1;  // -1 because of newline char
		const size_t len = (read_len_true/32)*32;
		//const bool need_terminal_kmer = (read_len_true != len);



		bool hit_flag = true;

		head:
		if (revcompl) {
			for (size_t i = 0; i < len /*read_len_true*/; i++) {
				char rev = '\0';
				switch (read[i]) {
				case 'a': case 'A': rev = 'T'; break;
				case 'c': case 'C': rev = 'G'; break;
				case 'g': case 'G': rev = 'C'; break;
				case 't': case 'T': rev = 'A'; break;
				default: hit_flag = false; break;
				}
				if (!hit_flag) break;
				read_revcompl[len /*read_len_true*/ - i - 1] = rev;
			}
			if (!hit_flag) {
				index_table.best = NULL;
				index_table.ambiguous = false;
				continue;
			}
			memcpy(read, read_revcompl, len /*read_len_true*/);  // newline and '\0' already in `read`
		}

		hit_flag = true;

		size_t kmer_count = 0;
		for (size_t i = 0; i < len; i += 32) {
			bool kmer_had_n;
			kmer_t kmer = encode_kmer(&read[i], &kmer_had_n);

			if (kmer_had_n)
			{
				hit_flag = false;
				break;
			}

			kmers[kmer_count++] = kmer;
		}

		if (!hit_flag) {
			index_table.best = NULL;
			index_table.ambiguous = false;
			continue;
		}

		// (possibly) one last k-mer to cover entire read
		/*
		if (need_terminal_kmer) {
			bool kmer_had_n;
			kmer_t kmer = encode_kmer(&read[read_len_true - 32], &kmer_had_n);

			if (kmer_had_n)
				goto nohit;

			kmers[kmer_count++] = kmer;
		}
		*/

		n_ref_hits = 0;
		n_snp_hits = 0;

		/* loop over k-mers, perform ref/SNP dict queries */
		for (size_t i = 0; i < kmer_count; i++) {
			const kmer_t kmer = kmers[i];
			const char qual_char = qual[i];
			//const uint32_t offset = (need_terminal_kmer && i == (kmer_count - 1)) ? (read_len_true - 32) : 32*i;
			const uint32_t offset = 32*i;

			struct kmer_entry *ref_hit = query_ref_dict(kmer, ref_jumpgate, ref_dict, ref_dict_size);
			struct snp_kmer_entry *snp_hit = query_snp_dict(kmer, snp_jumpgate, snp_dict, snp_dict_size);

			int block_size = check_block_size(kmer, ref_jumpgate, ref_dict_size);
			// [TODO] here we also need to check the block size of snp_jumpgate. 
			// But we can have an assumption that snp block is always smaller than ref block.

			const bool orig_ref_hit_not_null = (ref_hit != NULL);
			const bool orig_snp_hit_not_null = (snp_hit != NULL);

			if (orig_ref_hit_not_null && ref_hit->pos != POS_AMBIGUOUS) {
				if (ref_hit->ambig_flag == FLAG_UNAMBIGUOUS) {
					const uint32_t read_pos = ref_hit->pos - offset;
					ref_hit_contexts[n_ref_hits++] = (kmer_context){.kmer = kmer,
					                                                .position = read_pos,
					                                                .kmer_pos = ref_hit->pos,
																	.modified_pos = NO_MODIFICATION,
#if DEBUG
					                                                .is_neighbor = false
#endif
					                                            };
					index_table_add(&index_table, read_pos, ref_hit->pos, index_2_kmer_pos_set, false);
#if DEBUG
					++unambig_hits;
#endif
				} else if (ref_hit->ambig_flag == FLAG_AMBIGUOUS) {
					const struct aux_table *p = &ref_aux_table[ref_hit->pos];
					const uint32_t *pos_list = p->pos_list;

					for (int i = 0; i < AUX_TABLE_COLS; i++) {
						const uint32_t pos = pos_list[i];

						if (pos == 0) break;

						const uint32_t read_pos = pos - offset;
						ref_hit_contexts[n_ref_hits++] = (kmer_context){.kmer = kmer,
						                                                .position = read_pos,
						                                                .kmer_pos = pos,
																		.modified_pos = NO_MODIFICATION,
#if DEBUG
						                                                .is_neighbor = false
#endif
						                                            };
						index_table_add(&index_table, read_pos, pos, index_2_kmer_pos_set, false);
					}
				} else {
					assert(0);
				}
			}
#if DEBUG
			else if (orig_ref_hit_not_null && ref_hit->pos == POS_AMBIGUOUS) {
				++ambig_hits;
			}
#endif

			if (orig_snp_hit_not_null && snp_hit->pos != POS_AMBIGUOUS) {
				if (snp_hit->ambig_flag == FLAG_UNAMBIGUOUS) {
					const uint32_t read_pos = snp_hit->pos - offset;
					snp_hit_contexts[n_snp_hits++] = (kmer_context){.kmer = kmer,
					                                                .position = read_pos,
					                                                .kmer_pos = snp_hit->pos,
																	.modified_pos = NO_MODIFICATION,
#if DEBUG
					                                                .is_neighbor = false
#endif
					                                            };
					index_table_add(&index_table, read_pos, snp_hit->pos, index_2_kmer_pos_set, false);
#if DEBUG
					++unambig_hits;
#endif
				} else if (snp_hit->ambig_flag == FLAG_AMBIGUOUS) {
					const struct snp_aux_table *p = &snp_aux_table[snp_hit->pos];
					const uint32_t *pos_list = p->pos_list;

					for (int i = 0; i < AUX_TABLE_COLS; i++) {
						const uint32_t pos = pos_list[i];

						if (pos == 0) break;

						const uint32_t read_pos = pos - offset;
						snp_hit_contexts[n_snp_hits++] = (kmer_context){.kmer = kmer,
						                                                .position = read_pos,
						                                                .kmer_pos = pos,
																		.modified_pos = NO_MODIFICATION,
#if DEBUG
						                                                .is_neighbor = false
#endif
						                                            };
						index_table_add(&index_table, read_pos, pos, index_2_kmer_pos_set, false);
					}
				} else {
					assert(0);
				}
			}
#if DEBUG
			else if (orig_snp_hit_not_null && snp_hit->pos == POS_AMBIGUOUS) {
				++ambig_hits;
			}
#endif
			if (qual_char - QUALITY_SCORE >= 0) continue;
			
			/* loop over hamming neighbors of `kmer`, maybe */
			uint32_t ref_search_bound = 64;
			uint32_t snp_search_bound = 64;

#if REF_LITE
			uint64_t ref_kmer_lo = LO40(kmer);
			if (ref_bf->check_value(ref_kmer_lo) == 0) ref_search_bound = 40;
#else
			uint32_t ref_kmer_lo = LO(kmer);
			if (ref_bf->check_value(ref_kmer_lo) == 0) ref_search_bound = 32;
#endif
			uint64_t snp_kmer_lo = LO40(kmer);
			if (snp_bf->check_value(snp_kmer_lo) == 0) snp_search_bound = 40;

			// uint32_t max_search_bound = max(ref_search_bound, snp_search_bound);

			// this for loop can be dismissed if no more than 100 kmers in this block
			// instead, there will be another loop
			if (block_size >= BLOCK_SIZE_THRESHOLD)
			{
				/* loop over right half of int64 / left half of kmer */
				for (unsigned i = 0; i < 32; i += 2) {
					const unsigned diff_base_pos = i / 2;
					const uint64_t mask = 0x3UL << i; // 0x3UL == 0b11
					const uint64_t base = (kmer & mask) >> i;

					for (uint64_t j = 0; j < 0x4; j++) {
						if (j == base) continue;

						const kmer_t neighbor = (kmer & ~mask) | (j << i);

						struct kmer_entry *ref_hit = query_ref_dict(neighbor, ref_jumpgate, ref_dict, ref_dict_size);

						struct snp_kmer_entry *snp_hit = query_snp_dict(neighbor, snp_jumpgate, snp_dict, snp_dict_size);

						const size_t ref_hit_diff_loc = (ref_hit != NULL &&
							ref_hit->pos != POS_AMBIGUOUS &&
							ref_hit->ambig_flag == FLAG_UNAMBIGUOUS) ?
							(ref_hit->pos + diff_base_pos) :
							0;

						if (ref_hit != NULL && ref_hit->pos != POS_AMBIGUOUS) {
							if (ref_hit->ambig_flag == FLAG_UNAMBIGUOUS &&
#if PCOMPACT
								ptable_get(&ptable, ref_hit_diff_loc) == NULL
#else
								pileup_table[ref_hit_diff_loc].ref == 0 &&
								pileup_table[ref_hit_diff_loc].alt == 0
#endif
								) {

								const uint32_t read_pos = ref_hit->pos - offset;
								ref_hit_contexts[n_ref_hits++] = (kmer_context) {
									.kmer = neighbor,
										.position = read_pos,
										.kmer_pos = ref_hit->pos,
										.modified_pos = diff_base_pos,
#if DEBUG
										.is_neighbor = true
#endif
								};
								index_table_add(&index_table, read_pos, ref_hit->pos, index_2_kmer_pos_set);
								//benchmark_ref_neighbors.insert(neighbor);
#if DEBUG
								++unambig_hits;
#endif
							}
							else if (ref_hit->ambig_flag == FLAG_AMBIGUOUS) {
								const struct aux_table *p = &ref_aux_table[ref_hit->pos];
								const uint32_t *pos_list = p->pos_list;

								for (int i = 0; i < AUX_TABLE_COLS; i++) {
									const uint32_t pos = pos_list[i];

									if (pos == 0) break;

									const size_t ref_hit_diff_loc = pos + diff_base_pos;
									if (
#if PCOMPACT
										ptable_get(&ptable, ref_hit_diff_loc) == NULL
#else
										pileup_table[ref_hit_diff_loc].ref == 0 &&
										pileup_table[ref_hit_diff_loc].alt == 0
#endif
										) {
										const uint32_t read_pos = pos - offset;
										ref_hit_contexts[n_ref_hits++] = (kmer_context) {
											.kmer = neighbor,
												.position = read_pos,
												.kmer_pos = pos,
												.modified_pos = diff_base_pos,
#if DEBUG
												.is_neighbor = true
#endif
										};
										index_table_add(&index_table, read_pos, pos, index_2_kmer_pos_set);
										//benchmark_ref_neighbors.insert(neighbor);
									}
								}
							}
						}
#if DEBUG
						else if (ref_hit != NULL && ref_hit->pos == POS_AMBIGUOUS) {
							++ambig_hits;
						}
#endif

						if (snp_hit != NULL && snp_hit->pos != POS_AMBIGUOUS) {

							if (snp_hit->ambig_flag == FLAG_UNAMBIGUOUS && SNP_INFO_POS(snp_hit->snp) != diff_base_pos) {
								const uint32_t read_pos = snp_hit->pos - offset;
								snp_hit_contexts[n_snp_hits++] = (kmer_context) {
									.kmer = neighbor,
										.position = read_pos,
										.kmer_pos = snp_hit->pos,
										.modified_pos = diff_base_pos,
#if DEBUG
										.is_neighbor = true
#endif
								};
								index_table_add(&index_table, read_pos, snp_hit->pos, index_2_kmer_pos_set);
								//benchmark_snp_neighbors.insert(neighbor);
#if DEBUG
								++unambig_hits;
#endif
							}
							else if (snp_hit->ambig_flag == FLAG_AMBIGUOUS) {
								const struct snp_aux_table *p = &snp_aux_table[snp_hit->pos];
								const uint32_t *pos_list = p->pos_list;
								const uint8_t *snp_list = p->snp_list;

								for (int i = 0; i < AUX_TABLE_COLS; i++) {
									const uint32_t pos = pos_list[i];

									if (pos == 0) break;

									if (SNP_INFO_POS(snp_list[i]) != diff_base_pos) {
										const uint32_t read_pos = pos - offset;
										snp_hit_contexts[n_snp_hits++] = (kmer_context) {
											.kmer = neighbor,
												.position = read_pos,
												.kmer_pos = pos,
												.modified_pos = diff_base_pos,
#if DEBUG
												.is_neighbor = true
#endif
										};

										index_table_add(&index_table, read_pos, pos, index_2_kmer_pos_set);
										//benchmark_snp_neighbors.insert(neighbor);
									}
								}
							}
						}
#if DEBUG
						else if (snp_hit != NULL && snp_hit->pos == POS_AMBIGUOUS) {
							++ambig_hits;
						}
#endif
					} // end of for (uint64_t j = 0; j < 0x4; j++)
				}
			} // end of if (block_size >= BLOCK_SIZE_THRESHOLD)
			else // this else corresponds to if (block_size >= BLOCK_SIZE_THRESHOLD)
			{
				// check every hit
				struct kmer_entry *ref_hit_set[BLOCK_SIZE_THRESHOLD];
				struct snp_kmer_entry *snp_hit_set[BLOCK_SIZE_THRESHOLD];
				uint64_t ref_neighbors[BLOCK_SIZE_THRESHOLD];
				uint64_t snp_neighbors[BLOCK_SIZE_THRESHOLD];

				int ref_diff_base_pos_list[BLOCK_SIZE_THRESHOLD];
				int snp_diff_base_pos_list[BLOCK_SIZE_THRESHOLD];

				int ref_hit_set_size = 0;
				int snp_hit_set_size = 0;

				iterate_ref_dict(kmer, ref_jumpgate, ref_dict, ref_dict_size, ref_hit_set, ref_hit_set_size, ref_neighbors, ref_diff_base_pos_list);

				for (int h = 0; h < ref_hit_set_size; h++) {
					// get every ref_hit
					struct kmer_entry *ref_hit = ref_hit_set[h];
					int diff_base_pos = ref_diff_base_pos_list[h];

					const size_t ref_hit_diff_loc = (ref_hit != NULL && ref_hit->pos != POS_AMBIGUOUS && ref_hit->ambig_flag == FLAG_UNAMBIGUOUS) ? (ref_hit->pos + diff_base_pos) : 0;
					// get the neighbor
					auto neighbor = ref_neighbors[h];

					if (ref_hit == NULL || ref_hit->pos == POS_AMBIGUOUS) continue;
					if (ref_hit->ambig_flag == FLAG_UNAMBIGUOUS &&
#if PCOMPACT
						ptable_get(&ptable, ref_hit_diff_loc) == NULL
#else
						pileup_table[ref_hit_diff_loc].ref == 0 &&
						pileup_table[ref_hit_diff_loc].alt == 0
#endif
						) {
						const uint32_t read_pos = ref_hit->pos - offset;
						ref_hit_contexts[n_ref_hits++] = (kmer_context) { .kmer = neighbor, .position = read_pos, .kmer_pos = ref_hit->pos, .modified_pos = diff_base_pos};
						index_table_add(&index_table, read_pos, ref_hit->pos, index_2_kmer_pos_set);
						//query_ref_neighbors.insert(neighbor);
					}
					else if (ref_hit->ambig_flag == FLAG_AMBIGUOUS) {
						const struct aux_table *p = &ref_aux_table[ref_hit->pos];
						const uint32_t *pos_list = p->pos_list;
						for (int i = 0; i < AUX_TABLE_COLS; i++) {
							const uint32_t pos = pos_list[i];
							if (pos == 0) break;
							const size_t ref_hit_diff_loc = pos + diff_base_pos;
							if (
#if PCOMPACT
								ptable_get(&ptable, ref_hit_diff_loc) == NULL
#else
								pileup_table[ref_hit_diff_loc].ref == 0 &&
								pileup_table[ref_hit_diff_loc].alt == 0
#endif
								) {
								const uint32_t read_pos = pos - offset;
								ref_hit_contexts[n_ref_hits++] = (kmer_context) {.kmer = neighbor, .position = read_pos, .kmer_pos = pos, .modified_pos = diff_base_pos};
								index_table_add(&index_table, read_pos, pos, index_2_kmer_pos_set);
								//query_ref_neighbors.insert(neighbor);
							}
						}
					}
				}
				// end of for (int h = 0; h < ref_hit_set_size; h++)

				iterate_snp_dict(kmer, snp_jumpgate, snp_dict, snp_dict_size, snp_hit_set, snp_hit_set_size, snp_neighbors, snp_diff_base_pos_list);

				for (int h = 0; h < snp_hit_set_size; h++) {
					struct snp_kmer_entry *snp_hit = snp_hit_set[h];
					if (snp_hit == NULL || snp_hit->pos == POS_AMBIGUOUS) continue;
					auto neighbor = snp_neighbors[h];
					int diff_base_pos = snp_diff_base_pos_list[h];

					if (snp_hit->ambig_flag == FLAG_UNAMBIGUOUS && SNP_INFO_POS(snp_hit->snp) != diff_base_pos) {
						const uint32_t read_pos = snp_hit->pos - offset;
						snp_hit_contexts[n_snp_hits++] = (kmer_context) {.kmer = neighbor, .position = read_pos, .kmer_pos = snp_hit->pos, .modified_pos = diff_base_pos};
						index_table_add(&index_table, read_pos, snp_hit->pos, index_2_kmer_pos_set);
						//query_snp_neighbors.insert(neighbor);
					}
					else if (snp_hit->ambig_flag == FLAG_AMBIGUOUS) {
						const struct snp_aux_table *p = &snp_aux_table[snp_hit->pos];
						const uint32_t *pos_list = p->pos_list;
						const uint8_t *snp_list = p->snp_list;

						for (int i = 0; i < AUX_TABLE_COLS; i++) {
							const uint32_t pos = pos_list[i];
							if (pos == 0) break;
							if (SNP_INFO_POS(snp_list[i]) != diff_base_pos) {
								const uint32_t read_pos = pos - offset;
								snp_hit_contexts[n_snp_hits++] = (kmer_context) {.kmer = neighbor, .position = read_pos, .kmer_pos = pos, .modified_pos = diff_base_pos};
								index_table_add(&index_table, read_pos, pos, index_2_kmer_pos_set);
								//query_snp_neighbors.insert(neighbor);
							}
						}
					}
				}
			} // end of else corresponds to if (block_size >= BLOCK_SIZE_THRESHOLD)

			
			/* loop over left/high half of int64 / right half of kmer */
			for (unsigned i = 32; i < 64; i += 2) {

				if(i >= ref_search_bound && i >= snp_search_bound) break;

				const unsigned diff_base_pos = i / 2;
				const uint64_t mask = 0x3UL << i;
				const uint64_t base = (kmer & mask) >> i;

				for (uint64_t j = 0; j < 0x4; j++) {
					if (j == base) continue;

					const kmer_t neighbor = (kmer & ~mask) | (j << i);

					// search ref
					//[TODO] here need code review to review it again

					if ( (i < ref_search_bound) && (is_ref_lite == false || i>=40 || block_size >= BLOCK_SIZE_THRESHOLD) ) {
						
                        struct kmer_entry *ref_hit = query_ref_dict(neighbor, ref_jumpgate, ref_dict, ref_dict_size);
						const size_t ref_hit_diff_loc = (ref_hit != NULL &&
							ref_hit->pos != POS_AMBIGUOUS &&
							ref_hit->ambig_flag == FLAG_UNAMBIGUOUS) ?
							(ref_hit->pos + diff_base_pos) :
							0;

						if (ref_hit != NULL && ref_hit->pos != POS_AMBIGUOUS) {
							if (ref_hit->ambig_flag == FLAG_UNAMBIGUOUS &&
#if PCOMPACT
								ptable_get(&ptable, ref_hit_diff_loc) == NULL
#else
								pileup_table[ref_hit_diff_loc].ref == 0 &&
								pileup_table[ref_hit_diff_loc].alt == 0
#endif
								) {

								const uint32_t read_pos = ref_hit->pos - offset;
								ref_hit_contexts[n_ref_hits++] = (kmer_context) {
									.kmer = neighbor,
										.position = read_pos,
										.kmer_pos = ref_hit->pos,
										.modified_pos = diff_base_pos,
#if DEBUG
										.is_neighbor = true
#endif
								};
								index_table_add(&index_table, read_pos, ref_hit->pos, index_2_kmer_pos_set);
#if DEBUG
								++unambig_hits;
#endif
							}
							else if (ref_hit->ambig_flag == FLAG_AMBIGUOUS) {
								const struct aux_table *p = &ref_aux_table[ref_hit->pos];
								const uint32_t *pos_list = p->pos_list;

								for (int i = 0; i < AUX_TABLE_COLS; i++) {
									const uint32_t pos = pos_list[i];

									if (pos == 0) break;

									const size_t ref_hit_diff_loc = pos + diff_base_pos;
									if (
#if PCOMPACT
										ptable_get(&ptable, ref_hit_diff_loc) == NULL
#else
										pileup_table[ref_hit_diff_loc].ref == 0 &&
										pileup_table[ref_hit_diff_loc].alt == 0
#endif
										) {
										const uint32_t read_pos = pos - offset;
										ref_hit_contexts[n_ref_hits++] = (kmer_context) {
											.kmer = neighbor,
												.position = read_pos,
												.kmer_pos = pos,
												.modified_pos = diff_base_pos,
#if DEBUG
												.is_neighbor = true
#endif
										};
										index_table_add(&index_table, read_pos, pos, index_2_kmer_pos_set);
									}
								}
							}
						}
#if DEBUG
						else if (ref_hit != NULL && ref_hit->pos == POS_AMBIGUOUS) {
							++ambig_hits;
						}
#endif
					}

					// search snp

					/* if block_size >= BLOCK_SIZE_THRESHOLD, we need to iterate the rest 32
					* else, direct match applied, we only need to iterate the rest high 24, which is in range [40 : 64)
					* */
					if (block_size >= BLOCK_SIZE_THRESHOLD || i >= 40) {

						if (i >= snp_search_bound) continue;

						struct snp_kmer_entry *snp_hit = query_snp_dict(neighbor, snp_jumpgate, snp_dict, snp_dict_size);

						if (snp_hit != NULL && snp_hit->pos != POS_AMBIGUOUS) {

							if (snp_hit->ambig_flag == FLAG_UNAMBIGUOUS && SNP_INFO_POS(snp_hit->snp) != diff_base_pos) {
								const uint32_t read_pos = snp_hit->pos - offset;
								snp_hit_contexts[n_snp_hits++] = (kmer_context) {
									.kmer = neighbor,
										.position = read_pos,
										.kmer_pos = snp_hit->pos,
										.modified_pos = diff_base_pos,
#if DEBUG
										.is_neighbor = true
#endif
								};
								index_table_add(&index_table, read_pos, snp_hit->pos, index_2_kmer_pos_set);
#if DEBUG
								++unambig_hits;
#endif
							}
							else if (snp_hit->ambig_flag == FLAG_AMBIGUOUS) {
								const struct snp_aux_table *p = &snp_aux_table[snp_hit->pos];
								const uint32_t *pos_list = p->pos_list;
								const uint8_t *snp_list = p->snp_list;

								for (int i = 0; i < AUX_TABLE_COLS; i++) {
									const uint32_t pos = pos_list[i];

									if (pos == 0) break;

									if (SNP_INFO_POS(snp_list[i]) != diff_base_pos) {
										const uint32_t read_pos = pos - offset;
										snp_hit_contexts[n_snp_hits++] = (kmer_context) {
											.kmer = neighbor,
												.position = read_pos,
												.kmer_pos = pos,
												.modified_pos = diff_base_pos,
#if DEBUG
												.is_neighbor = true
#endif
										};

										index_table_add(&index_table, read_pos, pos, index_2_kmer_pos_set);
									}
							}
						}
					}
#if DEBUG
						else if (snp_hit != NULL && snp_hit->pos == POS_AMBIGUOUS) {
							++ambig_hits;
						}
#endif
				} // end of if(block_size >= BLOCK_SIZE_THRESHOLD || i >= 40){
			} // end of loop over all neighbor nt(A,T,G,C)
			}// end of loop nt positions

		} // end of loop over all kmers in single read

		/*
		 * Now we loop over our ref/SNP hits and find the ones that support the 'best' position
		 * according to our index table, and use those to update the pileup table. At the same
		 * time, we clear our index table for when we process the next read.
		 */

		const bool process_read = (index_table.best && (index_table.best->freq > 1) && !index_table.ambiguous);
		const uint32_t target_index = index_table.best ? index_table.best->index : 0;

#if DEBUG
		bool read_good = false;
#endif

		for (size_t i = 0; i < n_ref_hits; i++) {
			const uint32_t index = ref_hit_contexts[i].position;
			index_table_clear_index(&index_table, index);

			if (process_read && index == target_index) {
				const uint32_t kmer_pos = ref_hit_contexts[i].kmer_pos;
				const kmer_t kmer = ref_hit_contexts[i].kmer;
				const uint32_t modified_pos = ref_hit_contexts[i].modified_pos;
				for (unsigned i = 0; i < 32; i++) {
					if(i == modified_pos) continue;
					const unsigned base = kmer_get_base(kmer, i);

#if PCOMPACT
					struct pileup_entry *p = ptable_get(&ptable, kmer_pos + i);
#else
					struct packed_pileup_entry *p = &pileup_table[kmer_pos + i];
#endif

					if (
#if PCOMPACT
					    p != NULL
#else
					    p->ref != p->alt
#endif
					   ) {

						if (base == p->ref) {
#if DEBUG
							read_good = true;
#endif
							if (p->ref_cnt != MAX_COV)
								++p->ref_cnt;
#if DEBUG
							++ref_covs;
#endif
						}
						else if (base == p->alt) {
#if DEBUG
							read_good = true;
#endif
							if (p->alt_cnt != MAX_COV)
								++p->alt_cnt;
#if DEBUG
							++alt_covs;
#endif
						}
#if DEBUG
						else {
							++non_ref_or_alt_covs;
						}
#endif
					}
#if DEBUG && !PCOMPACT
					else {
						assert(p->ref == 0 && p->alt == 0);
					}
#endif
				}
			}
		}

		for (size_t i = 0; i < n_snp_hits; i++) {
			const uint32_t index = snp_hit_contexts[i].position;
			index_table_clear_index(&index_table, index);

			if (process_read && index == target_index) {
				const uint32_t kmer_pos = snp_hit_contexts[i].kmer_pos;
				const kmer_t kmer = snp_hit_contexts[i].kmer;
				const uint32_t modified_pos = snp_hit_contexts[i].modified_pos;
				for (unsigned i = 0; i < 32; i++) {
					if(i == modified_pos) continue;
					const unsigned base = kmer_get_base(kmer, i);

#if PCOMPACT
					struct pileup_entry *p = ptable_get(&ptable, kmer_pos + i);
#else
					struct packed_pileup_entry *p = &pileup_table[kmer_pos + i];
#endif

					if (
#if PCOMPACT
					    p != NULL
#else
					    p->ref != p->alt
#endif
					   ) {

						if (base == p->ref) {
#if DEBUG
							read_good = true;
#endif
							if (p->ref_cnt != MAX_COV)
								++p->ref_cnt;
#if DEBUG
							++ref_covs;
#endif
						}
						else if (base == p->alt) {
#if DEBUG
							read_good = true;
#endif
							if (p->alt_cnt != MAX_COV)
								++p->alt_cnt;
#if DEBUG
							++alt_covs;
#endif
						}
#if DEBUG
						else {
							++non_ref_or_alt_covs;
						}
#endif
					}
#if DEBUG && !PCOMPACT
					else {
						assert(p->ref == 0 && p->alt == 0);
					}
#endif
				}
			}
		}

		if (!process_read && !revcompl) {
			revcompl = true;
			index_table.best = NULL;
			index_table.ambiguous = false;
			index_2_kmer_pos_set.clear();
			goto head;
		}

#if DEBUG
		if (read_good)
			++good_reads;
		else
			++bad_reads;

		++total_count;

		if (index_table.best) {
			fprintf(read_data, "%s %d ", index_table.ambiguous ? "A" : "U", index_table.best->freq);

			for (size_t i = 0; i < n_ref_hits; i++) {
				const uint32_t index = ref_hit_contexts[i].position;

				if (index == target_index) {
					fprintf(read_data, "%u:%s ", index, ref_hit_contexts[i].is_neighbor ? "1" : "0");
				}
			}

			for (size_t i = 0; i < n_snp_hits; i++) {
				const uint32_t index = snp_hit_contexts[i].position;

				if (index == target_index) {
					fprintf(read_data, "%u:%s ", index, snp_hit_contexts[i].is_neighbor ? "1" : "0");
				}
			}

			fprintf(read_data, "\n");
		}

		if (index_table.best != NULL && index_table.best->freq > 1 && !index_table.ambiguous) {
			++match_count;
		} else {
			if (index_table.best != NULL && index_table.best->freq > 1 && index_table.ambiguous) {
				++multi_count;
			}
			else {
				++nohit_count;
			}
		}
#endif

		//nohit:
		index_table.best = NULL;
		index_table.ambiguous = false;
		index_2_kmer_pos_set.clear();
	}

#if DEBUG
	fclose(read_data);
#endif

	size_t ref_call_count = 0;
	size_t alt_call_count = 0;
	size_t het_call_count = 0;

#if PCOMPACT
	struct pileup_entry **table = ptable.table;
	const size_t pileup_size = ptable.size;
#endif

	for (size_t i = 0; i < pileup_size; i++) {
#if PCOMPACT
		for (struct pileup_entry *p = table[i]; p != NULL; p = p->next) {
#else
		struct packed_pileup_entry *p = &pileup_table[i];
#endif

			if (p->ref == p->alt) {
				continue;  // no SNP here
			}

#if PCOMPACT
			size_t index = p->key;
#else
			size_t index = i;
#endif

			/* index w.r.t. correct chromosome */
			size_t j;
			for (j = 0; j < num_chrs && index > chrlens[j].len; j++) {
				index -= chrlens[j].len;
			}

			const struct call call = choose_best_genotype(p->ref_cnt, p->alt_cnt, p->ref_freq, p->alt_freq);

			switch (call.genotype) {
			case GTYPE_NONE:
				break;
			case GTYPE_REF:
				++ref_call_count;
				fprintf(out, "%s %lu 0/0 %.15g\n", chrlens[j].name, index, call.confidence);
				break;
			case GTYPE_ALT:
				++alt_call_count;
				fprintf(out, "%s %lu 1/1 %.15g\n", chrlens[j].name, index, call.confidence);
				break;
			case GTYPE_HET:
				++het_call_count;
				fprintf(out, "%s %lu 0/1 %.15g\n", chrlens[j].name, index, call.confidence);
				break;
			}

#if PCOMPACT
		}
#endif
	}

	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Time: %f sec\n", time_spent);

#if DEBUG
	/*
	static const char bases[] = {'A', 'C', 'G', 'T'};
	FILE *counts = fopen("counts.txt", "w");
	FILE *all_snps = fopen("all_snps.txt", "w");
	for (size_t i = 0; i < pileup_size; i++) {
		struct pileup_entry *p = &pileup_table[i];
		if (p->ref == p->alt) {
			assert(p->ref == 0 && p->alt == 0);
			continue;  // no SNP here
		}

		size_t index = i;

		// index w.r.t. correct chromosome
		size_t j;
		for (j = 0; j < num_chrs && index > chrlens[j].len; j++) {
			index -= chrlens[j].len;
		}

		if (p->ref_cnt != 0 || p->alt_cnt != 0) {
			fprintf(counts, "%s %lu (%c:%f / %c:%f) : %u / %u\n",
			                chrlens[j].name,
			                index,
			                bases[p->ref],
			                p->ref_freq/255.0f,
			                bases[p->alt],
			                p->alt_freq/255.0f,
			                p->ref_cnt,
			                p->alt_cnt);
		}
		fprintf(all_snps, "%s %lu\n", chrlens[j].name, index);
	}

	fclose(all_snps);
	fclose(counts);
	*/

	printf("Total: %lu\n", total_count);
	printf("Match: %lu\n", match_count);
	printf("Multi: %lu\n", multi_count);
	printf("NoHit: %lu\n", nohit_count);
	printf("\n");
	printf("Unambig. hits: %lu\n", unambig_hits);
	printf("Ambig. hits:   %lu\n", ambig_hits);
	printf("\n");
	printf("Good reads: %lu\n", good_reads);
	printf("Bad reads: %lu\n", bad_reads);
	printf("\n");
	printf("Ref calls: %lu\n", ref_call_count);
	printf("Alt calls: %lu\n", alt_call_count);
	printf("Het calls: %lu\n", het_call_count);
	printf("\n");
	printf("Ref covs:         %lu\n", ref_covs);
	printf("Alt covs:         %lu\n", alt_covs);
	printf("Non ref/alt covs: %lu\n", non_ref_or_alt_covs);
#endif

	free(ref_jumpgate);
	free(ref_dict);
	free(ref_aux_table);
	free(snp_jumpgate);
	free(snp_dict);
	free(snp_aux_table);

#if PCOMPACT
	ptable_dealloc(&ptable);
#else
	free(pileup_table);
#endif
}

static inline struct call choose_best_genotype(const int ref_cnt,
                                               const int alt_cnt,
                                               const uint8_t ref_freq_enc,
                                               const uint8_t alt_freq_enc)
{
	static struct {
		double g0;  // P(counts|G0)/(ref_cnt + alt_cnt choose ref_cnt)
		double g1;  // P(counts|G1)/(ref_cnt + alt_cnt choose ref_cnt)
		double g2;  // P(counts|G2)/(ref_cnt + alt_cnt choose ref_cnt)
	} cache[MAX_COV + 1][MAX_COV + 1];

	static double poisson[2*MAX_COV + 1];

	static bool init = false;

	if (!init) {
		for (int ref_cnt = 0; ref_cnt <= MAX_COV; ref_cnt++) {
			for (int alt_cnt = 0; alt_cnt <= MAX_COV; alt_cnt++) {
				cache[ref_cnt][alt_cnt].g0 = pow(1.0 - ERR_RATE, ref_cnt) * pow(ERR_RATE, alt_cnt);
				cache[ref_cnt][alt_cnt].g1 = pow(0.5, ref_cnt + alt_cnt);
				cache[ref_cnt][alt_cnt].g2 = pow(ERR_RATE, ref_cnt) * pow(1.0 - ERR_RATE, alt_cnt);
			}
		}

		const double M = exp(-AVG_COV);
		for (int i = 0; i <= (2*MAX_COV); i++) {
			poisson[i] = (M*pow(AVG_COV, i))/exp(lgamma(i+1.0));
		}

		init = true;
	}

	if ((ref_cnt == 0 && alt_cnt == 0) || (ref_cnt == MAX_COV && alt_cnt == MAX_COV)) {
		return CALL(GTYPE_NONE, 0.0);
	}

	const double g0 = cache[ref_cnt][alt_cnt].g0;
	const double g1 = cache[ref_cnt][alt_cnt].g1;
	const double g2 = cache[ref_cnt][alt_cnt].g2;

	const double p = ref_freq_enc/255.0;
	const double q = alt_freq_enc/255.0;
	const double p2 = p*p;
	const double q2 = q*q;

	const double p_g0 = p2*g0;
	const double p_g1 = (1.0 - p2 - q2)*g1;
	const double p_g2 = q2*g2;
	const double total = p_g0 + p_g1 + p_g2;

	const int n = ref_cnt + alt_cnt;

	if (p_g0 > p_g1 && p_g0 > p_g2) {
		return CALL(GTYPE_REF, ((double)(p_g0/total))*poisson[n]);
	} else if (p_g1 > p_g0 && p_g1 > p_g2) {
		return CALL(GTYPE_HET, ((double)(p_g1/total))*poisson[n]);
	} else {
		return CALL(GTYPE_ALT, ((double)(p_g2/total))*poisson[n]);
	}
}


/* === Front-End === */

static void print_help(void)
{
	fprintf(stderr, "Usage: vargeno_lite <option> [option parameters ...]\n");
	fprintf(stderr, "Option  Description                   Parameters\n");
	fprintf(stderr, "------  -----------                   ----------\n");
	fprintf(stderr, "ucscd    Generate dictionary files, if known SNPs are in UCSC text file format     "
	                "<input FASTA> <input SNPs> <output ref dict> <output SNP dict>\n");
	fprintf(stderr, "vcfd    Generate dictionary files, if known SNPs are in VCF file format    "
	                "<input FASTA> <input SNPs> <output ref dict> <output SNP dict>\n");
	fprintf(stderr, "filt    Filter reference dictionary   "
		            "<ref dict> <snp_pos file> <output ref dict>\n");
	fprintf(stderr, "geno    Perform genotyping            "
	                "<input ref dict> <input SNP dict> <input FASTQ> <chrlens file> <ref Bloom filter> <snp Bloom filter> <output file>\n");
	fprintf(stderr, "------  -----------                   ----------\n");
	fprintf(stderr, "Note: to generate Bloom filters, please use command 'gbf_lite'\n");
}

static void arg_check(int argc, int expected)
{
	if (argc - 2 != expected) {  /* -2 because option params start at argv[2] */
		print_help();
		exit(EXIT_FAILURE);
	}
}

void printint(uint64_t x) {
	std::bitset<64> y(x);
	cout << y << endl;
}

int main(const int argc, const char *argv[])
{
	if (argc < 2) {
		print_help();
		exit(EXIT_FAILURE);
	}

	const char *opt = argv[1];

	if (STREQ(opt, "vcfd")) {
		cout << "VCF format support coming soon." << endl;
	}
	if (STREQ(opt, "ucscd")) {
		arg_check(argc, 4);
		const char *ref_filename = argv[2];
		const char *snp_filename = argv[3];
		const char *refdict_filename = argv[4];
		const char *snpdict_filename = argv[5];

		SeqVec ref = parse_fasta(ref_filename);

#define CHRLENS_EXT ".chrlens"
		char chrlens_filename[4096];
		assert(strlen(ref_filename) < (sizeof(chrlens_filename) - strlen(CHRLENS_EXT)));
		sprintf(chrlens_filename, "%s" CHRLENS_EXT, ref_filename);
		FILE *chrlens = fopen(chrlens_filename, "w");
		assert(chrlens);
#undef CHRLENS_EXT

		for (size_t i = 0; i < ref.size; i++) {
			fprintf(chrlens, "%s %lu\n", ref.seqs[i].name, ref.seqs[i].size);
		}

		fclose(chrlens);

		FILE *snp_file = fopen(snp_filename, "r");
		assert(snp_file);

		FILE *snpdict_file = fopen(snpdict_filename, "wb");
		assert(snpdict_file);

		bool *snp_locations;
		size_t snp_locs_size;
		make_snp_dict(ref, snp_file, snpdict_file, &snp_locations, &snp_locs_size);
		assert(snp_locations);

#if GEN_FLT_DATA
		FILE *snp_locs = fopen("snp_locs", "wb");
		assert(snp_locs);
		serialize_uint64(snp_locs, snp_locs_size);
		for (size_t i = 0; i < snp_locs_size; i++) {
			serialize_uint8(snp_locs, snp_locations[i]);
		}
		fclose(snp_locs);
#endif

		FILE *refdict_file = fopen(refdict_filename, "wb");
		assert(refdict_file);

		make_ref_dict(ref, refdict_file);

		fclose(refdict_file);
		free(snp_locations);

		seqvec_dealloc(&ref);
		fclose(snp_file);
		fclose(snpdict_file);
	} else if (STREQ(opt, "filt")) {
		arg_check(argc, 3);

		const char *refdict_filename = argv[2];
		const char *snp_pos_filename = argv[3];
		const char *out_filename = argv[4];

		FILE *refdict_file = fopen(refdict_filename, "rb");
		assert(refdict_file);

		FILE *snp_pos_file = fopen(snp_pos_filename, "rb");
		assert(snp_pos_file);

		FILE *out_file = fopen(out_filename, "wb");
		assert(out_file);

		dict_filt(refdict_file, snp_pos_file, out_file);
	} else if (STREQ(opt, "geno")) {

		clock_t start_time, end_time;
        
        start_time = clock();

        arg_check(argc, 7);
		const char *refdict_filename = argv[2];
		const char *snpdict_filename = argv[3];
		const char *fastq_filename = argv[4];
		const char *chrlens_filename = argv[5];
		string ref_bf_filename = argv[6];
		string snp_bf_filename = argv[7];
		const char *out_filename = argv[8];

#if REF_LITE
		ref_bf = new BloomFilter(BFGenerator::SNP_BF_RANGE);
#else
        ref_bf = new BloomFilter(BFGenerator::REF_BF_RANGE);
#endif
        snp_bf = new BloomFilter(BFGenerator::SNP_BF_RANGE);

		ref_bf->load(ref_bf_filename);
		snp_bf->load(snp_bf_filename);

// [TODO] here we'd better check if ref_bf contains 40 bits instead of 32 bits elements


		//FILE *out_file = fopen(out_filename, "w");
		//assert(out_file);

		//initialize mask
		for(int i = 0; i < 32; i++) {
    		odd_mask_64 <<= 2;
    		odd_mask_64 |= 2;
    		even_mask_64 <<= 2;
    		even_mask_64 |= 1;

		  if(i < 16){
			odd_mask_32 <<= 2;
			odd_mask_32 |= 2;
			even_mask_32 <<= 2;
			even_mask_32 |= 1;
		  }
		}

		//initialize diff base dictionary
		uint64_t base_mask_1 = 1;
		uint64_t base_mask_2 = 2;
		uint64_t base_mask_3 = 3;
		for (int left_shift_time = 0; left_shift_time < 32; left_shift_time++) {
		
			// if two kmers a and b has only one diff base, then a^b only have 01, 10, 11 for that base
			diff_base_dictionary[base_mask_1] = left_shift_time;
			diff_base_dictionary[base_mask_2] = left_shift_time;
			diff_base_dictionary[base_mask_3] = left_shift_time;
			//
			//cout << "position: " << left_shift_time << endl;
			//printint(base_mask_1);
			//printint(base_mask_2);
			//printint(base_mask_3);
			//
			base_mask_1 <<= 2;
			base_mask_2 <<= 2;
			base_mask_3 <<= 2;
		}


        FILE *refdict_file = fopen(refdict_filename, "rb");
	    assert(refdict_file);

	    FILE *snpdict_file = fopen(snpdict_filename, "rb");
	    assert(snpdict_file);

	    FILE *fastq_file = fopen(fastq_filename, "r");
	    assert(fastq_file);

	    FILE *chrlens_file = fopen(chrlens_filename, "r");
	    assert(chrlens_file);
		
        FILE *out_file = fopen(out_filename, "w");
		assert(out_file);
		
        genotype(refdict_file, snpdict_file, fastq_file, chrlens_file, out_file);
		    
        fclose(refdict_file);
	    fclose(snpdict_file);
	    fclose(fastq_file);
	    fclose(chrlens_file);
		fclose(out_file);

		if (ref_bf != NULL) {
			delete ref_bf;
			ref_bf = NULL;
		}

		if (snp_bf != NULL) {
			delete snp_bf;
			snp_bf = NULL;
		}

		//fclose(refdict_file);
		//fclose(snpdict_file);
		//fclose(fastq_file);
		//fclose(chrlens_file);
		//fclose(out_file);
	} else if (STREQ(opt, "help")) {
		print_help();
		exit(EXIT_SUCCESS);
	} else {
		print_help();
		exit(EXIT_FAILURE);
	}
}

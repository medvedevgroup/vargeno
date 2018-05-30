#ifndef LAVA_H
#define LAVA_H

#include <stdint.h>

/////////////////////////////
#define DEBUG        0
#define REF_LITE     1
#define PCOMPACT     0
#define GEN_FLT_DATA 0

#define READ_LEN       101
#define ERR_RATE       0.01
#define AVG_COV        7.1
#define MAX_MATES_DIST 2000
/////////////////////////////

#define BASE_A 0x0
#define BASE_C 0x1
#define BASE_G 0x2
#define BASE_T 0x3
#define BASE_N 0x4
#define BASE_X 0x7

#define MAX_COV ((1 << 6) - 1)

enum {
	GTYPE_NONE, GTYPE_REF, GTYPE_ALT, GTYPE_HET
};

#define POS_AMBIGUOUS ((uint32_t)(-1))

#define FLAG_UNAMBIGUOUS 0x00
#define FLAG_AMBIGUOUS   0x01

typedef uint64_t kmer_t;   /* k = 32 */

typedef uint8_t snp_info;  /* 5 bits: pos in kmer, 3 bits: ref base */

#define SNP_INFO_MAKE(pos, ref) ((((pos) & 0x1F) << 3) | ((ref) & 0x07))
#define SNP_INFO_POS(snp_info)  (((snp_info) & 0xF8) >> 3)
#define SNP_INFO_REF(snp_info)  ((snp_info) & 0x07)

/* Note:
 * Frequencies (between 0 and 1) are encoded as 8-bit uints as: freq*0xff.
 * They are decoded as: freq_enc/255.0f
 */

struct kmer_info {
	kmer_t kmer;
	uint32_t pos;
} __attribute__((packed));

struct snp_kmer_info {
	kmer_t kmer;
	uint32_t pos;
	snp_info snp;
	uint8_t ref_freq;
	uint8_t alt_freq;
} __attribute__((packed));

struct kmer_entry {
#if REF_LITE
	uint64_t kmer_lo40 : 40;
#else
	uint32_t kmer_lo;
#endif
	uint32_t pos;
	uint8_t ambig_flag;
} __attribute__((packed));

struct snp_kmer_entry {
	uint64_t kmer_lo40 : 40;
	snp_info snp;
	uint32_t pos;
	uint8_t ambig_flag;
} __attribute__((packed));

#if !PCOMPACT
struct packed_pileup_entry {
	unsigned ref : 2;
	unsigned alt : 2;
	unsigned ref_cnt : 6;
	unsigned alt_cnt : 6;
	uint8_t ref_freq;
	uint8_t alt_freq;
} __attribute__((packed));
#endif

/* table for storing multiple positions in case of ambiguous k-mers */
#define AUX_TABLE_COLS 10

#define AUX_TABLE_INIT_SIZE 75000000

struct aux_table {
	uint32_t pos_list[AUX_TABLE_COLS];
} __attribute__((packed));

#define SNP_AUX_TABLE_INIT_SIZE 10000000

#define BLOCK_SIZE_THRESHOLD 100

struct snp_aux_table_dictgen {
	kmer_t kmer;
	uint32_t pos_list[AUX_TABLE_COLS];
	snp_info snp_list[AUX_TABLE_COLS];
	uint8_t ref_freqs[AUX_TABLE_COLS];
	uint8_t alt_freqs[AUX_TABLE_COLS];
} __attribute__((packed));

struct snp_aux_table {
	uint32_t pos_list[AUX_TABLE_COLS];
	snp_info snp_list[AUX_TABLE_COLS];
} __attribute__((packed));

#endif /* LAVA_H */


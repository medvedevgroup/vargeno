#ifndef UTIL_H
#define UTIL_H

#include <stdbool.h>
#include <string.h>
#include "vartype.h"

#define NEG_ONE  ((uint32_t)(-1))
#define POW_2_32 (1UL << 32)
#define POW_2_24 (1UL << 24)

#define UNUSED(x) (void)(x);

#define SIZE(array) (sizeof(array)/sizeof((array)[0]))

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

#define STREQ(str1, str2) (strcmp(str1, str2) == 0)

#define HI(kmer)   (((kmer) & 0xFFFFFFFF00000000) >> 32)
#define LO(kmer)   ((kmer) & 0x00000000FFFFFFFF)

#define HI24(kmer) (((kmer) & 0xFFFFFF0000000000) >> 40)
#define LO40(kmer) ((kmer) & 0x000000FFFFFFFFFF)

void serialize_uint64(FILE *out, const uint64_t x);

void serialize_uint32(FILE *out, const uint32_t x);

void serialize_uint8(FILE *out, const uint8_t x);

uint64_t read_uint64(FILE *in);

uint32_t read_uint32(FILE *in);

uint8_t read_uint8(FILE *in);

int kmer_cmp(const void *p1, const void *p2);

int snp_kmer_cmp(const void *p1, const void *p2);

uint64_t encode_base(const char base);

kmer_t encode_kmer(const char *kmer, bool *kmer_had_n);

kmer_t shift_kmer(const kmer_t kmer, const char next_base);

unsigned kmer_get_base(const kmer_t kmer, unsigned base);

kmer_t rev_compl(const kmer_t orig);

void decode_kmer(const kmer_t kmer, char *buf);

void split_line(const char *str, char **out);

void copy_until_space(char *dest, const char *src);

bool equal_up_to_space(const char *a, const char *b);

#endif /* UTIL_H */


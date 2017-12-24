#ifndef PILEUP_H
#define PILEUP_H

#include <stdlib.h>
#include <stdint.h>

struct pileup_entry {
	unsigned ref : 2;
	unsigned alt : 2;
	unsigned ref_cnt : 6;
	unsigned alt_cnt : 6;
	uint8_t ref_freq;
	uint8_t alt_freq;
	uint32_t key;
	struct pileup_entry *next;
};

typedef struct {
	struct pileup_entry **table;
	size_t count;
	size_t size;
	size_t threshold;
} PileupTable;

void ptable_init(PileupTable *p, const size_t size);
void ptable_dealloc(PileupTable *p);
void ptable_add(PileupTable *p, const uint32_t key,
                unsigned ref, unsigned alt,
                uint8_t ref_freq, uint8_t alt_freq);

/*
 * Adapted from java.util.HashMap
 */
static inline uint32_t hash(uint32_t h)
{
	h ^= (h >> 20) ^ (h >> 12);
	return h ^ (h >> 7) ^ (h >> 4);
}

static inline struct pileup_entry *ptable_get(PileupTable *p, const uint32_t key)
{
	for (struct pileup_entry *e = p->table[hash(key) & (p->size - 1)];
	     e != NULL;
	     e = e->next) {
		if (e->key == key) {
			return e;
		}
	}
	return NULL;
}

#endif /* PILEUP_H */


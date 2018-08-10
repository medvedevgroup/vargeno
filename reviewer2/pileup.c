#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include "pileup.h"

#define LOAD_FACTOR 0.4

void ptable_init(PileupTable *p, const size_t size)
{
	assert((size & (size - 1)) == 0);  // `size` must be a power of 2
	p->table = (struct pileup_entry**)calloc(size, sizeof(*(p->table)));
	assert(p->table);
	p->count = 0;
	p->size = size;
	p->threshold = (size_t)(size * LOAD_FACTOR);
}

void ptable_dealloc(PileupTable *p)
{
	struct pileup_entry **table = p->table;
	const size_t size = p->size;
	for (size_t i = 0; i < size; i++) {
		struct pileup_entry *e = table[i];
		while (e != NULL) {
			struct pileup_entry *temp = e;
			e = e->next;
			free(temp);
		}
	}
	free(table);
}

static void grow(PileupTable *p)
{
	const size_t size = p->size;
	struct pileup_entry **table = p->table;

	const size_t new_size = 2*size;
	struct pileup_entry **new_table = (struct pileup_entry**)calloc(new_size, sizeof(*new_table));
	assert(new_table);

	for (size_t i = 0; i < size; i++) {
		struct pileup_entry *e = table[i];
		while (e != NULL) {
			struct pileup_entry *next = e->next;
			const uint32_t n = hash(e->key) & (new_size - 1);
			e->next = new_table[n];
			new_table[n] = e;
			e = next;
		}
	}

	p->table = new_table;
	p->size = new_size;
	p->threshold = (size_t)(new_size * LOAD_FACTOR);
	free(table);
}

void ptable_add(PileupTable *p, const uint32_t key,
                unsigned ref, unsigned alt,
                uint8_t ref_freq, uint8_t alt_freq)
{
	if (ptable_get(p, key) != NULL) {
		return;
	}

	const size_t size = p->size;
	struct pileup_entry **table = p->table;

	const uint32_t n = hash(key) & (size - 1);

	struct pileup_entry *e = (struct pileup_entry *)malloc(sizeof(*e));
	e->ref = ref;
	e->alt = alt;
	e->ref_cnt = 0;
	e->alt_cnt = 0;
	e->ref_freq = ref_freq;
	e->alt_freq = alt_freq;
	e->key = key;
	e->next = table[n];
	table[n] = e;

	if (p->count++ > p->threshold) {
		grow(p);
	}
}


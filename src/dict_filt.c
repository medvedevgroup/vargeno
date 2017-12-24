#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>
#include "vartype.h"
#include "util.h"

static bool ref_kmer_snp_proximity_check(const uint32_t pos, bool *snp_locations, size_t snp_locs_size)
{
	if (pos >= snp_locs_size)
		return false;

	const uint32_t lo = pos > (READ_LEN - 32) ? pos - (READ_LEN - 32) : 0;
	const uint32_t hi = (pos < snp_locs_size - (READ_LEN - 1)) ? pos + (READ_LEN - 1) : snp_locs_size - 1;
	for (uint32_t i = lo; i <= hi; i++) {
		if (snp_locations[i])
			return true;
	}
	return false;
}

void dict_filt(FILE *ref_dict, FILE *snp_pos, FILE *out)
{
	const uint64_t snp_locs_size = read_uint64(snp_pos);
	bool *snp_locations = (bool*)malloc(snp_locs_size);
	assert(snp_locations);
	for (uint64_t i = 0; i < snp_locs_size; i++) {
		snp_locations[i] = read_uint8(snp_pos);
	}
	fclose(snp_pos);

	const uint64_t ref_dict_size = read_uint64(ref_dict);
	const uint64_t ref_aux_table_size = read_uint64(ref_dict);

	serialize_uint64(out, 0UL);
	serialize_uint64(out, ref_aux_table_size);

	size_t removed = 0;
	size_t ref_dict_size_new = ref_dict_size;
	for (uint64_t i = 0; i < ref_dict_size; i++) {
		const kmer_t kmer = read_uint64(ref_dict);
		const uint32_t pos = read_uint32(ref_dict);
		const uint8_t ambig_flag = read_uint8(ref_dict);

		if (pos == POS_AMBIGUOUS ||
		    ambig_flag == FLAG_AMBIGUOUS ||
		    ref_kmer_snp_proximity_check(pos, snp_locations, snp_locs_size)) {

			serialize_uint64(out, kmer);
			serialize_uint32(out, pos);
			serialize_uint8(out, ambig_flag);
		} else {
			--ref_dict_size_new;
			++removed;
		}
	}
	printf("New size: %lu\n", ref_dict_size_new);
	printf("Removed:  %lu/%lu\n", removed, (size_t)ref_dict_size);

	for (uint64_t i = 0; i < ref_aux_table_size; i++) {
		uint32_t aux_row[AUX_TABLE_COLS] = { 0 };
		size_t next = 0;

		for (size_t j = 0; j < AUX_TABLE_COLS; j++) {
			const uint32_t pos = read_uint32(ref_dict);
			aux_row[next++] = pos;
		}

		for (size_t j = 0; j < AUX_TABLE_COLS; j++) {
			serialize_uint32(out, aux_row[j]);
		}
	}

	rewind(out);
	serialize_uint64(out, ref_dict_size_new);
	fclose(out);
	free(snp_locations);
}


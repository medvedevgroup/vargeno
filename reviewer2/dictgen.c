#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include "fasta_parser.h"
#include "util.h"
#include "dictgen.h"

static size_t ref_to_constituent_kmers(struct kmer_info *kmers,
                                       const char *ref,
                                       const size_t ref_len,
                                       uint32_t *index)
{
	assert(ref_len >= 32);
	const size_t kmers_len_max = ref_len - 32 + 1;
	size_t kmers_len_true = 0;

	kmer_t kmer;
	uint32_t index_true = *index;
	bool need_full_encode = true;
	bool kmer_had_n;

	for (size_t i = 0; i < kmers_len_max; i++) {
		const char next_base = ref[i + 32 - 1];

		if (need_full_encode) {
			kmer = encode_kmer(&ref[i], &kmer_had_n);
			need_full_encode = kmer_had_n;
		} else if (next_base != 'N') {
			kmer = shift_kmer(kmer, next_base);
			kmer_had_n = false;
		} else {
			need_full_encode = true;
			kmer_had_n = true;
		}

		if (!kmer_had_n) {
			kmers[kmers_len_true].kmer = kmer;
			kmers[kmers_len_true].pos = index_true;
			++kmers_len_true;
		}

		++index_true;
	}

	*index = index_true + 32 - 1;  // since last k-mer index != last base index
	return kmers_len_true;
}

static void sort_kmers(struct kmer_info *kmers, const size_t kmers_len)
{
	qsort(kmers, kmers_len, sizeof(*kmers), kmer_cmp);
}

static void sort_snp_kmers(struct snp_kmer_info *kmers, const size_t kmers_len)
{
	qsort(kmers, kmers_len, sizeof(*kmers), snp_kmer_cmp);
}

static void write_kmers(struct kmer_info *kmers, const size_t kmers_len, FILE *out)
{
	struct aux_table *aux_table = (struct aux_table*)malloc(AUX_TABLE_INIT_SIZE * sizeof(*aux_table));
	assert(aux_table);
	uint64_t aux_table_count = 0;
	uint64_t aux_table_cap = AUX_TABLE_INIT_SIZE;

	const uint64_t placeholder = 0UL;
	serialize_uint64(out, placeholder);  /* dict size (rows) placeholder */
	serialize_uint64(out, placeholder);  /* aux table size (rows) placeholder */

	uint64_t kmers_written = 0UL;
	size_t i = 0;

	/* keep track of a few statistics */
	const size_t total_kmers = kmers_len;
	size_t unambig_kmers = 0;
	size_t ambig_unique_kmers = 0;
	size_t ambig_total_kmers = 0;

	while (i < kmers_len) {
		const kmer_t kmer = kmers[i].kmer;
		const uint32_t pos = kmers[i].pos;
		serialize_uint64(out, kmer);

		++i;

		if (i < kmers_len && kmer == kmers[i].kmer) {
			++ambig_unique_kmers;
			++ambig_total_kmers;

			if (aux_table_count == aux_table_cap) {
				const size_t new_cap = (aux_table_cap * 3)/2 + 1;
				aux_table = (struct aux_table*)realloc(aux_table, new_cap * sizeof(*aux_table));
				assert(aux_table);
				aux_table_cap = new_cap;
			}

			aux_table[aux_table_count].pos_list[0] = pos;

			size_t k = 1;
			bool too_many_positions = false;
			do {
				if (k < AUX_TABLE_COLS) {
					aux_table[aux_table_count].pos_list[k++] = kmers[i].pos;
				} else {
					too_many_positions = true;
				}

				++ambig_total_kmers;
				++i;
			} while (i < kmers_len && kmer == kmers[i].kmer);

			if (too_many_positions) {  // these need not be included, but I'm including it anyway...
				serialize_uint32(out, POS_AMBIGUOUS);
				serialize_uint8(out, FLAG_AMBIGUOUS);
			} else {
				/* fill in remainder with 0s */
				while (k < AUX_TABLE_COLS) {
					aux_table[aux_table_count].pos_list[k++] = 0;
				}

				serialize_uint32(out, aux_table_count);
				serialize_uint8(out, FLAG_AMBIGUOUS);
				++aux_table_count;
			}
		} else {
			++unambig_kmers;
			serialize_uint32(out, pos);
			serialize_uint8(out, FLAG_UNAMBIGUOUS);
		}

		++kmers_written;
	}

	for (i = 0; i < aux_table_count; i++) {
		for (size_t j = 0; j < AUX_TABLE_COLS; j++) {
			serialize_uint32(out, aux_table[i].pos_list[j]);
		}
	}
	free(aux_table);

	rewind(out);
	serialize_uint64(out, kmers_written);
	serialize_uint64(out, aux_table_count);

	printf("Ref Dictionary\n");
	printf("Total k-mers:        %lu\n", total_kmers);
	printf("Unambig k-mers:      %lu\n", unambig_kmers);
	printf("Ambig unique k-mers: %lu\n", ambig_unique_kmers);
	printf("Ambig total k-mers:  %lu\n", ambig_total_kmers);
}

static void write_snp_kmers(struct snp_kmer_info *kmers, const size_t kmers_len, FILE *out)
{
	struct snp_aux_table_dictgen *aux_table = (struct snp_aux_table_dictgen *)malloc(SNP_AUX_TABLE_INIT_SIZE * sizeof(*aux_table));
	assert(aux_table);
	uint64_t aux_table_count = 0;
	uint64_t aux_table_cap = SNP_AUX_TABLE_INIT_SIZE;

	const uint64_t placeholder = 0UL;
	serialize_uint64(out, placeholder);  /* dict size (rows) placeholder */
	serialize_uint64(out, placeholder);  /* aux table size (rows) placeholder */

	uint64_t kmers_written = 0UL;
	size_t i = 0;

	/* keep track of a few statistics */
	const size_t total_kmers = kmers_len;
	size_t unambig_kmers = 0;
	size_t ambig_unique_kmers = 0;
	size_t ambig_total_kmers = 0;

	while (i < kmers_len) {
		const kmer_t kmer = kmers[i].kmer;
		const uint32_t pos = kmers[i].pos;
		const snp_info snp = kmers[i].snp;
		const uint8_t ref_freq = kmers[i].ref_freq;
		const uint8_t alt_freq = kmers[i].alt_freq;
		serialize_uint64(out, kmer);

		++i;

		if (i < kmers_len && kmer == kmers[i].kmer) {
			++ambig_unique_kmers;
			++ambig_total_kmers;

			if (aux_table_count == aux_table_cap) {
				const size_t new_cap = (aux_table_cap * 3)/2 + 1;
				aux_table = (struct snp_aux_table_dictgen*)realloc(aux_table, new_cap * sizeof(*aux_table));
				assert(aux_table);
				aux_table_cap = new_cap;
			}

			aux_table[aux_table_count].kmer = kmer;
			aux_table[aux_table_count].pos_list[0] = pos;
			aux_table[aux_table_count].snp_list[0] = snp;
			aux_table[aux_table_count].ref_freqs[0] = ref_freq;
			aux_table[aux_table_count].alt_freqs[0] = alt_freq;

			size_t k = 1;
			bool too_many_positions = false;
			do {
				if (k < AUX_TABLE_COLS) {
					aux_table[aux_table_count].pos_list[k] = kmers[i].pos;
					aux_table[aux_table_count].snp_list[k] = kmers[i].snp;
					aux_table[aux_table_count].ref_freqs[k] = kmers[i].ref_freq;
					aux_table[aux_table_count].alt_freqs[k] = kmers[i].alt_freq;
					++k;
				} else {
					too_many_positions = true;
				}

				++ambig_total_kmers;
				++i;
			} while (i < kmers_len && kmer == kmers[i].kmer);

			if (too_many_positions) {  // these need not be included, but I'm including it anyway...
				serialize_uint32(out, POS_AMBIGUOUS);
				serialize_uint8(out, 0);  // SNP info
				serialize_uint8(out, FLAG_AMBIGUOUS);
				serialize_uint8(out, 0);  // ref freq
				serialize_uint8(out, 0);  // alt freq
			} else {
				/* fill in remainder with 0s */
				while (k < AUX_TABLE_COLS) {
					aux_table[aux_table_count].pos_list[k] = 0;
					aux_table[aux_table_count].snp_list[k] = 0;
					aux_table[aux_table_count].ref_freqs[k] = 0;
					aux_table[aux_table_count].alt_freqs[k] = 0;
					++k;
				}

				serialize_uint32(out, aux_table_count);
				serialize_uint8(out, 0);  // SNP info
				serialize_uint8(out, FLAG_AMBIGUOUS);
				serialize_uint8(out, 0);  // ref freq
				serialize_uint8(out, 0);  // alt freq
				++aux_table_count;
			}
		} else {
			++unambig_kmers;
			serialize_uint32(out, pos);
			serialize_uint8(out, snp);
			serialize_uint8(out, FLAG_UNAMBIGUOUS);
			serialize_uint8(out, ref_freq);  // ref freq
			serialize_uint8(out, alt_freq);  // alt freq
		}

		++kmers_written;
	}

	for (i = 0; i < aux_table_count; i++) {
		serialize_uint64(out, aux_table[i].kmer);
		for (size_t j = 0; j < AUX_TABLE_COLS; j++) {
			serialize_uint32(out, aux_table[i].pos_list[j]);
			serialize_uint8(out, aux_table[i].snp_list[j]);
			serialize_uint8(out, aux_table[i].ref_freqs[j]);
			serialize_uint8(out, aux_table[i].alt_freqs[j]);
		}
	}
	free(aux_table);

	rewind(out);
	serialize_uint64(out, kmers_written);
	serialize_uint64(out, aux_table_count);

	printf("SNP Dictionary\n");
	printf("Total k-mers:        %lu\n", total_kmers);
	printf("Unambig k-mers:      %lu\n", unambig_kmers);
	printf("Ambig unique k-mers: %lu\n", ambig_unique_kmers);
	printf("Ambig total k-mers:  %lu\n", ambig_total_kmers);
}

void make_ref_dict(SeqVec ref, FILE *out)
{
	const size_t ref_len = ref.size;

	size_t total_kmers = 0;
	for (size_t i = 0; i < ref_len; i++) {
		total_kmers += ref.seqs[i].size - 32 + 1;
	}

	struct kmer_info *kmers = (struct kmer_info *)malloc(total_kmers * sizeof(*kmers));
	assert(kmers);
	size_t kmers_len = 0;
	uint32_t index = 1;

	for (size_t i = 0; i < ref_len; i++) {
		const char *seq = ref.seqs[i].seq;
		const size_t seq_len = ref.seqs[i].size;
		kmers_len += ref_to_constituent_kmers(&kmers[kmers_len], seq, seq_len, &index);
	}

	kmers = (struct kmer_info*)realloc(kmers, kmers_len * sizeof(*kmers));
	sort_kmers(kmers, kmers_len);
	write_kmers(kmers, kmers_len, out);
	free(kmers);
}

static Seq *find_seq_by_name(SeqVec ref, const char *name, unsigned int *start_index)
{
	unsigned int start_index_true = 1;

	for (size_t i = 0; i < ref.size; i++) {
		const char *seq_name = ref.seqs[i].name;

		if (strcmp(name, seq_name) == 0) {
			*start_index = start_index_true;
			return &ref.seqs[i];
		} else {
			start_index_true += ref.seqs[i].size;
		}
	}

	*start_index = 0;
	return NULL;
}

static char rev(const char c) {
	switch (c) {
	case 'A':
	case 'a':
		return 'T';
	case 'C':
	case 'c':
		return 'G';
	case 'G':
	case 'g':
		return 'C';
	case 'T':
	case 't':
		return 'A';
	default:
		return 'N';
	}
}

/*
 * `snp_file` in this case is a file in the UCSC SNP-txt format, consisting
 * of common SNP to be included in the dictionary.
 *
 * http://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=varRep&hgta_trac
 * k=snp141Common&hgta_table=snp141Common&hgta_doSchema=describe+table+schema
 *
 * Be sure to filter any SNPs with abnormal conditions (e.g. inconsistent alleles).
 */
void make_snp_dict(SeqVec ref, FILE *snp_file, FILE *out, bool **snp_locations, size_t *snp_locs_size)
{
#define CHROM_FIELD   1
#define INDEX_FIELD   2
#define STRAND_FIELD  6
#define REF1_FIELD    7
#define REF2_FIELD    8
#define ALT_FIELD     9
#define TYPE_FIELD    11
#define COUNT_FIELD   21
#define ALLELES_FIELD 22
#define FREQS_FIELD   24

	char line[6000];
	char chrom_name[50];
	char *line_split[30];
	Seq *chrom = NULL;

	size_t lines = 0;
	while (!feof(snp_file)) {
		if (fgetc(snp_file) == '\n')
			++lines;
	}
	rewind(snp_file);

	*snp_locs_size = 10;
	*snp_locations = (bool*)malloc(*snp_locs_size * sizeof(bool));
	assert(*snp_locations);
	memset(*snp_locations, false, *snp_locs_size);

	const size_t max_kmers_len = lines * 32;  /* 32 k-mers per line */
	struct snp_kmer_info *kmers = (struct snp_kmer_info*)malloc(max_kmers_len * sizeof(*kmers));
	assert(kmers);
	size_t kmers_len = 0;

	unsigned int start_index = 1;  // 1-based

	while (fgets(line, sizeof(line), snp_file)) {
		assert(!ferror(snp_file));

		if (line[0] == '#' || line[0] == '\n')
			continue;

		split_line(line, line_split);

		/* copy chromosome name into an independent buffer */
		size_t chrom_index;
		for (chrom_index = 0;
		     !isspace(line_split[CHROM_FIELD][chrom_index]) &&
		       (chrom_index < (sizeof(chrom_name) - 1));
		     chrom_index++) {

			chrom_name[chrom_index] = line_split[CHROM_FIELD][chrom_index];
		}
		chrom_name[chrom_index] = '\0';

		const char ref_base = toupper(line_split[REF1_FIELD][0]);
		const unsigned ref_base_u = encode_base(ref_base);

		if (ref_base_u == BASE_X ||
			strncmp(line_split[TYPE_FIELD], "single", strlen("single")) != 0 ||
		    ref_base != toupper(line_split[REF2_FIELD][0])) {

			continue;
		}

		/* check if reference sequences are 1 base long */
		if (!(isspace(line_split[REF1_FIELD][1]) && isspace(line_split[REF2_FIELD][1]))) {
			continue;
		}

		if (chrom == NULL || strcmp(chrom->name, chrom_name) != 0) {
			chrom = find_seq_by_name(ref, chrom_name, &start_index);

			if (chrom == NULL) {
				continue;
			}
		}

		const unsigned int index = atoi(line_split[INDEX_FIELD]);  // 0-based

		if (index >= chrom->size || toupper(chrom->seq[index]) != ref_base) {
			fprintf(stderr,
			        "Mismatch found between reference sequence and SNP file at 0-based index %u in %s.\n",
			        index,
			        chrom->name);
			exit(EXIT_FAILURE);
		}

		if (index < 32 || (index + 32) > chrom->size) {
			continue;
		}

		/* we should only process bi-allelic SNPs */
		if (*line_split[COUNT_FIELD] != '2') {
			continue;
		}

		const bool neg = (line_split[STRAND_FIELD][0] == '-');
		if (!neg)
			assert(line_split[STRAND_FIELD][0] == '+');

		const char a1 = neg ? rev(toupper(line_split[ALLELES_FIELD][0])) : toupper(line_split[ALLELES_FIELD][0]);
		const char a2 = neg ? rev(toupper(line_split[ALLELES_FIELD][2])) : toupper(line_split[ALLELES_FIELD][2]);

		assert((a1 == 'A' || a1 == 'C' || a1 == 'G' || a1 == 'T') &&
		       (a2 == 'A' || a2 == 'C' || a2 == 'G' || a2 == 'T'));

		if (a1 != ref_base && a2 != ref_base) {
			continue;
		}

		const size_t loc = start_index + index;
		if (loc >= *snp_locs_size) {
			*snp_locations = (bool*)realloc(*snp_locations, (loc + 1) * sizeof(bool));
			assert(*snp_locations);
			memset(*snp_locations + *snp_locs_size, false, loc - *snp_locs_size + 1);
			*snp_locs_size = loc + 1;
		}
		(*snp_locations)[loc] = true;

		char *p = line_split[FREQS_FIELD];
		float freq1 = atof(p);
		while (*p++ != ',');
		float freq2 = atof(p);

		if (a2 == ref_base) {
			const float tmp = freq1;
			freq1 = freq2;
			freq2 = tmp;
		}

		const uint8_t freq1_enc = (uint8_t)(freq1*0xff);
		const uint8_t freq2_enc = (uint8_t)(freq2*0xff);

		for (char *alt_p = line_split[ALT_FIELD]; !isspace(*alt_p); alt_p++) {
			struct snp_kmer_info snp_kmers[32];

			const char alt = neg ? rev(toupper(*alt_p)) : toupper(*alt_p);

			if (alt == ref_base || !(alt == 'A' || alt == 'C' || alt == 'G' || alt == 'T')) {
				continue;
			}

			assert(kmers_len + 32 <= max_kmers_len);

			const char *seq = chrom->seq;
			bool kmer_had_n;
			kmer_t kmer = encode_kmer(&seq[index - 32], &kmer_had_n);

			if (kmer_had_n)
				goto end;

			for (unsigned int i = 0; i < 32; i++) {
				const char next_base = (i ? seq[index + i] : alt);

				if (next_base == 'N' || next_base == 'n')
					goto end;

				kmer = shift_kmer(kmer, next_base);
				snp_kmers[i].kmer = kmer;
				snp_kmers[i].pos = start_index + index - 32 + 1 + i;
				snp_kmers[i].snp = SNP_INFO_MAKE(32 - 1 - i, ref_base_u);
				snp_kmers[i].ref_freq = freq1_enc;
				snp_kmers[i].alt_freq = freq2_enc;
			}

			memcpy(&kmers[kmers_len], snp_kmers, 32 * sizeof(*kmers));
			kmers_len += 32;

			end:
			break;
		}
	}

	kmers = (struct snp_kmer_info*)realloc(kmers, kmers_len * sizeof(*kmers));
	sort_snp_kmers(kmers, kmers_len);
	write_snp_kmers(kmers, kmers_len, out);
	free(kmers);

#undef CHROM_FIELD
#undef INDEX_FIELD
#undef STRAND_FIELD
#undef REF1_FIELD
#undef REF2_FIELD
#undef ALT_FIELD
#undef TYPE_FIELD
#undef COUNT_FIELD
#undef ALLELES_FIELD
#undef FREQS_FIELD
}


#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <ctype.h>
#include "fasta_parser.h"

static char norm(const char c)
{
	switch (c) {
	case 'A':
	case 'a':
		return 'A';
	case 'C':
	case 'c':
		return 'C';
	case 'G':
	case 'g':
		return 'G';
	case 'T':
	case 't':
		return 'T';
	default:
		return 'N';
	}
}

#define GENOME_VECTOR_INITIAL_SIZE    10
#define MAX_GENOME_NAME_LENGTH        64

/*
 * We perform the parsing with two passes of the file. The
 * first pass determines the size of each sequence, and the
 * second pass actually reads the sequences into memory.
 */
SeqVec parse_fasta(const char *filename)
{
	size_t seqs_capacity = GENOME_VECTOR_INITIAL_SIZE;
	Seq *seqs = (Seq*)malloc(seqs_capacity * sizeof(Seq));
	size_t seqs_size = 0;

	FILE *file = fopen(filename, "r");

	if (file == NULL) {
		fprintf(stderr,
		        "Error: Could not open file '%s' for writing.\n",
		        filename);

		exit(EXIT_FAILURE);
	}

	char c;
	while ((c = fgetc(file)) != EOF) {
		if (c == '>') {
			if (seqs_size == seqs_capacity) {
				seqs_capacity = (seqs_capacity * 3)/2 + 1;
				seqs = (Seq*)realloc(seqs, seqs_capacity * sizeof(Seq));
			}

			seqs[seqs_size].name = (char*)malloc(MAX_GENOME_NAME_LENGTH + 1);
			seqs[seqs_size].size = 0;

			int name_idx = 0;
			bool newline = false;
			while ((c = fgetc(file)) != EOF) {
				if (c == '|' ||
				    isspace(c) ||
				    name_idx == MAX_GENOME_NAME_LENGTH) {

					newline = (c == '\n');
					break;
				}

				seqs[seqs_size].name[name_idx++] = c;
			}
			seqs[seqs_size].name[name_idx] = '\0';

			if (!newline) {
				while (fgetc(file) != '\n');
			}

			while ((c = fgetc(file)) != EOF) {
				if (c == '>') {
					ungetc('>', file);
					break;
				}

				if (c == '\n') {
					continue;
				}

				++seqs[seqs_size].size;
			}

			seqs[seqs_size].seq = (char*)malloc(seqs[seqs_size].size + 1);
			++seqs_size;
		}
	}

	rewind(file);

	unsigned int idx = 0;
	while ((c = fgetc(file)) != EOF) {
		if (c == '>') {
			while (fgetc(file) != '\n');

			const size_t size = seqs[idx].size;

			size_t i = 0;
			while (i < size) {
				c = fgetc(file);

				if (c == '\n') {
					continue;
				}

				if (c == EOF || c == '>') {
					exit(EXIT_FAILURE);
				}

				seqs[idx].seq[i++] = norm(c);
			}

			seqs[idx].seq[size] = '\0';
			++idx;
		}
	}

	fclose(file);

	return (SeqVec){.seqs = seqs, .size = seqs_size};
}

void seq_dealloc(const Seq *seq)
{
	free(seq->seq);
	free(seq->name);
}

void seqvec_dealloc(const SeqVec *seqvec)
{
	const size_t size = seqvec->size;

	for (size_t i = 0; i < size; i++) {
		seq_dealloc(&seqvec->seqs[i]);
	}

	free(seqvec->seqs);
}

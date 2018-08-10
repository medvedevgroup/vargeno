#ifndef DICTGEN_H
#define DICTGEN_H

#include <stdio.h>
#include "fasta_parser.h"

void make_ref_dict(SeqVec ref, FILE *out);

void make_snp_dict(SeqVec ref, FILE *snp_file, FILE *out, bool **snp_locations, size_t *snp_locs_size);

#endif /* DICTGEN_H */


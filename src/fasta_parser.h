#ifndef FASTA_PARSER_H
#define FASTA_PARSER_H

typedef struct {
	char *name;
	char *seq;
	size_t size;
} Seq;

typedef struct {
	Seq *seqs;
	size_t size;
} SeqVec;

SeqVec parse_fasta(const char *filename);
void seq_dealloc(const Seq *seq);
void seqvec_dealloc(const SeqVec *seqvec);

#endif /* FASTA_PARSER_H */

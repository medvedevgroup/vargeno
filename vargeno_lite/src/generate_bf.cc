#include "generate_bf.h"
#include <cstdlib>

BFGenerator::BFGenerator() {

}

BFGenerator::~BFGenerator() {

}

/**
--------------
readFasta
	read genome sequence in FASTA file
--------------
*/
int BFGenerator::readFasta(string & fasta_filename) {

	std::ifstream input(fasta_filename);
	if (!input.good()) {
		std::cerr << "Error opening: " << fasta_filename << " . You have failed." << std::endl;
		return -1;
	}

	std::string line, id, DNA_sequence;
	// Don't loop on good(), it doesn't allow for EOF!!
	//    while (std::getline(input, line).good()) {
	while (std::getline(input, line)) {

		// line may be empty so you *must* ignore blank lines
		// or you have a crash waiting to happen with line[0]
		if (line.empty())
			continue;

		if (line[0] == '>') {
			// output previous line before overwriting id
			// but ONLY if id actually contains something
			if (!id.empty()) {
				GenomeSeq genome;
				genome.name = id;
				genome.seq = DNA_sequence;
				genome.size = DNA_sequence.size();
				genome_vector.emplace_back(genome);
			}

			id = line.substr(1);
			DNA_sequence.clear();
		}
		else {//  if (line[0] != '>'){ // not needed because implicit
			DNA_sequence += line;
		}
	}

	// output final entry
	// but ONLY if id actually contains something
	if (!id.empty()) {
		GenomeSeq genome;
		genome.name = id;
		genome.seq = DNA_sequence;
		genome.size = DNA_sequence.size();
		genome_vector.emplace_back(genome);
	}

#if DEBUG
	for (auto it = genome_vector.begin(); it != genome_vector.end(); ++it) {
		cout << it->name << ": " << it->seq.size() << endl;
	}
#endif // DEBUG

	input.close();
	return 0;
}

/**
------------
constructBfFromGenomeseq
	construct bloom filter 
	extract kmers from genome sequence
	encode kmers that do not contains 'N'
------------
Unit Tested on 3/24/2017 by Chen Sun
------------
Modified on 4/27/2017 by Chen Sun
new algorithm
------------
Code Reviewed on 4/28/2017 by Chen Sun
Unit Tested on 4/28/2017 by Chen Sun
*/
int BFGenerator::constructBfFromGenomeseq(string bf_filename, bool is_canonical = false) {

	/*set bloom filter size*/
	int64_t unique_kmers = REF_BF_BYTES;
	int64_t bf_size = unique_kmers * 8; // 8x the number of unique kmers, the resulting bloom filter is expected to be 3G

	/*initialize bloom filter*/
	BloomFilter* bf = new BloomFilter(bf_size, REF_BF_RANGE);
	/*iterate through reference sequences*/
	//unordered_set<uint32_t> lo32_set;
    //unordered_set<uint64_t> lo40_set;

    //string huge_bf_filename = bf_filename + ".huge.lite.bf";
    //int64_t huge_kmer_size = HUGE_LITE_BF_BYTES;
    //int64_t huge_bf_size = huge_kmer_size * 8;
    //BloomFilter* huge_bf = new BloomFilter(huge_bf_size, SNP_BF_RANGE);

	for (auto it = genome_vector.begin(); it != genome_vector.end(); ++it) {
		int ref_len = it->seq.size();
		assert(ref_len >= KMER_SIZE);
		const size_t kmers_len_max = ref_len - KMER_SIZE + 1;
		
		bool need_full_encode = true;
		bool kmer_had_n;
		kmer_t kmer;
		/*insert kmers into bloom filter*/
		for (size_t i = 0; i < kmers_len_max; i++) {
			const char next_base = it->seq[i + KMER_SIZE - 1];
			
			/*need to generate a new 32mer from pos i
				if the 32mer contains N
					then need_full_encode, set i to the next pos after N, and continue
				else
					set need_full_encode = false, encode kmer, add to bloom filter
			*/
			if (need_full_encode) {
				string kmer_string = it->seq.substr(i, KMER_SIZE);
				kmer = encode_kmer(kmer_string.c_str(), &kmer_had_n);
				need_full_encode = kmer_had_n;
			}
			/*else if need_full_code == false && next_base != 'N'*/
			else if (next_base != 'N' && next_base != 'n') {
				kmer = shift_kmer(kmer, next_base);
				kmer_had_n = false;
			}
			/*else if need_full_code == false && next_base == 'N', i should goto pos(N)*/
			else {
				need_full_encode = true;
				kmer_had_n = true;
			}

			if (!kmer_had_n) {
				/*
					get hi 32 of kmer
					insert into bloom filter
				*/
				uint32_t kmer_lo = LO(kmer);
				bf->set_value(kmer_lo);
                //kmer_t kmer_lo40 = LO40(kmer);
                //huge_bf->set_value(kmer_lo40);
				//lo32_set.insert(kmer_lo);
                //lo40_set.insert(LO40(kmer));
			}
		}
	}

	std::cout << "[BloomFilter constructBfFromGenomeseq] bit vector: " << bf->count_ones() << "/" << bf->num_filter_bits() << std::endl;
	//std::cout << "[BloomFilter constructBfFromGenomeseq] lite bit vector: " << huge_bf->count_ones() << "/" << huge_bf->num_filter_bits() << std::endl;
	//std::cout << "[BloomFilter constructBfFromGenomeseq] total ref lo_32: " << lo32_set.size() << std::endl;
	//std::cout << "[BloomFilter constructBfFromGenomeseq] total ref lo_40: " << lo40_set.size() << std::endl;
	//lo32_set.clear();
    //lo40_set.clear();
	bf->save(bf_filename);
    //huge_bf->save(huge_bf_filename);
	delete bf;
    //delete huge_bf;
	
	return 0;
}

/**
--------
updateBfFromVcf
	for each snp, create 32 kmers that cover this snp
	insert kmers into bloom filter
--------
Code Reviewed by Chen Sun 3/24/2017
--------
*/
int BFGenerator::constructBfFromVcf(const string & vcf_filename, string bf_filename, bool is_canonical = false) {
	string pre_chr_name = "XO";
	string seq;

#if DEBUG
	assert(bf != NULL);
#endif

	std::ifstream input(vcf_filename);
	if (!input.good()) {
		std::cerr << "Error opening: " << vcf_filename << " . You have failed." << std::endl;
		return -1;
	}
	string line;

	/*set bloom filter size*/
	int64_t unique_kmers = SNP_BF_BYTES;
	int64_t bf_size = unique_kmers * 8; // 4x the number of unique kmers, the resulting bloom filter is expected to be 3G

	/*initialize bloom filter*/
	BloomFilter* bf = new BloomFilter(bf_size, SNP_BF_RANGE);

	//string lite_bf_filename = bf_filename + ".ref.lite.bf";
    //int64_t lite_kmer_size = REF_LITE_BF_BYTES;
    //int64_t lite_bf_size = lite_kmer_size * 8;
    //BloomFilter* lite_bf = new BloomFilter(lite_bf_size, SNP_BF_RANGE);

	while (std::getline(input, line)) {
		if (line.empty()) continue;
		if (line[0] == '#') continue;
		vector<string> columns = split(line, '\t');
		string chr_name = columns[0];
		int pos = stoi(columns[1]) - 1; // 1-based to 0-based coordinate
		string ref_seq = columns[3];
		string alt_seq = columns[4];
		/*filter out SNPs*/
		if (ref_seq.size() > 1 || alt_seq.size() > 1) continue;

		/*check chr_name and find corresponding seq*/
		if (chr_name != pre_chr_name) {
			for (auto it = genome_vector.begin(); it != genome_vector.end(); ++it) {
				if (it->name == chr_name) {
					seq = it->seq;
					break;
				}
			}
			pre_chr_name = chr_name;
		}

		if (pos < KMER_SIZE || (pos + KMER_SIZE) > seq.size()) continue;

		char ref_nt = ref_seq[0];
		char alt_nt = alt_seq[0];

		/*check if ref_nt is right*/
		if (ref_nt != seq[pos] || ref_nt == alt_nt) continue;

		/*get the kmer before pos*/
		string kmer_string = seq.substr(pos - KMER_SIZE, KMER_SIZE);

		/*if the kmer before pos has N, then skip current SNP*/
		bool has_n = false;

		kmer_t kmer = encode_kmer(kmer_string.c_str(), &has_n);
		//kmer_t ref_kmer = encode_kmer(kmer_string.c_str(), &has_n);
		if (has_n) continue;
		/*generate all kmers cover current SNPs, and insert into bloom filter*/
		/*
		there is potential problem of this for loop:
			if a SNP is in the area with N inside, then some kmers might be added into bloom filter
			this will only slightly increase false positives, but that is fine, no big deal
		*/
		//vector<uint64_t> kmer_lo_list;
		for (unsigned int i = 0; i < KMER_SIZE; i++) {
			/* the first next base is alternative allele, the rest is reference nt */
			const char next_base = (i ? seq[pos + i] : alt_nt); // if i ==0, change it to alternative allele, if i > 0, then it is the original nt

			//const char ref_next_base = seq[pos+i];

			if (next_base == 'N' || next_base == 'n') {
				has_n = true;
				break;
			}

			/*shift kmer to contain the next base*/
			shift_kmer(kmer, next_base);
			//shift_kmer(ref_kmer, ref_next_base);

			uint64_t kmer_lo = LO40(kmer);
			bf->set_value(kmer_lo);

			//uint64_t ref_kmer_lo = LO40(ref_kmer);
			//lite_bf->set_value(ref_kmer_lo);
            //kmer_lo_list.emplace_back(kmer_lo);
		}
		if (has_n) continue;

		// add all kmer_hi into bloom filter
		//for (auto kmer_lo : kmer_lo_list) {
		//	bf->set_value(kmer_lo);
		//}
	}

	input.close();
	std::cout << "[BloomFilter constructBfFromGenomeseq] bit vector: " << bf->count_ones() << "/" << bf->num_filter_bits() << std::endl;
	//std::cout << "[BloomFilter constructBfFromGenomeseq] lite bit vector: " << lite_bf->count_ones() << "/" << lite_bf->num_filter_bits() << std::endl;
	bf->save(bf_filename);
	//lite_bf->save(lite_bf_filename);
	delete bf;
	//delete lite_bf;

	return 0;
}


// rewrite make_snp_dict in dectgen.c
/*
int BFGenerator::constructBfFromUcsc(const string & text_filename, string bf_filename, bool is_canonical = false) {
	
	const int CHROM_FIELD = 1;
	const int INDEX_FIELD = 2;
	const int STRAND_FIELD = 6;
	const int REF1_FIELD = 7;
	const int REF2_FIELD = 8;
	const int ALT_FIELD = 9;
	const int TYPE_FIELD = 11;
	const int COUNT_FIELD = 21;
	const int ALLELES_FIELD = 22;
	const int FREQS_FIELD = 24;

	string pre_chr_name = "XO";
	string seq;

#if DEBUG
	assert(bf != NULL);
#endif

	std::ifstream input(text_filename);
	if (!input.good()) {
		std::cerr << "Error opening: " << text_filename << " . You have failed." << std::endl;
		return -1;
	}
	string line;

	//unordered_set<uint64_t> kmer_set;


	int64_t unique_kmers = SNP_BF_BYTES;
	int64_t bf_size = unique_kmers * 8; // 4x the number of unique kmers, the resulting bloom filter is expected to be 3G


	BloomFilter* bf = new BloomFilter(bf_size, SNP_BF_RANGE);

	while (std::getline(input, line)) {
		if (line.empty()) continue;
		if (line[0] == '#') continue;
		vector<string> columns = split(line, '\t');
		string chr_name = columns[CHROM_FIELD];
		uint32_t pos = stoi(columns[INDEX_FIELD]); // 0-based coordinate
		
		const char ref_base = toupper(columns[REF1_FIELD][0]);
		const unsigned ref_base_u = encode_base(ref_base);
		if (ref_base_u == BASE_X || columns[TYPE_FIELD] != "single" || ref_base != toupper(columns[REF2_FIELD][0])) continue;

		string ref_seq = columns[REF2_FIELD];
		auto alleles = split(columns[ALT_FIELD], '/');
		if(alleles.size() != 2) continue;
		

		if (ref_seq.size() > 1) continue;


		if (chr_name != pre_chr_name) {
			bool find_chromosome = false;
			for (auto it = genome_vector.begin(); it != genome_vector.end(); ++it) {
				if (it->name == chr_name) {
					seq = it->seq;
					find_chromosome = true;
					break;
				}
			}
			if (!find_chromosome) continue;
			pre_chr_name = chr_name;
		}

		if (pos < KMER_SIZE || (pos + KMER_SIZE) > seq.size()) continue;

		if (pos >= seq.length() || toupper(seq[pos]) != ref_base) continue;

		if (pos < 32 || (pos + 32) > seq.length()) continue;

		const bool neg = (columns[STRAND_FIELD][0] == '-');
		if (!neg)
			assert(columns[STRAND_FIELD][0] == '+');

		const char a1 = neg ? rev(toupper(columns[ALLELES_FIELD][0])) : toupper(columns[ALLELES_FIELD][0]);
		const char a2 = neg ? rev(toupper(columns[ALLELES_FIELD][2])) : toupper(columns[ALLELES_FIELD][2]);

		//assert((a1 == 'A' || a1 == 'C' || a1 == 'G' || a1 == 'T') &&
		//	(a2 == 'A' || a2 == 'C' || a2 == 'G' || a2 == 'T'));
        
        if(a1 != 'A' && a1 != 'C' && a1 != 'G' && a1 != 'T'){
            std::cerr << "[BFGenerator] unrecognized allele: " << a1 << std::endl;
        }
        
        if(a2 != 'A' && a2 != 'C' && a2 != 'G' && a2 != 'T'){
            std::cerr << "[BFGenerator] unrecognized allele: " << a2 << std::endl;
        }

		if (a1 != ref_base && a2 != ref_base) {
			continue;
		}

		vector<string> p = split(columns[FREQS_FIELD], ',');
		float freq1 = stof(p[0]);
		float freq2 = stof(p[1]);

		char alt_base = a2;

		if (a2 == ref_base) {
			const float tmp = freq1;
			freq1 = freq2;
			freq2 = tmp;
			alt_base = a1;
		}


		string kmer_string = seq.substr(pos - KMER_SIZE, KMER_SIZE);


		bool has_n = false;

		kmer_t kmer = encode_kmer(kmer_string.c_str(), &has_n);
		if (has_n) continue;


		//vector<uint64_t> kmer_lo_list;
		for (unsigned int i = 0; i < KMER_SIZE; i++) {
			
			const char next_base = (i ? seq[pos + i] : alt_base); // if i ==0, change it to alternative allele, if i > 0, then it is the original nt

			if (next_base == 'N' || next_base == 'n') {
				has_n = true;
				break;
			}

			
			shift_kmer(kmer, next_base);

			uint64_t kmer_lo = LO40(kmer);
			bf->set_value(kmer_lo);
            //kmer_set.insert(kmer_lo);
            //kmer_lo_list.emplace_back(kmer_lo);
		}
		if (has_n) continue;

		// add all kmer_hi into bloom filter
		//for (auto kmer_lo : kmer_lo_list) {
		//	bf->set_value(kmer_lo);
			//kmer_set.insert(kmer_lo);
		//}
	}

	input.close();
	std::cout << "[BloomFilter constructBfFromUCSC] bit vector: " << bf->count_ones() << "/" << bf->num_filter_bits() << std::endl;
	//std::cout << "[BloomFilter constructBfFromUcsc] unique snp lo40: " << kmer_set.size() << std::endl;
	//kmer_set.clear();
	bf->save(bf_filename);
	delete bf;

	return 0;
}
*/

int BFGenerator::constructBfFromUcsc(const string & text_filename, string bf_filename, bool is_canonical = false) {
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

	string pre_chr_name = "XO";
	string seq_string;

	std::ifstream input(text_filename);
	if (!input.good()) {
		std::cerr << "Error opening: " << text_filename << " . You have failed." << std::endl;
		return -1;
	}
	/*set bloom filter size*/
	int64_t unique_kmers = SNP_BF_BYTES;
	int64_t bf_size = unique_kmers * 8; // 4x the number of unique kmers, the resulting bloom filter is expected to be 3G

	/*initialize bloom filter*/
	BloomFilter* bf = new BloomFilter(bf_size, SNP_BF_RANGE);
	//string lite_bf_filename = bf_filename + ".ref.lite.bf";
    //int64_t lite_kmer_size = REF_LITE_BF_BYTES;
    //int64_t lite_bf_size = lite_kmer_size * 8;
    //BloomFilter* lite_bf = new BloomFilter(lite_bf_size, SNP_BF_RANGE);

	size_t kmers_len = 0;

	string line_string;

	while (std::getline(input, line_string)) {

		if (line_string[0] == '#' || line_string[0] == '\n')
			continue;

		vector<string> columns = split(line_string, '\t');
		string chr_name = columns[CHROM_FIELD];

		const char ref_base = toupper(columns[REF1_FIELD][0]);
		const unsigned ref_base_u = encode_base(ref_base);

		if (ref_base_u == BASE_X ||
			columns[TYPE_FIELD] != "single" ||
		    ref_base != toupper(columns[REF2_FIELD][0])) {

			continue;
		}

		/* check if reference sequences are 1 base long */
		if (columns[REF1_FIELD].size() != 1 || columns[REF2_FIELD].size() != 1) {
			continue;
		}

		if (chr_name != pre_chr_name) {
			bool find_chromosome = false;
			for (auto it = genome_vector.begin(); it != genome_vector.end(); ++it) {
				if (it->name == chr_name) {
					seq_string = it->seq;
					find_chromosome = true;
					break;
				}
			}
			if (!find_chromosome) continue;
			pre_chr_name = chr_name;
		}

		const unsigned int index = stoi(columns[INDEX_FIELD]);  // 0-based

		if (index >= seq_string.size() || toupper(seq_string[index]) != ref_base) {
			cout << "Mismatch found between reference sequence and SNP file at 0-based index %u in %s.\n";
			exit(EXIT_FAILURE);
		}

		if (index < 32 || (index + 32) > seq_string.size()) {
			continue;
		}

		/* we should only process bi-allelic SNPs */
		if (columns[COUNT_FIELD] != "2") {
			continue;
		}

		const bool neg = (columns[STRAND_FIELD][0] == '-');
		if (!neg)
			assert(columns[STRAND_FIELD][0] == '+');

		const char a1 = neg ? rev(toupper(columns[ALLELES_FIELD][0])) : toupper(columns[ALLELES_FIELD][0]);
		const char a2 = neg ? rev(toupper(columns[ALLELES_FIELD][2])) : toupper(columns[ALLELES_FIELD][2]);

		assert((a1 == 'A' || a1 == 'C' || a1 == 'G' || a1 == 'T') &&
		       (a2 == 'A' || a2 == 'C' || a2 == 'G' || a2 == 'T'));

		if (a1 != ref_base && a2 != ref_base) {
			continue;
		}

		char *alt_p = const_cast<char*> (columns[ALT_FIELD].c_str());

		for (; !isspace(*alt_p); alt_p++) {

			const char alt = neg ? rev(toupper(*alt_p)) : toupper(*alt_p);

			if (alt == ref_base || !(alt == 'A' || alt == 'C' || alt == 'G' || alt == 'T')) {
				continue;
			}

			//assert(kmers_len + 32 <= max_kmers_len);

			const char *seq = seq_string.c_str();
			bool kmer_had_n;
			kmer_t kmer = encode_kmer(&seq[index - 32], &kmer_had_n);
			uint64_t kmer_lo = LO40(kmer);
			bf->set_value(kmer_lo);

			//kmer_t ref_kmer = encode_kmer(&seq[index - 32], &kmer_had_n);
			//uint64_t ref_kmer_lo = LO40(ref_kmer);
			//lite_bf->set_value(ref_kmer_lo);
			
            if (kmer_had_n)
				goto end;

			for (unsigned int i = 0; i < 32; i++) {
				const char next_base = (i ? seq[index + i] : alt);
				const char ref_next_base = seq[index + i];

				if (next_base == 'N' || next_base == 'n')
					goto end;

				kmer = shift_kmer(kmer, next_base);
				uint64_t kmer_lo = LO40(kmer);
				bf->set_value(kmer_lo);


				//ref_kmer = shift_kmer(ref_kmer, ref_next_base);
				//uint64_t ref_kmer_lo = LO40(ref_kmer);
				//lite_bf->set_value(ref_kmer_lo);
			}

			end:
			break;
		}
	}

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

	input.close();
	std::cout << "[BloomFilter constructBfFromUCSC] bit vector: " << bf->count_ones() << "/" << bf->num_filter_bits() << std::endl;
	bf->save(bf_filename);
	delete bf;

	//std::cout << "[BloomFilter constructBfFromUCSC] ref lite bit vector: " << lite_bf->count_ones() << "/" << lite_bf->num_filter_bits() << std::endl;
	//lite_bf->save(lite_bf_filename);
    //delete lite_bf;

	return 0;
}


char BFGenerator::rev(const char c) {
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


int BFGenerator::constructBfFromEncode(const string & encode_filename, string bf_filename, bool is_canonical = false) {
    
	std::ifstream input(encode_filename);
	if (!input.good()) {
		std::cerr << "Error opening: " << encode_filename << " . You have failed." << std::endl;
		return -1;
	}
	string line;

	//unordered_set<uint64_t> kmer_set;

	/*set bloom filter size*/
	int64_t unique_kmers = SNP_BF_BYTES;
	int64_t bf_size = unique_kmers * 8; // 4x the number of unique kmers, the resulting bloom filter is expected to be 3G

	/*initialize bloom filter*/
	BloomFilter* bf = new BloomFilter(bf_size, SNP_BF_RANGE);

	while (std::getline(input, line)) {
		if (line.empty()) continue;
        //if (line[0] == '#') continue;
        vector<string> columns = split(line, ' ');
		string num_string = columns[0];
        uint64_t encode = std::strtoull(num_string.c_str(),NULL,0);
        bf->set_value(encode);
        //kmer_set.insert(encode);
    }

    input.close();
    
	std::cout << "[BloomFilter constructBfFromEncode] bit vector: " << bf->count_ones() << "/" << bf->num_filter_bits() << std::endl;
	//std::cout << "unique snp kmer: " << kmer_set.size() << std::endl;
	//kmer_set.clear();
	bf->save(bf_filename);
	delete bf;
    
    return 0;
}



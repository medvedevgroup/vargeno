#include "generate_bf.h"

void print_help() {
	cout << "Usage: " << endl;
	cout << "\tvcf <input FASTA> <input VCF> <output ref BF> <output snp BF>" << endl;
	cout << "\tucsc <input FASTA> <input UCSC file> <output ref BF> <output snp BF>" << endl;
	cout << "\tencode <input FASTA> <input Encode file> <output ref BF> <output snp BF>" << endl;
    cout << "\tsnp <input FASTA> <input Encode file> <output snp BF>" << endl;
}

int main(const int argc, const char *argv[]) {

	// unit test
    //
	if (argc < 5) {
		print_help();
		exit(EXIT_FAILURE);
	}

	string snp_format = argv[1];
	if (snp_format != "vcf" && snp_format != "ucsc" && snp_format != "encode" && snp_format != "snp") {
		print_help();
		exit(EXIT_FAILURE);
	}
	
	string ref_filename = argv[2];
	string snp_filename = argv[3];
	string refbf_filename;
	string snpbf_filename;

    if(snp_format == "snp"){
        
	    if (argc != 5) {
		    print_help();
		    exit(EXIT_FAILURE);
	    }
        
        snpbf_filename = argv[4];
    }else{
	    
        if (argc < 6) {
		    print_help();
		    exit(EXIT_FAILURE);
	    }
        
        refbf_filename = argv[4];
        snpbf_filename = argv[5];
    }

	BFGenerator * bg = new BFGenerator();

	bg->readFasta(ref_filename);

    if(snp_format != "snp"){
        bg->constructBfFromGenomeseq(refbf_filename, false);
    }

	if (snp_format == "vcf") {
		bg->constructBfFromVcf(snp_filename, snpbf_filename, false);
	}
	else if (snp_format == "ucsc") {
		bg->constructBfFromUcsc(snp_filename, snpbf_filename, false);
    }
    else if (snp_format == "encode" || snp_format == "snp"){
        bg->constructBfFromEncode(snp_filename, snpbf_filename, false);
    }

	delete bg;

	return 0;
}

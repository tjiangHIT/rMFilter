#include "Ctrl_option.h"
#include "info.h"

#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>

const char *short_options = "hl:";
struct option long_options[] = {
	//{ "threads",     1,   NULL,    't'   },
	//{ "num",     1,   NULL,    'n' },
	{"help",	0,	NULL,	'h'},
	{"seed_length",	1,	NULL,	'l'},
	//{ "hit_max",		1,NULL,	'm'},
	//{ "auto_load", 0, NULL, 'a'},
	{0,	0,	0,	0   }
};

Ctrl_option::Ctrl_option(Options *opt)
{
	opt->len_kmer = 13;
}

int Ctrl_option::Usage()
{
	fprintf(stderr, "\n"); 
	fprintf(stderr, "Program:   Index\n"); 
	fprintf(stderr, "Version:   %s\n", PACKAGE_VERSION); 
	fprintf(stderr, "Contact:   %s\n\n", CONTACT); 
	fprintf(stderr, "Usage:     Index [Options] <HashIndexDir> <Reference>\n\n"); 
	fprintf(stderr, "Options:   -h, --help                   help\n"); 
 	//fprintf(stderr, "           -t, --threads       <int>    thread\n"); 
	//fprintf(stderr, "           -n, --num           <int>    candidate number [5]\n"); 
	//fprintf(stderr, "           -m, --hit_max       <int>    max hit times of a seed [1000]\n"); 
	fprintf(stderr, "           -l, --seed_length   <int>    seed length of hash index [13]\n"); 
	//fprintf(stderr, "           -a, --auto_load              load hash table from hash file without produce hash file");
	//fprintf(stderr, "           -c, --write_cigar            print cigar in XA fields [False]\n"); 
	fprintf(stderr, "\n"); 
	return 0;
}

int Ctrl_option::opt_parse(int argc, char *argv[], Options* opt)
{
	int c;
	int option_index = 0;
	if( argc == 1 )
		return Usage();
	while((c = getopt_long(argc, argv, short_options, long_options, &option_index))>=0)
	{
		switch(c)
		{
			case 'h':
				return Usage();
				break;
			case 'l':
				opt->len_kmer = atoi(optarg);
				break;
			default:
				fprintf(stderr,"Not proper parameters\n");
				return Usage();
		}
	}
	if(optind + 2 != argc)
	{
		fprintf(stderr, "[opt_parse]: index directory, reference file can't be omited!\n");
		return 0;
	}
	strncpy(opt->hash_dir, argv[optind++],sizeof(opt->hash_dir));
	strncpy(opt->ref_path, argv[optind++],sizeof(opt->ref_path));
	return 1;
}

void Ctrl_option::show_parameters(Options* opt)
{
	fprintf(stderr, "\n"); 
	fprintf(stderr, ":::: Simulation parameters :::\n"); 
	fprintf(stderr, "\n"); 
	fprintf(stderr, "seed_length:   %u\n", opt->len_kmer); 
	fprintf(stderr, "HashIndexDir:  %s\n", opt->hash_dir); 
	fprintf(stderr, "Reference:     %s\n", opt->ref_path); 
	fprintf(stderr, "\n"); 
	fprintf(stderr, ":::: Simulation parameters :::\n"); 
	fprintf(stderr, "\n"); 
}
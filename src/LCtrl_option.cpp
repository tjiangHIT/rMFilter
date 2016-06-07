/**  
 * All rights Reserved, Designed By HIT-Bioinformatics   
 * @Title:  LCtrl_option.cpp
 * @Package 
 * @Description:    Control the options of rMFilter alignment
 * @author: tjiang
 * @date:   June 7 2016
 * @version V1.0     
 */

#include "LCtrl_option.h"
#include "info.h"

#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>

const char *short_options = "r:hl:t:";
struct option long_options[] = {
	{ "ratio",     1,   NULL,    'r' },
	{"help",	0,	NULL,	'h'},
	{"kmerSize",	1,	NULL,	'k'},
	{ "threads",     1,   NULL,    't'   },
	//{ "hit_max",		1,NULL,	'm'},
	//{ "auto_load", 0, NULL, 'a'},
	{0,	0,	0,	0   }
};

LCtrl_option::LCtrl_option(Options *opt)
{
	opt->len_kmer = 15;
	opt->CandidateRatio = 0;
	opt->thread = 1;
}

int LCtrl_option::Usage()
{
	fprintf(stderr, "\n"); 
	fprintf(stderr, "Program:   Alignment\n"); 
	fprintf(stderr, "Version:   %s\n", PACKAGE_VERSION); 
	fprintf(stderr, "Contact:   %s\n\n", CONTACT); 
	fprintf(stderr, "Usage:     Alignment [Options] <HashIndexDir> <Reads>\n\n"); 
	fprintf(stderr, "Options:   -h, --help                   help\n"); 
 	fprintf(stderr, "           -t, --threads       <int>    number of threads [1]\n"); 
	fprintf(stderr, "           -r, --ratio         <int>    candidate ratio [0.00]\n"); 
	//fprintf(stderr, "           -m, --hit_max       <int>    max hit times of a seed [1000]\n"); 
	fprintf(stderr, "           -k, --kmerSize      <int>    kmer size of hash index [15]\n"); 
	//fprintf(stderr, "           -a, --auto_load              load hash table from hash file without produce hash file");
	//fprintf(stderr, "           -c, --write_cigar            print cigar in XA fields [False]\n"); 
	fprintf(stderr, "\n"); 
	return 0;
}

int LCtrl_option::opt_parse(int argc, char *argv[], Options* opt)
{
	int c;
	int option_index = 0;
	if( argc == 1 )
		return Usage();
	while((c = getopt_long(argc, argv, short_options, long_options, &option_index))>=0)
	{
		switch(c)
		{
			case 'r':
				opt->CandidateRatio = atof(optarg);
				break;
			case 'h':
				return Usage();
				break;
			case 'k':
				opt->len_kmer = atoi(optarg);
				break;
			case 't':
				opt->thread = atoi(optarg);
				break;
			default:
				fprintf(stderr,"Not proper parameters\n");
				return Usage();
		}
	}
	if(optind + 2 != argc)
	{
		fprintf(stderr, "[opt_parse]: index directory, read file can't be omited!\n");
		return 0;
	}
	strncpy(opt->hash_dir, argv[optind++],sizeof(opt->hash_dir));
	strncpy(opt->read_path, argv[optind++],sizeof(opt->read_path));
	return 1;
}

void LCtrl_option::show_parameters(Options* opt)
{
	fprintf(stderr, "\n"); 
	fprintf(stderr, ":::: Simulation parameters :::\n"); 
	fprintf(stderr, "\n"); 
	fprintf(stderr, "CandidateRatio: %.4f\n", opt->CandidateRatio); 
	fprintf(stderr, "threads:        %d\n", opt->thread); 
	fprintf(stderr, "kmerSize:       %u\n", opt->len_kmer); 
	fprintf(stderr, "HashIndexDir:   %s\n", opt->hash_dir); 
	fprintf(stderr, "Reads:          %s\n", opt->read_path); 
	fprintf(stderr, "\n"); 
	fprintf(stderr, ":::: Simulation parameters :::\n"); 
	fprintf(stderr, "\n"); 
}
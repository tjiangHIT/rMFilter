/**  
 * All rights Reserved, Designed By HIT-Bioinformatics   
 * @Title:  Index.cpp
 * @Package 
 * @Description:    rMFilter Index
 * @author: tjiang
 * @date:   June 7 2016
 * @version V1.0     
 */

#include <stdio.h>
#include <iostream>
#include <time.h>
#include <cstdlib>
#include "Ctrl_option.h"
#include "readfl.h"
#include "mk_hash.h"

using namespace std;

const std::string getCurrentDateTime() {
	time_t     now = time(0);
	struct tm  tstruct;
	char       buf[80];
	tstruct = *localtime(&now);
	// Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
	// for more information about date/time format
	strftime(buf, sizeof(buf), "[INFO] %Y-%m-%d,%X", &tstruct);

	return buf;
}

int main(int argc, char *argv[])
{
	Options *opt = new Options;
	Ctrl_option fm(opt);
	if( fm.opt_parse(argc,argv,opt) != 1 )
		exit(1);

	fm.show_parameters(opt);

	fprintf(stderr,"%s Index started\n",getCurrentDateTime().c_str());

	uint32_t len_genome = 0;
	read_file inseq;
	char *genome = inseq.read_ref( opt->ref_path, &len_genome );

	Hash hash;
	int hashtab;
	hashtab = hash.mk_hash(opt->hash_dir, genome, opt->len_kmer, len_genome);
	//cout << hashtab << endl;

	if( NULL != opt )
		delete opt;
	fprintf(stderr,"%s Index ended\n",getCurrentDateTime().c_str());
	return 0;
}
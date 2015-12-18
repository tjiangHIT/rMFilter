/*
	data 2015 6 24
	Alignment the read to find SV read
*/

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdint.h>
#include <sys/time.h> 
#include <stdio.h>
#include "Del_local_aln.h"

using namespace std;

const std::string getCurrentDateTime() {
	time_t     now = time(0);
	struct tm  tstruct;
	char       buf[80];
	tstruct = *localtime(&now);
	strftime(buf, sizeof(buf), "[INFO] %Y-%m-%d,%X", &tstruct);
	return buf;
}

int main(int argc, char *argv[])
{
	Options *opt = new Options;
	LCtrl_option fm(opt);
	if( fm.opt_parse(argc,argv,opt) != 1 )
		exit(1);

	fm.show_parameters(opt);

	fprintf(stderr,"%s Aln started\n",getCurrentDateTime().c_str());
	Local_aln Aln;
	Aln.process(opt);
	fprintf(stderr,"%s Aln ended\n",getCurrentDateTime().c_str());

	if( opt != NULL )
		delete[] opt;
	return 0;
}
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdint.h>
#include <sys/time.h> 
#include <stdio.h>
#include "local_aln.h"

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

int main(int argc, char const *argv[])
{
	fprintf(stderr,"%s Aln started\n",getCurrentDateTime().c_str());
	Local_aln Aln;
	Aln.Files_open("/home/tjiang/", 13);
	int h;
	h = Aln.local_aln("../_0001.fastq", 13);
	Aln.Cleaning();
	fprintf(stderr,"%s Aln ended\n",getCurrentDateTime().c_str());
	return 0;
}
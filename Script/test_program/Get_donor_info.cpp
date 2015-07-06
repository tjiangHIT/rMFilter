#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

using namespace std;

#define readLINE 30000

int main(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "Usage: %s <in.fastq> <read.L>\n", argv[0]);
		return 1;
	}

	gzFile fp;
	kseq_t *Seq;

	fp = gzopen(argv[1], "r");
	Seq = kseq_init(fp);

	ofstream out;
	out.open(argv[2]);

	while (kseq_read(Seq) >= 0) 
	{
		out << Seq->seq.l << endl;
		//cout << Seq->seq.s << endl;
		//break;
	}

	kseq_destroy(Seq);
	gzclose(fp);
	out.close();
	return 0;
}

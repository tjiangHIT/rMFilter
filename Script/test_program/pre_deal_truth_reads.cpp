#include <iostream>
#include <stdlib.h>
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

using namespace std;

#define MinL 1024

int main(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "Usage: %s <in.fastq> <out.fastq>\n", argv[0]);
		return 1;
	}
	/* code */
	gzFile fp;
	kseq_t *Seq;

	fp = gzopen(argv[1], "r");
	Seq = kseq_init(fp);

	FILE *out;
	out = fopen(argv[2], "w");
	int num = 0;

	while (kseq_read(Seq) >= 0) 
	{
		if( Seq->seq.l > MinL )
		{
			fprintf(out, "@%s %s\n", Seq->name.s, Seq->comment.s);
			fprintf(out, "%s\n", Seq->seq.s);
			fprintf(out, "+%s %s\n", Seq->name.s, Seq->comment.s);
			fprintf(out, "%s\n", Seq->qual.s);
			num++;
			//break;
		}
	}

	kseq_destroy(Seq);
	gzclose(fp);
	fclose(out);
	cout << num << endl;
	return 0;
}
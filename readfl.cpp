#include "readfl.h"
#include <iostream>
using namespace std;

#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread);

char *read_file::read_ref(char *path,unsigned *len_genome)
{
        char *genome = NULL;
        gzFile fp;
        kseq_t *trunk;
        fp = gzopen(path, "r");
        trunk = kseq_init(fp);
        if(kseq_read(trunk) >= 0)
        {
            genome = new char[trunk->seq.l+1];
            strncpy(genome,trunk->seq.s,trunk->seq.l);
        }
        *len_genome = trunk->seq.l;
       // printf("%d\t%d\t%d\n", n, slen, qlen);
        kseq_destroy(trunk);
        gzclose(fp);
        return genome;
}

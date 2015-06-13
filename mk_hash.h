#ifndef MK_HASH_
#define MK_HASH_
#include <stdint.h>
#define Hpart 32


typedef struct 
{
	uint32_t SR;
	//uint32_t freq;
}SRdata;

class Hash
{
public:
	int mk_hash(char *path, char *genome, uint8_t kmer, uint32_t len_genome);
};
	extern uint8_t trans[];
	uint32_t transfer(char *genome,uint32_t len_sed,uint32_t start);
#endif
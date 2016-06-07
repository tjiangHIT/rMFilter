/**  
 * All rights Reserved, Designed By HIT-Bioinformatics   
 * @Title:  mk_hash.h
 * @Package 
 * @Description:    make the hashtable
 * @author: tjiang
 * @date:   June 7 2016
 * @version V1.0     
 */

#ifndef MK_HASH_
#define MK_HASH_
#include <stdint.h>
#define Hpart 32

class Hash
{
public:
	int mk_hash(char *path, char *genome, uint32_t kmer, uint32_t len_genome);
};
	extern uint8_t trans[];
	//extern char random_N[];
	uint32_t transfer(char *genome,uint32_t len_sed,uint32_t start);
	char *get_path(char *path, const char *str2);
	int compare(const void *a, const void *b);
#endif
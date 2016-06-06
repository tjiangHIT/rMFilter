/**  
 * All rights Reserved, Designed By HIT-Bioinformatics   
 * @Title:  readfl.h
 * @Package 
 * @Description:    read files
 * @author: tjiang
 * @date:   June 7 2016
 * @version V1.0     
 */

#ifndef READ_H_
#define READ_H_
#include <stdint.h>
class read_file
{
	public:
		char *read_ref(char *path,uint32_t *len_genome, uint32_t *Start_pos, char **chrName,int *count_chrl);
		char *read_ref(char *path,uint32_t *len_genome);
	private:
		void preDealRef(char *ref, uint32_t len);

};

#endif
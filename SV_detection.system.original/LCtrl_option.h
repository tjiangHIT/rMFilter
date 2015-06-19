#include <stdint.h>

#define PATH_LEN 1024

typedef struct 
{
	uint8_t		len_kmer;
	char		read_path[PATH_LEN];
	char		hash_dir[PATH_LEN];
	uint32_t	CandidateNum;
	/*#define SEED_NUM 1000
	#define UpBound 1000
	#define TopChoice 5
	#define Hpart 32
	#define READ_MAX_LENGTH 100000
	#define Extend 1000
	#define COVER_SCORE 100
	#define DISTANCE_SCORE 750
	#define TransParameter 1.1*/
}Options;

class LCtrl_option
{
	Options	*opt;
public:
	LCtrl_option(Options *opt);
	int Usage();
	int opt_parse(int argc, char *argv[], Options* opt);
};
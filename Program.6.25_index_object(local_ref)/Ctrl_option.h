#include <stdint.h>

#define PATH_LEN 1024

typedef struct 
{
	uint8_t	len_kmer;
	char	ref_path[PATH_LEN];
	char	hash_dir[PATH_LEN];
}Options;

class Ctrl_option
{
	Options	*opt;
public:
	Ctrl_option(Options *opt);
	int Usage();
	int opt_parse(int argc, char *argv[], Options* opt);
	void show_parameters(Options* opt);
};
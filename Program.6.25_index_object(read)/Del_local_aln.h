/*
	Data 6 25
*/

#ifndef LOCAL_ALN_
#define LOCAL_ALN_
#include <stdint.h>

//#define K_MER 13
#define SEED_NUM 1000
#define UpBound 1000
//#define TopChoice 5
#define Hpart 32
#define READ_MAX_LENGTH 100000
#define Extend 1000
#define COVER_SCORE 65
#define DISTANCE_SCORE 750
#define TransParameter 1.1

typedef struct 
{
	uint32_t ps;
	uint32_t pe;
	uint32_t pf;
}CF;

typedef struct way_node
{
	uint32_t ID;	//which way
	uint32_t win;	//data
	bool operator<( const way_node & a )const
	{
		return win > a.win;
	}
}Node;

typedef struct node_read
{
	uint32_t	Win_Begin_start;
	uint32_t	Win_Begin_end;
	bool		direction;
	int		cover_score;
	//int		distance_score;
	bool operator<( const node_read & a )const
	{
		return cover_score < a.cover_score;
	}
}PreRead;

typedef struct node_read_
{
	uint32_t	Win_Begin_start;
	uint32_t	Win_Begin_end;
	bool		direction;
	int		cover_score;
	//int		distance_score;
	bool operator<( const node_read_ & a )const
	{
		return cover_score > a.cover_score;
	}
}PreRead_;

typedef struct
{
	uint32_t Kpointer1;
	uint32_t Kpointer2;
}ReadKmerHash;

typedef struct tuple
{
	uint32_t	read_begin;
	uint32_t	ref_begin;
	//uint32_t	tuplelength;
	bool operator<( const tuple & a )const
	{
		//return cover_score < a.cover_score;
		return ( read_begin == a.read_begin ) ? ref_begin > a.ref_begin : read_begin > a.read_begin;
	}
}Tuple;

class Local_aln
{
	char			*genome;
	uint32_t		*kmerHash;
	uint32_t		*SR_HASH;
	uint32_t		Hash2sum;
	uint32_t		len_genome;
	//uint32_t	kmer;
	ReadKmerHash		*readKmerhash;
	uint32_t		*readSRhash;
	uint32_t		*tempArray;
	uint32_t		*tempArrayNext;
	uint32_t 		TopChoice;
	uint32_t		**Result_way;
	char			*RCRead;
	
public:
	void Files_open(char *path, uint32_t kmer);
	int local_aln(char *path, uint32_t kmer);
	void ParaAssign( uint32_t num);
	void Cleaning();
	PreRead *FindStart(char *Seq, bool direction, uint32_t ReadLength, uint32_t kmer);
	int Judge_SV(char *Seq, uint32_t ReadLength, uint32_t kmer, uint32_t GU, uint32_t GD);
	int Tuple_link( Tuple t1, Tuple t2);
};

#endif
/*  Del_local_aln.h
 *  Author: tjiang (tjiang@hit.edu.cn)
 * -------------------------------------------------------------------
 * Description: 
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 18  2015 (rd)
 * Created:  Jan  30 2015 (rd)
 *-------------------------------------------------------------------
 */

#ifndef LOCAL_ALN_
#define LOCAL_ALN_
#include <stdint.h>
#include "ksw.h"
 #include "kseq.h"
 #include <zlib.h>
 #include "LCtrl_option.h"

#define SEED_NUM 1000
#define UpBound 1000
#define Hpart 32
#define READ_MAX_LENGTH 1000000
#define Extend 1000
#define DISTANCE_SCORE 448//{392,448,750}
#define FL_Dis 448//250
#define Bind_choice 50
// #define TransParameter 1//{1,1.1}

#define local_seed_length 11

 KSEQ_INIT(gzFile, gzread)

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
	uint32_t	cover_score;
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
	uint32_t	tuplelength;
	bool operator<( const tuple & a )const
	{
		return ( read_begin == a.read_begin ) ? ref_begin > a.ref_begin : read_begin > a.read_begin;
	}
}Tuple;

class Local_aln
{
	char		*genome;
	uint32_t		*kmerHash;
	uint32_t		*SR_HASH;
	uint32_t		Hash2sum;
	uint32_t		len_genome;
	uint8_t		*transu8;
	
public:
	void process(Options *opt);
	char local_aln(kseq_t *Seq, ReadKmerHash *readKmerhash,  uint32_t *readSRhash, uint32_t *tempArray, 
		uint32_t *tempArrayNext, uint32_t **Result_way, char *RCRead, int *Track, int *Score, Options *opt);
	// int local_aln(char *path, uint32_t kmer);
	PreRead FindStart(char *Seq, bool direction, uint32_t ReadLength, uint32_t kmer, uint32_t **Result_way);
	int Judge_SV(char *Seq, uint32_t ReadLength, uint32_t GU, uint32_t GD, ReadKmerHash *readKmerhash, uint32_t *readSRhash, 
		uint32_t *tempArray, uint32_t *tempArrayNext, int *Track, int *Score);
	int Tuple_link( Tuple t1, Tuple t2);
	int Tuple_link( Tuple t1, Tuple t2, int newWL );
	int transIntoDec(uint8_t *transtr,char *str, uint32_t length, uint32_t start);
	int n_seq_read(kstream_t *_fp, kseq_t *_Seq, int need);
};

typedef struct
{
	int tid;
	Local_aln *aln;
	int n_seq;
	kseq_t *seq;
	ReadKmerHash *readKmerhash;
	uint32_t *readSRhash;
	uint32_t *tempArray;
	uint32_t *tempArrayNext;
	uint32_t **Result_way;
	char *RCRead;
	int *Track;
	int *Score;
	Options *opt;
	char *out_flag;
}ParT;

#endif
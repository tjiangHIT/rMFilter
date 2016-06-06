/*  Del_local_aln.cpp
 *  Author: tjiang (tjiang@hit.edu.cn)
 * -------------------------------------------------------------------
 * Description: 
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 23  2015 (rd)
 * Created:  Jan  30 2015 (rd)
 *-------------------------------------------------------------------
 */

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdint.h>
#include <queue>
#include <vector>
#include <sys/time.h> 
// #include <zlib.h>
#include <stdio.h>
#include "Del_local_aln.h"
#include "mk_hash.h"
#include <pthread.h>

using namespace std;

int read_seq;
pthread_rwlock_t rwlock;

uint8_t rev[128]={
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
	4,84,4,71,4,4,4,67,4,4,4,4,4,4,67,4,
	4,4,4,4,65,4,4,4,4,4,4,4,4,4,4,4,4,
	84,4,71,4,4,4,67,4,4,4,4,4,4,67,4,4,
	4,4,4,65,4,4,4,4,4,4,4,4,4,4,4
};

long getCurrentTime()    
{
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return tv.tv_sec * 1000 + tv.tv_usec / 1000;    
} 

int RevComRead(char *rcread,char *read,int len_read)
{
	int i;
	for(i=0; i<len_read; ++i )
		rcread[i] = rev[read[len_read - 1 - i]];
	rcread[i] = '\0';
	return 0;
}


int Get_Max_range( uint32_t *array )
{
	int i = 0;
	int score = 0;
	int     MaxScore = 0;
	for( i = 0; i < SEED_NUM; ++i )
	{
		if( array[i] == 0 )
			score++;
		else
		{
			if( MaxScore < score )
				MaxScore = score;
			score = 0;
		}
	}
	if( array[i-1] == 0 )
		if( MaxScore < score )
			MaxScore = score;
	return MaxScore;
}

int Local_aln::Tuple_link( Tuple t1, Tuple t2)
{
	int distance_read_s;
	int distance_ref_s;

	distance_read_s = t2.read_begin - t1.read_begin - t1.tuplelength;
	distance_ref_s = t2.ref_begin - t1.ref_begin - t1.tuplelength;

	if( distance_read_s <= DISTANCE_SCORE && distance_read_s >= 0  && distance_ref_s >= 0 ) 
		return 1;
	else
		return 0;
}

int Local_aln::Tuple_link( Tuple t1, Tuple t2, int newWL)
{
	int distance_read_s;
	int distance_ref_s;

	distance_read_s = t2.read_begin - t1.read_begin - t1.tuplelength;
	distance_ref_s = t2.ref_begin - t1.ref_begin - t1.tuplelength;

	if( distance_read_s <= newWL && distance_read_s >= 0  && distance_ref_s >= 0 ) //No link_chioce 
		return 1;
	else
		return 0;
}

PreRead Local_aln::FindStart(char *Seq, bool Direction, uint32_t ReadLength, uint32_t kmer, uint32_t **Result_way)
{
	uint32_t seed_bet = (ReadLength - kmer) / SEED_NUM;
	uint32_t Each_way_count[SEED_NUM];

	for( int i = 0 ; i < SEED_NUM ; ++i )
	{
		uint32_t start = i * seed_bet;
		uint32_t l2r = transfer(Seq,kmer,start);
		
		if( kmerHash[l2r+1] <  kmerHash[l2r])
			Each_way_count[i] = 0;
		else if( (kmerHash[l2r+1] - kmerHash[l2r]) > UpBound )
			Each_way_count[i] = 0;
		else
		{
			uint32_t bias = i * seed_bet;
			bias = bias / Hpart;
			Each_way_count[i] = kmerHash[l2r+1] - kmerHash[l2r];
			for( uint32_t j = 0 ; j < Each_way_count[i] ; j++ )
			{
				Result_way[i][j] = SR_HASH[j+kmerHash[l2r]];
				if( Result_way[i][j] >= bias )
					Result_way[i][j] = Result_way[i][j] - bias;
				else    Result_way[i][j] = 0;
			}
		}
	}

	priority_queue<Node>        windowQueue;
	priority_queue<PreRead_>    pre_read_Queue;

	//uint32_t judge_element = -1, count_element = 0;
	CF region;
	region.ps = -1, region.pe = -1, region.pf = 0;

	uint32_t Each_way_flag[SEED_NUM];
	for( int i = 0 ; i < SEED_NUM ; i++ )
	{
		Each_way_flag[i] = 0;
		if( Each_way_count[i] != 0 )
		{
			Node node;
			node.ID = i;
			node.win = Result_way[i][0];
			windowQueue.push(node);
			Each_way_flag[i]++;
		}
	}

	while( !windowQueue.empty())
	{
		Node node;
		node.ID = windowQueue.top().ID;
		
		if( region.pe <= windowQueue.top().win && windowQueue.top().win <= region.pe + 2 )
		{
			region.pe = windowQueue.top().win;
			region.pf++;
			//con_array[node.ID]++;
		}
		else
		{
			PreRead_ temp_preRead;
			temp_preRead.Win_Begin_start = region.ps;
			temp_preRead.Win_Begin_end = region.pe;
			//temp_preRead.direction = Direction;
			temp_preRead.cover_score = region.pf;
			//cout << region.ps << "\t" << region.pe << "\t" << region.pf << "\n";
			//temp_preRead.distance_score = READ_MAX_LENGTH;
			pre_read_Queue.push(temp_preRead);
			if( pre_read_Queue.size() > 1 )
				pre_read_Queue.pop();

			region.ps = windowQueue.top().win;
			region.pe = windowQueue.top().win;
			region.pf = 1;
		}

		windowQueue.pop();
		if( Each_way_flag[node.ID] < Each_way_count[node.ID] )
		{
			node.win = Result_way[node.ID][Each_way_flag[node.ID]];
			Each_way_flag[node.ID]++;
			if( node.win != 0 )
			{
				windowQueue.push(node);
			}
		}
	}

	// cout << pre_read_Queue.size() << endl;
	//cout << pre_read_Queue.top().Win_Begin_start << "\t" << pre_read_Queue.top().Win_Begin_end << "\t" << pre_read_Queue.top().cover_score << "\n";

	PreRead Top_preRead;
	Top_preRead.Win_Begin_start = pre_read_Queue.top().Win_Begin_start;
	Top_preRead.Win_Begin_end = pre_read_Queue.top().Win_Begin_end;
	Top_preRead.direction = Direction;
	Top_preRead.cover_score = pre_read_Queue.top().cover_score;
	pre_read_Queue.pop();

	//cout << Top_preRead.Win_Begin_start << "\t" << Top_preRead.Win_Begin_end << "\t" << Top_preRead.cover_score << "\n";

	while( !pre_read_Queue.empty() )
		pre_read_Queue.pop();

	return Top_preRead;
}

int Local_aln::Judge_SV(char *Seq, uint32_t ReadLength, uint32_t GU, uint32_t GD, ReadKmerHash *readKmerhash, uint32_t *readSRhash, uint32_t *tempArray, uint32_t *tempArrayNext, int *Track, int *Score)
{
	// 0 stands for sv, 1 stands for normal.
	uint32_t kmer = local_seed_length;

	uint64_t SR_length = ReadLength + 1 - kmer;
	uint32_t pre = 0;
	uint32_t mask = 0xffffffff>>(32 - (kmer<<1));
	//produce all kmer value;
	pre = transfer(Seq, kmer, 0);
	tempArray[0] = pre; 
	tempArrayNext[0] = pre;

	// if you want to mask uncomment the next sentence
	for (uint32_t i=1; i< SR_length; ++i )
	{
		pre = ((pre<<2)&mask)|trans[Seq[i + kmer - 1 ]];
		tempArray[i] = pre;
		tempArrayNext[i] = pre;
	}

	qsort(tempArray, SR_length, sizeof(uint32_t), compare);

	uint32_t flag;
	for( uint32_t i = 0 ; i < SR_length ; ++i )
	{
		flag = tempArray[i];
		readKmerhash[flag].Kpointer1++;
	}

	flag = -1;
	uint64_t sum = 0;
	for( uint32_t i = 0 ; i < SR_length ; ++i )
	{
		if( flag != tempArray[i] )
		{
			flag = tempArray[i];
			uint32_t temp = readKmerhash[flag].Kpointer1;
			readKmerhash[flag].Kpointer1 = sum;
			readKmerhash[flag].Kpointer2 = sum;
			sum += temp;
		}
	}

	for( uint32_t i = 0 ; i < SR_length ; ++i )
	{
		flag = tempArrayNext[i];
		uint32_t add = readKmerhash[flag].Kpointer2;
		readSRhash[add] = i;
		readKmerhash[flag].Kpointer2++;
	}

	priority_queue<Tuple>   Tuple_Queue;    
	pre = transfer(genome, kmer, GU);
	pre = ((pre>>2)&mask);  
	for( uint32_t i = 0 ; i < GD - GU ; ++i )
	{
		Tuple temp_node;
		temp_node.ref_begin = i + GU;
		pre = ((pre<<2)&mask)|trans[genome[i + GU + kmer - 1 ]];
		uint32_t how_many_node = readKmerhash[pre].Kpointer2 - readKmerhash[pre].Kpointer1;
		for( uint32_t k = 0 ; k < how_many_node ; ++k )
		{
			temp_node.read_begin = readSRhash[k + readKmerhash[pre].Kpointer1];
			//maybe change the shorter substr
			Tuple_Queue.push(temp_node);
		}
	}

	for( uint32_t i = 0 ; i < SR_length ; i++ )
	{
		readKmerhash[tempArray[i]].Kpointer1 = 0;
		readKmerhash[tempArray[i]].Kpointer2 = 0;
	}

	if ( Tuple_Queue.size() > 100000 )
	{
		// cout << "Large queue." << endl;
		while( !Tuple_Queue.empty() )
			Tuple_Queue.pop();
		return 2;
	}

	
	if( Tuple_Queue.top().read_begin > FL_Dis )
	{
		// tempOutNum++;
		while( !Tuple_Queue.empty() )
			Tuple_Queue.pop();
		return 0;
	}
	else
	{
		vector<Tuple> Tuple_Vector;
		Tuple temp_node;
		temp_node.ref_begin = Tuple_Queue.top().ref_begin;
		temp_node.read_begin = Tuple_Queue.top().read_begin;
		temp_node.tuplelength = 1;
		Tuple_Queue.pop();
		Tuple_Vector.push_back(temp_node);

		while( !Tuple_Queue.empty() )
		{
			int flag = 0;
			for( int i = 0 ; i < Tuple_Vector.size() ; ++i )
			{
				uint32_t V_read = Tuple_Vector[i].read_begin + Tuple_Vector[i].tuplelength;
				uint32_t V_ref = Tuple_Vector[i].ref_begin + Tuple_Vector[i].tuplelength;
				uint32_t Q_read = Tuple_Queue.top().read_begin;
				uint32_t Q_ref = Tuple_Queue.top().ref_begin;
				if( Q_read -  V_read == 0 && Q_ref - V_ref == 0 )
				{
					Tuple_Vector[i].tuplelength++;
					flag = 1;
				}
			}
			if( flag == 0 )
			{
				temp_node.ref_begin = Tuple_Queue.top().ref_begin;
				temp_node.read_begin = Tuple_Queue.top().read_begin;
				temp_node.tuplelength = 1;
				Tuple_Vector.push_back(temp_node);
			}
			Tuple_Queue.pop();
		}

		for( int i = 0 ; i < Tuple_Vector.size() ; ++i )
		{
			Tuple_Vector[i].tuplelength = Tuple_Vector[i].tuplelength + kmer - 1;
			Track[i] = -1;
			Score[i] = 0;
		}

		if( ReadLength - Tuple_Vector[Tuple_Vector.size()-1].read_begin - Tuple_Vector[Tuple_Vector.size()-1].tuplelength > FL_Dis )
		{
			// tempOutNum++;
			if( !Tuple_Vector.empty() )
			{
				Tuple_Vector.clear();
			}
			return 0;
		}

		int maxScore = 0;
		int maxScore_pos = 0;

		for( int i = 1 ; i < Tuple_Vector.size() ; ++i )
		{
			for( int k =  i - 1 ; k >= 0 ; k-- )
			{
				int local_ans = 0;
				local_ans = Tuple_link( Tuple_Vector[k], Tuple_Vector[i]);
				if( local_ans == 1 )
				{
					if( Score[i] < Score[k] + Tuple_Vector[k].tuplelength )
					{
						Score[i] = Score[k] + Tuple_Vector[k].tuplelength;
						Track[i] = k;
						if( Score[i] >= maxScore )
						{
							maxScore = Score[i];
							maxScore_pos = i;
						}
					}
				}
			}
		}

		// for( int i = 0 ; i < Tuple_Vector.size() ; i++ )
		// {
		//  cout << i << "\t";
		//  cout << Tuple_Vector[i].read_begin << "\t";
		//  cout << Tuple_Vector[i].ref_begin << "\t";
		//  cout << Tuple_Vector[i].tuplelength << "\t";
		//  cout << Track[i] << "\t";
		//  cout << Score[i] << endl;
		// }
		// cout << maxScore << "\t" << maxScore_pos << endl;
		// cout << Tuple_Vector[maxScore_pos].read_begin << "\t" << Tuple_Vector[maxScore_pos].ref_begin << endl;

		if( maxScore == 0 && maxScore_pos == 0 )
		{
			if( !Tuple_Vector.empty() )
				Tuple_Vector.clear();
			return 0;
		}

		if( FL_Dis < ReadLength - Tuple_Vector[maxScore_pos].read_begin - Tuple_Vector[maxScore_pos].tuplelength )
		{
			// tempOutNum++;
			if( !Tuple_Vector.empty() )
				Tuple_Vector.clear();
			return 0;
		}
		else
		{
			//binary search
			int flag;
			for( int i = Track[maxScore_pos] ; i >= 0 ; i = Track[i] )
			// for( int i = maxScore_pos ; i > 0 ; i = Track[i] )
			{
				// cout << Tuple_Vector[i].read_begin << "\t";
				// cout << Tuple_Vector[i].ref_begin << "\t";
				flag = i;
				if( -1 == Track[i] )
					break;
				int link_chioce;
				uint32_t readDiff = Tuple_Vector[i].read_begin - Tuple_Vector[Track[i]].read_begin - Tuple_Vector[Track[i]].tuplelength;
				uint32_t refDiff = Tuple_Vector[i].ref_begin - Tuple_Vector[Track[i]].ref_begin - Tuple_Vector[Track[i]].tuplelength;
				link_chioce = readDiff - refDiff;
				// cout << link_chioce << "\t" << flag << endl;
				if(link_chioce <= Bind_choice && link_chioce >= -Bind_choice)
					continue;
				else
				{
					// cout << ReadLength << "\t" << link_chioce << endl;
					break;
				}
			}
			if( Tuple_Vector[flag].read_begin > FL_Dis )
			{
				if( !Tuple_Vector.empty() )
					Tuple_Vector.clear();
				return 0;
			}
		}

		if( !Tuple_Vector.empty() )
			Tuple_Vector.clear();
	}
	return 1;
}

PreRead Local_aln::acquirePre(kseq_t *Seq, uint32_t **Result_way, char *RCRead, Options *opt)
{
	PreRead     preread_1;
	PreRead     preread_2;

	RevComRead(RCRead, Seq->seq.s, Seq->seq.l);
	preread_1 = FindStart(Seq->seq.s, true, Seq->seq.l, opt->len_kmer, Result_way);
	preread_2 = FindStart(RCRead, false, Seq->seq.l, opt->len_kmer, Result_way);

	if( preread_1.cover_score > preread_2.cover_score ) return preread_1;
	else return preread_2;
}

void Local_aln::deal_preread(PreRead pre, uint32_t count)
{
	localscore[count] = pre.cover_score;
	localPreRead[count].Win_Begin_start = pre.Win_Begin_start;
	localPreRead[count].Win_Begin_end = pre.Win_Begin_end;
	localPreRead[count].direction = pre.direction;
	localPreRead[count].cover_score = pre.cover_score;
}

char Local_aln::local_aln(kseq_t *Seq, PreRead preread, ReadKmerHash *readKmerhash,  uint32_t *readSRhash, uint32_t *tempArray, uint32_t *tempArrayNext, char *RCRead, int *Track, int *Score, Options *opt)
{
	if( Threshold_below > preread.cover_score || Threshold_up < preread.cover_score )
	{
		// cout << Seq->name.s << "\tT\t" << preread.cover_score << endl;
		return 'T';
	}
	RevComRead(RCRead, Seq->seq.s, Seq->seq.l);
	int SV_flag = 1;
	uint32_t start = preread.Win_Begin_start;
	uint32_t end = preread.Win_Begin_end;
	start = start * Hpart;
	end = end * Hpart + Extend + Seq->seq.l;
	if( start < Extend ) start = 0;
	else start = start - Extend;
	if( end > len_genome - opt->len_kmer )
		end = len_genome - opt->len_kmer;
	int judge;
	if( preread.direction == true )
		judge = Judge_SV(Seq->seq.s, Seq->seq.l, start, end, readKmerhash, readSRhash, tempArray, tempArrayNext, Track, Score);
	else
		judge = Judge_SV(RCRead, Seq->seq.l, start, end, readKmerhash, readSRhash, tempArray, tempArrayNext, Track, Score);
	SV_flag = SV_flag * judge;
	if( SV_flag == 1 )
	{
		// cout << Seq->name.s << "\tN\t"  << preread.cover_score << endl;
		return 'N';
	}
	else if( SV_flag == 2 )
	{
		// cout << Seq->name.s << "\tJ\t" << preread.cover_score << endl;
		return 'J';
	}
	else
	{
		// cout << Seq->name.s << "\tT\t" << preread.cover_score << endl;
		return 'T';
	}
	return 'I';
}

int Local_aln::n_seq_read(kstream_t *_fp, kseq_t *_Seq, int n_needed)
{
	kseq_t *temp = _Seq;
	int i = 0;  
	while( i <n_needed && (temp[i].f = _fp) && kseq_read(temp+i)>=0 ) ++i;
	return i;
}

uint64_t transKmer2Length( uint32_t kmer )
{
	uint64_t length;
	length = 1 << ( kmer << 1 );
	length += 1;
	return length;
}

ParT *thread_initiate(ReadKmerHash **readKmerhash, uint32_t **readSRhash, uint32_t **tempArray, uint32_t **tempArrayNext, uint32_t ***Result_way, char **RCRead, int **Track, int **Score, Options *opt, Local_aln *aln)
{
	ParT *part_thread = new ParT[opt->thread];
	for ( int i = 0;i < opt->thread; ++i ) {
		part_thread[i].opt = opt;
		part_thread[i].aln = aln;
		part_thread[i].readKmerhash = readKmerhash[i];
		part_thread[i].readSRhash = readSRhash[i];
		part_thread[i].tempArray = tempArray[i];
		part_thread[i].tempArrayNext = tempArrayNext[i];
		part_thread[i].Result_way = Result_way[i];
		part_thread[i].RCRead = RCRead[i];
		part_thread[i].Track = Track[i];
		part_thread[i].Score = Score[i];
	}
	return part_thread;
}

void *thread_worker(void *data)
{
	ParT *part_thread = (ParT *)data;
	int _read_seq;
	while (1) {
		pthread_rwlock_wrlock(&rwlock);
		_read_seq = read_seq++;
		pthread_rwlock_unlock(&rwlock);
		if (_read_seq < part_thread->n_seq)
		{
			char flag;
			flag = part_thread->aln->local_aln(part_thread->seq + _read_seq, part_thread->aln->localPreRead[part_thread->round_n + _read_seq], part_thread->readKmerhash,  part_thread->readSRhash, part_thread->tempArray, part_thread->tempArrayNext, part_thread->RCRead, part_thread->Track, part_thread->Score, part_thread->opt);
			part_thread->out_flag[_read_seq] = flag;
		}
		else
			break;
	}
}

void *thread_worker_score(void *data)
{
	ParT *part_thread = (ParT *)data;
	int _read_seq;
	while (1) {
		pthread_rwlock_wrlock(&rwlock);
		_read_seq = read_seq++;
		pthread_rwlock_unlock(&rwlock);
		if (_read_seq < part_thread->n_seq)
		{
			PreRead thispreread;
			thispreread = part_thread->aln->acquirePre(part_thread->seq + _read_seq, part_thread->Result_way, part_thread->RCRead, part_thread->opt);
			part_thread->aln->deal_preread(thispreread, part_thread->round_n + _read_seq);
		}
		else
			break;
	}
}

void output( kseq_t *Seq, char flag )
{
	if( flag == 'T' || flag == 'J' )
	{
		if(Seq->qual.l == 0)
		{
			cout << ">" << Seq->name.s << "\n";
			int pointer = 0;
			while( pointer < Seq->seq.l )
			{
				cout << Seq->seq.s[pointer];
				pointer++;
				if( pointer % 80 == 0 ) cout << endl;
			}
			if( Seq->seq.l % 80 != 0 ) cout << endl;
		}
		else cout << "@" << Seq->name.s << "\n" << Seq->seq.s << "\n+\n" << Seq->qual.s << endl;
	}
}

void Local_aln::process(Options *opt)
{
	fstream LengthFile;

	char *LF_path = NULL;
	LF_path = get_path(opt->hash_dir, "LengthFile");
	LengthFile.open(LF_path);
	delete[] LF_path;
	char temp[1000];
	string str;
	LengthFile.getline(temp, 1000, '\n');
	str = temp;
	Hash2sum = atoi(str.c_str());
	LengthFile.getline(temp, 1000, '\n');
	str = temp;
	len_genome = atoi(str.c_str());
	LengthFile.close();

	genome = new char[len_genome];
	if (NULL == genome) { fprintf(stderr,"Fail to alloc genome, now exit"); exit(1); } 
	FILE * file_ref;
	char *ref_path = NULL;
	ref_path = get_path(opt->hash_dir, "Reference");
	file_ref = fopen(ref_path, "rb");
	delete[] ref_path;
	fread( genome, sizeof(char), len_genome, file_ref );
	fclose(file_ref);
	if (NULL == genome) { fprintf(stderr,"Fail to load Genome reference, now exit"); exit(1); } 

	// uint32_t movement = 2 * len_kmer;
	uint64_t Kmer_code_length = transKmer2Length( opt->len_kmer );
	// Kmer_code_length += 1;
	kmerHash = new uint32_t[Kmer_code_length];
	if (NULL == kmerHash) { fprintf(stderr,"Fail to alloc kmerHash, now exit");     exit(1); } 

	FILE * file_hash_1;
	char *hash1_path = NULL;
	hash1_path = get_path(opt->hash_dir, "Hash_1");
	file_hash_1 = fopen(hash1_path, "rb");
	delete[] hash1_path;
	fread( kmerHash, sizeof(uint32_t), Kmer_code_length, file_hash_1 );
	fclose(file_hash_1);
	if (NULL == kmerHash) { fprintf(stderr,"Fail to load kmerHash, now exit");  exit(1); }

	SR_HASH = new uint32_t[Hash2sum];
	if (NULL == SR_HASH) { fprintf(stderr,"Fail to alloc SR_HASH, now exit"); exit(1); } 

	FILE * file_hash_2;
	char *hash2_path = NULL;
	hash2_path = get_path(opt->hash_dir, "Hash_2");
	file_hash_2 = fopen(hash2_path, "rb");
	delete[] hash2_path;
	fread( SR_HASH, sizeof(uint32_t), Hash2sum, file_hash_2 );
	fclose(file_hash_2); 
	if (NULL == SR_HASH) { fprintf(stderr,"Fail to load SR_HASH, now exit"); exit(1); }

	gzFile fp_cal;
	read_cal = 0;
	fp_cal = gzopen(opt->read_path, "r");
	kseq_t *seq_cal = kseq_init(fp_cal);
	while (kseq_read(seq_cal) >= 0) read_cal++;
	kseq_destroy(seq_cal);
	gzclose(fp_cal);

	localscore = new uint32_t[read_cal];
	if (NULL == localscore) { fprintf(stderr,"Fail to alloc localscore, now exit"); exit(1); } 

	localPreRead = new PreRead[read_cal];
	if (NULL == localPreRead) { fprintf(stderr,"Fail to alloc localPreRead, now exit"); exit(1); } 

	gzFile fp_score;
	fp_score = gzopen(opt->read_path, "r");

	gzFile fp;
	fp = gzopen(opt->read_path, "r");

	if( opt->thread <= 1)
	{
		uint64_t Kmer_read_length = transKmer2Length( local_seed_length );

		ReadKmerHash *readKmerhash = new ReadKmerHash[Kmer_read_length];
		if (NULL == readKmerhash) { fprintf(stderr,"Fail to alloc readKmerhash, now exit"); exit(1); } 
		for( uint32_t i = 0 ; i < Kmer_read_length ; ++i ) readKmerhash[i].Kpointer1 = 0;

		uint32_t* readSRhash = new uint32_t[READ_MAX_LENGTH];
		if (NULL == readSRhash) { fprintf(stderr,"Fail to alloc readSRhash, now exit"); exit(1); }

		uint32_t* tempArray = new uint32_t[READ_MAX_LENGTH];
		if (NULL == tempArray) { fprintf(stderr,"Fail to alloc tempArray, now exit"); exit(1); } 

		uint32_t* tempArrayNext = new uint32_t[READ_MAX_LENGTH];
		if (NULL == tempArrayNext) { fprintf(stderr,"Fail to alloc tempArrayNext, now exit"); exit(1); }

		// transu8 = new uint8_t[READ_MAX_LENGTH];
		// if (NULL == transu8) { fprintf(stderr,"Fail to alloc transu8, now exit"); exit(1); }

		uint32_t** Result_way = new uint32_t*[SEED_NUM]; 
		for( int i = 0 ; i < SEED_NUM ; ++i ) Result_way[i] = new uint32_t[UpBound];

		char* RCRead = new char[READ_MAX_LENGTH];
		if (NULL == RCRead) { fprintf(stderr,"Fail to alloc RCRead, now exit"); exit(1); }

		int* Track = new int[READ_MAX_LENGTH];
		if (NULL == Track) { fprintf(stderr,"Fail to alloc Track, now exit"); exit(1); } 
		int* Score = new int[READ_MAX_LENGTH];
		if (NULL == Score) { fprintf(stderr,"Fail to alloc Score, now exit"); exit(1); } 

		kseq_t *Seq_score = kseq_init(fp_score);
		uint32_t cycletime = 0;
		while(kseq_read(Seq_score) >= 0)
		{
			PreRead thispreread;
			thispreread = acquirePre(Seq_score, Result_way, RCRead, opt);
			deal_preread(thispreread, cycletime);
			cycletime++;
		}
		kseq_destroy(Seq_score);

		//sort
		qsort(localscore, read_cal, sizeof(uint32_t), compare);
		uint32_t Threshold_pos_below = read_cal * opt->CandidateRatio;
		uint32_t Threshold_pos_up = read_cal * ( 1.0 - opt->CandidateRatio );
		if( Threshold_pos_up == read_cal ) Threshold_pos_up = Threshold_pos_up - 1;
		Threshold_below = localscore[Threshold_pos_below];
		Threshold_up = localscore[Threshold_pos_up];

		// cout << Threshold << endl;
		// cout << Threshold_below << "\t" << Threshold_up << endl;

		kseq_t *Seq = kseq_init(fp);
		cycletime = 0;
		while (kseq_read(Seq) >= 0) 
		{
			char flag;
			flag = local_aln(Seq, localPreRead[cycletime], readKmerhash,  readSRhash, tempArray, tempArrayNext, RCRead, Track, Score, opt);
			cycletime++;
			output( Seq, flag);
		}
		kseq_destroy(Seq);

		if (NULL != readKmerhash) delete[] readKmerhash;
		if (NULL != readSRhash) delete[] readSRhash;
		if (NULL != tempArray) delete[] tempArray;
		if (NULL != tempArrayNext) delete[] tempArrayNext;
		// if (NULL != transu8) delete[] transu8;
		if (NULL != RCRead) delete[] RCRead;
		if (NULL != Result_way) {
			for (int i=0;i<SEED_NUM;++i) { delete[] Result_way[i]; }
			delete[] Result_way;
		}
		if( Score != NULL ) delete[] Score;
		if( Track != NULL ) delete[] Track;
	}
	else
	{
		uint64_t Kmer_read_length = transKmer2Length( local_seed_length );
		ReadKmerHash **readKmerhash = new ReadKmerHash*[opt->thread];
		if (NULL == readKmerhash) { fprintf(stderr,"Fail to alloc readKmerhash, now exit"); exit(1); } 
		for( int i = 0; i< opt->thread; ++i )
			readKmerhash[i] = new ReadKmerHash[Kmer_read_length];
		for( int i = 0; i< opt->thread; ++i )
		{
			for( int j = 0; j < Kmer_read_length; ++j ) readKmerhash[i][j].Kpointer1 = 0;
		}
		uint32_t** readSRhash = new uint32_t*[opt->thread];
		if (NULL == readSRhash) { fprintf(stderr,"Fail to alloc readSRhash, now exit"); exit(1); } 
		for( int i = 0; i< opt->thread; ++i )
			readSRhash[i] = new uint32_t[READ_MAX_LENGTH];
		uint32_t** tempArray = new uint32_t*[opt->thread];
		if (NULL == tempArray) { fprintf(stderr,"Fail to alloc tempArray, now exit"); exit(1); } 
		for( int i = 0; i< opt->thread; ++i )
			tempArray[i] = new uint32_t[READ_MAX_LENGTH];
		uint32_t** tempArrayNext = new uint32_t*[opt->thread];
		if (NULL == tempArrayNext) { fprintf(stderr,"Fail to alloc tempArrayNext, now exit"); exit(1); } 
		for( int i = 0; i< opt->thread; ++i )
			tempArrayNext[i] = new uint32_t[READ_MAX_LENGTH];
		uint32_t*** Result_way = new uint32_t**[opt->thread];
		if (NULL == Result_way) { fprintf(stderr,"Fail to alloc Result_way, now exit"); exit(1); } 
		for( int i = 0; i< opt->thread; ++i )
		{
			Result_way[i] = new uint32_t*[SEED_NUM];
			for( int j = 0; j < SEED_NUM; ++j )
				Result_way[i][j] = new uint32_t[UpBound];
		}
		char** RCRead = new char*[opt->thread];
		if (NULL == RCRead) { fprintf(stderr,"Fail to alloc RCRead, now exit"); exit(1); } 
		for( int i = 0; i< opt->thread; ++i )
			RCRead[i] = new char[READ_MAX_LENGTH];
		int** Track = new int*[opt->thread];
		if (NULL == Track) { fprintf(stderr,"Fail to alloc Track, now exit"); exit(1); } 
		for( int i = 0; i< opt->thread; ++i )
			Track[i] = new int[READ_MAX_LENGTH];
		int** Score = new int*[opt->thread];
		if (NULL == Score) { fprintf(stderr,"Fail to alloc Score, now exit"); exit(1); } 
		for( int i = 0; i< opt->thread; ++i )
			Score[i] = new int [READ_MAX_LENGTH];


		//step2
		kstream_t *_fp_score = ks_init(fp_score);
		int seq_num;
		int n_needed = 5000;
		kseq_t *seq_socre = (kseq_t *)calloc(n_needed, sizeof(kseq_t));
		if( NULL == seq_socre ) { fprintf(stderr, "Failed when applying for new space! now exit"); exit(1);}

		ParT *part_thread = thread_initiate(readKmerhash, readSRhash, tempArray, tempArrayNext, Result_way, RCRead, Track, Score, opt, this);
		uint32_t cycletime = 0;
		pthread_rwlock_init(&rwlock, NULL);
		while((seq_num = n_seq_read(_fp_score, seq_socre, n_needed))>0)
		{
			read_seq = 0;
			pthread_t *tid;
			tid = (pthread_t*)calloc(opt->thread, sizeof(pthread_t));
			for (int j = 0; j < opt->thread; ++j)
			{
				part_thread[j].tid = j;
				part_thread[j].n_seq = seq_num;
				part_thread[j].seq = seq_socre;
				part_thread[j].round_n = cycletime * n_needed;
				pthread_create(&tid[j], NULL, thread_worker_score, part_thread + j);
			}
			for (int j = 0; j < opt->thread; ++j) pthread_join(tid[j], 0);
			free(tid);
			cycletime++;
		}
		pthread_rwlock_destroy(&rwlock);
		if (NULL != seq_socre) free(seq_socre);
		if (NULL != _fp_score)  free(_fp_score);

		qsort(localscore, read_cal, sizeof(uint32_t), compare);
		uint32_t Threshold_pos_below = read_cal * opt->CandidateRatio;
		uint32_t Threshold_pos_up = read_cal * ( 1.0 - opt->CandidateRatio );
		if( Threshold_pos_up == read_cal ) Threshold_pos_up = Threshold_pos_up - 1;
		Threshold_below = localscore[Threshold_pos_below];
		Threshold_up = localscore[Threshold_pos_up];

		// cout << read_cal << endl;
		// cout << opt->CandidateRatio << "\t" << 1.0 - opt->CandidateRatio << endl;
		// cout << Threshold_below << "\t" << Threshold_up << endl;

		// cout << Threshold << endl;

		//deal
		cycletime = 0;
		kstream_t *_fp = ks_init(fp);
		int n_Seq;
		kseq_t *Seq = (kseq_t *)calloc(n_needed, sizeof(kseq_t));
		if (Seq == NULL) { fprintf(stderr, "Failed when applying for new space! now exit"); exit(1);}
		char *out_flag = new char[n_needed];
		if (NULL == out_flag) { fprintf(stderr,"Fail to alloc out_flag, now exit"); exit(1); } 

		// ParT *part_thread = thread_initiate(readKmerhash, readSRhash, tempArray, tempArrayNext, Result_way, RCRead, Track, Score, opt, this);
		pthread_rwlock_init(&rwlock, NULL);

		while((n_Seq = n_seq_read(_fp, Seq, n_needed))>0)
		{
			read_seq = 0;
			pthread_t *tid;
			tid = (pthread_t*)calloc(opt->thread, sizeof(pthread_t));
			for (int j = 0; j < opt->thread; ++j)
			{
				part_thread[j].tid = j;
				part_thread[j].n_seq = n_Seq;
				part_thread[j].seq = Seq;
				part_thread[j].out_flag = out_flag;
				part_thread[j].round_n = cycletime * n_needed;
				pthread_create(&tid[j], NULL, thread_worker, part_thread + j);
			}
			for (int j = 0; j < opt->thread; ++j) pthread_join(tid[j], 0);
			free(tid);
			for( int i = 0; i < n_Seq; ++i ) output(Seq+i, out_flag[i]);
			cycletime++;
		}
		//free
		pthread_rwlock_destroy(&rwlock);
		if (NULL != Seq) free(Seq);
		if (NULL != _fp)  free(_fp);
		if (NULL != out_flag) delete[] out_flag;


		if (NULL != readKmerhash) {
			for (int i=0;i<opt->thread;++i) { delete[] readKmerhash[i]; }
			delete[] readKmerhash;
		}
		if (NULL != readSRhash) {
			for (int i=0;i<opt->thread;++i) { delete[] readSRhash[i]; }
			delete[] readSRhash;
		}
		if (NULL != tempArray) {
			for (int i=0;i<opt->thread;++i) { delete[] tempArray[i]; }
			delete[] tempArray;
		}
		if (NULL != tempArrayNext) {
			for (int i=0;i<opt->thread;++i) { delete[] tempArrayNext[i]; }
			delete[] tempArrayNext;
		}
		if (NULL != Result_way) {
			for (int i=0;i<opt->thread;++i) 
			{
				for( int j = 0; j < SEED_NUM; ++j )
					delete[] Result_way[i][j]; 
				delete[] Result_way[i]; 
			}
			delete[] Result_way;
		}
		if (NULL != RCRead) {
			for (int i=0;i<opt->thread;++i) { delete[] RCRead[i]; }
			delete[] RCRead;
		}
		if (NULL != Track) {
			for (int i=0;i<opt->thread;++i) { delete[] Track[i]; }
			delete[] Track;
		}
		if (NULL != Score) {
			for (int i=0;i<opt->thread;++i) { delete[] Score[i]; }
			delete[] Score;
		}
	}
	gzclose(fp_score);
	gzclose(fp);

	if (NULL != genome) delete[] genome;
	if (NULL != kmerHash) delete[] kmerHash;
	if (NULL != SR_HASH) delete[] SR_HASH;
	if (NULL != localscore) delete[] localscore;
	if (NULL != localPreRead) delete[] localPreRead;
	// return Threshold;
}

/*
	Data 2015/6/26
	local aln change the index object
*/

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdint.h>
#include <queue>
#include <sys/time.h> 
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include "Del_local_aln.h"
#include "mk_hash.h"

using namespace std;

KSEQ_INIT(gzFile, gzread)

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
	for(i=0; i<len_read; i++)
		rcread[i] = rev[read[len_read - 1 - i]];
	rcread[i] = '\0';
	return 0;
}


int Get_Max_range( uint32_t *array )
{
	int	i = 0;
	int	score = 0;
	int 	MaxScore = 0;
	for( i = 0; i < SEED_NUM; i++ )
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

void Local_aln::ParaAssign(uint32_t num)
{
	TopChoice = num;
}

PreRead *Local_aln::FindStart(char *Seq, bool Direction, uint32_t ReadLength, uint32_t kmer)
{
	uint32_t	seed_bet = (ReadLength - kmer) / SEED_NUM;
	//uint32_t 	Result_way[SEED_NUM][UpBound];
	uint32_t 	Each_way_count[SEED_NUM];

	for( int i = 0 ; i < SEED_NUM ; i++ )
	{
		uint32_t start = i * seed_bet;
		uint32_t l2r = transfer(Seq,kmer,start);
		// should be better
		
		if( kmerHash[l2r+1] <  kmerHash[l2r])
			Each_way_count[i] = 0;
		else if( (kmerHash[l2r+1] - kmerHash[l2r]) > UpBound )
			Each_way_count[i] = 0;
		else
		{
			int bias = i * seed_bet / TransParameter;
			bias = bias / Hpart;
			Each_way_count[i] = kmerHash[l2r+1] - kmerHash[l2r];
			for( int j = 0 ; j < Each_way_count[i] ; j++ )
			{
				Result_way[i][j] = SR_HASH[j+kmerHash[l2r]];
				if( Result_way[i][j] >= bias )
					Result_way[i][j] = Result_way[i][j] - bias;
				else	Result_way[i][j] = 0;
			}
		}
	}

	priority_queue<Node>		windowQueue;
	priority_queue<PreRead_>	pre_read_Queue;

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
			//temp_preRead.distance_score = READ_MAX_LENGTH;
			pre_read_Queue.push(temp_preRead);
			if( pre_read_Queue.size() > TopChoice )
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

	PreRead *Top_preRead  = new PreRead[TopChoice];
	if( NULL == Top_preRead )
	{
		fprintf(stderr,"Fail to alloc Top_preRead, now exit");
		exit(1);
	}

	//cout << WinQueue_L << endl;

	for( int i = 0 ; i < TopChoice ; i++ )
	{
		//uint32_t con_array [1000] = {0};
		Top_preRead[i].Win_Begin_start = pre_read_Queue.top().Win_Begin_start;
		Top_preRead[i].Win_Begin_end = pre_read_Queue.top().Win_Begin_end;
		Top_preRead[i].direction = Direction;
		Top_preRead[i].cover_score = pre_read_Queue.top().cover_score;
		// for( int k = 0 ; k < SEED_NUM ; k++ )
		// {
		// 	// Judge the bonder maybe useful.

		// 	for( int l = 0 ; l < Each_way_count[k] ; l++ )
		// 	{
		// 		if( Result_way[k][l] >= Top_preRead[i].Win_Begin_start && Result_way[k][l] <= Top_preRead[i].Win_Begin_end )
		// 		{
		// 			con_array[k]++;
		// 			break;
		// 		}
		// 	}

		// 	/*if( WinQueue[k].win >= Top_preRead[i].Win_Begin_start &&  WinQueue[k].win <= Top_preRead[i].Win_Begin_end )
		// 		con_array[WinQueue[k].ID]++;*/
		// }
		
		/*for(int k = 0 ; k < 1000 ; k++)
			cout  << k << " " << con_array[k] << endl;*/

		//Top_preRead[i].distance_score = Get_Max_range(con_array) * seed_bet;
		pre_read_Queue.pop();
		//break;
	}
	while( !pre_read_Queue.empty() )
		pre_read_Queue.pop();

	return Top_preRead;
}


int Local_aln::local_aln(char *path, uint32_t kmer)
{

	//long start_time = getCurrentTime();

	ofstream Answer;
	char *Final_path = NULL;
	Final_path = get_path(path, ".FinalAnswer_local_change");
	Answer.open(Final_path);
	delete[] Final_path;
	//Answer.open("../Final_Answer");

	gzFile fp;
	kseq_t *Seq;

	fp = gzopen(path, "r");
	Seq = kseq_init(fp);

	int		cycle_time = 0;
	uint32_t	ReadLength;
	char		*PositiveRead = NULL;
	//char		*RCRead = NULL;
	//char		*RCRead = new char[READ_MAX_LENGTH];
	if (NULL == RCRead) {
		fprintf(stderr,"Fail to allocate RCRead, now exit");
		exit(1);
	} 

	while (kseq_read(Seq) >= 0) 
	{
		ReadLength = Seq->seq.l;
		PositiveRead = Seq->seq.s;
		RevComRead(RCRead, Seq->seq.s, ReadLength);

		PreRead *preread_1 = NULL;
		PreRead *preread_2 = NULL;
		preread_1 = FindStart(PositiveRead, true, ReadLength, kmer);
		preread_2 = FindStart(RCRead, false, ReadLength, kmer);

		/*for( int i = 0 ; i < 5 ; i++ )
		{
			cout << preread_1[i].Win_Begin_start << "\t";
			cout << preread_1[i].Win_Begin_end << "\t";
			cout << preread_1[i].direction << "\t";
			cout << preread_1[i].distance_score << "\t";
			cout << preread_1[i].cover_score << "\n";

			cout << preread_2[i].Win_Begin_start << "\t";
			cout << preread_2[i].Win_Begin_end << "\t";
			cout << preread_2[i].direction << "\t";
			cout << preread_2[i].distance_score << "\t";
			cout << preread_2[i].cover_score << "\n";
		}
		return 0;*/

		priority_queue<PreRead>	final_read;

		for( int i = 0 ; i < TopChoice ; i++ )
		{
			final_read.push(preread_1[i]);
			final_read.push(preread_2[i]);
		}
		if( preread_1 != NULL )
			delete[] preread_1;
		if( preread_2 != NULL )
			delete[] preread_2;

		int SV_flag = 0;
		for( int i = 0 ; i < TopChoice ; i++ )
		{
			int Score_A = final_read.top().cover_score;
			//int Score_B = final_read.top().distance_score;
			//if( Score_A >= COVER_SCORE && Score_B <= DISTANCE_SCORE )
			if( Score_A >= COVER_SCORE )
			{
				/*Answer << cycle_time + 1 << "\t";
				Answer << ReadLength << "\n";
				Answer << final_read.top().Win_Begin_start << "\t";
				Answer << final_read.top().Win_Begin_end << "\t";
				Answer << final_read.top().cover_score << "\n";*/

				SV_flag = 1;
				uint32_t start = final_read.top().Win_Begin_start;
				uint32_t end = final_read.top().Win_Begin_end;
				start = start * Hpart;
				end = end * Hpart + Extend + ReadLength;
				if( start < Extend )
					start = 0;
				else
					start = start - Extend;
				if( end > len_genome - kmer )
					end = len_genome - kmer;
				int judge;
				if( final_read.top().direction == true )
					judge = Judge_SV(PositiveRead, ReadLength, kmer, start, end);// 0 stands for sv, 1 stands for normal.
				if( final_read.top().direction == false )
					judge = Judge_SV(RCRead, ReadLength, kmer, start, end);
				SV_flag = SV_flag * judge;
				if( SV_flag == 1 )
					break;
			}
			final_read.pop();
		}

		while( !final_read.empty() )
			final_read.pop();

		if( SV_flag == 1 )
			Answer << cycle_time + 1 << "\tN" << endl;
		else
			Answer << cycle_time + 1 << "\tT" << endl;
		cycle_time++;

		//break;
		//if(cycle_time==1000)break;
	}
	kseq_destroy(Seq);
	gzclose(fp);
	Answer.close();

	//long end_time = getCurrentTime();
	//cout << "Times:" << end_time - start_time << "msec." << endl;

	return 0;
}


int Local_aln::Judge_SV(char *Seq, uint32_t ReadLength, uint32_t kmer, uint32_t GU, uint32_t GD)
{
	// 0 stands for sv, 1 stands for normal.

	uint64_t SR_length = ReadLength + 1 - kmer;
	uint32_t pre = 0;
	uint32_t mask = 0xffffffff>>(32 - (kmer<<1));
	//produce all kmer value;
	pre = transfer(Seq, kmer, 0);
	tempArray[0] = pre; 
	tempArrayNext[0] = pre;

	// if you want to mask uncomment the next sentence
	for (uint32_t i=1; i< SR_length; i++ )
	{
		pre = ((pre<<2)&mask)|trans[Seq[i + kmer - 1 ]];
		tempArray[i] = pre;
		tempArrayNext[i] = pre;
	}

	qsort(tempArray, SR_length, sizeof(uint32_t), compare);

	uint32_t flag;
	for( int i = 0 ; i < SR_length ; i++ )
	{
		flag = tempArray[i];
		readKmerhash[flag].Kpointer1++;
	}

	flag = -1;
	uint64_t sum = 0;
	for( int i = 0 ; i < SR_length ; i++ )
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

	for( int i = 0 ; i < SR_length ; i++ )
	{
		flag = tempArrayNext[i];
		uint32_t add = readKmerhash[flag].Kpointer2;
		//cout << SR_length << "\t" << flag << "\t" << add << endl;
 		readSRhash[add] = i;
		readKmerhash[flag].Kpointer2++;
	}

	/*ofstream out;
	out.open("../test");*/

	priority_queue<Tuple>	Tuple_Queue;	
	pre = transfer(genome, kmer, GU);
	pre = ((pre>>2)&mask);	
	for( int i = 0 ; i < GD - GU ; i++ )
	{
		Tuple temp_node;
		temp_node.ref_begin = i + GU;
		//temp_node.tuplelength = kmer;
		//uint32_t l2r = transfer(genome, kmer, i + GU);
		pre = ((pre<<2)&mask)|trans[genome[i + GU + kmer - 1 ]];
		uint32_t how_many_node = readKmerhash[pre].Kpointer2 - readKmerhash[pre].Kpointer1;
		for( int k = 0 ; k < how_many_node ; k++ )
		{
			temp_node.read_begin = readSRhash[k + readKmerhash[pre].Kpointer1];
			//maybe change the shorter substr
			Tuple_Queue.push(temp_node);
			/*out << temp_node.ref_begin << "\t";
			out << temp_node.read_begin << "\t";
			out << temp_node.tuplelength << endl;*/
		}
	}

	//out.close();

	for( int i = 0 ; i < SR_length ; i++ )
	{
		readKmerhash[tempArray[i]].Kpointer1 = 0;
		readKmerhash[tempArray[i]].Kpointer2 = 0;
	}

	//cout << Tuple_Queue.size() << endl;

	if( Tuple_Queue.top().read_begin > 250 )
	{
		while( !Tuple_Queue.empty() )
			Tuple_Queue.pop();
		return 0;
	}
	else
	{
		Tuple temp_node;
		temp_node.ref_begin = Tuple_Queue.top().ref_begin;
		temp_node.read_begin = Tuple_Queue.top().read_begin;
		//temp_node.tuplelength = Tuple_Queue.top().tuplelength;
		Tuple_Queue.pop();
		int local_ans = 0;
			/*out << temp_node.ref_begin << "\t";
			out << temp_node.read_begin << "\n";*/
			//out << temp_node.tuplelength << endl;
		while( !Tuple_Queue.empty() )
		{
			local_ans = Tuple_link( temp_node, Tuple_Queue.top());
			if( local_ans == 1 )
			{
				temp_node.ref_begin = Tuple_Queue.top().ref_begin;
				temp_node.read_begin = Tuple_Queue.top().read_begin;
				//temp_node.tuplelength = Tuple_Queue.top().tuplelength;
				/*out << temp_node.ref_begin << "\t";
				out << temp_node.read_begin << "\n";*/
				//out << temp_node.tuplelength << endl;
			}
			Tuple_Queue.pop();
		}

		if( 250 < ReadLength - temp_node.read_begin )
		{
			return 0;
		}
	}
	return 1;
}

int Local_aln::Tuple_link( Tuple t1, Tuple t2)
{
	int	distance_read_s;
	int	distance_ref_s;
	double	link_chioce;

	distance_read_s = t1.read_begin - t2.read_begin;
	distance_ref_s = t1.ref_begin - t2.ref_begin;

	link_chioce = distance_read_s - TransParameter * distance_ref_s;
	//Maybe useful
	//distance_read_s = distance_read_s + t1.tuplelength;

	/*cout << t1.read_begin << "\t" << t2.read_begin << "\t";
	cout << t1.ref_begin << "\t" << t2.ref_begin << "\t";
	cout << distance_read_s << "\t" << link_chioce << endl;*/

	if( distance_read_s <= DISTANCE_SCORE && link_chioce <= 50 && link_chioce >= -50 )
		return 1;
	else
		return 0;
}

/*char *get_path(char *path, const char *str2)
{
	const char *str = path;
	const size_t len = strlen(str)+strlen(str2);
	char *n_str = new char[len+1];
	strcpy(n_str , str);
	strcat(n_str ,str2);
	return n_str;
}*/

void Local_aln::Files_open(char* path, uint32_t kmer)
{
	fstream LengthFile;
	char *LF_path = NULL;
	LF_path = get_path(path, "LengthFile");
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
	/*fscanf(LengthFile, "%u\n", Hash2sum);
	fscanf(LengthFile, "%u\n", len_genome);*/
	LengthFile.close();

	genome = new char[len_genome];
	if (NULL == genome) {
		fprintf(stderr,"Fail to alloc genome, now exit");
		exit(1);
	} 
	FILE * file_ref;
	char *ref_path = NULL;
	ref_path = get_path(path, "Reference");
	file_ref = fopen(ref_path, "rb");
	delete[] ref_path;
	fread( genome, sizeof(char), len_genome, file_ref );
	fclose(file_ref);
	if (NULL == genome) {
		fprintf(stderr,"Fail to load Genome reference, now exit");
		exit(1);
	} 

	int movement = 2 * kmer;
	uint64_t Kmer_code_length = 1 << movement;

	kmerHash = new uint32_t[Kmer_code_length];
	if (NULL == kmerHash) {
		fprintf(stderr,"Fail to alloc kmerHash, now exit");
		exit(1);
	} 
	FILE * file_hash_1;
	char *hash1_path = NULL;
	hash1_path = get_path(path, "Hash_1");
	file_hash_1 = fopen(hash1_path, "rb");
	delete[] hash1_path;
	fread( kmerHash, sizeof(uint32_t), Kmer_code_length, file_hash_1 );
	fclose(file_hash_1);
	if (NULL == kmerHash) {
		fprintf(stderr,"Fail to load kmerHash, now exit");
		exit(1);
	}

	SR_HASH = new uint32_t[Hash2sum];
	if (NULL == SR_HASH) {
		fprintf(stderr,"Fail to alloc SR_HASH, now exit");
		exit(1);
	} 
	FILE * file_hash_2;
	char *hash2_path = NULL;
	hash2_path = get_path(path, "Hash_2");
	file_hash_2 = fopen(hash2_path, "rb");
	delete[] hash2_path;
	fread( SR_HASH, sizeof(uint32_t), Hash2sum, file_hash_2 );
	fclose(file_hash_2); 
	if (NULL == SR_HASH) {
		fprintf(stderr,"Fail to load SR_HASH, now exit");
		exit(1);
	}

	readKmerhash = new ReadKmerHash[Kmer_code_length];
	if (NULL == readKmerhash) {
		fprintf(stderr,"Fail to alloc readKmerhash, now exit");
		exit(1);
	} 
	for( int i = 0 ; i < Kmer_code_length ; i++ )
		readKmerhash[i].Kpointer1 = 0;

	readSRhash = new uint32_t[READ_MAX_LENGTH];
	if (NULL == readSRhash) {
		fprintf(stderr,"Fail to alloc readSRhash, now exit");
		exit(1);
	} 

	tempArray = new uint32_t[READ_MAX_LENGTH];
	if (NULL == tempArray) {
		fprintf(stderr,"Fail to alloc tempArray, now exit");
		exit(1);
	} 

	tempArrayNext = new uint32_t[READ_MAX_LENGTH];
	if (NULL == tempArrayNext) {
		fprintf(stderr,"Fail to alloc tempArrayNext, now exit");
		exit(1);
	}

	Result_way = new uint32_t*[SEED_NUM];
	for( int i = 0 ; i < SEED_NUM ; i++ )
		Result_way[i] = new uint32_t[UpBound];

	RCRead = new char[READ_MAX_LENGTH];
	if (NULL == RCRead) {
		fprintf(stderr,"Fail to alloc RCRead, now exit");
		exit(1);
	}
}


void Local_aln::Cleaning()
{
	if (NULL != genome)
		delete[] genome;

	if (NULL != kmerHash)
		delete[] kmerHash;

	if (NULL != SR_HASH)
		delete[] SR_HASH;

	if (NULL != readKmerhash)
		delete[] readKmerhash;

	if (NULL != readSRhash)
		delete[] readSRhash;

	if (NULL != tempArray)
		delete[] tempArray;

	if (NULL != tempArrayNext)
		delete[] tempArrayNext;

	if (NULL != RCRead)
		delete[] RCRead;

	if (NULL != Result_way) {
		for (int i=0;i<SEED_NUM;++i) { delete[] Result_way[i]; }
		delete[] Result_way;
	}
}

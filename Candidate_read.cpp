#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <queue>
#include <algorithm>
#include <sys/time.h> 
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"

using namespace std;

#define readLINE 30000
#define K_MER 13
#define SEED_NUM 1000
#define UpBound 1000
#define TopChoice 5
#define MAX_ARRAY 30
#define Hpart 32

KSEQ_INIT(gzFile, gzread)

typedef struct 
{
	uint32_t SR;
	//uint32_t freq;
}SRdata;

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
	//uint32_t freq;
	bool operator<( const way_node & a )const
	{
		return win > a.win;
	}
}Node;

typedef struct cal
{
	uint32_t win;
	uint32_t freq;
	int flag;
	int array[MAX_ARRAY];
	bool operator<( const cal & a )const
	{
		return ( flag == a.flag ) ? freq < a.freq : flag > a.flag;
		//return (freq == a.freq) ? flag > a.flag : freq < a.freq;
		/*if( freq == a.freq )
		{
			return flag > a.flag;
		}
		else
			return freq < a.freq;*/
	}
}Freq;

uint8_t trans[128]= {
		4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
		4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
		4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
		4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
		4,0,4,1,4,4,4,2,4,4,4,4,4,4,2,4,
		4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,
		4,0,4,1,4,4,4,2,4,4,4,4,4,4,2,4,
		4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4
	};

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

uint32_t transfer(char *genome, uint32_t len_sed, uint32_t start)
{

	uint32_t l2r = 0;
	uint32_t temp;
	for(uint32_t i = 0;i<len_sed;i++) {
		temp = trans[genome[start+i]];
		l2r = l2r << 2;
		l2r = l2r | temp;
	}
	return l2r;
}

/*char *RevComRead(char *read, uint32_t ReadLength)
{
	int i;
	char rcread[ReadLength];
	for( i = 0 ; i < ReadLength ; i++ )
		rcread[i] = rev[read[ReadLength - 1 - i]];
	rcread[ReadLength] = '\0';
	char* str = rcread;
	return str;
}*/

int RevComRead(char *rcread,char *read,int len_read)
{
	int i;
	for(i=0; i<len_read; i++)
		rcread[i] = rev[read[len_read - 1 - i]];
	rcread[i] = '\0';
	return 0;
}

long getCurrentTime()    
{
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return tv.tv_sec * 1000 + tv.tv_usec / 1000;    
} 

int Get_Max_range( uint32_t *array )
{
	int	i = 0;
	int	score = 0;
	int 	MaxScore = 0;
	for( i = 0; i < 1000; i++ )
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
	if( array[i] == 0 )
		if( MaxScore < score )
			MaxScore = score;
	return MaxScore;
}

int main(int argc, char const *argv[])
{
	//create the hsah table

	long start_time = getCurrentTime();

	int movement = 2 * K_MER;
	uint64_t Kmer_code_length = 1 << movement;
	uint32_t * kmerHash = (uint32_t*)malloc(sizeof(uint32_t)*Kmer_code_length);
	if( NULL == kmerHash ) {
		cout<<"Fail to allocate new space to hash table"<<endl;
		exit(1);
	}
	FILE* inHash_1;
	inHash_1 = fopen("../Hash_1","rb");
	fread(kmerHash,sizeof(uint32_t),Kmer_code_length,inHash_1);
	fclose(inHash_1);

	uint64_t Hash2sum = 225371321;//223830061;
	/*
	file NOTE
	maybe read from the file !
	*/
	SRdata * SR_HASH = (SRdata*)malloc(Hash2sum*sizeof(SRdata));
	if( NULL == SR_HASH) {
		cout<<"Fail to allocate new space to SR_HASH"<<endl;
		exit(1);
	}
	FILE* inHash_2;
	inHash_2 = fopen("../Hash_2","rb");
	fread(SR_HASH,sizeof(SRdata),Hash2sum,inHash_2);
	fclose(inHash_2);

	ofstream fout;
	fout.open("../LowHITread.newtry.2.Version_kmer..");

	gzFile fp;
	kseq_t *Seq;

	fp = gzopen("../_0001.fastq", "r");
	Seq = kseq_init(fp);

	int		cycle_time = 0;
	uint32_t	ReadLength;
	char		*PositiveRead = NULL;
	//char		*RCRead = NULL;
	uint32_t	seed_bet;
	char *RCRead = new char[100000];

	while (kseq_read(Seq) >= 0) {

		/*if(cycle_time < 3)
		{
			cycle_time++;
			continue;
		}*/

		ReadLength = Seq->seq.l;
		/*
			This model need to be change if given a new reads data.
			Due to the distribution is not like "+/-/+/-...".
			The model should be finding the max score of the "+" and "-" read.
		*/
		if( cycle_time % 2 != 0 )
		{
			RevComRead(RCRead, Seq->seq.s, ReadLength);
			PositiveRead = RCRead;
		}
		else
			PositiveRead = Seq->seq.s;
		
		seed_bet = (ReadLength - K_MER) / SEED_NUM;
		uint32_t seed_found = 0;
		uint32_t window_sum = ReadLength / 1.1 / Hpart;
		uint32_t window_hash[window_sum];
		for( int k = 0 ; k < window_sum ; k++ )
			window_hash[k] = 0;

		uint32_t Result_way[SEED_NUM][UpBound];
		uint32_t Result_way_freq[SEED_NUM][UpBound];
		uint32_t Each_way_count[SEED_NUM];

		
		for( int i = 0 ; i < SEED_NUM ; i++ )
		{
			uint32_t start = i * seed_bet;
			uint32_t l2r = transfer(PositiveRead,K_MER,start);
			/*for( int fuck = 0; fuck <K_MER ;fuck++)
			{
				cout << PositiveRead[start+fuck];
			}
			cout << endl;
			cout << i << "\t" << l2r << endl;
			cout << kmerHash[l2r] << endl;*/
			
			if( kmerHash[l2r+1] <  kmerHash[l2r])
				Each_way_count[i] = 0;
			else if( (kmerHash[l2r+1] - kmerHash[l2r]) > UpBound )
				Each_way_count[i] = 0;
			else
			{
				int bias = i * seed_bet / 1.1;
				bias = bias / Hpart;
				Each_way_count[i] = kmerHash[l2r+1] - kmerHash[l2r];
				for( int j = 0 ; j < Each_way_count[i] ; j++ )
				{
					Result_way[i][j] = SR_HASH[j+kmerHash[l2r]].SR;
					if( Result_way[i][j] >= bias )
						Result_way[i][j] = Result_way[i][j] - bias;
					else	Result_way[i][j] = 0;
				}
				//break;
			}
		}

		//cout << ReadLength << endl;

		priority_queue<Node>	windowQueue;
		priority_queue<Node>	KmerQueue;
		//priority_queue<Freq>	windowCount;

		uint32_t judge_element = -1, count_element = 0;
		CF region, MaxRegion;
		MaxRegion.pf = 0;
		region.ps = -1, region.pe = -1, region.pf = 0;
		uint32_t con_array [1000] = {0};
		//uint32_t Max_con_array[1000];

		uint32_t Each_way_flag[SEED_NUM];
		for( int i = 0 ; i < SEED_NUM ; i++ )
		{
			Each_way_flag[i] = 0;
			if( Each_way_count[i] != 0 )
			{
				Node node;
				node.ID = i;
				node.win = Result_way[i][0];
				//node.freq = Result_way_freq[i][0];
				windowQueue.push(node);
				KmerQueue.push(node);
				Each_way_flag[i]++;
			}
		}

		while( !windowQueue.empty())
		{
			// if( judge_element != windowQueue.top().win )
			// {
			// 	Freq ele;
			// 	ele.win = judge_element;
			// 	ele.freq = count_element;
			// 	int this_score = 0;
			// 	for( int k = 0 ; k < window_sum ; k++ )
			// 		this_score += window_hash[k];
			// 	if( this_score == window_sum )
			// 		ele.flag = 0;
			// 	else
			// 	{
			// 		ele.flag = window_sum - this_score;
			// 		int tempL = 0;
			// 		for( int k = 0 ; k < window_sum ; k++ )
			// 		{
			// 			if( window_hash[k] == 0 )
			// 			{
			// 				ele.array[tempL] = k;
			// 				tempL++;
			// 			}
			// 		}
			// 	}	
			// 	windowCount.push(ele);

			// 	for( int k = 0 ; k < window_sum ; k++ )
			// 		window_hash[k] = 0;
			// 	judge_element = windowQueue.top().win;
			// 	//count_element = windowQueue.top().freq;
			// 	count_element = 1;
			// 	int which_win = windowQueue.top().ID * seed_bet / 1.1 / Hpart;
			// 	window_hash[which_win] = 1;
			// }
			// else
			// {
			// 	//count_element += windowQueue.top().freq;
			// 	count_element++;
			// 	int which_win = windowQueue.top().ID * seed_bet / 1.1 / Hpart;
			// 	window_hash[which_win] = 1;
			// }

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
				if( region.pf > MaxRegion.pf )
				{
					MaxRegion.pf = region.pf;
					MaxRegion.ps = region.ps;
					MaxRegion.pe = region.pe;
					/*for( int lm = 0 ; lm < 1000 ; lm++ )
						Max_con_array[lm] = con_array[lm];*/
				}
				region.ps = windowQueue.top().win;
				region.pe = windowQueue.top().win;
				region.pf = 1;
				/*for( int lm = 0 ; lm < 1000 ; lm++ )
				{
					con_array[lm] = 0;
				}
				con_array[node.ID]++;*/
			}


			windowQueue.pop();
			if( Each_way_flag[node.ID] < Each_way_count[node.ID] )
			{
				node.win = Result_way[node.ID][Each_way_flag[node.ID]];
				Each_way_flag[node.ID]++;
				if( node.win != 0 )
				{
					windowQueue.push(node);
					KmerQueue.push(node);
				}
			}
		}

		while( !KmerQueue.empty() )
		{
			if( KmerQueue.top().win >= MaxRegion.ps && KmerQueue.top().win <= MaxRegion.pe )
			{
				con_array[KmerQueue.top().ID]++;
			}
			KmerQueue.pop();
		}

		int MaxScore = Get_Max_range(con_array);

		cycle_time++;
		fout <<"Read:" << cycle_time << "\t" << ReadLength << "\t";
		//fout << MaxScore * seed_bet / Hpart;
		fout << MaxScore;
		fout << endl;
		fout << MaxRegion.ps << "\t" << MaxRegion.pe << "\t" << MaxRegion.pf << endl;

		if( cycle_time == 1000 ) break;
		
	}
	kseq_destroy(Seq);
	gzclose(fp);
	fout.close();
	free(kmerHash);
	free(SR_HASH);
	delete[] RCRead;

	long end_time = getCurrentTime();
	cout << "Times:" << end_time - start_time << "msec." << endl;

	return 0;
}
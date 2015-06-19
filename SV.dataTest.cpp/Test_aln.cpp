#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdint.h>
#include <queue>
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

long getCurrentTime()    
{
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return tv.tv_sec * 1000 + tv.tv_usec / 1000;    
} 

typedef struct
{
	uint32_t winS;
	uint32_t winE;
}Read;

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

int RevComRead(char *rcread,char *read,int len_read)
{
	int i;
	for(i=0; i<len_read; i++)
		rcread[i] = rev[read[len_read - 1 - i]];
	rcread[i] = '\0';
	return 0;
}

int compare(const void *a, const void *b)
{
	uint32_t *pa = (uint32_t*)a;
	uint32_t *pb = (uint32_t*)b;
	return *pa > *pb;  //从小到大排序
}

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
		//return cover_score < a.cover_score;
		return ( read_begin == a.read_begin ) ? ref_begin > a.ref_begin : read_begin > a.read_begin;
	}
}Tuple;

int Tuple_link( Tuple t1, Tuple t2)
{
	int	distance_read_s;
	int	distance_ref_s;
	double	link_chioce;

	distance_read_s = t1.read_begin - t2.read_begin;
	distance_ref_s = t1.ref_begin - t2.ref_begin;

	link_chioce = distance_read_s - 1.1 * distance_ref_s;
	//Maybe useful
	//distance_read_s = distance_read_s + t1.tuplelength;

	// cout << t1.read_begin << "\t" << t2.read_begin << "\t";
	// cout << t1.ref_begin << "\t" << t2.ref_begin << "\t";
	// cout << distance_read_s << "\t" << link_chioce << endl;

	if( distance_read_s <= 750 && link_chioce <= 50 && link_chioce >= -50 )
	{
		return 1;
	}
	else
		return 0;
}

int main()
{
	long start_time = getCurrentTime();

	int movement = 2 * K_MER;
	uint64_t Kmer_code_length = 1 << movement;

	ReadKmerHash	*readKmerhash = new ReadKmerHash[Kmer_code_length];
	uint32_t *readSRhash = new uint32_t[100000];
	uint32_t *tempArray = new uint32_t[100000];
	uint32_t *tempArrayNext = new uint32_t[100000];

	char *FLAG = new char[622656];
	Read *read = new Read[622656];
	ifstream flagIN;
	flagIN.open("../Nextctrl_data");
	char temp[readLINE];
	string str;
	for( int i = 0; i < 622656; i++ )
	{
		flagIN.getline(temp, readLINE, '\n');
		FLAG[i] = temp[0];
	}
	flagIN.close();

	ifstream StartIN;
	StartIN.open("../LowHITread.newtry.2.Version_bp");
	for( int i = 0; i < 622656; i++ )
	{
		StartIN.getline(temp, readLINE, '\n');
		StartIN.getline(temp, readLINE, '\t');
		str = temp;
		read[i].winS = atoi(str.c_str());
		StartIN.getline(temp, readLINE, '\t');
		str = temp;
		read[i].winE = atoi(str.c_str());
		StartIN.getline(temp, readLINE, '\n');
		//break;
	}
	StartIN.close();

	char *genome = new char[249250621];
	FILE * file_ref;
	file_ref = fopen("/home/tjiang/Reference", "rb");
	fread( genome, sizeof(char), 249250621, file_ref );
	fclose(file_ref);

/*	cout << FLAG[9] << endl;
	cout << read[0].winS << " " <<  read[0].winE << endl;
	cout << read[1].winS << " " <<  read[1].winE << endl;*/

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

	for( int i = 0 ; i < Kmer_code_length ; i++ )
		readKmerhash[i].Kpointer1 = 0;

	ofstream out;
	out.open("../local_index_results.150");

	while (kseq_read(Seq) >= 0) 
	{
		/*if(cycle_time < 288)
		{
			cycle_time++;
			continue;
		}*/

		if( FLAG[cycle_time] != 'T' )
		{
			ReadLength = Seq->seq.l;
			if( cycle_time % 2 != 0 )
			{
				RevComRead(RCRead, Seq->seq.s, ReadLength);
				PositiveRead = RCRead;
			}
			else
				PositiveRead = Seq->seq.s;


			uint64_t SR_length = ReadLength + 1 - K_MER;

			uint32_t pre = 0;
			uint32_t mask = 0xffffffff>>(32 - (K_MER<<1));

			//produce all kmer value;
			pre = transfer(PositiveRead, K_MER, 0);
			tempArray[0] = pre; 
			tempArrayNext[0] = pre;

			// if you want to mask uncomment the next sentence
			for (uint32_t i=1; i< SR_length; i++ )
			{
				pre = ((pre<<2)&mask)|trans[PositiveRead[i + K_MER - 1 ]];
				//prseq(i+seq, 13,false);
				//cout<<"\t"<<i<<endl;
				tempArray[i] = pre;
				tempArrayNext[i] = pre;
				// seq_num[i] = i;
				// int count_kmer

				//order[i] = i;
			}

			
			/*for( int i = 0 ; i < SR_length ; i++ )
			{
				tempArray[i] = transfer(PositiveRead, K_MER, i);
				tempArrayNext[i] = tempArray[i];
			}*/

			qsort(tempArray, SR_length, sizeof(uint32_t), compare);

			uint32_t flag = -1;
			for( int i = 0 ; i < SR_length ; i++ )
			{
				//if( flag != tempArray[i] )
				//{
					flag = tempArray[i];
					readKmerhash[flag].Kpointer1++;
					//readKmerhash[flag].Kpointer2++;
					//kmerHash[flag]++;
				//}
			}

			/*uint64_t sum = 0;
			for( int i = 0 ; i < Kmer_code_length ; i++ )
			{
				uint32_t temp = readKmerhash[i].Kpointer1;
				readKmerhash[i].Kpointer1 = sum;
				readKmerhash[i].Kpointer2 = sum;
				sum = sum + temp;
			}
*/
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

			/*uint64_t sum = readKmerhash[tempArray[0]].Kpointer1;
			for( int i = 0 ; i <= tempArray[0] ; i++ )
				readKmerhash[i].Kpointer1 = 0;*/
			//readKmerhash[tempArray[0]].Kpointer1 = 0;

			/*flag = -1;
			for( int i = 1 ; i < SR_length ; i++ )
			{
				flag = tempArray[i];
				if( flag != tempArray[i] )
				{
					uint32_t temp = readKmerhash[flag].Kpointer1;
					for( int k = tempArray[i-1]+1 ; k < flag ; k++ )
					{
						readKmerhash[k].Kpointer1 = 0;
					}
					readKmerhash[flag].Kpointer1 = sum;
					sum+=temp;
				}
				
			}*/

			/*for( int i = 0 ; i < Kmer_code_length ; i++ )
				readKmerhash[i].Kpointer2 = readKmerhash[i].Kpointer1;*/

			//flag = -1;
			for( int i = 0 ; i < SR_length ; i++ )
			{
				//if( flag != tempArrayNext[i] )
				//{
					flag = tempArrayNext[i];
					uint32_t add = readKmerhash[flag].Kpointer2;
					//cout << SR_length << "\t" << flag << "\t" << add << endl;
 					readSRhash[add] = i;
					readKmerhash[flag].Kpointer2++;
				//}
			}

			uint32_t GU = read[cycle_time].winS * Hpart - 1000;
			uint32_t GD = read[cycle_time].winE *Hpart + ReadLength + 1000;

			priority_queue<Tuple>	Tuple_Queue;		
			for( int i = 0 ; i < GD - GU ; i++ )
			{
				Tuple temp_node;
				temp_node.ref_begin = i + GU;
				temp_node.tuplelength = K_MER;
				uint32_t l2r = transfer(genome, K_MER, i + GU);
				uint32_t how_many_node = readKmerhash[l2r].Kpointer2 - readKmerhash[l2r].Kpointer1;
				for( int k = 0 ; k < how_many_node ; k++ )
				{
					temp_node.read_begin = readSRhash[k + readKmerhash[l2r].Kpointer1];
					Tuple_Queue.push(temp_node);
				}

			}

			for( int i = 0 ; i < SR_length ; i++ )
			{
				readKmerhash[tempArray[i]].Kpointer1 = 0;
				readKmerhash[tempArray[i]].Kpointer2 = 0;
			}

			if( Tuple_Queue.top().read_begin > 150 )
			{
				while( !Tuple_Queue.empty() )
					Tuple_Queue.pop();
				/*delete[] readKmerhash;
				delete[] readSRhash;*/
				FLAG[cycle_time] = 'T';
			}
			else
			{
				Tuple temp_node;
				temp_node.ref_begin = Tuple_Queue.top().ref_begin;
				temp_node.read_begin = Tuple_Queue.top().read_begin;
				temp_node.tuplelength = Tuple_Queue.top().tuplelength;
				Tuple_Queue.pop();
				int local_ans = 0;
				while( !Tuple_Queue.empty() )
				{
					local_ans = Tuple_link( temp_node, Tuple_Queue.top());
					if( local_ans == 1 )
					{
						temp_node.ref_begin = Tuple_Queue.top().ref_begin;
						temp_node.read_begin = Tuple_Queue.top().read_begin;
						temp_node.tuplelength = Tuple_Queue.top().tuplelength;
					}
					Tuple_Queue.pop();
				}

				if( 150 < ReadLength - temp_node.read_begin )
				{
					FLAG[cycle_time] = 'T';
				}
			}

			

			/*ofstream out;
			out.open("../local_index_results");
			while( !Tuple_Queue.empty() )
			{
				out << Tuple_Queue.top().read_begin << "\t";
				out << Tuple_Queue.top().ref_begin << "\t";
				out << Tuple_Queue.top().tuplelength << "\n";
				Tuple_Queue.pop();
			}
			out.close();*/


			// cout << ReadLength << endl;
			//break;
		}
		out << cycle_time+1 << "\t"<< FLAG[cycle_time] << endl;
		cycle_time ++;
		//break;
		//if( cycle_time == 1000 ) break;
	}

	kseq_destroy(Seq);
	gzclose(fp);

	
	/*for( int i = 0; i < 622656; i++ )
	{
		out << FLAG[i] << "\n";
	}*/
	out.close();

	delete[] FLAG;
	delete[] read;
	delete[] readKmerhash;
	delete[] readSRhash;
	delete[] tempArrayNext;
	delete[] tempArray;

	long end_time = getCurrentTime();
	cout << "Times:" << end_time - start_time << "msec." << endl;
	return 0;
}
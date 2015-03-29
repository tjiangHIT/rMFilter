#include <iostream>
#include <stdlib.h>
#include "readfl.h"
#include <stdint.h>
#include <stdio.h>
#include <fstream>
#include <string>

#define K_MER 13
#define Hpart 512
#define Fpart 1024

using namespace std;

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

uint32_t transfer(char *genome,uint32_t len_sed,uint32_t start)
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

int compare(const void *a, const void *b)
{
	uint32_t *pa = (uint32_t*)a;
	uint32_t *pb = (uint32_t*)b;
	return *pa > *pb;  //从小到大排序
}

typedef struct 
{
	uint32_t SR;
	uint32_t freq;
}SRdata;

int main()
{
	char *path = "../hg19_chr1.fa";
	unsigned len_genome;
	read_file inseq;
	char *genome = inseq.read_ref( path, &len_genome );



	//uint64_t Kmer_code_length = 1 << ((2*K_MER)-1);
	int movement = 2*K_MER;
	uint64_t Kmer_code_length = 1 << movement;
	//cout << a << endl;
	uint32_t * kmerHash = (uint32_t*)malloc(sizeof(uint32_t)*Kmer_code_length);
	if( NULL == kmerHash) {
		cout<<"Fail to allocate new space to hash table"<<endl;
		exit(1);
	}

	for( int i = 0 ; i < Kmer_code_length ; i++ )
		kmerHash[i] = 0;

	/*cout << transfer( genome, 13, 0) << endl;
	string str = genome;
	cout << str.substr(0,13) << endl;
	cout << Kmer_code_length << endl;*/


	uint64_t SR_length = len_genome / Hpart;
	//cout << SR_length << endl;
	for( int i = 0 ; i < SR_length ; i++ )
//		for( int i = 0 ; i < 50 ; i++ )
	{
		uint32_t SQ_start = i * Hpart;
		uint32_t *tempArray = (uint32_t*)malloc((Hpart - K_MER + 1)*sizeof(uint32_t));
		if( NULL == tempArray) {
			cout<<"Fail to allocate new space to tempArray"<<endl;
			exit(1);
		}
		for( int j = 0 ; j < Hpart - K_MER + 1 ; j++ )
		{
			uint32_t l2r = transfer(genome,K_MER,j+SQ_start);
			tempArray[j] = l2r;
			//Qsort
			//delete duplication
			//caculation
			//cout << j << " " << tempArray[j] << endl;
		}
		qsort(tempArray, Hpart - K_MER + 1, sizeof(uint32_t), compare);
		uint32_t flag = -1;
		for( int j = 0 ; j < Hpart - K_MER + 1 ; j++ )
		{
			if( flag != tempArray[j] )
			{
				flag = tempArray[j];
				kmerHash[flag]++;
				//cout << i << " " << j << " " << tempArray[j] << " " << flag << " " << kmerHash[flag] << endl;
			}
			//cout << j << " " << tempArray[j] << endl;
		}
		//cout << i << endl;
		//break;
	}
	uint64_t Hash2sum = 0;
	for( int i = 0 ; i < Kmer_code_length ; i++ )
	{
		uint32_t temp = kmerHash[i];
		kmerHash[i] = Hash2sum;
		Hash2sum = Hash2sum + temp;
		//kmerHash[i] = Hash2sum - temp;
	}
	cout << Hash2sum << endl;
	//cout << kmerHash[SR_length-1] << endl;
	//cout << kmerHash[1377648] << endl;

	SRdata * SR_HASH = (SRdata*)malloc(Hash2sum*sizeof(SRdata));
	if( NULL == SR_HASH) {
		cout<<"Fail to allocate new space to SR_HASH"<<endl;
		exit(1);
	}

	/*ofstream out;
	out.open("../fuck.txt");
	for( int i = 0 ; i < SR_length ; i++ )
		out << kmerHash[i] << endl;*/
	FILE * outfile;
	outfile = fopen("../hash_1", "wb");
	fwrite( kmerHash, sizeof(uint32_t), Kmer_code_length,outfile );
	fclose(outfile);


	for( int i = 0 ; i < SR_length ; i++ )
//		for( int i = 19 ; i < 21 ; i++ )
	{
		uint32_t SQ_start = i * Hpart;
		uint32_t *tempArray = (uint32_t*)malloc((Hpart - K_MER + 1)*sizeof(uint32_t));
		if( NULL == tempArray) {
			cout<<"Fail to allocate new space to tempArray"<<endl;
			exit(1);
		}
		for( int j = 0 ; j < Hpart - K_MER + 1 ; j++ )
		{
			uint32_t l2r = transfer(genome,K_MER,j+SQ_start);
			tempArray[j] = l2r;
		}
		qsort(tempArray, Hpart - K_MER + 1, sizeof(uint32_t), compare);
		uint32_t flag = tempArray[0];
		uint32_t flag_n = 0;
		for( int j = 0 ; j < Hpart - K_MER + 1 ; j++ )
		{
			if( flag == tempArray[j] )
			{
				flag_n++;
			}
			else
			{
				uint32_t add = kmerHash[flag];
				SR_HASH[add].SR = i;
				SR_HASH[add].freq = flag_n;
				//cout << flag << " " << kmerHash[flag] << " " << SR_HASH[add].SR << " " << SR_HASH[add].freq << endl;
				kmerHash[flag]++;
				flag = tempArray[j];
				flag_n = 1;
			}
		}
		if( tempArray[0] == tempArray[Hpart - K_MER] )
		{
			uint32_t add = kmerHash[flag];
			SR_HASH[add].SR = i;
			SR_HASH[add].freq = flag_n;
			//cout << flag << " " << kmerHash[flag] << " " << SR_HASH[add].SR << " " << SR_HASH[add].freq << " " << kmerHash[flag+1] << endl;
			kmerHash[flag]++;
		}
		//cout << i << endl;
		//break;
	}
	/*for( int i = 0 ; i < 100 ; i++ )
		cout << SR_HASH[i].SR << " " << SR_HASH[i].freq << endl;*/
	FILE * outfile2;
	outfile2 = fopen("../hash_2", "wb");
	fwrite( SR_HASH, sizeof(SRdata), Hash2sum,outfile2 );
	fclose(outfile2);
	return 0;
}
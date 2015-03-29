#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <map>
#include <vector>
#include <algorithm>

using namespace std;

#define readLINE 30000
#define K_MER 13
#define Hpart 512
#define SEED_NUM 1000
typedef struct 
{
	uint32_t SR;
	uint32_t freq;
}SRdata;

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

uint32_t transfer(string genome,uint32_t len_sed,uint32_t start)
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

string RevComRead(string read)
{
	int i;
	char rcread[read.length()];
	for( i = 0 ; i < read.length() ; i++ )
		rcread[i] = rev[read[read.length() - 1 - i]];
	rcread[read.length()] = '\0';
	string str = rcread;
	return str;
}

int compare(const void *a, const void *b)
{
	SRdata *pa = (SRdata*)a;
	SRdata *pb = (SRdata*)b;
	return pa->freq < pb->freq;  //从小到大排序
}

int cmp(const void *a, const void *b)
{
	uint32_t *pa = (uint32_t*)a;
	uint32_t *pb = (uint32_t*)b;
	return pa < pb;  //从小到大排序
}

int main(int argc, char const *argv[])
{
	//create the hsah table

	int movement = 2 * K_MER;
	uint64_t Kmer_code_length = 1 << movement;
	uint32_t * kmerHash = (uint32_t*)malloc(sizeof(uint32_t)*Kmer_code_length);
	if( NULL == kmerHash ) {
		cout<<"Fail to allocate new space to hash table"<<endl;
		exit(1);
	}
	FILE* inHash_1;
	inHash_1 = fopen("../hash_1","rb");
	fread(kmerHash,sizeof(uint32_t),Kmer_code_length,inHash_1);
	fclose(inHash_1);

	uint64_t Hash2sum = 218595088; // read from file
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
	inHash_2 = fopen("../hash_2","rb");
	fread(SR_HASH,sizeof(SRdata),Hash2sum,inHash_2);
	fclose(inHash_2);

	ofstream fout;
	fout.open("../NEW_results.version(1.1,100)");

	ifstream filein;
	filein.open("../_0001.fastq");

	ifstream filein_sup;
	filein_sup.open("../read_ori");

	ofstream low_cover_read;
	low_cover_read.open("../low_cover_read.version(1.1,100)");

	char temp[readLINE];
	int fuck = 0;
	int seed_found_L = 0 ;
	uint64_t sum_cover = 0;
	while( !filein.eof() )
	{
		string strtemp;
		filein_sup.getline(temp,readLINE,'\t');
		filein_sup.getline(temp,readLINE,'\n');
		strtemp = temp;
		uint32_t Each_Read_start;
		Each_Read_start = atoi(strtemp.c_str());

		string seq;
		filein.getline(temp,readLINE,'\n');
		filein.getline(temp,readLINE,'\n');
		seq = temp;
		if( seq == "" ) break;
		if( fuck % 2 != 0 ) seq = RevComRead(seq);

		//*********************change********************
		uint32_t seed_bet;
		seed_bet = (seq.length() - K_MER) / SEED_NUM;
		uint32_t seed_found = 0;

		for( int i = 0 ; i < SEED_NUM ; i++ )
		{
			uint32_t start = i * seed_bet;
			uint32_t l2r = transfer(seq,K_MER,start);
			//cout << l2r << " " << kmerHash[l2r] << " " << kmerHash[l2r+1] << endl;
			uint32_t True_window = Each_Read_start + start / 1.1;
			// uint32_t True_window = Each_Read_start + start / 1.1;

			if( kmerHash[l2r+1] <  kmerHash[l2r])
				continue;
			//else if( (kmerHash[l2r+1] - kmerHash[l2r]) > 1000 )
			//	continue;
			else
			{
				for( int j = kmerHash[l2r] ; j < kmerHash[l2r+1] ; j++ )
				{
					//binary_search
					if( SR_HASH[j].SR * Hpart <= True_window  && (SR_HASH[j].SR + 1) * Hpart >= True_window )
					{
						seed_found = seed_found + SR_HASH[j].freq;
					//	cout << i << " " << SR_HASH[j].SR << endl;
						//seed_found++;
					}
					//cout << SR_HASH[j].SR << endl;
				}
				//break;
			}
		}
		//cout << seed_found <<endl;
		//if( fuck == 10 )
		//	break;

		filein.getline(temp,readLINE,'\n');
		filein.getline(temp,readLINE,'\n');
		//break;
		fuck++;
		fout << "Read:" << fuck << "\t" << seq.length() << "\t" << seed_bet << "\t" << seed_found << "\t" << Each_Read_start << endl;
		//if( fuck == 100 ) break;
		if( seed_found < 100 )
		{
			seed_found_L++;
			low_cover_read << "Read:" << fuck << "\t" << seq.length() << "\t" << seed_bet << "\t" << seed_found << "\t" << Each_Read_start << endl;
		}
		sum_cover = sum_cover + seed_found;
	}
	filein.close();
	filein_sup.close();
	fout.close();
	low_cover_read.close();
	cout << seed_found_L << endl;
	cout << sum_cover << endl;
	cout << sum_cover / fuck << endl; //average coverage
	free(kmerHash);
	free(SR_HASH);
	return 0;
}
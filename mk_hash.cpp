#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include "mk_hash.h"

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

int compare(const void *a, const void *b)
{
	uint32_t *pa = (uint32_t*)a;
	uint32_t *pb = (uint32_t*)b;
	return *pa > *pb;  //从小到大排序
}

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

int Hash::mk_hash(char *path, char *genome, uint8_t kmer, uint32_t len_genome)
{
	int movement = 2 * kmer;
	uint64_t Kmer_code_length = 1 << movement;

	uint32_t *kmerHash = (uint32_t*)malloc(sizeof(uint32_t)*Kmer_code_length);
	if( NULL == kmerHash) {
		cout<<"Fail to allocate new space to hash table"<<endl;
		exit(1);
	}

	for( int i = 0 ; i < Kmer_code_length ; i++ )
		kmerHash[i] = 0;

	uint64_t SR_length = len_genome / Hpart;

	//cout << SR_length << endl;

	for( int i = 0 ; i < SR_length ; i++ )
	{
		uint32_t SQ_start = i * Hpart;
		uint32_t *tempArray = (uint32_t*)malloc(Hpart*sizeof(uint32_t));
		if( NULL == tempArray) {
			cout<<"Fail to allocate new space to tempArray"<<endl;
			exit(1);
		}
		for( int j = 0 ; j < Hpart ; j++ )
		{
			uint32_t l2r = transfer(genome,kmer,j+SQ_start);
			tempArray[j] = l2r;
			//Qsort
			//delete duplication
			//caculation
			//cout << j << " " << tempArray[j] << endl;
		}
		qsort(tempArray, Hpart, sizeof(uint32_t), compare);
		uint32_t flag = -1;
		for( int j = 0 ; j < Hpart ; j++ )
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
		free(tempArray);
	}
	uint64_t Hash2sum = 0;
	for( int i = 0 ; i < Kmer_code_length ; i++ )
	{
		uint32_t temp = kmerHash[i];
		kmerHash[i] = Hash2sum;
		Hash2sum = Hash2sum + temp;
		//kmerHash[i] = Hash2sum - temp;
	}
	
	//cout << Hash2sum << endl;

	SRdata * SR_HASH = (SRdata*)malloc(Hash2sum*sizeof(SRdata));
	if( NULL == SR_HASH) {
		cout<<"Fail to allocate new space to SR_HASH"<<endl;
		exit(1);
	}

	FILE * outfile;
	//outfile = fopen("../hash_1", "wb");
	const char *str = path;
 	const char *str2 = "Hash_1";
 	const size_t len = strlen(str)+strlen(str2);
	char *n_str = new char[len+1];
	strcpy(n_str , str);
	strcat(n_str ,str2);
	//cout<<n_str<<endl;
	outfile = fopen(n_str, "wb");
	delete [] n_str;
	fwrite( kmerHash, sizeof(uint32_t), Kmer_code_length,outfile );
	fclose(outfile);

	for( int i = 0 ; i < SR_length ; i++ )
//		for( int i = 19 ; i < 21 ; i++ )
	{
		uint32_t SQ_start = i * Hpart;
		uint32_t *tempArray = (uint32_t*)malloc(Hpart*sizeof(uint32_t));
		if( NULL == tempArray) {
			cout<<"Fail to allocate new space to tempArray"<<endl;
			exit(1);
		}
		for( int j = 0 ; j < Hpart ; j++ )
		{
			uint32_t l2r = transfer(genome,kmer,j+SQ_start);
			tempArray[j] = l2r;
		}
		qsort(tempArray, Hpart, sizeof(uint32_t), compare);

		uint32_t flag = -1;
		for( int j = 0 ; j < Hpart ; j++ )
		{
			if( flag != tempArray[j] )
			{
				flag = tempArray[j];
				uint32_t add = kmerHash[flag];
				SR_HASH[add].SR = i;
				kmerHash[flag]++;
				//cout << i << " " << j << " " << tempArray[j] << " " << flag << " " << kmerHash[flag] << endl;
			}
			//cout << j << " " << tempArray[j] << endl;
		}

		/*uint32_t flag = tempArray[0];
		uint32_t flag_n = 0;
		for( int j = 0 ; j < Hpart ; j++ )
		{
			if( flag == tempArray[j] )
			{
				flag_n++;
			}
			else
			{
				uint32_t add = kmerHash[flag];
				SR_HASH[add].SR = i;
				//SR_HASH[add].freq = flag_n;
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
			//SR_HASH[add].freq = flag_n;
			//cout << flag << " " << kmerHash[flag] << " " << SR_HASH[add].SR << " " << SR_HASH[add].freq << " " << kmerHash[flag+1] << endl;
			kmerHash[flag]++;
		}*/
		//cout << i << endl;
		//break;
		free(tempArray);
	}
	/*for( int i = 0 ; i < 100 ; i++ )
		cout << SR_HASH[i].SR << " " << SR_HASH[i].freq << endl;*/
	FILE * outfile2;
	//outfile2 = fopen("../hash_2", "wb");

	const char *str_ = path;
 	const char *str2_ = "Hash_2";
 	const size_t len_ = strlen(str_)+strlen(str2_);
	char *n_str_ = new char[len_+1];
	strcpy(n_str_, str_);
	strcat(n_str_ ,str2_);
	//cout<<n_str<<endl;
	outfile2 = fopen(n_str_, "wb");
	delete [] n_str_;

	fwrite( SR_HASH, sizeof(SRdata), Hash2sum,outfile2 );
	fclose(outfile2);
	free(kmerHash);
	free(SR_HASH);
	return Hash2sum;
}

/**  
 * All rights Reserved, Designed By HIT-Bioinformatics   
 * @Title:  mk_hash.cpp
 * @Package 
 * @Description:    make the hashtable
 * @author: tjiang
 * @date:   June 7 2016
 * @version V1.0     
 */

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <time.h>
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

//char random_N[4] = {'A', 'C', 'G', 'T'};

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

char *get_path(char *path, const char *str2)
{
	const char *str = path;
	const size_t len = strlen(str)+strlen(str2);
	char *n_str = new char[len+1];
	strcpy(n_str , str);
	strcat(n_str ,str2);
	return n_str;
}

int Hash::mk_hash(char *path, char *genome, uint32_t kmer, uint32_t len_genome)
{
	// change the reference
	// N and n changing the random character{A, C ,G ,T }
	/*srand((int)time(NULL));
	for( int i = 0; i < len_genome; i++ )
	{
		if( genome[i] == 'N')//maybe 'n' is the wrong.
			genome[i] = random_N[rand()%4];
	}*/
	FILE * file_ref;
	char *ref_path = NULL;
	ref_path = get_path(path, "Reference");
	file_ref = fopen(ref_path, "wb");
	delete[] ref_path;
	fwrite( genome, sizeof(char), len_genome, file_ref );
	fclose(file_ref);




	//int movement = 2 * kmer;
	uint32_t movement = kmer << 1;
	uint64_t Kmer_code_length = 1 << movement;

	//uint32_t *kmerHash = (uint32_t*)malloc(sizeof(uint32_t)*Kmer_code_length);
	uint32_t *kmerHash = new uint32_t[Kmer_code_length];
	if( NULL == kmerHash) {
		cout<<"Fail to allocate new space to hash table"<<endl;
		exit(1);
	}

	uint32_t *tempArray = new uint32_t[Hpart];
	if( NULL == tempArray )
	{
		cout<<"Fail to allocate new space to tempArray"<<endl;
		exit(1);
	}

	for( uint32_t i = 0 ; i < Kmer_code_length ; i++ )
		kmerHash[i] = 0;

	uint64_t SR_length = len_genome / Hpart;

	for( uint32_t i = 0 ; i < SR_length ; i++ )
	{
		uint32_t SQ_start = i * Hpart;
		/*uint32_t *tempArray = (uint32_t*)malloc(Hpart*sizeof(uint32_t));
		if( NULL == tempArray) {
			cout<<"Fail to allocate new space to tempArray"<<endl;
			exit(1);
		}*/

		uint32_t pre = 0;
		uint32_t mask = 0xffffffff >> ( 32 - ( kmer << 1 ));
		pre = transfer(genome, kmer, SQ_start);
		tempArray[0] = pre;

		for( int j = 1 ; j < Hpart ; j++ )
		{
			pre = ((pre<<2)&mask)|trans[genome[j+SQ_start + kmer - 1 ]];
			//uint32_t l2r = transfer(genome,kmer,j+SQ_start);
			tempArray[j] = pre;
			//Qsort
			//delete duplication
			//caculation
		}
		qsort(tempArray, Hpart, sizeof(uint32_t), compare);
		uint32_t flag = -1;
		for( int j = 0 ; j < Hpart ; j++ )
		{
			if( flag != tempArray[j] )
			{
				flag = tempArray[j];
				kmerHash[flag]++;
			}
		}
		//free(tempArray);
	}
	uint64_t Hash2sum = 0;
	for( uint32_t i = 0 ; i < Kmer_code_length ; i++ )
	{
		uint32_t temp = kmerHash[i];
		kmerHash[i] = Hash2sum;
		Hash2sum = Hash2sum + temp;
	}
	
	//cout << Hash2sum << endl;

	uint32_t *SR_HASH = new uint32_t[Hash2sum];
	if( NULL == SR_HASH) {
		cout<<"Fail to allocate new space to SR_HASH"<<endl;
		exit(1);
	}

	FILE * outfile;
	char *hash1_path = NULL;
	hash1_path = get_path(path, "Hash_1");
	outfile = fopen(hash1_path, "wb");
	delete[] hash1_path;
	fwrite( kmerHash, sizeof(uint32_t), Kmer_code_length, outfile );
	fclose(outfile);

	for( uint32_t i = 0 ; i < SR_length ; i++ )
	{
		uint32_t SQ_start = i * Hpart;
		/*uint32_t *tempArray = (uint32_t*)malloc(Hpart*sizeof(uint32_t));
		if( NULL == tempArray) {
			cout<<"Fail to allocate new space to tempArray"<<endl;
			exit(1);
		}*/
		uint32_t pre = 0;
		uint32_t mask = 0xffffffff >> ( 32 - ( kmer << 1 ));
		pre = transfer(genome, kmer, SQ_start);
		tempArray[0] = pre;
		for( int j = 1 ; j < Hpart ; j++ )
		{
			pre = ((pre<<2)&mask)|trans[genome[j+SQ_start + kmer - 1 ]];
			tempArray[j] = pre;
		}
		/*for( int j = 0 ; j < Hpart ; j++ )
		{
			uint32_t l2r = transfer(genome,kmer,j+SQ_start);
			tempArray[j] = l2r;
		}*/
		qsort(tempArray, Hpart, sizeof(uint32_t), compare);

		uint32_t flag = -1;
		for( int j = 0 ; j < Hpart ; j++ )
		{
			if( flag != tempArray[j] )
			{
				flag = tempArray[j];
				uint32_t add = kmerHash[flag];
				SR_HASH[add] = i;
				kmerHash[flag]++;
			}
		}
		//free(tempArray);
	}

	FILE * outfile2;
	char *hash2_path = NULL;
	hash2_path = get_path(path, "Hash_2");
	outfile2 = fopen(hash2_path, "wb");
	delete[] hash2_path;
	fwrite( SR_HASH, sizeof(uint32_t), Hash2sum, outfile2 );
	fclose(outfile2);

	delete[] tempArray;
	delete[] kmerHash;
	delete[] SR_HASH;
	//may be change to outside

	FILE *LengthFile;
	char *LF_path = NULL;
	LF_path = get_path(path, "LengthFile");
	LengthFile = fopen(LF_path,"w");
	delete[] LF_path;
	fprintf(LengthFile, "%lu\n", Hash2sum);
	fprintf(LengthFile, "%u\n", len_genome);
	fclose(LengthFile);
	/*
		Maybe mote the number of Hash2sum.
	*/
	return Hash2sum;
}

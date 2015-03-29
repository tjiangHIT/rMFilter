#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

#define readLINE 30000
#define SV_RESULTS 322
#define LOW_READ 622656

typedef struct 
{
	string chr;
	int start;
	int end;
	int length;
	string flag;
}SV_DASTA;

typedef struct
{	
	string name;
	int length;
	int start;
}POS;

int compare(const void *a, const void *b)
{
	POS *pa = (POS*)a;
	POS *pb = (POS*)b;
	return pa->start - pb->start;  //从小到大排序
}

int cmp(const void *a, const void *b)
{
	int *pa = (int *)a;
	int *pb = (int *)b;
	return *pa > *pb;
}

int main()
{
	ifstream filein_svresults;
	filein_svresults.open("../chr1_sort.RSV");
	char temp[readLINE];
	//SV_DASTA * data_sv = (SV_DASTA*)malloc(sizeof(SV_DASTA)*SV_RESULTS);
	SV_DASTA * data_sv = new SV_DASTA[SV_RESULTS];
	if( NULL == data_sv) {
		cout<<"Fail to allocate new space to data_sv"<<endl;
		exit(1);
	}
	//SV_DASTA data_sv[SV_RESULTS];
	int sum = 0;
	int INS_L = 0;
	for( int i = 0 ; i < SV_RESULTS ; i++ )
	{
		string str;
		filein_svresults.getline( temp, readLINE, '\t' );
		data_sv[i].chr = temp;
		//cout << data_sv[i].chr << "\t";
		filein_svresults.getline( temp, readLINE, '\t' );
		str = temp;
		data_sv[i].start = atoi(str.c_str());
		//data_sv[i].start = atoi(str.c_str()) + INS_L;
		//cout << data_sv[i].start <<  "\t";
		filein_svresults.getline( temp, readLINE, '\t' );
		str = temp;
		data_sv[i].end = atoi(str.c_str());
		//data_sv[i].end = atoi(str.c_str()) + INS_L;
		//cout << data_sv[i].end <<  "\t";
		filein_svresults.getline( temp, readLINE, '\t' );
		str = temp;
		data_sv[i].length = atoi(str.c_str());
		//cout << data_sv[i].length <<  "\t";
		filein_svresults.getline( temp, readLINE, '\n' );
		data_sv[i].flag = temp;
		/*if( data_sv[i].flag == "INS" )
		{
			data_sv[i].end = data_sv[i].end + data_sv[i].length;
			INS_L = INS_L + data_sv[i].length;
		}*/
		//cout << data_sv[i].flag << endl;
		//sum++;
		//break;
	}
	filein_svresults.close();
	/*for( int i = 0 ; i < SV_RESULTS ; i++ )
	{
		cout  << i << "\t" << data_sv[i].start << "\t" << data_sv[i].end << "\t" << data_sv[i].length << "\t" << data_sv[i].flag << endl;
	}*/

	ofstream fout;
	fout.open("../SV_read.version(1.1,100)");

	POS *read = new POS[LOW_READ];
	if( NULL == read) {
		cout<<"Fail to allocate new space to read"<<endl;
		exit(1);
	}

	ifstream filein_length;
	filein_length.open("../low_cover_read.version(1.1,100)");
	//filein_length.open("../NEW_results.version(1.1,50)");

	//************************
	int SV_note[readLINE];
	int SV_note_L = 0;
	//************************

	for( int i = 0 ; i < LOW_READ ; i++ )
	{
		string str;
		filein_length.getline(temp,readLINE,'\t');
		read[sum].name = temp;
		//cout << read[sum].name << "\t";
		filein_length.getline(temp,readLINE,'\t');
		str = temp;
		read[sum].length = atoi(str.c_str());
		//cout << read[sum].length << "\t";
		filein_length.getline(temp,readLINE,'\t');
		filein_length.getline(temp,readLINE,'\t');
		filein_length.getline(temp,readLINE,'\n');
		str = temp;
		read[sum].start = atoi(str.c_str());
		//cout << read[sum].start << endl;
		//break;
		int read_length = read[sum].start + read[sum].length / 1.1;
		for( int i = 0 ; i < SV_RESULTS ; i++ )
		{
			if( (read[sum].start >= data_sv[i].start && read[sum].start <= data_sv[i].end) ||
				(read_length >= data_sv[i].start && read_length <= data_sv[i].end) || (read[sum].start <= data_sv[i].start && read_length >=data_sv[i].end) )
			{
				fout << read[sum].name << "\t" << read[sum].start << "\t" << read_length << "\t" << i+1 << "\t" << data_sv[i].flag << endl;
				SV_note[SV_note_L] = i+1;
				SV_note_L++;
			}
		}
		sum++;
	}
	//filein_real_start.close();
	fout.close();
	filein_length.close();

	//test SV_not_found

	ofstream NoFound_note;
	NoFound_note.open("../NoFound_note");
	char note_flag[SV_RESULTS+1];
	for( int i = 0 ; i < SV_RESULTS+1 ; i++ )
	{
		note_flag[i] = 'N';
	}
	for( int i = 0 ; i < SV_note_L ; i++ )
	{
		note_flag[SV_note[i]] = 'Y';
	}
	for( int i = 1 ; i < SV_RESULTS+1 ; i++ )
	{
		if( note_flag[i] == 'N' )
		{
			NoFound_note << i << endl;
		}
	}
	NoFound_note.close();
	delete[] data_sv;
	delete[] read;
	return 0;
}
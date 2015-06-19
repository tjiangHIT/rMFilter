#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

using namespace std;

#define readLINE 30000
#define SV_RESULTS 322
#define LOW_READ 622656
#define Standard 100
#define KMER_STANDARD 750

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
	int kmer_size;
	int start;
	int end;
	int score;
}POS;

int main()
{
	ifstream filein_svresults;
	filein_svresults.open("../chr1_sort.RSV");
	char 		temp[readLINE];
	string 		str;
	//SV_DASTA * data_sv = (SV_DASTA*)malloc(sizeof(SV_DASTA)*SV_RESULTS);
	SV_DASTA * data_sv = new SV_DASTA[SV_RESULTS];
	if( NULL == data_sv) {
		cout<<"Fail to allocate new space to data_sv"<<endl;
		exit(1);
	}
	//SV_DASTA data_sv[SV_RESULTS];
	for( int i = 0 ; i < SV_RESULTS ; i++ )
	{
		filein_svresults.getline( temp, readLINE, '\t' );
		data_sv[i].chr = temp;
		filein_svresults.getline( temp, readLINE, '\t' );
		str = temp;
		data_sv[i].start = atoi(str.c_str());
		filein_svresults.getline( temp, readLINE, '\t' );
		str = temp;
		data_sv[i].end = atoi(str.c_str());
		filein_svresults.getline( temp, readLINE, '\t' );
		str = temp;
		data_sv[i].length = atoi(str.c_str());
		filein_svresults.getline( temp, readLINE, '\n' );
		data_sv[i].flag = temp;
	}
	filein_svresults.close();
	/*for( int i = 0 ; i < SV_RESULTS ; i++ )
	{
		cout  << i << "\t" << data_sv[i].start << "\t" << data_sv[i].end << "\t" << data_sv[i].length << "\t" << data_sv[i].flag << endl;
	}*/

	POS *read = new POS[LOW_READ];
	if( NULL == read) {
		cout<<"Fail to allocate new space to read"<<endl;
		exit(1);
	}

	ifstream Filein;
	Filein.open("../read_start.test");
	for(int i = 0; i < LOW_READ; i++ )
	{
		Filein.getline( temp, readLINE, '\t');
		str = temp;
		read[i].start = atoi(str.c_str());
		Filein.getline( temp, readLINE, '\n');
		str = temp;
		read[i].end = atoi(str.c_str());
	}
	Filein.close();

	ifstream filein;
	filein.open("../LowHITread.newtry.2.Version_bp");
	ofstream Judge;
	Judge.open("../Answer.test");
	ifstream Nextctrl;
	Nextctrl.open("../local_index_results.150");
	
	int		SV_note[readLINE];
	int		SV_note_L = 0;
	int		sumLowScoreRead = 0;
	int		LowScoreSVRead = 0;
	
	for(int i = 0; i < LOW_READ; i++ )
	{
		Nextctrl.getline( temp, readLINE, '\t');
		Nextctrl.getline( temp, readLINE, '\n');
		char finalFlag = temp[0];

		filein.getline( temp, readLINE, '\t');
		read[i].name = temp;
		filein.getline( temp, readLINE, '\t');
		str = temp;
		read[i].length = atoi(str.c_str());
		filein.getline( temp, readLINE, '\n');
		str = temp;
		read[i].kmer_size = atoi(str.c_str());

		filein.getline( temp, readLINE, '\t');
		// str = temp;
		// read[i].start_1 = atoi(str.c_str());
		filein.getline( temp, readLINE, '\t');
		// str = temp;
		// read[i].start_2 = atoi(str.c_str());
		filein.getline( temp, readLINE, '\n');
		str = temp;
		read[i].score = atoi(str.c_str());
		uint32_t START = read[i].start;
		uint32_t END = read[i].end;
		//if( read[i].score <= Standard && read[i].score > StandardLow )
		//if( read[i].score <= Standard || read[i].kmer_size >= KMER_STANDARD )
		/*{

		}
		else*/
		if( finalFlag == 'T' )
		{
			//Nextctrl << "T" << endl;
			sumLowScoreRead++;
			for( int k = 0 ; k < SV_RESULTS ; k++ )
			{
				if( (START >= data_sv[k].start && START <= data_sv[k].end) || (START >= data_sv[k].start && END <= data_sv[k].end) ||
					(END >= data_sv[k].start && END <= data_sv[k].end) || (START <= data_sv[k].start && END >= data_sv[k].end) )
				{
					LowScoreSVRead++;
					SV_note[SV_note_L] = k+1;
					SV_note_L++;
					Judge << k+1 << "\t" << data_sv[k].start << "\t" << data_sv[k].end << "\t" << data_sv[k].flag << "\t";
					Judge << i+1 << "\t" << START << "\t" << END << endl;
				}
				/*if( (START > data_sv[k].start && END < data_sv[k].end) )
				{
					cout << i << endl;
					cout << read[i].length << endl;
					cout << START << " " << END << endl;
					cout << data_sv[k].start << " " << data_sv[k].end << endl;
					return 0;
				}*/
			}
		}
		/*else
			Nextctrl << "N" << endl;*/
	}
	filein.close();
	Judge.close();
	Nextctrl.close();

	//test SV_not_found

	ofstream NoFound_note;
	NoFound_note.open("../NoFound_note");
	char		note_flag[SV_RESULTS+1];
	int		NotFoundNum = 0;
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
			NotFoundNum++;
		}
	}
	NoFound_note.close();
	delete[] data_sv;
	delete[] read;
	double percent_below = sumLowScoreRead * 1.0 / LOW_READ;
	double percent_SV = LowScoreSVRead * 1.0 / sumLowScoreRead;
	cout << "Below " << Standard  << " & " << KMER_STANDARD << " is " << sumLowScoreRead << "(" << percent_below << ")." << endl;
	cout << "Including SV is "<< LowScoreSVRead <<"(" << percent_SV << ")." << endl;
	cout << "No finding SV is "<< NotFoundNum << endl;
	return 0;
}
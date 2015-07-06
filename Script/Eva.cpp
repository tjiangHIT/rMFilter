/*
# Date: Jul 6, 2015 $
# Author: tjiang $
# Purpose: Evaluation results #
*/

#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

typedef struct
{
	char flag;
	int Length;
	int clippingL;
}SamNode;

typedef struct
{
	char flag;
}ReadNode;

int main( int argc, char *argv[] )
{
	if (argc < 3) {
		fprintf(stderr, "Usage: %s <EvaSam> <EvaRead>\n", argv[0]);
		return 1;
	}

	ifstream file1, file2;
	file1.open(argv[1]);
	file2.open(argv[2]);
	char temp[1000];
	string str;
	int TrueAns = 0;
	int SUM = 0;
	int MaybeTrue = 0;
	while( !file1.eof() )
	//for( int i = 0 ; i < 10 ; i++ )
	{
		ReadNode noderead;
		SamNode nodesam;

		file1.getline(temp, 1000, '\t');
		file1.getline(temp, 1000, '\t');
		nodesam.flag = temp[0];
		file1.getline(temp, 1000, '\t');
		str = temp;
		nodesam.Length = atoi(str.c_str());
		file1.getline(temp, 1000, '\t');
		file1.getline(temp, 1000, '\n');
		str = temp;
		nodesam.clippingL = atoi(str.c_str());

		file2.getline(temp, 1000, '\t');
		file2.getline(temp, 1000, '\n');
		noderead.flag = temp[0];

		if( nodesam.flag == 'Y' )
		{
			float percentage = nodesam.clippingL * 1.0 / nodesam.Length;
			if( percentage > 0.2 && noderead.flag == 'T' )
				TrueAns++;
			if( percentage > 0.2 )
				MaybeTrue++;
		}
		else
		{
			MaybeTrue++;
			if( noderead.flag == 'T' )
				TrueAns++;
		}

		if( noderead.flag == 'T' )
			SUM++;
	}

	file1.close();
	file2.close();
	cout << "TrueAns: " << TrueAns << endl;
	cout << "SUM: " << SUM << endl; 
	cout << "MaybeTrue: " << MaybeTrue << endl; 
	return 0;
}
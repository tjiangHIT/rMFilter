#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>

using namespace std;

#define SV_NUM 368
#define READ_NUM 623212

typedef struct{
	int start;
	int end;
	int length;
	int bia;
	string sv;
}AnswerWin;

int main(int argc, char *argv[])
{
	if (argc < 5) {
		fprintf(stderr, "Usage: %s <in.RSV> <read.L> <read.donor> <read.test>\n", argv[0]);
		return 1;
	}

	AnswerWin		a_win[SV_NUM];
	char			temp[1000];
	string			str;
	int			tempbia = 0;
	ifstream answer_in;
	answer_in.open(argv[1]);

	for( int k = 0; k < SV_NUM; k++ )
	{
		answer_in.getline(temp, 1000, '\t');
		answer_in.getline(temp, 1000, '\t');
		str = temp;
		int sW = atoi(str.c_str());
		answer_in.getline(temp, 1000, '\t');
		str = temp;
		int eW = atoi(str.c_str());
		answer_in.getline(temp, 1000, '\t');
		str = temp;
		a_win[k].length = atoi(str.c_str());

		a_win[k].start = sW + tempbia;

		answer_in.getline(temp, 1000, '\n');
		str = temp;
		a_win[k].sv = str;
		if( str == "DEL" )
		{
			tempbia -= a_win[k].length;
			a_win[k].end = eW + tempbia;
		}
		if( str == "INS" )
		{
			tempbia += a_win[k].length;
			a_win[k].end = eW + tempbia;
		}
		a_win[k].bia = tempbia;
	}
	answer_in.close();

	/*ofstream answer_out;
	answer_out.open("../chr1_sort.RSV.change");
	for( int k = 0; k < SV_NUM; k++ )
	{
		answer_out << a_win[k].start << "\t" << a_win[k].end << "\t" << a_win[k].bia << endl;
	}
	answer_out.close();*/

	ifstream		startIN;
	int			startOri[READ_NUM];
	int			endOri[READ_NUM];
	int			Length[READ_NUM];

	startIN.open(argv[3]);
	for( int i = 0 ; i < READ_NUM ; i++ )
	{
		startIN.getline(temp, 1000, '\t');
		startIN.getline(temp, 1000, '\n');
		str = temp;
		startOri[i] = atoi(str.c_str());
	}
	startIN.close();

	ifstream LengthIN;
	LengthIN.open(argv[2]);
	for( int i = 0 ; i < READ_NUM ; i++ )
	{
		LengthIN.getline(temp, 1000, '\n');
		str = temp;
		Length[i] = atoi(str.c_str());
		//if( i % 2 == 0 )
			endOri[i] = startOri[i] + Length[i];
		/*else
		{
			if( startOri[i] > Length[i] )
				endOri[i] = startOri[i] - Length[i];
			else
				endOri[i] = 0;
		}*/
	}
	LengthIN.close();

	for( int i = 0 ; i < READ_NUM ; i++ )
	{
		int Bia = 0;
		int flag = 0;
		for( int k = 0; k < SV_NUM; k++ )
		{
			if( startOri[i] > a_win[k].start )
			{
				Bia = a_win[k].bia;
				if( a_win[k].sv == "INS" )
				{
					if( startOri[i] <= a_win[k].end )
					{
						startOri[i] = a_win[k].end - Bia;
						flag = 1;
						break;
					}
				}
			}
			else
			{
				//startOri[i] -= Bia;
				break;
			}
		}
		if( flag == 0 )
			startOri[i] -= Bia;
	}

	for( int i = 0 ; i < READ_NUM ; i++ )
	{
		int Bia = 0;
		int flag = 0;
		for( int k = 0; k < SV_NUM; k++ )
		{
			if( endOri[i] > a_win[k].start )
			{
				Bia = a_win[k].bia;
				if( a_win[k].sv == "INS" )
				{
					if( endOri[i] <= a_win[k].end )
					{
						flag = 1;
						endOri[i] = a_win[k].end - Bia;
						break;
					}
				}
			}
			else
			{
				//endOri[i] -= Bia;
				break;
			}
		}
		if( flag == 0 )
			endOri[i] -= Bia;
	}

	ofstream Out;
	Out.open(argv[4]);
	int score = 0;
	for( int i = 0 ; i < READ_NUM ; i++ )
	{
/*		if( i % 2 == 0 )
		{*/
			Out << startOri[i] << "\t" << endOri[i] << endl;
			int LS = endOri[i] - startOri[i];
			if( LS < 0 )
				score++;
/*		}
		else
		{
			Out << endOri[i] << "\t" << startOri[i] << endl;
			int LS = startOri[i] -  endOri[i];
			if( LS != Length[i] )
				score++;
		}*/
	} 
	Out.close();
	cout << score << endl;
	return 0;
}
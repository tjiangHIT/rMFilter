#include <iostream>
#include <fstream>
#include <stdint.h>
#include <stdlib.h>

using namespace std;

#define readLINE 30000

int main(int argc, char *argv[])
{
	if (argc < 4) {
		fprintf(stderr, "Usage: %s <in.maf> <start.donor> <Length(int)>\t\n", argv[0]);
		return 1;
	}

	int num = 0;
	ofstream file_out;
	file_out.open(argv[2]);
	ifstream filein_sup;
	filein_sup.open(argv[1]);
	char temp[readLINE];
	string str = argv[3];
	int L = atoi(str.c_str());
	for( int i = 0 ; i < L ; i++ )
	{
		string strtemp;
		filein_sup.getline(temp,readLINE,'\n');
		filein_sup.getline(temp,readLINE,' ');
		filein_sup.getline(temp,readLINE,' ');
		filein_sup.getline(temp,readLINE,' ');
		strtemp = temp;
		while( strtemp.length() == 0 )
		{
			filein_sup.getline(temp,readLINE,' ');
			strtemp = temp;
		}
		//uint32_t Each_Read_start;
		//Each_Read_start = atoi(strtemp.c_str());
		file_out << num++ << "\t" << strtemp << endl;
		//cout << num++ << "\t" << strtemp << endl;
		//cout << fuck << " " << Each_Read_start << " " << strtemp.length() << endl;
		filein_sup.getline(temp,readLINE,'\n');
		filein_sup.getline(temp,readLINE,'\n');
		filein_sup.getline(temp,readLINE,'\n');
	}
	file_out.close();
	filein_sup.close();
	cout << num << endl;
	return 0;
}
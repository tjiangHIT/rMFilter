#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdint.h>

#define readLINE 30000
#define readSum 622656
#define betWeen 10

using namespace std;

int main()
{
	ifstream filein;
	filein.open("../LowHITread.newtry.2");
	
	char 		temp[readLINE];
	string 		str;
	uint64_t	sumScore = 0;
	uint32_t	numRead[betWeen] = {0};
	double		percent[betWeen];
	int		arrayLen[readSum];
	
	for(int i = 0; i < readSum; i++ )
	{
		filein.getline( temp, readLINE, '\t');
		filein.getline( temp, readLINE, '\n');
		str = temp;
		arrayLen[i] = atoi(str.c_str());
		filein.getline( temp, readLINE, '\t');
		filein.getline( temp, readLINE, '\t');
		filein.getline( temp, readLINE, '\n');
		str = temp;
		uint32_t score = atoi(str.c_str());
		sumScore += score;
		for( int k = 0 ; k < betWeen ; k++ )
			if( score < betWeen * k + betWeen )
				numRead[k]++;
	}

	filein.close();
	cout  << "Average Score is " << sumScore / readSum << endl;
	for( int k = 0 ; k < betWeen ; k++ )
	{
		percent[k] = numRead[k] * 1.0 / readSum;
		cout << "Score below " << k*betWeen+betWeen << " is " << percent[k] << endl;
	}
	// ofstream out;
	// out.open("../read_L");
	// for(int i = 0; i < readSum; i++ )
	// 	out << arrayLen[i] << endl;
	// out.close();
	return 0;
}
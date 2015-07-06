# Date: Jul 6, 2015 $
# Author: tjiang $
# Purpose: Screening longer reads #
# Usage: python(at least v 3.0) change_reads.py reads.fastq > somefile 

import sys,re

def main():
	ori_readsFile = open(sys.argv[1])
	line = ori_readsFile.readline()
	while line:
		Group_A = line
		line = ori_readsFile.readline()
		Group_B = line
		line = ori_readsFile.readline()
		Group_C = line
		line = ori_readsFile.readline()
		Group_D = line
		line = ori_readsFile.readline()

		if len(Group_B) > 1025:
			print Group_A,
			print Group_B,
			print Group_C,
			print Group_D,


if __name__ == '__main__':
	main()
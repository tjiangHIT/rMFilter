# Date: Jul 6, 2015 $
# Author: tjiang $
# Purpose: Evaluation results #
# Usage: python(at least v 3.0) EvaSam EvaReads

from __future__ import division  
import sys,re

def main():
	File_1 = open(sys.argv[1])
	File_2 = open(sys.argv[2])

	RightAns = 0
	MaybeAns = 0
	SumAns = 0
	RightAns_outN = 0
	MaybeAns_outN = 0
	SumAns_outN = 0

	for line in File_1:
		seq_1 = line.strip('\n').split('\t')
		seq_2 = File_2.readline().strip('\n').split('\t')

		if seq_2[1] == 'T':
			SumAns = SumAns + 1

		if seq_1[1] == 'N':
			MaybeAns = MaybeAns + 1
			if seq_2[1] == 'T':
				RightAns = RightAns + 1

		if seq_1[1] == 'Y':
			if seq_2[1] == 'T':
				SumAns_outN = SumAns_outN + 1

			if int(seq_1[4]) / int(seq_1[2]) >= 0.2:
				MaybeAns = MaybeAns + 1
				MaybeAns_outN = MaybeAns_outN + 1
				if seq_2[1] == 'T':
					RightAns = RightAns + 1
					RightAns_outN = RightAns_outN + 1


	print "Containing erroneous data."
	print ( "  RightAns:\t%d" % RightAns )
	print ( "  MaybeAns:\t%d" % MaybeAns )
	print ( "  SumAns:\t\t%d" % SumAns )
	print ( "  The recognition rate is %.6f." % round( RightAns / SumAns , 6))
	print ( "  The recall rate is %.6f." % round( RightAns / MaybeAns , 6))

	print "Removing erroneous data."
	print ( "  RightAns_outN:\t%d" % RightAns_outN )
	print ( "  MaybeAns_outN:\t%d" % MaybeAns_outN )
	print ( "  SumAns_outN:\t\t%d" % SumAns_outN )
	print ( "  The recognition rate is %.6f." % round( RightAns_outN / SumAns_outN , 6))
	print ( "  The recall rate is %.6f." % round( RightAns_outN / MaybeAns_outN , 6))


if __name__ == '__main__':
	main()
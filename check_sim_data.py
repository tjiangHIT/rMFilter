#version Python 2.7
#usage: python simdata.clipping simdata.Finalanswer simdata.originalgap
#check simdata extra_filter_percentage
#simdata.clipping come from checksam.py
#simdata.Finalanswer come from VRDT-alignment
#simdata.originalgap come from chm1-PrintGaps.py
#setup the per by handinput(0.2)
#data:2015.12.1

from __future__ import division
import sys,re

def readReadName(dicRead):
	File = open(sys.argv[2])
	for line in File:
		seq = line.strip('\n').split('\t')
		if( seq[1] != 'N' ):
			if seq[0] not in dicRead:
				dicRead[seq[0]] = 0
			else:
				dicRead[seq[0]] += 1

def chm1Name(dicchm1):
	File = open(sys.argv[3])
	for line in File:
		seq = line.strip('\n').split('\t')
		if seq[7].split('/')[0] not in dicchm1:
			dicchm1[seq[7].split('/')[0]] = 0
		else:
			dicchm1[seq[7].split('/')[0]] += 1

def judgeSAM(dicSAM, per):
	File_1 = open(sys.argv[1])
	count = 0
	for line in File_1:
		count += 1
		seq = line.strip('\n').split('\t')
		if( seq[1] == 'N' ):
			# count += 1
			if seq[0] not in dicSAM:
				dicSAM[seq[0]] = 0
			else:
				dicSAM[seq[0]] += 1
			# print seq[0]
		else:
			RL = int(seq[2])
			CL = int(seq[4])
			percent = CL / RL
			if( percent >= per ):
				# count += 1
				if seq[0] not in dicSAM:
					dicSAM[seq[0]] = 0
				else:
					dicSAM[seq[0]] += 1
				# print seq[0]
	return count


if __name__ == '__main__':
	dicSAM = {}
	dicRead = {}
	dicchm1 = {}
	per = 0.2
	print "  *percentage",per
	count = judgeSAM(dicSAM,per)
	readReadName(dicRead)
	print " Filter reads",len(dicRead)
	chm1Name(dicchm1)
	print " chm1 number",len(dicchm1)
	# print count
	print " clipping count",len(dicSAM)
	outlier = 0
	for key in dicSAM:
		if key not in dicchm1:
			if key in dicRead:
				outlier += 1
	print " Wrong number",outlier
	extra = len(dicRead) - outlier
	outlier = 0
	for key in dicchm1:
		if key in dicRead:
			outlier += 1
	print " True number",outlier
	extra = extra - outlier
	print " Extra",extra, '%.4f'%(extra / count)

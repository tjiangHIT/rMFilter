# rMFilter
Acceleration of long read-based structure variation detection with chimeric read filtering

---
###Getting Start
	$ git clone https://github.com/tjiangHIT/rMFilter.git
	$ cd rMFilter
	$ make
	$ ./rMFilter-indexer indexDir ref.fa
	$ ./rMFilter-aligner indexDir read.fq > read.filter.fq
---	
###Introduction

---
###Memory usage

The memory usage of rMFilter can fit the configurations of most modern servers and workstations.
Its peak memory footprint is about 18.00 Gigabytes (default setting), on a server with Intel Xeon CPU at 2.00 GHz, 1 Terabytes RAM running Linux Ubuntu 14.04. These reads were aligned to human reference genome GRCh37/hg19.

---
###Installation

Current version of rMFilter needs to be run on Linux operating system.
The source code is written in C++, and can be directly download from: https://github.com/hitbc/rMFilter 
A mirror is also in: https://github.com/tjiangHIT/rMFilter
The makefile is attached. Use the make command for generating the executable file.

---
###Synopsis
Reference genome indexing
	
	rMFilter-indexer [-k kmerSize] <HashIndexDir> <Reference>
Read alignment & filtering
	
	rMFilter-aligner [-k kmerSize] [-t threadNumber] <HashIndexDir> <ReadFile>

---
###Commands and options

	rMFilter-indexer:
	-k, --kmerSize       [INT]    The size of the k-mers extracted from the reference genome for indexing. [15]

	rMFilter-aligner:
	-k, --kmerSize       [INT]    The size of the k-mers extracted from the reference genome for reading index. [15] 

	-r, --ratio          [INT]    The candidate ratio of the filtering. [0.05]

	-t, --threads        [INT]    The number of threads. [1]

---
###Reference

---
###Contact
For advising, bug reporting and requiring help, please contact ydwang@hit.edu.cn or tjiang@hit.edu.cn


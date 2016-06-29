# rMFilter: 
Acceleration of long read-based structure variation calling by chimeric read filtering

---
###Getting Start
	$ git clone https://github.com/hitbc/rMFilter.git (git clone https://github.com/tjiangHIT/rMFilter.git)
	$ cd rMFilter
	$ make clean && make
	$ ./rMFilter-indexer indexDir ref.fa
	$ ./rMFilter-aligner indexDir read.fq > read.filter.fq
---	
###Introduction
rMFilter is an efficient tool to filter chimeric noisy long reads produced by 3rd generation sequencing platform, such as PacBio SMRT sequencing, to accelerate long read-based detection of genome structural variations (SVs). It improves the overall efficiency of SV calling pipeline by directly filtering potential SV spanning reads. The filtration is based on the analysis of the short token matches between the reads and local regions in reference genome. With the filtration, the numbers of the reads input into SV calling pipelines can be greatly reduced; meanwhile, most of the SV-spanning reads can be retained.

rMFilter has been tested on real and simulated SMRT datasets from human genome, the results demonstrate that the tool can fast filter the reads to substantially improve the overall speed of long read-based SV calling pipelines. Moreover, rMFilter can also correctly handle most of the reads to retain the effectiveness of SV calling pipelines.

---
###Memory usage

The memory usage of rMFilter can fit the configurations of most modern servers and workstations.
Its peak memory footprint is about 18.50 Gigabytes (default setting), on a server with Intel Xeon CPU at 2.00 GHz, 1 Terabytes RAM running Linux Ubuntu 14.04. These reads were aligned to human reference genome GRCh37/hg19.

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
rMFilter: acceleration of long read-based structure variation calling by chimeric read filtering. *Under review*

---
###Contact
For advising, bug reporting and requiring help, please contact ydwang@hit.edu.cn or tjiang@hit.edu.cn


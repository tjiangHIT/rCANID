# rCANID: 
rCANID - read Clustering and Assembly-based Novel insertion Detection tool

---
### Getting Start
	           _____    ____    __    _   _____   _____
	    _ __  / __  \  / __ \  |  \  | | |_   _| |  __ \
	   | ^__| | | |_| / /  \ \ |   \ | |   | |   | |  \ \
	   | |    | |  _  | |__| | | |\ \| |   | |   | |  | |
	   | |    | |_| | |  __  | | | \   |  _| |_  | |__/ /
	   |_|    \_____/ |_|  |_| |_|  \__| |_____| |_____/
     
	
	$ git clone https://github.com/hitbc/rCANID.git (git clone https://github.com/tjiangHIT/rCANID.git)
	$ cd rCANID/
	$ bash INSTALL.sh
	$ ./rCANID

---

### Introduction

Novel sequence insertion (NSI) is a class of genome structural variations (SVs) with potentially important biological functions and correlations with phenotypes and diseases. The rapid development of long read sequencing technologies provide the opportunity to more comprehensively study NSIs, since the much longer reads are helpful to the assembly and location of novel sequences. However, state-of-the-art long read-based generic SV detection approaches either only use the signals of chimerically aligned reads or the contigs of de novo assembly, which are not good at NSI detection and/or computationally expensive. Herein, we proposed Read Clustering and Assembly-based Novel Insertion Detection tool (rCANID), a novel long read-based NSI detection approach. rCANID fully takes advantages of chimerically aligned and unaligned reads by its specifically designed read clustering and lightweight local read assembly methods to effectively reconstruct inserted sequences with relatively low computational cost. The benchmarking on both of simulated and real datasets demonstrate that rCANID can sensitively discover NSIs, especially for those having large inserted novel sequences, which could be hard to state-of-the-art approaches. rCANID is suited to be integrated into many computational pipelines to play important roles in many genomic studies.

---

### Simulated datasets

The simulated datasets use for benchmarking are available at: https://drive.google.com/open?id=1TyA-fz7BBk-d2VOlpF5nvWRdXVdPcKej

---
### Dependences
	
	1. pysam
	2. Biopython
	3. ngmlr
	4. samtools
	5. cigar
    6. MECAT
    7. Minimap2

	Python version 2.7

---
### Installation

Current version of rCANID needs to be run on Linux operating system.
The source code is written in python, and can be directly download from: https://github.com/hitbc/rCANID 
A mirror is also in: https://github.com/tjiangHIT/rCANID
The INSTALL.sh is attached. Use the bash command for generating the executable file.

---
### Synopsis
Cluster all of signal reads and unmapped reads respectively.

	rCANID Cluster <alignments> <temp_dir> <output>

Generate high-quality contigs for each cluster.

	rCANID Assemble <Folder>

Detect novel sequence insertions.

	rCANID calling <FASTA/FA> <reference> <temp_dir> <out_type> <signals>

---
### Reference
rCANID: read Clustering and Assembly-based Novel Insertion Detection tool. *(BIBM 2018 Under Review)*

---
### Contact
For advising, bug reporting and requiring help, please contact ydwang@hit.edu.cn or tjiang@hit.edu.cn

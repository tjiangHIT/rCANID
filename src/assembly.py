#!/usr/bin/env python

''' 
 * All rights Reserved, Designed By HIT-Bioinformatics   
 * @Title:  assembly.py
 * @Package: 
 * @Description: Control rCANID pipeline
 * @author: Fu Yilei
 * @date: Augst 17 2018
 * @version V1.0.1     
'''

import sys
import argparse
import time
import os
import logging
from CommandRunner import *

MEMORY = 40
Threads = 4
# coverage = 10


USAGE="""\
    rCANID assembly step.
    A specific folder needs to be given.
"""

def parseArgs(argv):
	parser = argparse.ArgumentParser(prog="process.py Assemble", description=USAGE, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-f', '--folder', help = "Folder path.",  type = str)
	parser.add_argument('-o', '--output', help = "The prefix of output.",  type = str)
	parser.add_argument('-c', '--coverage', help = "The coverage of longest corrected reads to be extracted.", type = int)
	args = parser.parse_args(argv)
	return args


def filecombination(filepath, output):
    with open(output+"merged.contigs.fasta","a+") as k:
            for i in os.listdir(filepath):
                    if i[0] != ".":
                            f = open(filepath+i)
                            k.write(f.read()+"\n")

def run_mecat(file, output, coverage, file_size, Folder):
	logging.info("Runing mecat2pw.")
	cmd = ("cd %s && mecat2pw -j 0 -d %s.fa -o %s.pm.can -t %d -w wrk_dir -n 200 -a 50 -k 2" % (Folder, file, file, Threads))
	r, o, e = exe(cmd)
	logging.info("Runing mecat2cns.")
	cmd = ("cd %s && mecat2cns -i 0 -t %d %s.pm.can %s.fa %s.collect.fa -r 0.2 -a 50 -c 2 -l 50" % (Folder, Threads, file, file, file))
	r, o, e = exe(cmd)
	logging.info("Runing extract_sequences.")
	cmd = ("cd %s && extract_sequences %s.collect.fa %s.%dx %s %d" % (Folder, file, file, coverage, file_size, coverage))
	r, o, e = exe(cmd)
	logging.info("Runing mecat2canu.")
	cmd = ("cd %s && mecat2canu -trim-assemble -p %s -d %s corOutCoverage=%d minOverlapLength=50 minReadLength=50 genomeSize=%s ErrorRate=0.2 maxMemory=40 maxThreads=%d useGrid=0 Overlapper=mecat2asmpw -pacbio-corrected %s.%dx.fasta" % (Folder, file, file, coverage, file_size, Threads, file, coverage))
	r, o, e = exe(cmd)
	cmd = ("cd %s && cat %s/%s.contigs.fasta %s/%s.unassembled.fasta > %s%s.contigs.fasta" % (Folder, file, file, file, file, output, file))
	r, o, e = exe(cmd)
	cmd = ("cd %s && rm  *.fa.fastq *.gkpStore* *.qual *.qv *.frg *.collect.fa *.can *.can.part0 *.partition_files -r wrk_dir %s *.%dx.fasta" % (Folder, file, coverage))
	r, o, e = exe(cmd)

def Merge_data(output):
    logging.info("Merging contigs.")
    cmd = ("cd %s && cat *.contigs.fasta > %sContigs.fa" % (output, output))
    r, o, e = exe(cmd)

def run(argv):
    args = parseArgs(argv)
    setupLogging(False)
    starttime = time.time()

    if args.folder[-1] == '/':
        Folder = args.folder
    else:
        Folder = args.folder + '/'

    files = os.listdir(Folder)
    flist = list()
    for file in files:
        if len(file) == 0:
            continue
        if file.split('.')[-1] == 'fa':
            flist.append(file[:-3])

    for i in flist:
        fsize = os.path.getsize(Folder + i + ".fa")
        try:
            run_mecat(i, args.output, args.coverage, fsize, Folder)
        except:
            pass
    Merge_data(args.output)

    logging.info("Finished in %0.2f seconds."%(time.time() - starttime))
    # filecombination("contigs/", args.prefix_output)
    # os.system("rm -r contigs")

if __name__ == '__main__':
    # os.system("mkdir contigs")
    run(sys.argv[1:])

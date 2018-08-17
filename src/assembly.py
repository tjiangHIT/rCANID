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

USAGE="""\
    This is the reads assembly tools for rCANID. Please enter the folder path.
"""

def parseArgs(argv):
	parser = argparse.ArgumentParser(prog="process.py assembly", description=USAGE, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-f', '--folder', help = "Folder path contains",  type = str)
	args = parser.parse_args(argv)
	return args


def filecombination(filepath):
    with open("merged.contigs.fasta","a+") as k:
            for i in os.listdir(filepath):
                    if i[0] != ".":
                            f = open(filepath+i)
                            k.write(f.read()+"\n")


def run(argv):
    args = parseArgs(argv)
    starttime = time.time()
    filetype = -1
    file = os.listdir(args.folder)
    flist = []
    for i in file:
        if len(i) == 0:
            continue
        if i.split('.')[1] == 'fq':
            flist.append(i.split(".")[0])
            filetype = 0
        elif i.split('.')[1] == 'fastq':
            flist.append(i.split(".")[0])
            filetype = 1
        elif i.split('.')[1] == 'fa':
            flist.append(i.split(".")[0])
            filetype = 2
        elif i.split('.')[1] == 'fasta':
            flist.append(i.split(".")[0])
            filetype = 3

    if filetype == 0:
        for i in flist:
            fsize = os.path.getsize(args.folder + i + ".fq")
            os.system("bash run_MECAT.sh " + args.folder + " " + i + " " + str(fsize) + " fq")
    elif filetype == 1:
        for i in flist:
            fsize = os.path.getsize(args.folder + i +".fastq")
            os.system("bash run_MECAT.sh " + args.folder + " " + i + " " + str(fsize) + " fastq")
    elif filetype == 2:
        for i in flist:
            fsize = os.path.getsize(args.folder + i +".fa")
            os.system("bash run_MECAT.sh " + args.folder + " " + i + " " + str(fsize) + " fa")
    elif filetype == 3:
        for i in flist:
            fsize = os.path.getsize(args.folder + i +".fasta")
            os.system("bash run_MECAT.sh " + args.folder + " " + i + " " + str(fsize) + " fasta")

    logging.info("Finished in %0.2f seconds."%(time.time() - starttime))
    filecombination("contigs/")
    os.system("rm -r contigs")

if __name__ == '__main__':
    os.system("mkdir contigs")
    run(sys.argv[1:])
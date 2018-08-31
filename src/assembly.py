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
from transfer_contig_name import * 
from multiprocessing import Pool

MEMORY = 40
Threads = 4
# coverage = 10


USAGE="""\
    rCANID assembly step.
    A specific folder needs to be given.
"""

def parseArgs(argv):
    parser = argparse.ArgumentParser(prog="process.py Assemble", description=USAGE, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-f', '--folder', help = "Folder path of signal reads.",  type = str)
    parser.add_argument('-F', '--Folder', help = "Folder path of unmapped reads.",  type = str)
    parser.add_argument('-o', '--output', help = "The prefix of output.",  type = str)
    parser.add_argument('-c', '--coverage', help = "The coverage of longest corrected reads to be extracted.", type = int)
    parser.add_argument('-t', '--threads', help = "Number of threads to use.[%(default)s]", default = 16, type = int)

    args = parser.parse_args(argv)
    return args


def filecombination(filepath, output):
    with open(output+"merged.contigs.fasta","a+") as k:
            for i in os.listdir(filepath):
                    if i[0] != ".":
                            f = open(filepath+i)
                            k.write(f.read()+"\n")

def run_mecat(file, output, coverage, file_size, Folder, flag, id):
    logging.info("Assemble No.%d contigs." %(id))
    # logging.info("Runing mecat2pw.")
    cmd = ("cd %s && mecat2pw -j 0 -d %s.fa -o %s.pm.can -t %d -w wrk_dir_%d -n 200 -a 50 -k 2" % (Folder, file, file, Threads, id))
    r, o, e = exe(cmd)
    # logging.info("Runing mecat2cns.")
    cmd = ("cd %s && mecat2cns -i 0 -t %d %s.pm.can %s.fa %s.collect.fa -r 0.2 -a 50 -c 2 -l 50" % (Folder, Threads, file, file, file))
    r, o, e = exe(cmd)
    # logging.info("Runing extract_sequences.")
    cmd = ("cd %s && extract_sequences %s.collect.fa %s.%dx %s %d" % (Folder, file, file, coverage, file_size, coverage))
    r, o, e = exe(cmd)
    # logging.info("Runing mecat2canu.")
    cmd = ("cd %s && mecat2canu -trim-assemble -p %s -d %s corOutCoverage=%d minOverlapLength=50 minReadLength=50 genomeSize=%s ErrorRate=0.2 maxMemory=40 maxThreads=%d useGrid=0 Overlapper=mecat2asmpw -pacbio-corrected %s.%dx.fasta" % (Folder, file, file, coverage, file_size, Threads, file, coverage))
    r, o, e = exe(cmd)
    if flag == 0:
        cmd = ("cd %s && cat %s/%s.contigs.fasta %s/%s.unassembled.fasta > %s%s.contigs.fasta" % (Folder, file, file, file, file, output, file))
    else:
        cmd = ("cd %s && cat %s/%s.contigs.fasta > %s%s.contigs.fasta" % (Folder, file, file, output, file))
    r, o, e = exe(cmd)
    cmd = ("cd %s && rm  *.fa.fastq *.gkpStore* *.qual *.qv *.frg *.collect.fa *.can *.can.part0 *.partition_files -r wrk_dir_%d %s *.%dx.fasta" % (Folder, id, file, coverage))
    r, o, e = exe(cmd)

def Merge_data(Folder):
	logging.info("Merging contigs.")
	cmd = ("cd %s && cat *.contigs.fasta > preContigs.fa" % (Folder))
	r, o, e = exe(cmd)

def single_pipe(i, output, coverage, fsize, Folder, flag, round_n):
    try:
        run_mecat(i, output, coverage, fsize, Folder, flag, round_n)
    except:
        pass

def multi_run_wrapper(args):
   return single_pipe(*args)

def local_assembly(Folder, output, coverage, flag, threads):
    files = os.listdir(Folder)
    flist = list()
    for file in files:
        if len(file) == 0:
            continue
        if file.split('.')[-1] == 'fa':
            flist.append(file[:-3])

    if flag == 1:
        tag = "Unmapped"
    else:
        tag = "Signaling"
    logging.info("%d cluster reads of %s need to be assembled."%(len(flist), tag))
    round_n = 0
    analysis_pools = Pool(processes=int(threads))

    for i in flist:
        round_n += 1
        fsize = os.path.getsize(Folder + i + ".fa")
        para = [(i, output, coverage, fsize, Folder, flag, round_n)]
        analysis_pools.map_async(multi_run_wrapper, para)
        # try:
        #     run_mecat(i, output, coverage, fsize, Folder, flag)
        # except:
        #     pass
        # if round_n % 10 == 0:
        #     logging.info("%d assembly have been constructed." % (round_n))
    analysis_pools.close()
    analysis_pools.join()

def run_minimap2(output, threads):
    logging.info("Runing minimap2.")
    cmd = ("minimap2 -x ava-pb %s %s -t %d > %s" % (output + "Contigs.fa", output + "unmapped_Contigs.fa",threads, output+"overlap_m_um.paf"))
    r, o, e = exe(cmd)

def run(argv):
    args = parseArgs(argv)
    setupLogging(False)
    starttime = time.time()

    if args.folder[-1] == '/':
        folder = args.folder
    else:
        folder = args.folder + '/'

    if args.Folder[-1] == '/':
        Folder = args.Folder
    else:
        Folder = args.Folder + '/'

    # unmapped reads assembly
    local_assembly(Folder, args.output, args.coverage, 1, args.threads)
    Merge_data(Folder)
    trans_contig_name(Folder+"preContigs.fa", args.output+"unmapped_Contigs.fa")

    # signaling reads assembly
    local_assembly(folder, args.output, args.coverage, 0, args.threads)            
    Merge_data(folder)
    trans_contig_name(folder+"preContigs.fa", args.output+"Contigs.fa")


    # # generate super-contigs from signal & unmapped contigs
    # run_minimap2(args.output, args.threads)
    # # output super-contigs
    

    logging.info("Finished in %0.2f seconds."%(time.time() - starttime))
    # filecombination("contigs/", args.prefix_output)
    # os.system("rm -r contigs")

if __name__ == '__main__':
    # os.system("mkdir contigs")
    run(sys.argv[1:])

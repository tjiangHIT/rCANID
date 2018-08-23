#! /usr/bin/env python

"""
extract_reads.py
edit by tjiang
2018.8
"""

import pysam, sys
from string import maketrans
import logging

revComp = maketrans("ATCGNatcgn","TAGCNtagcn")

def get_names(names):
    with open(names, 'r') as infile:
        n = infile.read().splitlines()
    if '' in n:
        n.remove('')
    return n


def extract_reads(options):
    n = get_names(options.names)

    bamfile = pysam.AlignmentFile(options.bam, 'rb')
    name_indexed = pysam.IndexedReads(bamfile)
    name_indexed.build()
    header = bamfile.header.copy()

    # out = pysam.Samfile(options.out, 'wb', header=header)

    for name in n:
        try:
            name_indexed.find(name)
        except KeyError:
            pass
        else:
            iterator = name_indexed.find(name)
            for read in iterator:
                # out.write(x)
                if read.is_reverse:
                    if read.qual:
                        sys.stdout.write("@{0}\n{1}\n+\n{2}\n".format(read.qname, read.seq.translate(revComp)[::-1], read.qual[::-1]))
                    else:
                        sys.stdout.write(">{0}\n{1}\n".format(read.qname, read.seq.translate(revComp)[::-1]))
                else:
                    if read.qual:
                        sys.stdout.write("@{0}\n{1}\n+\n{2}\n".format(read.qname, read.seq, read.qual))
                    else:
                        sys.stdout.write(">{0}\n{1}\n".format(read.qname, read.seq))

    # sys.stdout.write("@{0}\n{1}\n+\n{2}\n".format(read.qname, read.seq, read.qual))

# def extract_reads_to_file(id_list, path):


def parse_cluster(path, bam_path, out_path, opt):
    # load the bam file
    logging.info("Loading the Bam file.")
    bamfile = pysam.AlignmentFile(bam_path, 'rb')
    name_indexed = pysam.IndexedReads(bamfile)
    name_indexed.build()

    # load cluster info
    logging.info("Loading the cluster file.")
    num = 0
    file = open(path, 'r')
    for line in file:
        num += 1
        if num % 100 == 0:
            logging.info("Finished %d clusters."%(num))
        seq = line.strip('\n').split('\t')
        chr = seq[0]
        breakpoint = seq[1]+'_'+seq[2]+'_'+seq[3]+'_'+str(len(seq[4:]))
        id_list = seq[4:]

        if len(id_list) < 5:
            continue

        if opt == "fq":
            file_path = "%s%s_%s.fq"%(out_path, chr, breakpoint)
        else:
            file_path = "%s%s_%s.fa"%(out_path, chr, breakpoint)
        out_file = open(file_path, 'w')

        for name in id_list:
            try:
                name_indexed.find(name)
            except KeyError:
                pass
            else:
                iterator = name_indexed.find(name)
                for read in iterator:
                    if read.is_reverse:
                        if opt == 'fq':
                            out_file.write("@{0}\n{1}\n+\n{2}\n".format(read.qname, read.seq.translate(revComp)[::-1], read.qual[::-1]))
                        else:
                            out_file.write(">{0}\n{1}\n".format(read.qname, read.seq.translate(revComp)[::-1]))
                            # out_file.write("@{0}\n{1}\n+\n{2}\n".format(read.qname, read.seq.translate(revComp)[::-1], read.qual[::-1]))
                    else:
                        if opt == 'fq':
                            out_file.write("@{0}\n{1}\n+\n{2}\n".format(read.qname, read.seq, read.qual))
                        else:
                            out_file.write(">{0}\n{1}\n".format(read.qname, read.seq))
        out_file.close()

    file.close()

if __name__ == "__main__":
    # from argparse import ArgumentParser

    # parser = ArgumentParser(description='Extract reads by read name from bam file')
    # parser.add_argument('-b', '--bam', help='bam file', required=True)
    # parser.add_argument('-n', '--names', help='list of read names to extract', required=True)
    # parser.add_argument('-o', '--out', help='file name for extracted alignments', required=True)   
    # options = parser.parse_args()
    # extract_reads(options)
    parse_cluster(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    
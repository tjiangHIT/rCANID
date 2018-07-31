#!/usr/bin/env python

''' 
 * All rights Reserved, Designed By HIT-Bioinformatics   
 * @Title:  simulate donor genome with novel sequence
 * @Package: sys, random, Bio
 * @Description: 
 * @author: tjiang
 * @date: Jul 31 2018
 * @version V1.0     
'''
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import random

def load_ref(ref_g):
	return SeqIO.to_dict(SeqIO.parse(ref_g, "fasta"))

def establish_fake_genome(ref_path, data, out_path):
	print("[INFO]: Loading reference genome ...")
	hg38 = load_ref(ref_path)

	Locus_pos = random.sample(xrange(len(hg38['chr1'])), 150)
	Locus_pos.sort()

	fake_genome = str()
	Answer = list()

	count = 0
	pointer = 0
	offect = 0

	for pos in Locus_pos:
		fake_genome += str(hg38['chr1'].seq[pointer:pos])
		pointer = pos
		fake_genome += data[count]
		# Answer.append(["chr1", pos+offect, len(data[count]), data[count]])
		Answer.append("%s\t%d\t%d\t%s\n"%("chr1", pos+offect, len(data[count]), data[count]))
		offect += len(data[count])
		count += 1
	fake_genome += str(hg38['chr1'].seq[pointer:])

	print("[INFO]: Construct fake genome.")
	_convert_local_fake_genome_ = SeqIO.SeqRecord(seq = Seq(fake_genome), id = '', name = 'chr1', description = 'chr1')
	# # _convert_local_fake_genome_.seq = Seq("".join(_local_fake_genome_))
	# # SeqIO.write(_convert_local_fake_genome_, pre_out + key + ".fasta", "fasta")
	# fake_genome.append(_convert_local_fake_genome_)
	print("[INFO]: Write fake genome on disk.")
	SeqIO.write([_convert_local_fake_genome_], out_path + "simulation.fa", "fasta")

	# return Answer
	file = open(out_path+"answer.txt", 'w')
	for line in Answer:
		file.write(line)
	file.close()

N_num = 150

def parse_sequence(seq):
	return seq.find('N')
	
def parse_alignment(path):
	AlignmentFile = open(path, 'r')
	A_id_list = list()
	sta_num = 0
	for line in AlignmentFile:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '@':
			continue
		read_id = seq[0]
		flag = seq[1]
		sequence = seq[9]
		# print sequence
		# break
		if flag == '4' and parse_sequence(sequence) == -1:
			# if read_id not in A_id_list:
			# 	A_id_list[read_id] = sequence
			# 	sta_num += 1
			# 	if sta_num == N_num:
			# 		break
			A_id_list.append(sequence)
			sta_num += 1
			if sta_num == N_num:
				break
	AlignmentFile.close()
	return	A_id_list

def run():
    # call(sys.argv[1])
 	A_id_list = parse_alignment(sys.argv[1])
 	# for id in A_id_list:
 	# 	# print id, len(A_id_list[id])
 	# 	print A_id_list[id]
 	establish_fake_genome(sys.argv[2], A_id_list, sys.argv[3])

if __name__ == '__main__':
	run()
#!/usr/bin/env python

''' 
 * All rights Reserved, Designed By HIT-Bioinformatics   
 * @Title:  process.py
 * @Package: 
 * @Description: Control the deNSF pipeline
 * @author: tjiang
 * @date: June 11 2018
 * @version V1.0.1     
'''

import argparse
import sys
import logging
import pysam
import cigar

INS_flag = {1:'I'}
# DEL_flag = {2:'D'}
clip_flag = {4:'S', 5:'H'}
# CLIP_note = dict()

def acquire_clip_pos(deal_cigar):
	seq = list(cigar.Cigar(deal_cigar).items())
	if seq[0][1] == 'S':
		first_pos = seq[0][0]
	else:
		first_pos = 0
	if seq[-1][1] == 'S':
		last_pos = seq[-1][0]
	else:
		last_pos = 0
	bias = 0
	for i in seq:
		if i[1] == 'M' or i[1] == 'D':
			bias += i[0]
	return [first_pos, last_pos, bias]

def detect_flag(Flag):
	# Signal
	Normal_foward = 1 >> 1
	Abnormal = 1 << 2
	Reverse_complement = 1 << 4
	Supplementary_map = 1 << 11

	signal = {Abnormal: 0, Normal_foward: 1, Reverse_complement: 2, Supplementary_map:3, Reverse_complement | Supplementary_map:4}
	if Flag in signal:
		return signal[Flag]
	else:
		return 0

def organize_split_signal(chr, primary_info, Supplementary_info, total_L, low_bandary):
	overlap = list()
	for i in Supplementary_info:
		seq = i.split(',')
		local_chr = seq[0]
		local_start = int(seq[1])
		local_cigar = seq[3]
		dic_starnd = {1:'+', 2:'-', 3:'+', 4:'-'}

		# print primary_info[4]
		if dic_starnd[primary_info[4]] != seq[2]:
			continue
		if chr != local_chr:
			continue
		local_set = acquire_clip_pos(local_cigar)
		if primary_info[0] < local_start:
			if primary_info[3]+local_set[0]-total_L > low_bandary:
				overlap.append([total_L - primary_info[3], local_set[0], primary_info[1]])
		else:
			if local_set[1]+primary_info[2]-total_L > low_bandary:
				overlap.append([total_L - local_set[1], primary_info[2], local_start+local_set[2]-1])
	return overlap

def parse_read_final(read, low_bandary):
	INS_pos = list()
	process_signal = detect_flag(read.flag)
	if process_signal == 0:
		return INS_pos

	# read.query_name
	qname = read.query_name
	# if int(qname.split('_')[2].split('=')[1]) < 10:
	# if qname.split('_')[5][-1] == 'm':
	# 	# return INS_pos
	# 	evidence_tag = 'w'
	# else:
	# 	evidence_tag = 's'

	evidence_tag = ''

	pos_start = read.reference_start
	shift = 0
	_shift_read_ = 0
	pos_end = read.reference_end
	primary_clip_0 = 0
	primary_clip_1 = 0

	for element in read.cigar:
		if element[0] == 0 or element[0] == 2:
			shift += element[1]
		if element[0] != 2:
			_shift_read_ += element[1]
		if element[0] in INS_flag and element[1] > low_bandary:
			shift += 1
			NSI_contig = read.query_sequence[_shift_read_ - element[1]:_shift_read_]
			# chr_name, breakpoint, insert_len, insert_seq
			INS_pos.append([read.reference_name, pos_start + shift, element[1], NSI_contig, evidence_tag])

		if element[0] in clip_flag:
			if shift == 0:
				primary_clip_0 = element[1]
			else:
				primary_clip_1 = element[1]

	if process_signal != 0:
		Tags = read.get_tags()
		chr = read.reference_name
		primary_info = [pos_start, pos_end, primary_clip_0, primary_clip_1, process_signal]

		for i in Tags:
			if i[0] == 'SA':
				Supplementary_info = i[1].split(';')[:-1]
				overlap = organize_split_signal(chr, primary_info, Supplementary_info, read.query_length, low_bandary)
				for k in overlap:
					NSI_contig = read.query_sequence[k[0]:k[1]]
					# [the breakpoint on reference]
					# [the insertion size]
					# [the NSI] 
					INS_pos.append([chr, k[2], k[1] - k[0], NSI_contig, evidence_tag])

	# uniq_list = dict()
	# for i in INS_pos:
	# 	key = "%s%d%d"%(i[0], i[1], i[2])
	# 	if key not in uniq_list:
	# 		uniq_list[key] = 0
	# 	else:
	# 		i[4] = 0

	return INS_pos

def load_final_alignment(sam_path, out_path):
	samfile = pysam.AlignmentFile(sam_path)
	outfile = open(out_path, 'w')

	for read in samfile.fetch():
		feed_back = parse_read_final(read, 50)
		if len(feed_back) > 0:
			for i in feed_back:
				if i[4] == 0:
					continue
				# outfile.write("%s\t%d\t%d\t%s\n"%(i[0], i[1], i[2], i[3])) 
				outfile.write(">%s_%d_%d_%s\n"%(i[0], i[1], i[2], i[4]))
				outfile.write("%s\n"%(i[3]))
	samfile.close()
	outfile.close()

if __name__ == '__main__':
	load_final_alignment(sys.argv[1], sys.argv[2])

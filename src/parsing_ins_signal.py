#!/usr/bin/env python

''' 
 * All rights Reserved, Designed By HIT-Bioinformatics   
 * @Title:  process.py
 * @Package: 
 * @Description: parse signaling read
 * @author: tjiang
 * @date: June 11 2018
 * @version V1.0.1     
'''

import argparse
import sys
import logging
import pysam
import cigar
from multiprocessing import Pool
from CommandRunner import *


INS_flag = {1:'I'}
DEL_flag = {2:'D'}
clip_flag = {4:'S', 5:'H'}
CLIP_note = dict()

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

def organize_split_signal(chr, primary_info, Supplementary_info, total_L, low_bandary):
	overlap = list()
	for i in Supplementary_info:
		seq = i.split(',')
		local_chr = seq[0]
		local_start = int(seq[1])
		local_cigar = seq[3]

		dic_starnd = {1:'+', 2: '-'}

		if dic_starnd[primary_info[4]] != seq[2]:
			continue
		if chr != local_chr:
			# return overlap
			continue
		local_set = acquire_clip_pos(local_cigar)
		# if len(local_set) == 0:
		# 	continue
		if primary_info[0] < local_start:
			if primary_info[3]+local_set[0]-total_L > low_bandary:
				overlap.append([total_L - primary_info[3], local_set[0], primary_info[1]])
		else:
			if local_set[1]+primary_info[2]-total_L > low_bandary:
				overlap.append([total_L - local_set[1], primary_info[2], local_start+local_set[2]-1])
			# exist some bugs

		# if local_start <= primary_end + 50 and local_start >= primary_end - 50:
		# 	local_set = acquire_clip_pos(local_cigar)
		# 	if len(local_set) != 0:
		# 		overlap.append([primary_clip, sum(local_set) - primary_clip])
	return overlap

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

def store_clip_pos(locus, chr, flag, name):
	# about collecting breakpoint from clipping 
	hash_1 = int(locus /10000)
	mod = locus % 10000
	hash_2 = int(mod / 50)
	element = [locus, flag, name]

	if hash_1 not in CLIP_note[chr]:
		CLIP_note[chr][hash_1] = dict()
		CLIP_note[chr][hash_1][hash_2] = list()
		CLIP_note[chr][hash_1][hash_2].append(element)
	else:
		if hash_2 not in CLIP_note[chr][hash_1]:
			CLIP_note[chr][hash_1][hash_2] = list()
			CLIP_note[chr][hash_1][hash_2].append(element)
		else:
			CLIP_note[chr][hash_1][hash_2].append(element)

def parse_read(read, Chr_name, low_bandary):
	'''
	Check:	1.Flag
			2.Supplementary mapping
			3.Seq
	'''

	# print read.query_name

	# INS_ME_pos = list()

	# Nonmapped_reads

	process_signal = detect_flag(read.flag)
	# if process_signal == 0:
		# return INS_ME_pos

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
			# MEI_contig = read.query_sequence[_shift_read_ - element[1]:_shift_read_]

			# INS_ME_pos.append([pos_start + shift, element[1], MEI_contig])

			# [the breakpoint on reference]
			# [the insertion size]
			# [the read name] 
			# INS_ME_pos.append([pos_start + shift, element[1], read.query_name])
			store_clip_pos(pos_start + shift, Chr_name, 2, read.query_name)

			# print read.query_name, "I", pos_start + shift
			# print MEI_contig

		if element[0] in clip_flag:

			if shift == 0:
				primary_clip_0 = element[1]
			else:
				primary_clip_1 = element[1]

			if element[1] > low_bandary:
				if shift == 0:
					clip_pos = pos_start - 1
					# clip_contig = read.query_sequence[:element[1]]
					store_clip_pos(clip_pos, Chr_name, 0, read.query_name)

					# primary_clip_0 = element[1]
					# left clip size

				else:
					clip_pos = pos_start + shift - 1
					# primary_clip = read.query_length - element[1]
					# clip_contig = read.query_sequence[read.query_length - element[1]:]
					store_clip_pos(clip_pos, Chr_name, 1, read.query_name)

					# primary_clip_1 = read.query_length - element[1]
					# right clip size

	if process_signal == 1 or process_signal == 2:
		Tags = read.get_tags()
		chr = Chr_name
		# primary_clip = pos_start
		primary_info = [pos_start, pos_end, primary_clip_0, primary_clip_1, process_signal]

		for i in Tags:
			if i[0] == 'SA':
				Supplementary_info = i[1].split(';')[:-1]
				# print process_signal
				# print chr, primary_info, read.query_length
				# print read.cigar
				# print i[1].split(';')[-1]
				overlap = organize_split_signal(chr, primary_info, Supplementary_info, read.query_length, low_bandary)
				for k in overlap:
					# print k
					# MEI_contig = read.query_sequence[k[0]:k[1]]
					# [the breakpoint on reference]
					# [the insertion size]
					# [the read name] 
					# INS_ME_pos.append([k[2], k[1] - k[0], read.query_name])
					store_clip_pos(k[2], Chr_name, 2, read.query_name)

	# return INS_ME_pos

def acquire_clip_locus(down, up, chr):
	list_clip = list()
	if int(up/10000) == int(down/10000):
		key_1 = int(down/10000)
		if key_1 not in CLIP_note[chr]:
			return list_clip
		for i in xrange(int((up%10000)/50)-int((down%10000)/50)+1):
			# exist a bug ***********************************
			key_2 = int((down%10000)/50)+i
			if key_2 not in CLIP_note[chr][key_1]:
				continue
			for ele in CLIP_note[chr][key_1][key_2]:
				if ele[0] >= down and ele[0] <= up:
					list_clip.append(ele)
	else:
		key_1 = int(down/10000)
		if key_1 in CLIP_note[chr]:
			for i in xrange(200-int((down%10000)/50)):
				# exist a bug ***********************************
				key_2 = int((down%10000)/50)+i
				if key_2 not in CLIP_note[chr][key_1]:
					continue
				for ele in CLIP_note[chr][key_1][key_2]:
					if ele[0] >= down and ele[0] <= up:
						list_clip.append(ele)
		key_1 += 1
		if key_1 not in CLIP_note[chr]:
			return list_clip
		for i in xrange(int((up%10000)/50)+1):
			# exist a bug ***********************************
			key_2 = i
			if key_2 not in CLIP_note[chr][key_1]:
				continue
			for ele in CLIP_note[chr][key_1][key_2]:
				if ele[0] >= down and ele[0] <= up:
					list_clip.append(ele)
	return list_clip

def construct_concensus_info(Ins_list, Clip_list, evidence_read, SV_size):
	# total_count = len(Ins_list) + len(Clip_list)
	unique_read = dict()
	for i in Ins_list:
		if i[2] not in unique_read:
			unique_read[i[2]] = 0
	for i in Clip_list:
		if i[2] not in unique_read:
			unique_read[i[2]] = 0
	total_count = len(unique_read)

	if total_count < evidence_read:
		return 0
	breakpoint = list()
	insert_size = list()

	boundary = list()

	for i in Ins_list:
		breakpoint.append(i[0])
		insert_size.append(i[1])
		# boundary.append(i[0])
	for i in Clip_list:
		# boundary.append(i[0])
		if i[1] == 1:
			breakpoint.append(i[0])

	# ==============method_1=====================
	# Prob_pos_1 = Counter(breakpoint).most_common(1)[0][0]
	# ==============method_2=====================
	Prob_pos_2 = sum(breakpoint)/len(breakpoint)
	Average_size = int(sum(insert_size)/len(insert_size))

	if Average_size < SV_size:
		return 0

	# # print Average_size
	# local_info = list()
	# # local_name = "_%d_%d_"%(Prob_pos_2, Average_size)
	# local_name = [Prob_pos_2, Average_size]
	# local_id = 0
	# for i in Ins_list:
	# 	# info = local_name + str(local_id) + '\n' + i[2] + '\n'
	# 	info = local_name + [str(local_id), i[2]]
	# 	# info = local_name + ["%d:%d:%d"%(local_id, min(boundary), max(boundary)), i[2]]
	# 	local_id += 1
	# 	local_info.append(info)
	# for i in Clip_list:
	# 	# info = local_name + str(local_id) + '\n' + i[1] + '\n'
	# 	info = local_name + [str(local_id), i[1]]
	# 	# info = local_name + ["%d:%d:%d"%(local_id, min(boundary), max(boundary)), i[1]]
	# 	local_id += 1
	# 	local_info.append(info)
	# 		# print(">%d\n%s"%(i[0], i[2])) 

	# # for i in local_info:
	# # 	print i
	read_list = list()
	for key in unique_read:
		read_list.append(key)
	return [Prob_pos_2, Average_size, total_count, read_list]

def merge_pos(pos_list, chr, evidence_read, SV_size):
	start = list()
	end = list()
	for ele in pos_list:
		start.append(ele[0])
		end.append(ele[0] + ele[1])

	search_down = min(start) - 10
	search_up = max(start) + 10

	temp_clip = acquire_clip_locus(search_down, search_up, chr)
	# temp_clip = list()
	# gc.collect()

	# concensus, ref_pos = construct_concensus_seq(pos_list, temp_clip)
	result = construct_concensus_info(pos_list, temp_clip, evidence_read, SV_size)

	if result != 0:
		# for i in xrange(len(result)):
		# 	# result[i] = ">INS_" + chr + result[i]
		# 	result[i] = ["INS", chr] + result[i] + [len(result)]
		return result
	else:
		return 0

def cluster(pos_list, chr, evidence_read, SV_size, low_bandary):
	_cluster_ = list()
	temp = list()
	temp.append(pos_list[0])

	result = single_clip(0, temp[0][0], chr)

	for pos in pos_list[1:]:
		# if temp[-1][0] + temp[-1][1] < pos[0]:
		if temp[-1][0] + low_bandary < pos[0]:
			result = merge_pos(temp, chr, evidence_read, SV_size)
			if result != 0:
				_cluster_.append(result)
			temp = list()
			temp.append(pos)
		else:
			temp.append(pos)
	result = merge_pos(temp, chr, evidence_read, SV_size)
	if result != 0:
		_cluster_.append(result)
	# _cluster_.append(merge_pos(temp, chr))
	return _cluster_

def single_clip(chr, Chr_length, read_evidence):
	Clip_list = list()
	# temp_breakpoint_s = 0
	# temp_breakpoint_e = 0
	# read_list = list()
	read_list = dict()
	start_list = list()
	strength = 0
	for i in xrange(int(Chr_length/50)):
		temp_clip = acquire_clip_locus(50*i, 50*(i+1), chr)
		if len(temp_clip) > 0:
			# read_list = list()
			for ele in temp_clip:
				start_list.append(ele[0])
				if ele[1] == 2:
					strength += 1
				# read_list.append(ele[2])
				if ele[2] not in read_list:
					read_list[ele[2]] = 0
			# if temp_breakpoint_s != 0:
			# 	temp_breakpoint_s = min(start_list)
			# temp_breakpoint_e = max(start_list)
		else:
			if len(start_list) == 0:
				continue
			else:
				if len(read_list) >= read_evidence:
					unique_read = list()
					for key in read_list:
						unique_read.append(key)
					Clip_list.append([chr, min(start_list), max(start_list), strength, unique_read])
				start_list = list()
				read_list = dict()
				strength = 0
	return Clip_list

def single_pipe(sam_path, out_path, Chr_name):
	samfile = pysam.AlignmentFile(sam_path)
	Chr_length = samfile.get_reference_length(Chr_name)
	# print Chr_name, Chr_length
	logging.info("Resolving the chromsome %s."%(Chr_name))
	if Chr_name not in CLIP_note:
		CLIP_note[Chr_name] = dict()
	cluster_pos_INS = list()
	for read in samfile.fetch(Chr_name):
		parse_read(read, Chr_name, 50)
	Cluster_INS = single_clip(Chr_name, Chr_length, 10)
	logging.info("%d INS signal loci in the chromsome %s."%(len(Cluster_INS), Chr_name))

	local_path = "%s_%s.txt"%(out_path, Chr_name)
	file = open(local_path, 'w')
	for line in Cluster_INS:
		# breakpoint = line[0]
		# sv_size = line[1]
		# evidence_num = line[2]
		read_list = "\t".join(line[4])
		file.write("%s\t%d\t%d\t%d\t%s\n"%(line[0], line[1], line[2], line[3],read_list))
	file.close()
	samfile.close()

def multi_run_wrapper(args):
   return single_pipe(*args)

def load_sam(sam_path, out_path, threads):
	'''
	Load_BAM_File
	library:	pysam.AlignmentFile

	load_Ref_Genome
	library:	Bio
	'''
	# p1 = args.AlignmentFile
	# p2 = args.Output_prefix
	# p3 = args.Reference
	# Ref = load_ref(p3)
	setupLogging(False)
	samfile = pysam.AlignmentFile(sam_path)
	# print(samfile.get_index_statistics())
	contig_num = len(samfile.get_index_statistics())
	# print contig_num
	logging.info("The total number of chromsomes: %d"%(contig_num))
	# print("The total number of chromsomes: %d"%(contig_num))
	# Acquire_Chr_name

	process_list = list()
	for i in samfile.get_index_statistics():
		process_list.append([i[0], i[3]])
		# #chr #read
	process_list = sorted(process_list, key = lambda x:-x[1])
	analysis_pools = Pool(processes=int(threads))

	for i in process_list:
		para = [(sam_path, out_path, i[0])]
		analysis_pools.map_async(multi_run_wrapper, para)
	analysis_pools.close()
	analysis_pools.join()

	# for _num_ in xrange(contig_num):
	# 	Chr_name = samfile.get_reference_name(_num_)
	# 	logging.info("Resolving the chromsome %s."%(Chr_name))
	# 	# print("Resolving the chromsome %s."%(Chr_name))
	# 	Chr_length = samfile.lengths[_num_]
	# 	# print Chr_length
	# 	if Chr_name not in CLIP_note:
	# 		# CLIP_note[Chr_name] = [0] * Chr_length
	# 		# CLIP_note[Chr_name] = Q.PriorityQueue()
	# 		CLIP_note[Chr_name] = dict()

	# 	cluster_pos_INS = list()
	# 	# cluster_pos_DEL = list()
	# 	for read in samfile.fetch(Chr_name):
	# 		# feed_back = parse_read(read, Chr_name, 50)
	# 		parse_read(read, Chr_name, 50)

	# 	# 	if len(feed_back) > 0:
	# 	# 		# print read.query_name
	# 	# 		for i in feed_back:
	# 	# 			# print i
	# 	# 			cluster_pos_INS.append(i)
	# 	# 		# break
	# 	# 	# if len(feed_back_del) > 0:
	# 	# 	# 	for i in feed_back_del:
	# 	# 	# 		cluster_pos_DEL.append(i)
	# 	# # while not CLIP_note[Chr_name].empty():
	# 	# # 	print Chr_name, CLIP_note[Chr_name].get()
	# 	# # print CLIP_note[Chr_name][6]
	# 	# cluster_pos_INS = sorted(cluster_pos_INS, key = lambda x:x[0])
	# 	# # cluster_pos_DEL = sorted(cluster_pos_DEL, key = lambda x:x[0])
	# 	# if len(cluster_pos_INS) == 0:
	# 	# 	Cluster_INS = list()
	# 	# else:
	# 	# 	# for i in cluster_pos_INS:
	# 	# 	# 	print i
	# 	# 	Cluster_INS = cluster(cluster_pos_INS, Chr_name, 5, 50, 30)

	# 	Cluster_INS = single_clip(Chr_name, Chr_length, 10)

	# 	logging.info("%d INS signal loci in the chromsome %s."%(len(Cluster_INS), Chr_name))		

	# 	local_path = "%s_%s.txt"%(out_path, Chr_name)
	# 	file = open(local_path, 'w')
	# 	for line in Cluster_INS:
	# 		# breakpoint = line[0]
	# 		# sv_size = line[1]
	# 		# evidence_num = line[2]
	# 		read_list = "\t".join(line[4])
	# 		file.write("%s\t%d\t%d\t%d\t%s\n"%(line[0], line[1], line[2], line[3],read_list))
	# 	file.close()

	# 	# logging.info("%d INS signal loci in the chromsome %s."%(len(Cluster_INS), Chr_name))
	samfile.close()

if __name__ == '__main__':
	load_sam(sys.argv[1], sys.argv[2], sys.argv[3])

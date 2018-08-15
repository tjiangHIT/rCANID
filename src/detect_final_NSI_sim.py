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

import sys
import cigar

def analysis_cigar(deal_cigar, ins_l):
	seq = list(cigar.Cigar(deal_cigar).items())
	SoftClip_len = 0

	if seq[0][1] == 'S':
		SoftClip_len += seq[0][0]

	if seq[-1][1] == 'S':
		SoftClip_len += seq[-1][0]

	if SoftClip_len * 2 > ins_l:
		return 0
	else:
		return 1

def judgement(flag, MAPQ, cigar, NSI_len, evidence_tag):
	if flag == '4':
		return 1
		# NSI
	else:
		if MAPQ < 10:
			return 0
			# invaild insertion
		else:
			info = analysis_cigar(cigar, NSI_len)
			if info	== 0:
				return 1
				# NSI
			else:
				return 2 
				# normal insertion

def final_call(sam_path, out_path):
	callset = dict()
	file = open(sam_path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '@':
			continue

		chr = seq[0].split('_')[0]
		breakpoint = int(seq[0].split('_')[1])
		NSI_len = int(seq[0].split('_')[2])
		flag = seq[1]
		NSI = seq[9]
		MAPQ = int(seq[4])
		cigar = seq[5]
		evidence_tag = seq[0].split('_')[3]
		local_key = "%d_%d"%(breakpoint, NSI_len)

		if chr not in callset:
			callset[chr] = dict()
			# callset[chr][breakpoint] = list()
			info = judgement(flag, MAPQ, cigar, NSI_len, evidence_tag)
			if info	!= 0:
				callset[chr][local_key] = list()
				callset[chr][local_key].append([NSI_len, NSI, info, evidence_tag])
		else:
			# if breakpoint not in callset[chr]:
				# callset[chr][breakpoint] = list()
			info = judgement(flag, MAPQ, cigar, NSI_len, evidence_tag)
			if info	!= 0:
				if local_key not in callset[chr]:
					callset[chr][local_key] = list()
				callset[chr][local_key].append([NSI_len, NSI, info, evidence_tag])
		# if flag != '4':
		# 	# continue
		# 	if analysis_cigar(cigar, NSI_len) == 1:
		# 		continue
		# if chr not in callset:
		# 	callset[chr] = list()
		# callset[chr].append([breakpoint, NSI_len, NSI, evidence_tag])
		# # if chr not in callset:
		# # 	callset[chr] = dict()
		# # 	callset[chr][breakpoint] = [NSI_len, NSI]
		# # else:
		# # 	if breakpoint not in callset[chr]:
		# # 		callset[chr][breakpoint] = [NSI_len, NSI]
		# # 	else:
		# # 		if NSI_len < callset[chr][breakpoint][0]:
		# # 			callset[chr][breakpoint] = [NSI_len, NSI]
	file.close()

	temp_calling = dict()
	for chr in callset:
		for breakpoint in callset[chr]:
			signal = 0
			max_ins = 0
			min_ins = sys.maxint
			for ele in callset[chr][breakpoint]:
				# print breakpoint, ele[0]
				if ele[0] > max_ins:
					max_ins = ele[0]
				if ele[0] < min_ins:
					min_ins = ele[0]
				if ele[2] == 2:
					signal = 1
			if signal != 1:
				if (max_ins - min_ins)*20 > (max_ins + min_ins):
					continue
				for ele in callset[chr][breakpoint]:
					if chr not in  temp_calling:
						temp_calling[chr] = list()
					temp_calling[chr].append([int(breakpoint.split('_')[0]), ele[0], ele[1], ele[3]])


	file = open(out_path, 'w')
	for chr in temp_calling:
		temp_calling[chr] = sorted(temp_calling[chr], key=lambda x:x[0])

		# for i in temp_calling[chr]:
		# 	print i[0],i[1],i[3]

		temp = list()
		temp.append(temp_calling[chr][0])
		for ele in temp_calling[chr][1:]:
			if temp[-1][0] + 20 < ele[0]:
				temp = sorted(temp, key=lambda x:x[1])
				if len(temp) > 3:
					# file.write("%s\t%d\t%d\t%s\n"%(chr, temp[len(temp)/2][0], temp[len(temp)/2][1], temp[len(temp)/2][2]))
					file.write("%s\t%d\t%d\t%s\n"%(chr, temp[0][0], temp[0][1], temp[0][2]))
				else:
					for i in temp:
						if i[3] == 's':
							# file.write("%s\t%d\t%d\t%s\n"%(chr, temp[len(temp)/2][0], temp[len(temp)/2][1], temp[len(temp)/2][2]))
							file.write("%s\t%d\t%d\t%s\n"%(chr, temp[0][0], temp[0][1], temp[0][2]))
							break
				# temp = sorted(temp, key=lambda x:x[1])
				# file.write("%s\t%d\t%d\t%s\n"%(chr, temp[len(temp)/2][0], temp[len(temp)/2][1], temp[len(temp)/2][2]))		
				temp = list()
				temp.append(ele)
			else:
				temp.append(ele)
		temp = sorted(temp, key=lambda x:x[1])
		if len(temp) > 3:
			# file.write("%s\t%d\t%d\t%s\n"%(chr, temp[len(temp)/2][0], temp[len(temp)/2][1], temp[len(temp)/2][2]))
			file.write("%s\t%d\t%d\t%s\n"%(chr, temp[0][0], temp[0][1], temp[0][2]))
		else:
			for i in temp:
				if i[3] == 's':
					# file.write("%s\t%d\t%d\t%s\n"%(chr, temp[len(temp)/2][0], temp[len(temp)/2][1], temp[len(temp)/2][2]))
					file.write("%s\t%d\t%d\t%s\n"%(chr, temp[0][0], temp[0][1], temp[0][2]))
					break
		# temp = sorted(temp, key=lambda x:x[1])
		# file.write("%s\t%d\t%d\t%s\n"%(chr, temp[len(temp)/2][0], temp[len(temp)/2][1], temp[len(temp)/2][2]))
	# file = open(out_path, 'w')
	# for chr in callset:
	# 	for key in callset[chr]:
	# 		file.write("%s\t%s\t%d\t%s\n"%(chr, key, callset[chr][key][0], callset[chr][key][1]))
	file.close()

if __name__ == '__main__':
	final_call(sys.argv[1], sys.argv[2])

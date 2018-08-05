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

def final_call(sam_path, out_path):
	callset = dict()
	file = open(sam_path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '@':
			continue

		chr = seq[0].split('_')[0]
		breakpoint = seq[0].split('_')[1]
		NSI_len = int(seq[0].split('_')[2])
		flag = seq[1]
		NSI = seq[9]

		if flag != '4':
			continue

		if chr not in callset:
			callset[chr] = dict()
			callset[chr][breakpoint] = [NSI_len, NSI]
		else:
			if breakpoint not in callset[chr]:
				callset[chr][breakpoint] = [NSI_len, NSI]
			else:
				if NSI_len < callset[chr][breakpoint][0]:
					callset[chr][breakpoint] = [NSI_len, NSI]
	file.close()

	file = open(out_path, 'w')
	for chr in callset:
		for key in callset[chr]:
			file.write("%s\t%s\t%d\t%s\n"%(chr, key, callset[chr][key][0], callset[chr][key][1]))
	file.close()

if __name__ == '__main__':
	final_call(sys.argv[1], sys.argv[2])

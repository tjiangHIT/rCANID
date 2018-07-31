import sys
import copy
key_list = dict()

def collect_sublist(key):
	return key_list[key]


def cluster_reads(path):

	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		pair_1 = seq[0]
		pair_2 = seq[5]

		if pair_1 not in key_list:
			key_list[pair_1] = list()
		key_list[pair_1].append(pair_2)
	file.close()
	print len(key_list)

	cluster_dic = copy.deepcopy(key_list)
	for key in key_list:
		for ele in key_list[key]:
			if ele in key_list:



if __name__ == '__main__':
	cluster_reads(sys.argv[1])
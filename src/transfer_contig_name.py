import sys

def trans_contig_name(path):
	file = open(path, 'r')
	_id_ = 0
	for line in file:
		seq = line.strip('\n').split(' ')
		if seq[0][0] == '>':
			tag = seq[2].split('=')[1].split('.')[0]
			print(">%d_%s"%(_id_, tag))
			_id_ += 1
		# else:
		# 	print line.strip('\n')

	

if __name__ == "__main__":
	trans_contig_name(sys.argv[1])
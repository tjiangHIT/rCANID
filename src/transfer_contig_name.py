import sys

def trans_contig_name(path, opath):
	file = open(path, 'r')
	outfile = open(opath, 'w')
	_id_ = 0
	for line in file:
		seq = line.strip('\n').split(' ')
		if seq[0][0] == '>':
			# tag = seq[2].split('=')[1].split('.')[0]
			# tag = seq[1].split('=')[1]
			tag = '_'.join(seq[1:])
			# print(">%d_%s"%(_id_, tag))
			outfile.write(">%d_%s\n"%(_id_, tag))
			_id_ += 1
		else:
			# print line.strip('\n')
			outfile.write(line)
	file.close()
	outfile.close()	

def filter_contigs(path):
	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split(' ')
		if seq[0][0] == '>':
			print seq[5]

if __name__ == "__main__":
	trans_contig_name(sys.argv[1], sys.argv[2])

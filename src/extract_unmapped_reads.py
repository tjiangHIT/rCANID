import pysam, sys
import time

def load_alignment(path):
	# print path
	samfile = pysam.AlignmentFile(path)
	for read in samfile.fetch():
		# print read.query_name, read.flag
		# print read.query_sequence
		# print read.qual
		# break
		if read.flag == 4:
			print read.query_name
			print read.query_sequence
			print "+"
			print read.qual

def parse_alignment(path):
	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '@':
			continue

		# print seq[0]		

		read_name = seq[0]
		read_flag = seq[1]
		read_seq = seq[9]
		read_qual = seq[10]

		if read_flag == '4':
			print '@'+read_name
			print read_seq
			print '+'+read_name
			print read_qual
			# break
	file.close()

if __name__ == '__main__':
	# starttime = time.time()
	# load_alignment(sys.argv[1])
	parse_alignment(sys.argv[1])
	# print("[INFO]: Finished in %0.2f seconds."%(time.time() - starttime))

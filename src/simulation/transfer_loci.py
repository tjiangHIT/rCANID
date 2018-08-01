import sys

def run(path):
	file = open(path, 'r')
	offect = 0
	for line in file:
		seq = line.strip('\n').split('\t')
		print("%s\t%d\t%s\t%s"%(seq[0], int(seq[1])-offect, seq[2], seq[3]))
		offect += int(seq[2])

if __name__ == '__main__':
	run(sys.argv[1])
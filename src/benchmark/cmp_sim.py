import sys

def run(p1, p2):
	answer = list()
	right = 0
	flag = 0
	file = open(p1, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		answer.append([int(seq[1]), 0])
	file.close()

	sta = 0

	file = open(p2, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		# breakpoint = int(seq[0])
		low = int(seq[1])
		up = int(seq[2])

		if len(seq[3:]) >= 10:
			sta +=1

		for pos in answer:
			# if pos[0] - 50 <= breakpoint and breakpoint <= pos[0] + 50:
			if low <= pos[0] and pos[0] <= up and len(seq[3:]) >= 10:
			# if low <= pos[0] and pos[0] <= up:	
				right += 1
				pos[1] = 1
				print "[INFO] correct " + str(len(seq[3:]))
				# print("[INFO] correct %d\t%d\t%0.3f"%(len(seq[3:]), up-low, float(len(seq[3:])/(up-low))))
				flag = 1
			# else:
			# 	print "[INFO] wrong " + str(len(seq[3:]))
		if flag == 0:
			print "[INFO] wrong " + str(len(seq[3:]))
			# print("[INFO] wrong %d\t%d\t%0.3f"%(len(seq[3:]), up-low, float(len(seq[3:])/(up-low))))
		flag = 0
	file.close()
	print right
	print sta
	# for i in answer:
	# 	if i[1] == 0:
	# 		print i

def run_list(p1, p2):
	answer = list()
	right = 0
	flag = 0
	file = open(p1, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		answer.append([int(seq[1]), 0])
	file.close()

	file = open(p2, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		breakpoint = int(seq[1])
		for pos in answer:
			if pos[0] - 50 <= breakpoint and breakpoint <= pos[0] + 50:
				right += 1
				pos[1] = 1
		# 		flag = 1
		# if flag == 0:
	file.close()
	print right
	for i in answer:
		if i[1] == 0:
			print i[0]

def run_list_real(p1, p2):
	# answer = list()
	answer = dict()
	right = 0
	flag = 0
	file = open(p1, 'r')
	for line in file:
		seq = line.strip('\n').split(',')
		chr = seq[13][1:-1]
		local_list = list()
		for k in xrange(6):
			num = seq[14+k][1:-1]
			if len(num) > 0:
				local_list.append(int(num))
		low = min(local_list)
		up = max(local_list)
		if chr not in answer:
			answer[chr] = list()
		# print chr, low, up
		answer[chr].append([low, up, 0])
	file.close()

	file = open(p2, 'r')
	for line in file:
		seq = line.strip('\n').split('_')
		# breakpoint = int(seq[1])
		# print seq
		chr = seq[0]
		low = int(seq[1])
		up = int(seq[2])
		read_count = int(seq[3].split('.')[0])

		if chr in answer:
			for pos in answer[chr]:
				if pos[0] <= up and pos[1] >= low:
					right += 1
					pos[2] = 1

	file.close()
	print right
	for key in answer:
		for i in answer[key]:
			if i[2] == 0:
				print key, i[0], i[1]
	# for i in answer:
	# 	if i[1] == 0:
	# 		print i[0]

if __name__ == '__main__':
	run_list_real(sys.argv[1], sys.argv[2])

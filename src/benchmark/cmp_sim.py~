import sys

def load_sim_answer(path):
	answer = list()
	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		answer.append([int(seq[1]), 0, int(seq[2])])
	file.close()
	return answer

def load_real_answer(path):
	answer = dict()
	file = open(path, 'r')
	for line in file:
		# seq = line.strip('\n').split('\t')
		# answer.append([int(seq[1]), 0, int(seq[2])])
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
	return answer

def run(p1, p2):
	answer = load_sim_answer(p1)
	right = 0
	flag = 0
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
	answer = load_sim_answer(p1)
	right = 0
	flag = 0
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
	answer = load_real_answer(p1)
	right = 0
	flag = 0
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

def run_real_rCANID(p1, p2):
	# answer = list()
	answer = load_real_answer(p1)
	right = 0
	flag = 0
	totally_call = 0
	file = open(p2, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		# breakpoint = int(seq[1])
		# print seq
		chr = seq[0]
		if chr != '1':
			continue
		breakpoint = int(seq[1])
		LEN = seq[2]
		# read_count = int(seq[3].split('.')[0])
		# print chr, breakpoint
		totally_call += 1
		if chr in answer:
			for pos in answer[chr]:
				# print chr, pos
				if pos[0] - 50 <= breakpoint and breakpoint <= pos[1]+50:
				# if pos[0] - 50  <= breakpoint and pos[1] + 50 >= breakpoint:
					# right += 1
					print chr, pos[0], pos[1], LEN
					pos[2] = 1
		# break
	file.close()
	totally_answer = 0
	for key in answer:
		totally_answer += len(answer[key])
		for i in answer[key]:
			if i[2] == 1:
				right += 1
			else:
				print key, i
	print("[INFO] The total nember of answers is %d."%(totally_answer))
	print("[INFO] rCANID calls %d SVs."%(totally_call))
	print("[INFO] Correct calling is %d."%(right))

def run_sim_rCANID(p1, p2):
	answer = load_sim_answer(p1)
	right = 0
	file = open(p2, 'r')
	totally_call = 0
	for line in file:
		seq = line.strip('\n').split('\t')
		# if seq[0][0] == '#':
		# 	continue
		flag = 0
		totally_call += 1
		breakpoint = int(seq[1])
		svlength = int(seq[2])

		for pos in answer:
			if pos[0] - 50 <= breakpoint and breakpoint <= pos[0] + 50:
				# right += 1
				# pos[1] = 1
				# flag = 1
				if abs(svlength - pos[2])*10 <= pos[2]:
					pos[1] = 1
					flag = 1
		if flag == 0:
			print breakpoint, svlength
	file.close()
	# print right
	for i in answer:
		if i[1] == 1:
			right += 1
		else:
			print i
	print("[INFO] The total nember of answers is %d."%(len(answer)))
	print("[INFO] rCANID calls %d SVs."%(totally_call))
	print("[INFO] Correct calling is %d."%(right))


def run_sim_sniffles(p1, p2):
	answer = load_sim_answer(p1)
	right = 0
	flag = 0
	file = open(p2, 'r')
	totally_call = 0
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue
		totally_call += 1
		breakpoint = int(seq[13])
		svlength = int(seq[16])
		svtype = seq[10]
		if svtype != "INS":
			continue
		for pos in answer:
			if pos[0] - 50 <= breakpoint and breakpoint <= pos[0] + 50:
				# right += 1
				# pos[1] = 1
		# 		flag = 1
				if abs(svlength - pos[2])*10 <= pos[2]:
					pos[1] = 1
		# if flag == 0:
	file.close()
	# print right
	for i in answer:
		if i[1] == 1:
			right += 1
		else:
			print i
	print("[INFO] The total nember of answers is %d."%(len(answer)))
	print("[INFO] Sniffles calls %d SVs."%(totally_call))
	print("[INFO] Correct calling is %d."%(right))

def run_real_sniffles(p1, p2):
	answer = load_real_answer(p1)
	right = 0
	flag = 0
	file = open(p2, 'r')
	totally_call = 0
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue
		chr = seq[0]
		totally_call += 1
		breakpoint = int(seq[13])
		# svlength = int(seq[16])
		svtype = seq[10]
		LEN = seq[16]
		if svtype != "INS":
			continue

		if chr in answer:
			for pos in answer[chr]:
				if pos[0] - 50 <= breakpoint and breakpoint <= pos[1] + 50:
					# right += 1
					pos[2] = 1
					print chr, pos[0], pos[1], LEN
		# 		flag = 1
		# if flag == 0:
	file.close()
	# print right
	totally_answer = 0
	for key in answer:
		totally_answer += len(answer[key])
		for i in answer[key]:
			if i[2] == 1:
				right += 1
	# print("[INFO] The total nember of answers is %d."%(totally_answer))
	# print("[INFO] Sniffles calls %d SVs."%(totally_call))
	# print("[INFO] Correct calling is %d."%(right))

def cmp_sniffles_rCANID(p1, p2):
	file = open(p1, 'r')
	sniffles = dict()
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue
		chr = seq[0]
		breakpoint = int(seq[13])
		svtype = seq[10]
		if svtype != 'INS':
			continue
		if chr not in sniffles:
			sniffles[chr] = list()
		sniffles[chr].append(breakpoint)
	file.close()

	file = open(p2, 'r')
	num = 0
	for line in file:
		seq = line.strip('\n').split('\t')
		chr = seq[0]
		breakpoint = int(seq[1])
		# print chr, breakpoint
		flag = 0
		if chr in sniffles:
			for ele in sniffles[chr]:
				if ele - 50 <= breakpoint and breakpoint <= ele + 50:
					num += 1
					flag = 1
			# if flag == 0:
			# 	print chr, breakpoint
	file.close()
	print("[INFO] total number is %d"%num)

if __name__ == '__main__':
	run_real_rCANID(sys.argv[1], sys.argv[2])

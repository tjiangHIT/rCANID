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

if __name__ == '__main__':
	run(sys.argv[1], sys.argv[2])
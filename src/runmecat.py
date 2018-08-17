# coding=utf-8   
import os          

filetype = -1
file = os.listdir(".")
flist = []
for i in file:
	if i.split('.')[1] == 'fq':
		flist.append(i.split(".")[0])
		filetype = 0
	elif i.split('.')[1] == 'fastq':
		flist.append(i.split(".")[0])
		filetype = 1
	elif i.split('.')[1] == 'fa':
		flist.append(i.split(".")[0])
		filetype = 2
	elif i.split('.')[1] == 'fasta':
		flist.append(i.split(".")[0])
		filetype = 3

if filetype == 0:
	for i in flist:
		fsize = os.path.getsize(i+".fq")
		os.system("bash run_MECAT.sh " + i + " " + str(fsize) + " fq")
elif filetype == 0:
	for i in flist:
		fsize = os.path.getsize(i+".fastq")
		os.system("bash run_MECAT.sh " + i + " " + str(fsize) + " fastq")
elif filetype == 0:
	for i in flist:
		fsize = os.path.getsize(i+".fa")
		os.system("bash run_MECAT.sh " + i + " " + str(fsize) + " fa")
elif filetype == 0:
	for i in flist:
		fsize = os.path.getsize(i+".fasta")
		os.system("bash run_MECAT.sh " + i + " " + str(fsize) + " fasta")

# os.system("rm *.fa.fastq *.gkpStore* *.collect.fa -r *.gkpStore")
# os.system("mv *.fasta ../contigs")

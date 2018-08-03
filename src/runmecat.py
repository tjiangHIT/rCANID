# coding=utf-8   
import os          


file = os.listdir(".")
flist = []
for i in file:
	if i.split('.')[1] == 'fq':
		flist.append(i.split(".")[0])
		
for i in flist:
	fsize = os.path.getsize(i+".fq")
	os.system("bash run_MECAT.sh " + i + " " + str(fsize))

# os.system("rm *.fa.fastq *.gkpStore* *.collect.fa -r *.gkpStore")
# os.system("mv *.fasta ../contigs")

#!/bin/bash

FASTQ=$1
genome_size=$2
#filedic=$3
memory=40
threads=4
coverage=30
prefix="nonmapped"
directory="nonmapped"
#setformat=$7


echo "[INFO] Detecting overlapping candidates with MECAT"
echo "[CODES] mecat2pw -j 0 -d ${FASTQ}.fq -o ${FASTQ}.pm.can -t ${threads} -w wrk_dir"
mecat2pw -j 0 -d ${FASTQ}.fq -o ${FASTQ}.pm.can -t ${threads} -w wrk_dir

echo "[INFO] Correcting noisy reads with MECAT"
echo "[CODE] mecat2cns -i 0 -t ${threads} ${FASTQ}.pm.can ${FASTQ}.fq ${FASTQ}.collect.fa"
mecat2cns -i 0 -t ${threads} ${FASTQ}.pm.can ${FASTQ}.fq ${FASTQ}.collect.fa

echo "[INFO] Extracting the longest ${coverage}X corrected reads with MECAT"
echo "[CODE] extract_sequences ${FASTQ}.collect.fa ${FASTQ}.${coverage}x.fasta ${genome_size} ${coverage}"
extract_sequences ${FASTQ}.collect.fa ${FASTQ}.${coverage}x ${genome_size} ${coverage}

# echo "[INFO] Assemble the longest ${coverage}X corrected reads using mecat2cacu"
# echo "[CODE] mecat2canu -trim-assemble -p ${prefix} -d ${directory} genomeSize=${genome_size} ErrorRate=0.02 maxMemory=${memory} maxThreads=${threads} useGrid=0 Overlapper=mecat2asmpw -pacbio-corrected ${FASTQ}.${coverage}x.fasta"
# mecat2canu -trim-assemble -p ${prefix} -d ${directory} genomeSize=${genome_size} ErrorRate=0.02 maxMemory=${memory} maxThreads=${threads} useGrid=0 Overlapper=mecat2asmpw -pacbio-corrected ${FASTQ}.${coverage}x.fasta

echo "[INFO] Deleting all useless files..."
echo "[CODE] rm *.qual *.qv *.frg *.fa *.can *.can.part0 *.partition_files -r wrk_dir nonmapped"
rm *.qual *.qv *.frg *.fa *.can *.can.part0 *.partition_files -r wrk_dir nonmapped



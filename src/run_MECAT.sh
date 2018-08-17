#!/bin/bash

FASTQ=$1
genome_size=$2
filetype = $3
memory=40
threads=4
coverage=30
prefix=${FASTQ}
directory=${FASTQ}
#setformat=$7

LIMIT=10m

echo "[INFO] Detecting overlapping candidates with MECAT"
echo "[CODE] mecat2pw -j 0 -d ${FASTQ}.${filetype} -o ${FASTQ}.pm.can -t ${threads} -w wrk_dir -n 200 -a 50 -k 2"
mecat2pw -j 0 -d ${FASTQ}.${filetype} -o ${FASTQ}.pm.can -t ${threads} -w wrk_dir -n 200 -a 50 -k 2

echo "[INFO] Correcting noisy reads with MECAT"
echo "[CODE] mecat2cns -i 0 -t ${threads} ${FASTQ}.pm.can ${FASTQ}.${filetype} ${FASTQ}.collect.fa -r 0.2 -a 50 -c 2 -l 50"
mecat2cns -i 0 -t ${threads} ${FASTQ}.pm.can ${FASTQ}.${filetype} ${FASTQ}.collect.fa -r 0.2 -a 50 -c 2 -l 50

echo "[INFO] Extracting the longest ${coverage}X corrected reads with MECAT"
echo "[CODE] extract_sequences ${FASTQ}.collect.fa ${FASTQ}.${coverage}x.fasta ${genome_size} ${coverage}"
extract_sequences ${FASTQ}.collect.fa ${FASTQ}.${coverage}x.fasta ${genome_size} ${coverage}

echo "[INFO] Assemble the longest ${coverage}X corrected reads using mecat2cacu"
echo "[CODE] mecat2canu -trim-assemble -p ${prefix} -d ${directory} corOutCoverage=${coverage} minOverlapLength=50 minReadLength=50 genomeSize=${genome_size} ErrorRate=0.2 maxMemory=${memory} maxThreads=${threads} useGrid=0 Overlapper=mecat2asmpw -pacbio-corrected ${FASTQ}.${coverage}x.fasta"
mecat2canu -trim-assemble -p ${prefix} -d ${directory} corOutCoverage=${coverage} minOverlapLength=50 minReadLength=50 genomeSize=${genome_size} ErrorRate=0.2 maxMemory=${memory} maxThreads=${threads} useGrid=0 Overlapper=mecat2asmpw -pacbio-corrected ${FASTQ}.${coverage}x.fasta


echo "[INFO] Deleting all useless files..."
echo "[CODE] rm *.qual *.qv *.frg *.fa *.can *.can.part0 *.partition_files -r wrk_dir ${directory}"
#filesize=ls -l ${directory}/${prefix}.contigs.fasta | awk '{print $5}'
if -s ${directory}/${prefix}.contigs.fasta ; then
    cat ${directory}/${prefix}.contigs.fasta > ${prefix}.contigs.fasta
else
    cat ${directory}/${prefix}.unassembled.fasta > ${prefix}.unassembled.fasta
fi
rm  *.fa.fastq *.gkpStore* *.qual *.qv *.frg *.collect.fa *.can *.can.part0 *.partition_files -r wrk_dir ${directory} *.${coverage}x.fasta


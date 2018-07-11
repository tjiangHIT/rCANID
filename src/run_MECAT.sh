# #!/bin/bash

# FASTQ=$1
# genome_size=$2
# threads=$3
# coverage=$4
# prefix=$5
# directory=$6
# setformat=$7


# echo "[INFO] Detecting overlapping candidates with MECAT"
# mecat2pw -j 0 -d ${FASTQ} -o ${FASTQ}.pm.can -t ${threads} -w wrk_dir

# echo "[INFO] Correcting the noisy reads with MECAT"
# mecat2cns -i 0 -t ${threads} ${FASTQ}.pm.can ${FASTQ} ${FASTQ}.collect.fa

# echo "[INFO] Extracting the longest ${4}X corrected reads with MECAT"
# extract_sequences ${FASTQ}.collect.fa ${FASTQ}.${coverage}x.fasta ${genome_size} ${coverage}

# echo "[INFO] Assemble the longest ${4}X corrected reads using mecat2cacu"
# mecat2canu -trim-assemble -p ${prefix} -d ${directory} genomeSize=${genome_size} ErrorRate=0.02 maxMemory=40 maxThreads=${threads} useGrid=0 Overlapper=mecat2asmpw -pacbio-corrected ${FASTQ}.${coverage}x.fasta


#!/bin/bash

FASTQ=$1
genome_size=4800000
memory=40
threads=16
coverage=25
prefix="nonmapped"
directory="nonmapped"
#setformat=$7


echo "[INFO] Detecting overlapping candidates with MECAT"
echo "[CODES] mecat2pw -j 0 -d ${FASTQ}.fastq -o ${FASTQ}.pm.can -t ${threads} -w wrk_dir"
mecat2pw -j 0 -d ${FASTQ}.fastq -o ${FASTQ}.pm.can -t ${threads} -w wrk_dir

 echo "[INFO] Correcting the noisy reads with MECAT"
echo "[CODE] mecat2cns -i 0 -t ${threads} ${FASTQ}.pm.can ${FASTQ}.fastq ${FASTQ}.collect.fa"
mecat2cns -i 0 -t ${threads} ${FASTQ}.pm.can ${FASTQ}.fastq ${FASTQ}.collect.fa

echo "[INFO] Extracting the longest ${coverage}X corrected reads with MECAT"
echo "[CODE] extract_sequences ${FASTQ}.collect.fa ${FASTQ}.${coverage}x.fasta ${genome_size} ${coverage}"
extract_sequences ${FASTQ}.collect.fa ${FASTQ}.${coverage}x ${genome_size} ${coverage}

echo "[INFO] Assemble the longest ${coverage}X corrected reads using mecat2cacu"
echo "[CODE] mecat2canu -trim-assemble -p ${prefix} -d ${directory} genomeSize=${genome_size} ErrorRate=0.02 maxMemory=${memory} maxThreads=${threads} useGrid=0 Overlapper=mecat2asmpw -pacbio-corrected ${FASTQ}.${coverage}x.fasta"
mecat2canu -trim-assemble -p ${prefix} -d ${directory} genomeSize=${genome_size} ErrorRate=0.02 maxMemory=${memory} maxThreads=${threads} useGrid=0 Overlapper=mecat2asmpw -pacbio-corrected ${FASTQ}.${coverage}x.fasta

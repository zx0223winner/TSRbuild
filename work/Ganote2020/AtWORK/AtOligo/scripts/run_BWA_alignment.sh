#!/bin/bash

#Setting variables:
#
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${DIR}/configfile

BWAdir=${BASEDIR}/${EXPERIMENT}/BWAdir

mkdir ${BWAdir}
mkdir ${BWAdir}/${GENOME_DIR}
ln -s ${BASEDIR}/${GENOME_DIR}/${GENOME_FILE} ${BWAdir}/${GENOME_DIR}/${GENOME_FILE}


cd ${BWAdir}

FQ_DIR=../fastq/
READS=${EXPERIMENT}_trimmed_tagdusted.fq


echo "Indexing the Arabidopsis genome (TAIR10) using bwa ..."
${BWA} index ${GENOME_DIR}/${GENOME_FILE}

echo "Aligning the A_thaliana reads to the TAIR10 index using BWA ..."
echo "bwa mem -t ${THREADS} ${GENOME_DIR}/${GENOME_FILE} ${FQ_DIR}/${READS} | ${SAMTOOLS} view -b -h -F 4 | samtools sort -o ${EXPERIMENT}_bwa.bam"
${BWA} mem -t ${THREADS} ${GENOME_DIR}/${GENOME_FILE} ${FQ_DIR}/${READS} | \
	${SAMTOOLS} view -b -h -F 4 | ${SAMTOOLS} sort -o ${EXPERIMENT}_bwa.bam
echo "... done"

${SAMTOOLS} index ${EXPERIMENT}_bwa.bam 
${SAMTOOLS} view -b -h ${EXPERIMENT}_bwa.bam 'Chr1' 'Chr2' 'Chr3' 'Chr4' 'Chr5' > ${EXPERIMENT}.genome.bam

${SAMTOOLS} view -h ${EXPERIMENT}.genome.bam > genome.sam
${SAMTOOLS} view -h ${EXPERIMENT}_bwa.bam > incorg.sam

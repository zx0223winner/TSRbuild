#!/bin/bash

#Setting variables:
#
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${DIR}/configfile

STARdir=${BASEDIR}/${EXPERIMENT}/STARdir

mkdir ${STARdir}
mkdir ${STARdir}/${GENOME_DIR}
ln -s ${BASEDIR}/${GENOME_DIR}/${GENOME_FILE} ${STARdir}/${GENOME_DIR}/${GENOME_FILE}


cd ${STARdir}

FQ_DIR=../fastq/
READ1=${EXPERIMENT}_trimmed_tagdusted_READ1.fq
READ2=${EXPERIMENT}_trimmed_tagdusted_READ2.fq


echo "Indexing the Arabidopsis genome (TAIR10) using STAR ..."
${STAR} --runMode genomeGenerate --runThreadN ${THREADS} --genomeDir ${GENOME_DIR} --genomeFastaFiles ${GENOME_DIR}/${GENOME_FILE}

echo "Aligning the A_thaliana reads to the TAIR10 index using STAR ..."
${STAR} --runMode alignReads --limitBAMsortRAM 19114876484 --runThreadN ${THREADS} --outSAMtype BAM SortedByCoordinate --outSAMorder Paired  --outFileNamePrefix ${EXPERIMENT}_STAR_ --genomeDir ${GENOME_DIR}  --readFilesIn ${FQ_DIR}/${READ1} ${FQ_DIR}/${READ2}
echo "... done"

ln -s ${STARdir}/${EXPERIMENT}_STAR_Aligned.sortedByCoord.out.bam ${EXPERIMENT}_STAR.bam
time ${SAMTOOLS} view -b -f 2 -@ 8 ${EXPERIMENT}_STAR.bam > ${EXPERIMENT}_STAR_proper.bam 

${SAMTOOLS} index ${EXPERIMENT}_STAR_proper.bam 
${SAMTOOLS} view -b -h ${EXPERIMENT}_STAR_proper.bam 'Chr1' 'Chr2' 'Chr3' 'Chr4' 'Chr5' > ${EXPERIMENT}.genome.bam

${SAMTOOLS} view -h ${EXPERIMENT}.genome.bam > genome.sam
${SAMTOOLS} view -h ${EXPERIMENT}_STAR_proper.bam > incorg.sam

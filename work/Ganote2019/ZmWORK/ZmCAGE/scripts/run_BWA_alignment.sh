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

echo "Indexing the Zea mays genome using bwa ..."
echo "bwa index ${GENOME_DIR}/${GENOME_FILE}"
${BWA} index ${GENOME_DIR}/${GENOME_FILE}

echo ${fastqDIR}

echo "Starting alignments ..."
for fq in ${fastqDIR}/*_demultiplexed_trimmed_filtered_tagdusted.fq;
do
	SAMPLE=$fq
	SORTED_BAM=$(basename $fq .fq)_sorted.bam
	FILTERED_BAM=$(basename $fq .fq)_sorted_mapq_20.bam

	echo "bwa aln -t ${THREADS} -n 3 ${GENOME_DIR}/${GENOME_FILE} -f $(basename ${SAMPLE} .fq).sai ${SAMPLE}"
	${BWA} aln -t ${THREADS} -n 3 ${GENOME_DIR}/${GENOME_FILE} -f $(basename ${SAMPLE} .fq).sai ${SAMPLE}

	echo "bwa samse ${GENOME_DIR}/${GENOME_FILE} $(basename ${SAMPLE} .fq).sai ${SAMPLE} | \
		samtools view -uS - | \
		samtools sort -O BAM - > ${SORTED_BAM}"

	${BWA} samse ${GENOME_DIR}/${GENOME_FILE} $(basename ${SAMPLE} .fq).sai ${SAMPLE} | \
		${SAMTOOLS} view -uS - | \
		${SAMTOOLS} sort -O BAM - > ${SORTED_BAM}

	echo "samtools index -b ${SORTED_BAM}"
	${SAMTOOLS} index -b ${SORTED_BAM}

	#.. post-alignment filtering for proper alignments and MAPQ >= 20:
	#
	echo "samtools view -F 4 -q 20 -u ${SORTED_BAM} | samtools sort -O BAM -@ 10 - > ${FILTERED_BAM}"
	${SAMTOOLS} view -F 4 -q 20 -u ${SORTED_BAM} | ${SAMTOOLS} sort -O BAM -@ 10 - > ${FILTERED_BAM}

	echo "samtools index -b ${FILTERED_BAM}"
	${SAMTOOLS} index -b ${FILTERED_BAM}
done

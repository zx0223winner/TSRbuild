#!/bin/bash

#Setting variables:
#
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${DIR}/configfile


Reads=${EXPERIMENT}
suffix=fq


ReadsIN="${Reads}_R1.${suffix} ${Reads}_R2.${suffix}"
ReadsOUT="${Reads}_trimmed_R1.${suffix} ${Reads}_unpaired_R1.${suffix} ${Reads}_trimmed_R2.${suffix} ${Reads}_unpaired_R2.${suffix}"

echo "Clipping adapaters from the reads using Trimmomatic ..."

java -classpath ${TRIMMOMATIC} \
	org.usadellab.trimmomatic.TrimmomaticPE \
	-threads ${THREADS} \
	${ReadsIN} \
	${ReadsOUT} \
	ILLUMINACLIP:${TrueSeq2PE}:2:30:10 TRAILING:20 MINLEN:25

#!/bin/bash

#Setting variables:
#
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${DIR}/configfile


Reads=${EXPERIMENT}
suffix=fq


ReadsIN=${Reads}.${suffix}
ReadsOUT=${Reads}_trimmed.${suffix}

echo "Clipping adapaters from the reads using Trimmomatic ..."

java -classpath ${TRIMMOMATIC} \
	org.usadellab.trimmomatic.TrimmomaticSE \
	-threads ${THREADS} \
	${ReadsIN} \
	${ReadsOUT} \
	ILLUMINACLIP:${TrueSeq2SE}:2:30:10 TRAILING:20 MINLEN:25

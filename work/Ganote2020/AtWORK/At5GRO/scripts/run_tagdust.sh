#!/bin/bash

#Setting variables:
#
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${DIR}/configfile


echo "Running tagdust on ${EXPERIMENT} samples"

R1=${EXPERIMENT}_trimmed.fastq
OP=${EXPERIMENT}_trimmed_tagdusted

${TAGDUST} -ref ${RNAfile} -dust 97 -t ${THREADS} -fe 3 -1 R:N ${R1} -o ${OP}

echo "Job Complete."

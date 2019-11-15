#!/bin/bash

#Setting variables:
#
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${DIR}/configfile


echo "Running tagdust on ${EXPERIMENT} samples"

R1=${EXPERIMENT}_trimmed_R1.fq
R2=${EXPERIMENT}_trimmed_R2.fq
OP=${EXPERIMENT}_trimmed_tagdusted

${TAGDUST} -ref ${RNAfile} -dust 97 -t ${THREADS} -fe 3 -1 R:N ${R1} ${R2} -o ${OP}

echo "Job Complete."

#!/bin/bash

echo "Running fastqc on trimmed ZmCAGE samples"

cd ../fastqc

for FQ in ../fastq/*_demultiplexed_trimmed_filtered.fastq; do

    fastqc $FQ

done

echo "Job Complete."

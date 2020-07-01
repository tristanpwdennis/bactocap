#!/bin/bash -ue
fastp -i C0101_S58_R1_001.fastq.gz -I C0101_S58_R2_001.fastq.gz -o C0101_S58_trim_R1.fastq.gz -O C0101_S58_trim_R2.fastq.gz

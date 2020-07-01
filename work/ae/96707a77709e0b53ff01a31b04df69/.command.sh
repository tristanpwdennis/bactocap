#!/bin/bash -ue
fastp -i BZA-1125_S66_R1_001.fastq.gz -I BZA-1125_S66_R2_001.fastq.gz -o BZA-1125_S66_trim_R1.fastq.gz -O BZA-1125_S66_trim_R2.fastq.gz

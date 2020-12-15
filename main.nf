#!/usr/bin/env nextflow
NXF_OPTS='-Xms512M -Xmx2G'
VERSION="0.1"

log.info "===================================================================="
log.info "                B A C T O C A P       (v${VERSION})                 "
log.info "===================================================================="

params.help = ""
if (params.help) {
  log.info " "
  log.info "The BACTOCAP workflow will run on whichever dataset is passed as an argument as shown below. You will be able to specify your own read data and reference genome at some point."
  log.info " "
  log.info "USAGE: "
  log.info " "
  log.info "nextflow run main.nf --dataset <dataset>"
  log.info " "
  log.info "Arguments:"
  log.info "    --dataset  STRING: anthrax, mlst, mycoplasma  (e.g. --dataset anthrax)  Pick whether to run BACTOCAP on anthrax, mlst, or mycoplasma datasets"
  log.info " "
  log.info "Optional arguments:"
  log.info "    --mappingonly     Will not run variant calling"
  log.info "===================================================================="
  exit 1
}

/*
 * Defines some parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */

if (params.dataset == "anthrax") {

  dirpath = "$baseDir/datasets/anthrax"

} else if (params.dataset == "mlst") {

  dirpath = "$baseDir/datasets/mlst"

} else if (params.dataset == "mycoplasma") {

  dirpath = "$baseDir/datasets/mycoplasma"

} else {

  println "Organism not specified, please specify organism as shown in --help"

}

params.reads = "${dirpath}/raw_reads/*{R1_001,R2_001}.fastq.gz"
params.annotation = "${dirpath}/ancillary/*.bed"
params.fasta = "${dirpath}/ref/*.fasta"
params.dict = "${dirpath}/ref/*.dict"
params.fai = "${dirpath}/ref/*.fasta.fai"
params.bwt = "${dirpath}/ref/*.fasta.bwt"
params.ann = "${dirpath}/ref/*.fasta.ann"
params.pac = "${dirpath}/ref/*.fasta.pac"
params.sa = "${dirpath}/ref/*.fasta.sa"
params.amb = "${dirpath}/ref/*.fasta.amb"
params.results = "${dirpath}/results"


/*
 * Take params and turn into channels
 */


if (params.reads) {
Channel
    .fromFilePairs( params.reads, flat:true )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs_ch } 
}
if (params.annotation) {
Channel
    .fromPath( params.annotation)
    .ifEmpty { error "Cannot find any annotation matching: ${params.annotation}" }
    .set { bed } 
}
if (params.fasta) {
Channel
    .fromPath( params.fasta )
    .ifEmpty { error "Cannot find any FASTA file matching: ${params.fasta}" }
    .into { fasta1; fasta2; fasta3; fasta4; fasta5 } 
}
if (params.dict) {
Channel
    .fromPath( params.dict )
    .ifEmpty { error "Cannot find any GATK sequence dictionary matching: ${params.dict}" }
    .into { dict1; dict2; dict3; dict4} 
}
if (params.fai) {
Channel
    .fromPath( params.fai )
    .ifEmpty { error "Cannot find any samtools .fai index matching: ${params.fai}" }
    .into { fai1; fai2; fai3; fai4} 
}
if (params.bwt) {
Channel
    .fromPath( params.bwt )
    .ifEmpty { error "Cannot find any bwa .bwt file matching: ${params.bwt}" }
    .set { bwt } 
}
if (params.ann) {
Channel
    .fromPath( params.ann )
    .ifEmpty { error "Cannot find any bwa .ann file matching: ${params.ann}" }
    .set { ann } 
}
if (params.pac) {
Channel
    .fromPath( params.pac )
    .ifEmpty { error "Cannot find any bwa .pac file matching: ${params.pac}" }
    .set { pac } 
}
if (params.sa) {
Channel
    .fromPath( params.sa )
    .ifEmpty { error "Cannot find any bwa .sa file matching: ${params.sa}" }
    .set { sa } 
}
if (params.amb) {
Channel
    .fromPath( params.amb )
    .ifEmpty { error "Cannot find any bwa .amb file matching: ${params.amb}" }
    .set { amb } 
}

//Allocate memory
int threads    = Runtime.getRuntime().availableProcessors()
threadmem      = (((Runtime.getRuntime().maxMemory() * 4) / threads) as nextflow.util.MemoryUnit)
threadmem_more = 4 * threadmem


process ReadTrimming {
  tag "Trimming ${pair_id}"
  memory threadmem_more
  cpus 4

  publishDir "$params.results/${pair_id}"

  input:
  tuple val(pair_id), file(mate1), file(mate2) from read_pairs_ch

  output:
  tuple val(pair_id), file("${pair_id}_trim_R1.fastq.gz"), file("${pair_id}_trim_R2.fastq.gz") into trimmed_reads, readsforqc

  script:

  """
  fastp -i $mate1 -I $mate2 -o ${pair_id}_trim_R1.fastq.gz -O ${pair_id}_trim_R2.fastq.gz
  """

}

bwa_index = ann.merge(bwt, pac, sa, amb, fasta1)
ref1 = trimmed_reads.combine(bwa_index)

//align with bwa
process bwaAlign{
  tag "Aligning ${pair_id} to reference: ${fasta1}"
  memory threadmem_more
  cpus 4

  input:
  tuple val(pair_id), path(read1), path(read2), file(bwt), file(ann), file(pac), file(sa), file(amb), file(fasta1) from ref1

  output:
  tuple val(pair_id), file("${pair_id}.sam") into bwa_bam

  script:

  """
  bwa mem -t ${task.cpus} ${fasta1} $read1 $read2 > ${pair_id}.sam
  """

}

//sambam conversion and sort
process bwaSort{
  tag "Samtools sorting ${pair_id}"
  memory '4 GB'
  input:
  tuple val(pair_id), file(sam) from bwa_bam

  output:
  tuple val(pair_id), file("${pair_id}.srt.bam") into bwa_sorted

  script:

  """
 samtools view -bS $sam | samtools sort - -o ${pair_id}.srt.bam
  """

}

//Mark duplicate reads
process MarkDuplicates {
  tag "Marking duplicate reads for: ${pair_id}"

  input:
  tuple val(pair_id), path(bam_file) from bwa_sorted

  output:
  tuple val(pair_id), file("${pair_id}.dupmarked.bam") into dupmarked_ch

  script:
  """
  gatk MarkDuplicates \
  -I $bam_file \
  -M ${pair_id}.metrics \
  -O ${pair_id}.dupmarked.bam 
  """
}

//add read groups
process AddOrReplaceReadGroups {
  tag "Adding read groups to: ${pair_id}"
  publishDir "$params.results/${pair_id}", mode: 'copy'
  input:
  tuple val(pair_id), file(bam_file) from dupmarked_ch

  output:
  tuple val(pair_id), file("${pair_id}.rg.bam"), file("${pair_id}.rg.bai") into rg_bam, bamqc1, bamqc2, bamqc3

  script:
  """
  gatk AddOrReplaceReadGroups \
  -I ${bam_file} \
  -O ${pair_id}.rg.bam \
  -LB ${pair_id} \
  -PL ILLUMINA \
  -PU NA \
  -SM ${pair_id} 

  samtools index ${pair_id}.rg.bam ${pair_id}.rg.bai
  """
}


bits = bed.merge(fasta5, fai4, dict4)
docchannel = bamqc3.combine(bits)

process CollectDoC {
  tag "Collecting coverage stats from ${pair_id}"
  publishDir "$params.results/${pair_id}"
  input:
  tuple val(pair_id), file("${pair_id}.rg.bam"), file("${pair_id}.rg.bai"), file(bed), file(fasta5), file(fai4), file(dict4) from docchannel
  output:
  file ("${pair_id}.sample_cumulative_coverage_counts") 
  file ("${pair_id}.sample_cumulative_coverage_proportions")
  file ("${pair_id}.sample_interval_statistics")
  file ("${pair_id}.sample_interval_summary")
  file ("${pair_id}.sample_statistics")
  file ("${pair_id}.sample_summary")

script:
 """
gatk DepthOfCoverage \
-R ${fasta5} \
-O ${pair_id} \
-I ${pair_id}.rg.bam \
-L ${bed}
  """  
}

process FlagstatRun {
  tag "Collecting read metadata from ${pair_id}"
  publishDir "$baseDir/results/individual_reports"

  input:
  tuple val(pair_id), file(bam), file(bai) from bamqc2

  output:
  file ("${pair_id}.stats.txt") into flagstat_results
  file ("${pair_id}.stats.txt") into flagstat_collect


  script:
  """
  samtools flagstat ${pair_id}.rg.bam > ${pair_id}.stats.txt
  """
}

process FlagstatCollect  {
  publishDir "$params.results/${pair_id}"
  input:
  file(flagfile) from flagstat_collect.toList()


  output:
  file("mapping_stats.csv") into flagtextfile_ch
  script:
  """
  printf sample_id,total,secondary,supplementary,duplicates,mapped,paired,read1,read2,properly_paired,with_itself_and_mate_mapped,singletons,mate_on_diff_chr,mate_on_diff_chrover5,"\n" > mapping_stats.csv

  for f in *txt
  do
   flag=`< \$f cut -d \\+ -f 1 | tr -s '[:blank:]' ','`
   echo "\${f%.stats.txt}",\$flag | tr -d ' ' >> mapping_stats.csv
  done
  """
}


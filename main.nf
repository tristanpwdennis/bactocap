#!/usr/bin/env nextflow

VERSION="0.1"

log.info "===================================================================="
log.info "This is the BACTOCAP pipeline (v${VERSION})                        "
log.info "===================================================================="

params.help = ""
if (params.help) {
  log.info " "
  log.info "The BACTOCAP workflow will run on whichever dataset is passed as an argument as shown below. You will be able to specify your own read data and reference genome at some point."
  log.info " "
  log.info "USAGE: "
  log.info " "
  log.info "nextflow run main.nf --organism <organism>"
  log.info " "
  log.info "Arguments:"
  log.info "    --dataset  STRING: anthrax, mlst, mycoplasma  (e.g. --organism anthrax)  Pick whether to run BACTOCAP on anthrax, mlst, or mycoplasma datasets"
  log.info " "
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
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs_ch } 
}
if (params.fasta) {
Channel
    .fromPath( params.fasta )
    .ifEmpty { error "Cannot find any FASTA file matching: ${params.fasta}" }
    .into { fasta1; fasta2; fasta3; fasta4 } 
}
if (params.dict) {
Channel
    .fromPath( params.dict )
    .ifEmpty { error "Cannot find any GATK sequence dictionary matching: ${params.dict}" }
    .into { dict1; dict2; dict3 } 
}
if (params.fai) {
Channel
    .fromPath( params.fai )
    .ifEmpty { error "Cannot find any samtools .fai index matching: ${params.fai}" }
    .into { fai1; fai2; fai3 } 
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
  tuple val(pair_id), path(reads) from read_pairs_ch

  output:
  tuple val(pair_id), file("*.fq.gz") into trimmed_reads, readsforqc

  script:

  """
  trim_galore --paired -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT $reads
  """

}

process FastQC {
  tag "Performing fastqc on ${pair_id}"
  input:
  tuple val(pair_id), path(reads) from readsforqc
  output:
  file("fastqc_${pair_id}_logs") into fastqc_ch

  script:
    """
    mkdir fastqc_${pair_id}_logs

    fastqc -o fastqc_${pair_id}_logs -f fastq -q ${reads}
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
  tuple val(pair_id), path(reads), file(bwt), file(ann), file(pac), file(sa), file(amb), file(fasta1) from ref1

  output:
  tuple val(pair_id), file("${pair_id}.sam") into bwa_bam

  script:

  """
  bwa mem -t ${task.cpus} ${fasta1} $reads > ${pair_id}.sam
  """

}

//sambam conversion and sort
process bwaSort{
  tag "Samtools sorting ${pair_id}"
  memory threadmem_more
  cpus 4

  input:
  tuple val(pair_id), file(sam) from bwa_bam

  output:
  tuple val(pair_id), file("${pair_id}.srt.bam") into bwa_sorted

  script:

  """
  samtools sort -o ${pair_id}.srt.bam -O BAM $sam
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
  publishDir "$params.results/${pair_id}"
  input:
  tuple val(pair_id), file(bam_file) from dupmarked_ch

  output:
  tuple val(pair_id), file("${pair_id}.rg.bam"), file("${pair_id}.rg.bai") into rg_bam, bamqc1, bamqc2

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

process Qualimap {
  tag "${pair_id}"
  publishDir "$baseDir/results/individual_reports"
  memory threadmem_more

  input:
  tuple val(pair_id), file(bam), file(bai) from bamqc1

  output:
  file ("${pair_id}") into qualimap_results

  script:
  """
  qualimap bamqc -bam ${pair_id}.rg.bam -outdir ${pair_id}
  """


}

process FlagstatRun {
  tag "${pair_id}"
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
  publishDir "$baseDir/results/"
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

process MultiQC {
  tag "Performing multiqc on fastqc output"
  publishDir "$params.results/"
  input:
  file('*')  from fastqc_ch.collect()
  file('*') from qualimap_results.collect()
  file('*')from flagstat_results.collect()

  output:
  file('multiqc_report.html')

  script:
  """
  export LC_ALL=C.UTF-8
  export LANG=C.UTF-8

  multiqc .
  """  
}


//merge reference and bam data into one channel
gatk_dict = fasta2.merge(dict1, fai1)
haplotypecaller = rg_bam.combine(gatk_dict)

process HaplotypeCaller {
  tag "Calling variants for: ${pair_id}"
  publishDir "$baseDir/results/${pair_id}"
  input:

  tuple val(pair_id), file(bam), file(bai), file(fasta2), file(dict1), file(fai1) from haplotypecaller

  output:
  val(pair_id) into pair_id
  file "${pair_id}.g.vcf" into gvcf_channel
  file "${pair_id}.g.vcf.idx" into gvcf_index

  script:
   """
   gatk HaplotypeCaller \
    -R $fasta2 \
    -O ${pair_id}.g.vcf \
    -I $bam \
    -ERC GVCF \
    -isr INTERSECTION \
    --native-pair-hmm-threads 1 \
    --max-alternate-alleles 3 \
    -contamination 0 \
    --QUIET
    """
}

//collect fasta, dict and fai data into a channel
dbdict = fasta3.merge(dict2, fai2)


process GenomicsDBImport {
  tag "Collecting gvcf into genomicsDB"
  input:
  file(gvcf) from gvcf_channel.collect()
  file(gvcf_idx) from gvcf_index.collect()
  val(pair_id) from pair_id.collect()
  tuple file(fasta3), file(dict2), file(fai2) from dbdict

  output: path("input_variant_files.list") into listch
          file("genomics_db") into genodb_ch

  script:
"""
  #make list of vcf files for import to genomicsdb

  for vcf in \$(ls *.vcf); do
    echo \$vcf >> input_variant_files.list
  done

  sed 's/.vcf//g' input_variant_files.list | awk '{print\$0"\t"\$0".vcf"}' > vcflist.list


  gatk GenomicsDBImport \
  --genomicsdb-workspace-path genomics_db \
  -L NC_007530_Bacillus_anthracis_Ames_Ancestor \
    ${gvcf.collect { "-V $it " }.join()} \
  --reader-threads 3 

"""

}

genodict = fasta4.merge(dict3, fai3)

process GenotypeGVCF {
  tag "Generating final vcf from gvcfs and genomicdb"
  publishDir "$baseDir/results/"
  input:
  tuple file(fasta4), file(dict3), file(fai3) from genodict
  file(genomics_db) from genodb_ch

  output:
  file("final_anthrax.vcf") into finalvcf_ch

  script:
  """
  gatk GenotypeGVCFs \
  -R ${fasta4} \
  -V gendb://genomics_db \
  -O final_anthrax.vcf
  """


}



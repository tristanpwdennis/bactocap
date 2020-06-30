#!/usr/bin/env nextflow

VERSION="0.1"

log.info "===================================================================="
log.info "This is the BACTOCAP pipeline (v${VERSION})                        "
log.info "===================================================================="

params.help = ""
if (params.help) {
  log.info " "
  log.info "The BACTOCAP workflow will run on whichever organism is passed as an argument as shown below. You will be able to specify your own read data and reference genome at some point."
  log.info " "
  log.info "USAGE: "
  log.info " "
  log.info "nextflow run main.nf --organism <organism>"
  log.info " "
  log.info "Arguments:"
  log.info "    --organism  STRING: anthrax, mlst, mycoplasma  (e.g. --organism anthrax)  Pick whether to run BACTOCAP on anthrax, mlst, or mycoplasma datasets"
  log.info " "
  log.info "===================================================================="
  exit 1
}

if (params.organism == "anthrax") {

  dirpath = "$baseDir/datasets/anthrax"

} else if (params.organism == "mlst") {

  dirpath = "$baseDir/datasets/mlst"

} else if (params.organism == "mycoplasma") {

  dirpath = "$baseDir/datasets/mycoplasma"

} else {

  println "Organism not specified, please specify organism as shown in --help"

}


params.reads = "${dirpath}/raw_reads/*{_uniq1,_uniq2}.fastq.gz"
params.fasta = "${dirpath}/ref/NC_007530.fasta"
params.dict = "${dirpath}/ref/NC_007530.dict"
params.fai = "${dirpath}/ref/NC_007530.fasta.fai"
params.bwt = "${dirpath}/ref/NC_007530.fasta.bwt"
params.ann = "${dirpath}/ref/NC_007530.fasta.ann"
params.pac = "${dirpath}/ref/NC_007530.fasta.pac"
params.sa = "${dirpath}/ref/NC_007530.fasta.sa"
params.amb = "${dirpath}/ref/NC_007530.fasta.amb"
params.results = "${dirpath}/results"



params.reads.println()


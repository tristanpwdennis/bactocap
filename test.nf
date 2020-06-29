#!/usr/bin/env nextflow

VERSION="0.1"

log.info "===================================================================="
log.info "This is the BACTOCAP pipeline (v${VERSION})                        "
log.info "===================================================================="

params.help = ""
if (params.help) {
  log.info " "
  log.info "The BACTOCAP workflow will run automatically on any paired short sequencing reads placed in the raw_reads directory"
  log.info " "
  log.info " "
  log.info "USAGE: "
  log.info " "
  log.info "nextflow run main.nf"
  log.info " "
  log.info "Optional arguments:"
  log.info "    --organism         STRING: anthrax, mlst, mycoplasma  (e.g. --organism anthrax)  Pick whether to run BACTOCAP on anthrax, mlst, or mycoplasma datasets"
  log.info " "
  log.info "===================================================================="
  exit 1
}
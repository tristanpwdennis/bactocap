# bactocap: Anthrax and Mycoplasma
### Target-enrichment sequencing yields valuable genomic data for difficult-to-culture bacteria of public health importance

This repo contains all the materials required to reproduce the analysis and workflow from the targeted sequence capture project in *Dennis et al, 2022*, see: LINK

### Getting started

This project uses Docker to manage all the dependencies, and nextflow to run the analysis. To get started, make sure you have docker installed. Installation instructions by platform are here:
https://docs.docker.com/engine/install/ . Once you're finished, fire up the terminal and doublecheck with ```docker -v```
Then, install nextflow: https://www.nextflow.io/docs/latest/getstarted.html . And check in the terminal again to make sure it runs ok.

### Setup and Docker installations

Clone the repo
```
git clone https://github.com/tristanpwdennis/bactocap.git
```
Enter the repo
```cd bactocap```

Now we need to build the custom Docker image for this project, and also download the GATK Docker image.
This command will build the dennistpw/align Docker image. This will take a few minutes.
```
DOCKER_BUILDKIT=1 docker build -t dennistpw/align --no-cache . 
```
Next we need to pull the GATK docker image
```
docker pull broadinstitute/gatk
```
Now let's check to make sure both of the images are ok
```
docker images
```
You should see the gatk and align repos are in the list.

### Data
The raw data are located at PRJEB46822 (_B. anthracis_) and PRJEB50216 (_M. amphoriforme_). Download the raw reads into the corresponding dataset/organism/raw_read directories.

### Running the workflow
It's as simple as running 
```
nextflow run main.nf --help
```
This will prompt the USAGE statement and some brief pointers.
```
===================================================================
This is the BACTOCAP pipeline (VERSION)                        
===================================================================
The BACTOCAP workflow will run on whichever dataset is passed as an argument as shown below. 

USAGE: 

nextflow run main.nf --dataset <dataset>

Arguments:
   --dataset  STRING: anthrax, mycoplasma  (e.g. --organism anthrax)  Pick whether to run BACTOCAP on anthrax, or mycoplasma datasets

====================================================================
```
Nextflow caches all the steps, so you don't have to go back to square one with each reanalysis. Just add more data to the raw_reads directory, or restart if you accidentally shut off your machine with
```
nextflow run main.nf -resume
```

Note, I quite like running these scripts in screen sessions: https://linuxize.com/post/how-to-use-linux-screen/
This allows me to run the workflow, check on it periodically as it runs on the other screen, whilst I tool about doing other stuff. It also reduces the likelihood of that scenario where you accidentally close your laptop when you have a terminal session running and halt your analysis - TD


### Output

The final bam files and mapping/sequencing stats will be published in the ```results``` directory in each dataset directory according to sample name
The individual fastqc and bamqc data will be published in the ``` individual_reports``` subdirectory and agglomerated in the ```multiqc_report.html``` document.
A tab delimited text file ```mapping_stats.csv``` contains the flagstat data for analysis etc.

### Analysis

Running the Rscript ```generate_bactocap_metadata.R``` will take the mapping output and parse it into a CSV containing sample metadata, mapping, duplicate and coverage information for anthrax and mycoplasma. Running ```bcanalysisfull.R``` in RStudio will generate plots and tables in the *output* directory. Code and plots can be examined more closely in the ```anthrax_mycoplasma_bactocap.Rmd``` notebook.

### Pointers

Annotations and metadata are located in the `ancillary` directory. Reference genomes are contained in organism-specific directories under `datasets`






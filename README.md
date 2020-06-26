# bactocap
This repo contains all the materials required to reproduce the analysis and workflow from the bactocap project

### Getting started

This project uses Docker to manage all the dependencies. To get started, make sure you have docker installed. Installation instructions by platform are here:
https://docs.docker.com/engine/install/ . Once you're finished, fire up the terminal and doublecheck with ```docker -v```

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
When the project is further advanced, the read data for the bactocap project will be available on sra, and I will include fastq-dump commands in the workflow that will enable download of the data. Right now, however, we will have to make do the cheapo way, so please put some trimmed read files of your own into the raw_reads directory.

### Running the workflow
It's as simple as running 
```
nextflow run main.nf
```
Nextflow caches all the steps, so you don't have to go back to square one with each reanalysis. Just add more data to the raw_reads directory, or restart if you accidentally shut off your machine with
```
nextflow run main.nf -resume
```

The final vcfs and bam files will be published in the ```results``` directory. 







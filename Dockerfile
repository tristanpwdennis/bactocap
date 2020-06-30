# Set the base image to ubuntu 18.04
FROM ubuntu:18.04

# File Author / Maintainer
MAINTAINER Tristan Dennis <tristanpwdennis@gmail.com>

#configure time zone and install tzdata for base -r installation

RUN export DEBIAN_FRONTEND=noninteractive && ln -fs /usr/share/zoneinfo/Europe/London /etc/localtime


#get bits and pieces
RUN apt-get update && apt-get install --yes --no-install-recommends \
    wget \
    locales \
    git \
    cmake \
    build-essential \
    gcc-multilib \
    python3 \
    openjdk-8-jre \
    python3-pip \
    libpython2.7-dev \
    autoconf \
    automake \
    make \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libncurses5-dev \
    ant \
    software-properties-common \
    gnupg2 \
    datamash \
    bwa \
    git-lfs \
    curl \
    unzip \
    python3-setuptools \
    tzdata \
    r-base 

#get fastqc
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
RUN unzip fastqc_v0.11.9.zip && \
chmod 755 FastQC/fastqc && \
ln -s $PWD/FastQC/fastqc /usr/local/bin/

#get fastp
RUN git clone https://github.com/OpenGene/fastp.git && cd fastp && make 
RUN ln -s $PWD/fastp/fastp /usr/local/bin

#install multiqc
RUN pip3 install multiqc
RUN export LC_ALL=C.UTF-8 && export LANG=C.UTF-8

#get and install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
RUN tar -xvf samtools-1.10.tar.bz2 && rm samtools-1.10.tar.bz2
WORKDIR samtools-1.10
RUN ./configure && make && make install
WORKDIR /

RUN Rscript -e "install.packages('optparse')"
RUN Rscript -e "install.packages('BiocManager', repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "BiocManager::install(c('NOISeq', 'Repitools'))"

RUN cd /opt && wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.2.1.zip
RUN cd /opt && unzip qualimap_v2.2.1.zip && rm qualimap_v2.2.1.zip
RUN Rscript /opt/qualimap_v2.2.1/scripts/installDependencies.r

ENV PATH="/opt/qualimap_v2.2.1:$PATH"




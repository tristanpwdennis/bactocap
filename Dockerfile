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


#install cutadapt
RUN pip3 install cutadapt

#get fastqc
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
RUN unzip fastqc_v0.11.9.zip && \
chmod 755 FastQC/fastqc && \
ln -s $PWD/FastQC/fastqc /usr/local/bin/

#install trim galore!
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.5.tar.gz -o trim_galore.tar.gz && tar xvzf trim_galore.tar.gz
RUN mv TrimGalore-0.6.5/trim_galore /usr/local/bin/ && rm -rf TrimGalore-0.6.5

#install multiqc
RUN pip3 install multiqc
RUN export LC_ALL=C.UTF-8 && export LANG=C.UTF-8


#get and install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
RUN tar -xvf samtools-1.10.tar.bz2 && rm samtools-1.10.tar.bz2
WORKDIR samtools-1.10
RUN ./configure && make && make install
WORKDIR /








# Set the base image to ubuntu 18.04
FROM ubuntu:18.04

# File Author / Maintainer
MAINTAINER Tristan Dennis <tristanpwdennis@gmail.com>

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
    git \
    ant \
    software-properties-common \
    gnupg2 \
    datamash \
    bwa \
    git-lfs


#install cutadapt
RUN python3 -m pip install --user --upgrade cutadapt

#get and install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
RUN tar -xvf samtools-1.10.tar.bz2 && rm samtools-1.10.tar.bz2
WORKDIR samtools-1.10
RUN ./configure && make && make install
WORKDIR /

RUN git clone --recursive git://github.com/ekg/freebayes.git
RUN cd freebayes && make -j4






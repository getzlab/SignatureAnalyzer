FROM nvidia/cuda:10.0-cudnn7-runtime-ubuntu18.04
MAINTAINER Shankara Anand

# -----------------------------
# Install Basics
# -----------------------------
RUN apt-get update && apt-get install -y software-properties-common && \
    apt-get update && apt-get install -y \
        apt-transport-https \
        build-essential \
        cmake \
        curl \
        libboost-all-dev \
        libbz2-dev \
        libcurl3-dev \
        liblzma-dev \
        libncurses5-dev \
        libssl-dev \
        python3 \
        python3-pip \
        sudo \
        unzip \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/

RUN python3 -m pip install --upgrade setuptools

# -----------------------------
# Install Signature Analyzer
# -----------------------------
RUN mkdir signatureanalyzer
COPY . /signatureanalyzer/
RUN python3 -m pip install -e ./signatureanalyzer/.

# Test
RUN signatureanalyzer -h

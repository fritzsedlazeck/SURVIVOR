# use the ubuntu base image
FROM ubuntu:18.04

# install required packages
RUN apt-get update && apt-get install -y \
    autoconf \
    build-essential \
    cmake \
    g++ \
    gfortran \
    git \
    libcurl4-gnutls-dev \
    hdf5-tools \
    libboost-date-time-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev \
    libbz2-dev \
    libhdf5-dev \
    libncurses-dev \
    liblzma-dev \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# set environment
ENV BOOST_ROOT /usr

# install alfred
RUN cd /opt \
    && git clone --recursive https://github.com/fritzsedlazeck/SURVIVOR.git \
    && cd /opt/SURVIVOR/Debug/ \
    && make

# Workdir
WORKDIR /root/

# Add Alfred to PATH
ENV PATH="/opt/SURVIVOR/Debug:${PATH}"

# by default /bin/sh is executed
CMD ["/bin/bash"]

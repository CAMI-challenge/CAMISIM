# Use Ubuntu 20.04 as base image
FROM ubuntu:20.04

# Set non-interactive mode
ENV DEBIAN_FRONTEND=non-interactive

# Update and install required packages
RUN apt-get update \
    && apt-get install -y \
    python \
    perl \
    pip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN pip install biopython==1.81
RUN pip install ete3==3.1.2

RUN apt-get update \
    && apt-get install -y samtools

CMD ["/bin/bash"]

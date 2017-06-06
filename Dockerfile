FROM bioboxes/base
# cami/binning
MAINTAINER Peter Hofmann, peter.hofmann@hhu.de

# package installation
RUN ${BBX_BINDIR}/dockerfile-install-packages \
python \
perl \
less \
nano \
bzip2 \
wget \
cmake


## Copy all files
ADD ./ /opt

## Set up Miniconda environment for python2
# COPY docker/snake/Miniconda-3.9.1-Linux-x86_64.sh /opt/miniconda.sh
ENV MINICONDA /opt/docker/snake/Miniconda-3.9.1-Linux-x86_64.sh
RUN chmod +x ${MINICONDA}
RUN ${MINICONDA} -p /opt/miniconda -b
ENV PATH /opt/miniconda/bin:$PATH
RUN conda update --yes conda && conda install --yes python=2.7

# Install Snakemake within a conda environment
RUN conda create --yes -n snakemake python=3.4 pip pyyaml && /opt/miniconda/envs/snakemake/bin/pip install snakemake
#RUN mkdir /opt/snake
#ADD docker/snake/Snakefile /opt/snake/Snakefile
#ADD docker/snake/config.json /opt/snake/config.json


## Install Python requirements
RUN ${BBX_BINDIR}/dockerfile-install-packages \
python2.7-numpy \
python-matplotlib \
python-biopython

## Install required modules of perl
RUN ${BBX_BINDIR}/dockerfile-install-packages libxml-simple-perl


# Get and extract all files
#ADD data /opt/data
#ADD bin /opt/tools
RUN cp /opt/tools/hmmer-2.3.2/hmmsearch /opt/tools/rnammer-1.2/hmmsearch


# Switch to bash as default shell
ENV SHELL /bin/bash

# add tasks
COPY docker/tasks $DCKR_TASKDIR

# change working dir
WORKDIR /opt

ADD docker/run.sh /usr/local/bin/run
RUN chmod a+x /usr/local/bin/run

#ENTRYPOINT ["/usr/local/bin/run"]

docker-MSP
==========

This is a docker container for the annotation of genomes and simulation of metagenome datasets.

Usage
-----

1. Get Docker build files

        git clone https://github.com/CAMI-challenge/MetagenomeSimulationPipeline.git
        git checkout docker

2. Build the Docker image

        docker build -t="pehofmann/emsep" ./

3. Create the folders to mount


4. Run Docker container

        sudo docker run \
        -v /data/cami/test/docker/unittest_MSP/nobackup/tmp:/bbx/mnt/tmp/ \
        -v /data/cami/test/docker/unittest_MSP/nobackup/output:/bbx/mnt/output/ \
        -v /data/cami/test/docker/unittest_MSP/unittest_ga/input:/bbx/mnt/input:ro \
        -v /data/cami/test/docker/unittest_MSP/unittest_ga/ref:/bbx/mnt/ref/:ro \
        -i -t --rm pehofmann/emsep --shell

        sudo docker run \
        -v /data/cami/test/docker/unittest_MSP/nobackup/tmp:/bbx/mnt/tmp/ \
        -v /data/cami/test/docker/unittest_MSP/nobackup/output:/bbx/mnt/output/ \
        -v /data/cami/test/docker/unittest_MSP/unittest_ms/input:/bbx/mnt/input:ro \
        -v /data/cami/test/docker/unittest_MSP/unittest_ga/ref:/bbx/mnt/ref/:ro \
        -i -t --rm pehofmann/emsep --shell

python /opt/tools/MetagenomeSimulationPipeline/genomeannotation.py /bbx/mnt/input/config.cfg
python /opt/tools/MetagenomeSimulationPipeline/metagenomesimulation.py /bbx/mnt/input/config.cfg

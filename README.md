[![Smoke Test](https://github.com/CAMI-challenge/CAMISIM/actions/workflows/ci.yml/badge.svg)](https://github.com/CAMI-challenge/CAMISIM/actions/workflows/ci.yml)
[![Unit Tests](https://github.com/CAMI-challenge/CAMISIM/actions/workflows/tests.yml/badge.svg)](https://github.com/CAMI-challenge/CAMISIM/actions/workflows/tests.yml)

# CAMISIM

CAMISIM is a software to model abundance distributions of microbial communities and to simulate corresponding shotgun metagenome datasets.\
It was mainly developed for the [Critical Assessment of Metagenome Annotation (CAMI)](http://microbiome-cosi.org/cami) challenge, but should be suitable for general use. Please don't hesitate to [open a new issue](https://github.com/CAMI-challenge/CAMISIM/issues) if you run into problems or need help.

### CAMISIM 2.0
CAMISIM received a major update to version 2.0, using nextflow. You can still use the old python standalone using the tag 1.31-final, it will not receive updates however.\
You can use the script "convert_config.py" to convert your CAMISIM1 config file to CAMISIM2 (no guarantees for correctness, please check yourself).\
This version has been tested, but if you encounter any unforeseen difficulties or differences from what you expect your simulation to look like, please raise an Issue.

### Documentation 
* [User manual](https://github.com/CAMI-challenge/CAMISIM/wiki/User-manual)
* [Configuration File Options](https://github.com/CAMI-challenge/CAMISIM/wiki/Configuration-File-Options)
* [File Formats](https://github.com/CAMI-challenge/CAMISIM/wiki/File-Formats)

### Citation

If you use CAMISIM, please cite the publication at *Microbiome*:
* Fritz*, Hofmann*, et al. (2019). **CAMISIM: Simulating metagenomes and microbial communities.** *Microbiome*, 2019, 7:17. doi:[10.1186/s40168-019-0633-6](https://doi.org/10.1186/s40168-019-0633-6)

A part of CAMISIM's functionality was also described in the CAMI manuscript, thus you may also cite:
* Sczyrba*, Hofmann*, Belmann*, et al. (2017). **Critical Assessment of Metagenome Interpretation—a benchmark of metagenomics software.** *Nature Methods*, 14, 11:1063–1071. doi:[10.1038/nmeth.4458](https://doi.org/10.1038/nmeth.4458)

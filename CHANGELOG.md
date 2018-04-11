# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/) and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.1.0]

### Added
- Added script for creating joint/hybrid gold standards
- Added script for creating json metadata of data sets

### Changed
- Output folder will now be created if it is not present

### Removed
- (Temporarily?) removed the possibility of using readsimulationwrapper as standalone

## [1.0.0] - 2018-03-16

### Added
- Added default config for long read simulation
- Added docker support
- Added no-replace option for genomes in from_profile (default True)
- Genomes are per default split at N positions
- Added new 2018 taxonomy
- Added changelog!

### Changed
- Loads of Bugfixes
- Changed libraries for downloading genome (now uses https)
- from_profile now uses ete2 API for NCBI

## [0.2.2] - 2017-09-12

### Changed
- Argument parsing now has more helptext
- CAMISIM from_profile prints help if called with no options

## [0.2.1] - 2017-09-08

### Added
- Added new ART error profile (MBARC-26)

### Changed
- Improved error handling for ftp-retrieval (re-try after timeouts)

## [0.2] - 2017-08-22

### Added
- from_profile script to create data sets from BIOM profile
- mapping files for downloading genomes from NCBI
- first merge of the original pipeline, the new read simulators and the from_profile branches

### Changed
- Reworked config and config sections
- i.e. config file reflects added read simulators
- Zero-abundance genomes are now allowed
- wgsim error rates are not fixed anymore
- Updated readme to reflect current changes

### Removed
- Removed some assertions which caused CAMISIM to terminate prematurely

## [0.1] - 2017-08-08

### Added
- Pipeline for the creation of metagenome data sets
- Version used for the creation of the [CAMI 1 data sets](https://data.cami-challenge.org/participate)
- Added defaults, default config files and files to execute a test run
- Implemented PBsim, NanoSim, wgsim read simulators
- Option to start from distribution files instead of completely *de novo*

### Changed
- Folder structure now separates external tools and internal scripts

[Unreleased](https://github.com/CAMI-challenge/CAMISIM/compare/1.1.0...HEAD)
[1.1.0](https://github.com/CAMI-challenge/CAMISIM/compare/1.0.0...1.1.0)
[1.0.0](https://github.com/CAMI-challenge/CAMISIM/compare/0.2.2...1.0.0)
[0.2.2](https://github.com/CAMI-challenge/CAMISIM/compare/0.2.1...0.2.2)
[0.2.1](https://github.com/CAMI-challenge/CAMISIM/compare/0.2...0.2.1)
[0.2](https://github.com/CAMI-challenge/CAMISIM/compare/0.1...0.2)
[0.1](https://github.com/CAMI-challenge/CAMISIM/releases/tag/0.1)

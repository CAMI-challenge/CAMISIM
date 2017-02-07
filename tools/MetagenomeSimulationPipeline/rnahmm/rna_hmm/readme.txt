This software is freely available at
http://tools.camera.calit2.net/camera/meta_rna/index.htm

If you use the software please cite:

1. Ying Huang, Paul Gilna and Weizhong Li
Identification of ribosomal RNA genes in metagenomic fragments
Bioinformatic Vol. 25 no. 10 2009, pages 1338-1340




Installation
=============

1) Checking system requirements

- Python (version 2.3 or higher)
- HMMER (tested on version 3.0b3, http://hmmer.janelia.org/)

2) Download rna_hmm_fs.tar.gz from our website, save it
  to your desired working directory
  
  tar -zxvf rna_hmm_fs.tar.gz 
  cd rna_hmm_fs
  ./rna_hmm.py -h, for help message
  ./rna_hmm.py -i, test.fna -o test.gff for a test example
  
3) Meaning of input parameters

  -i: input nucleotide sequence file in fasta format
  -o: output annotation file in gff format, check header for information
  -L: directory to store *.hmm files, default "HMMs/"
  -E: e-value cut-off for hmmsearch, default 0.01
  -k: kingdom for search, can be "arc","bac" or "arc,bac", if "arc,bac"
      is used, the program will select best kingdom for the sequence
      according to evalue
  -m: interested molecule type, can be
      "lsu": 23S rRNA; "ssu": 16S rRNA; "tsu": 5S rRNA
      or their combination joined by ","
      default value is "lsu,ssu,tsu"
      

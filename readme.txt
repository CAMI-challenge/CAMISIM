run.py
  main script that calls all others

gather_16s.sh
  searches through lists of genomes and writes all 16S genes found into a single file

cluster.sh
  aligns and clusters the 16S genes
  
otu_prediction_to_meta_table.py
  using the previous clustering, taxonomic classifications and otus are predicted

ani_prediction_to_meta_table.py
  if reference genome(s) available, the smallest avarage nucleotide identity is calculated
  TODO: calculate ANI between all unpublished genomes

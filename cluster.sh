#!/bin/bash
#HOME_DIR="/home/user03/pyprojects/taxonomic_classification_plus_otu/"
#t1="/home/user03/output/nobackup/2014_09_24_random1000_ref_silva_hmmer3/16S_rRNA.fna"
#t2="/home/user03/output/nobackup/2014_09_24_random1000_ref_silva_hmmer3/mothur_cluster_list.txt"
#t3="/home/user03/pyprojects/taxonomic_classification_plus_otu/references/"
#t4="0.05"
#t5="no longer used"
#t6="45"

HOME_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#getting passed arguments
INPUT_FASTA_FILE=$(basename "$1")
OUTPUT_FILE="$2"
OUTPUT_FOLDER=$(dirname "$OUTPUT_FILE")
extension="${INPUT_FASTA_FILE##*.}"
filename="${INPUT_FASTA_FILE%.*}"

INPUT_FASTA="$1"
SILVA_REFERENCE_DIRECTORY="$3"
CLUSTER_CUTOFF="$4"
FAST_TREE="$5"
PROCESSORS="$6"

CURRENT_TIME=$(date +"%Y_%m_%d_%T")
TOOL_DIR="$HOME_DIR/tools"

######################
#TMP_DIR="$HOME_DIR/2_tmp_($EUID)$CURRENT_TIME"
TMP_DIR=`mktemp -d`
######################


#set path to programs bases on assumed tool folder
MOTHUR="$TOOL_DIR/mothur"
FASTA_TO_PHYLIP="$TOOL_DIR/Fasta2Phylip.pl"
R_DIST_SCRIPT="$TOOL_DIR/dist_from_tree.R"

#list of possible output file names from mothur
#MERGED_ALIGNMENT="$filename.merged.align"
BAD_SEQ_FILE="$filename.unique.flip.accnos"


######################
#mkdir "$TMP_DIR"
#TMP_DIR="$TMP_DIR/_$filename"
#rm -rf "$TMP_DIR"
#mkdir "$TMP_DIR"
######################

SYMBOLIC_LINK="$TMP_DIR/$INPUT_FASTA_FILE"
ln -s "$INPUT_FASTA" "$SYMBOLIC_LINK"
cd "$TMP_DIR"

LOCAL_ALIGN_DIST="ref.align.dist"
REF_ALIGN_DIST="$SILVA_REFERENCE_DIRECTORY/mothur_ref_distances"
REF_UNIQUE_NAMES="$SILVA_REFERENCE_DIRECTORY/mothur_ref_names"
REF_FASTA="$SILVA_REFERENCE_DIRECTORY/mothur_alignment_ref.fasta"
# ref data can be huge 53 GB, but needs copied since it will be merged with the new distances
cp "$REF_ALIGN_DIST" "$LOCAL_ALIGN_DIST"

INPUT_NAMES_FILE="$filename.names"
INPUT_PICK_NAMES_FILE="$filename.pick.names"
MERGED_NAMES="$filename.merged.names"

echo "align.seqs"
mothur_commands="unique.seqs(fasta=$INPUT_FASTA_FILE)
align.seqs(candidate=current, template=$REF_FASTA, align=gotoh, flip=t, processors=$PROCESSORS)
remove.seqs(accnos=$BAD_SEQ_FILE, fasta=current, name=current)
dist.seqs(oldfasta=$REF_FASTA, column=$LOCAL_ALIGN_DIST, cutoff=$CLUSTER_CUTOFF, processors=$PROCESSORS, calc=onegap, countends=F)
merge.files(input=$INPUT_NAMES_FILE-$REF_UNIQUE_NAMES, output=$MERGED_NAMES)
merge.files(input=$INPUT_PICK_NAMES_FILE-$REF_UNIQUE_NAMES, output=$MERGED_NAMES)
set.current(name=$MERGED_NAMES, column=$LOCAL_ALIGN_DIST)
cluster(cutoff=$CLUSTER_CUTOFF, method=furthest, precision=100)"
echo -e "$mothur_commands" | "$MOTHUR"

if mv *.list "$OUTPUT_FOLDER/mothur_otu.txt"; then
# 2>/dev/null
	#echo "debug, clean up yourself!"
	echo "clean up folders"
	rm -rf "$TMP_DIR"
else
	echo "clustering failed"
	echo "intermediate data: $TMP_DIR"
fi
exit;
#!/bin/bash
HOME_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#if [ "$#" -eq 4 ]; then
#	FAST_TREE=0
#elif [ "$#" -eq 5 ]; then
#	FAST_TREE="$5"
#elif [ "$#" -ne 4 ]; then
#	echo "Illegal number of parameters: $#"
#	echo "$0 marker_gene.fasta output_otu.txt reference_db.align cluster_cutoff fast_tree"
#	exit;
#fi

#getting passed arguments
INPUT_FASTA_FILE=$(basename "$1")
OUTPUT_FILE="$2"
OUTPUT_FOLDER=$(dirname "$2")
extension="${INPUT_FASTA_FILE##*.}"
filename="${INPUT_FASTA_FILE%.*}"

INPUT_FASTA="$1"
REFERENCE_DB_ALIGNMENTS="$3"
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
MERGED_ALIGNMENT="$filename.merged.align"
BAD_SEQ_FILE="$filename.flip.accnos"


######################
#mkdir "$TMP_DIR"
#TMP_DIR="$TMP_DIR/_$filename"
#rm -rf "$TMP_DIR"
#mkdir "$TMP_DIR"
######################

SYMBOLIC_LINK="$TMP_DIR/$INPUT_FASTA_FILE"
ln -s "$INPUT_FASTA" "$SYMBOLIC_LINK"
cd "$TMP_DIR"

echo "align.seqs"
mothur_commands_00="align.seqs(candidate=$INPUT_FASTA_FILE, template=$REFERENCE_DB_ALIGNMENTS, align=gotoh, flip=t, processors=$PROCESSORS)
remove.seqs(accnos=$BAD_SEQ_FILE, fasta=current)"
echo -e "$mothur_commands_00" | "$MOTHUR"
#exit;

if [[ -f "$BAD_SEQ_FILE" ]]; then
	mv "$BAD_SEQ_FILE" "$OUTPUT_FOLDER/" 2>/dev/null
fi

INPUT_ALIGN="$filename.align"
INPUT_PICK_ALIGN="$filename.pick.align"
INPUT_INPUT_ALIGN="$filename.input.align"

if [[ -f "$INPUT_PICK_ALIGN" ]]; then
	mv "$INPUT_PICK_ALIGN" "$INPUT_INPUT_ALIGN"
elif [[ -f "$INPUT_ALIGN" ]]; then
	mv "$INPUT_ALIGN" "$INPUT_INPUT_ALIGN"
fi

echo "clustering"
mothur_commands_01="merge.files(input=$INPUT_INPUT_ALIGN-$REFERENCE_DB_ALIGNMENTS, output=$MERGED_ALIGNMENT)
set.current(fasta=$MERGED_ALIGNMENT)
unique.seqs()
filter.seqs()
dist.seqs(cutoff=$CLUSTER_CUTOFF, processors=$PROCESSORS, calc=onegap, countends=F)
cluster(cutoff=$CLUSTER_CUTOFF, method=furthest, precision=1000)"
#cluster.split(cutoff=$CLUSTER_CUTOFF, method=furthest, precision=1000, processors=$PROCESSORS)"
#cluster(cutoff=$CLUSTER_CUTOFF, method=nearest, precision=10000)"
echo -e "$mothur_commands_01" | "$MOTHUR"

#RESULT_FILE="*.list"
#if [[ ! -f "$RESULT_FILE" ]]; then
#	echo -e "Clustering failed, no result!"
#fi

if mv *.list "$OUTPUT_FOLDER/mothur_otu.txt"; then
# 2>/dev/null
#	#clean up folders, removed for testing
	echo "clean up folders"
	rm -rf "$TMP_DIR"
else
	echo "clustering failed"
	echo "intermediate data: $TMP_DIR"
fi
exit;
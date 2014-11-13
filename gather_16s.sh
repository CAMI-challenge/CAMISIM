#!/bin/bash
HOME_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#cd /net/metagenomics/projects/cami_2014/02_otu_clustering
#SOFTWARE_MANAGEMENT_REPOSITORY=/net/programs
#source /net/programs/software-management.source
importpackage parallel

#getting passed arguments
INPUT_GENOME_LIST_FILE_PATH="$1"
INPUT_REFERENCE_GENOME_LIST_FILE_PATH="$2"
OUTPUT_FOLDER="$3"
PROCESSORS="$4"
REFERENCE_AVAILABLE=0

# if referece "path_to_ref_genomefile" is used (boolean)
if [[ "$#" -eq 5 ]]; then
	REFERENCE_AVAILABLE=1
fi

# if test genome file "path_to_genomefile" does not exist
if [[ ! -f "$INPUT_GENOME_LIST_FILE_PATH" ]]; then
	echo -e "File does not exist: $INPUT_GENOME_LIST_FILE_PATH"
	exit;
fi

if [[ ! -f "$INPUT_REFERENCE_GENOME_LIST_FILE_PATH" ]]; then
	echo -e "File does not exist: $INPUT_REFERENCE_GENOME_LIST_FILE_PATH"
	exit;
fi

CURRENT_TIME=$(date +"%Y_%m_%d__%H.%M.%S")
TOOL_DIR="$HOME_DIR/tools"
#TMP_DIR="$HOME_DIR/1_tmp_$EUID$CURRENT_TIME"
#mkdir "$TMP_DIR"
TMP_DIR=`mktemp -d`
SYM_LINK_FOLDER="$TMP_DIR/sym_links"
WORKING_DIR="$TMP_DIR/working"
OUTPUT_DIR="$TMP_DIR/output"
#MARKER_GENES_DIR="$TMP_DIR/marker_genes"
DETECTOR16S_PATH="$TOOL_DIR/16SDetector"

# better as argument
SUFFIX28S="28S_rRNA.fna"
SUFFIX18S="18S_rRNA.fna"
SUFFIX08S="8S_rRNA.fna"
SUFFIX23S="23S_rRNA.fna"
SUFFIX16S="16S_rRNA.fna"
SUFFIX05S="5S_rRNA.fna"

# log file: marker gene found in a file?
MG_LOG_FILE="$OUTPUT_FOLDER/mg_log.txt"
# log file: stdout from 16S analysis
LOG_FILE="$OUTPUT_FOLDER/log.txt"
# log file: stderr from 16S analysis
ERROR_LOG_FILE="$OUTPUT_FOLDER/error_logfile.txt"

# list of fasta files
FASTA_LIST_FILE="$OUTPUT_FOLDER/fasta_list.txt"
# list of commands for "parallel" processing
COMMANDS_FILE="$OUTPUT_FOLDER/commands.txt"

if [[ -f "$COMMANDS_FILE" ]]; then
	rm "$COMMANDS_FILE"
fi
if [[ -f "$FASTA_LIST_FILE" ]]; then
	rm "$FASTA_LIST_FILE"
fi


## create temporary folders
#rm -rf "$TMP_DIR"
mkdir "$SYM_LINK_FOLDER"
#mkdir "$MARKER_GENES_DIR"
mkdir "$OUTPUT_FOLDER"
cd "$TMP_DIR"

# establish symbolic link to fasta files
function establish_symbolic
{
	echo "establish symbolic link to fasta files: $1"
	while read GENOME_ID FILE_PATH; do
		#echo "ID=$GENOME_ID PATH=$FILE_PATH"
		FILE_PATH=$([[ "$FILE_PATH" =~ [[:space:]]*([^[:space:]]|[^[:space:]].*[^[:space:]])[[:space:]]* ]]; echo -n "${BASH_REMATCH[1]}")
		if [[ ! -f "$FILE_PATH" ]]; then
			echo "File not found! ($FILE_PATH)"
			#cd "$HOME_DIR"
			#exit;
		else
			INPUT_FASTA_FILE=$(basename "$FILE_PATH")
			#INPUT_FILE_NAME="${INPUT_FASTA_FILE%.*}"
			SYMBOLIC_LINK="$SYM_LINK_FOLDER/$INPUT_FASTA_FILE"
			ln -s "$FILE_PATH" "$SYMBOLIC_LINK"
			echo -e "$GENOME_ID\t$INPUT_FASTA_FILE" >> "$FASTA_LIST_FILE"

			echo -e "$DETECTOR16S_PATH/run.py" -c "$TOOL_DIR/config.cfg" -nn -hmm 3 -i "$SYMBOLIC_LINK" -out "$TMP_DIR" >> "$COMMANDS_FILE"
		fi
	done < "$1"
}
establish_symbolic "$INPUT_GENOME_LIST_FILE_PATH"

if [[ "$REFERENCE_AVAILABLE" -eq 0 ]]; then
	establish_symbolic "$INPUT_REFERENCE_GENOME_LIST_FILE_PATH"
fi


if [[ ! -f "$COMMANDS_FILE" ]]; then
	echo -e "Could not create command-file for hmmsearch!"
	exit;
fi



echo -e "finding markergenes with hmmsearch\n"
cat "$COMMANDS_FILE" | parallel -j "$PROCESSORS" --no-notice 2> "$ERROR_LOG_FILE" > "$LOG_FILE"

# gather all marker genes into single files
function gather_genes
{
	SUFFIX="$1"
	OUT_FILE="$2"
	CUT_OFF="$3"
	GENE_NAME="$4"
	echo -e "gather $GENE_NAME genes into single file\n"
	while read GENOME_ID INPUT_FASTA_FILE; do
		GENE_FASTA_FILE="$TMP_DIR/working/$INPUT_FASTA_FILE.ids.$SUFFIX"
		if [[ -f "$GENE_FASTA_FILE" ]]; then
			"$TOOL_DIR/concat_fasta_on_fasta.py" -i "$GENE_FASTA_FILE" -o "$OUT_FILE" -id "$GENOME_ID" -c "$CUT_OFF" >> "$MG_LOG_FILE"
		else
			echo -e "$GENE_NAME not found in: $INPUT_FASTA_FILE" >> "$MG_LOG_FILE"
		fi
	done < "$FASTA_LIST_FILE"
}

#gather_genes "$SUFFIX05S" "$OUTPUT_FOLDER/$SUFFIX05S" "100" "5S"
gather_genes "$SUFFIX16S" "$OUTPUT_FOLDER/$SUFFIX16S" "900" "16S"
#gather_genes "$SUFFIX23S" "$OUTPUT_FOLDER/$SUFFIX23S" "1000" "23S"

#gather_genes "$SUFFIX08S" "$OUTPUT_FOLDER/$SUFFIX08S" "100" "8S"
#gather_genes "$SUFFIX18S" "$OUTPUT_FOLDER/$SUFFIX18S" "900" "18S"
#gather_genes "$SUFFIX28S" "$OUTPUT_FOLDER/$SUFFIX28S" "1000" "28S"

#add reference marker genes, if gathered in previous run
if [[ "$REFERENCE_AVAILABLE" -eq 1 ]]; then
	cp "$OUTPUT_FOLDER/$SUFFIX16S" "$OUTPUT_FOLDER/no_ref_$SUFFIX16S"
	cat "$INPUT_REFERENCE_GENOME_LIST_FILE_PATH" >> "$OUTPUT_FOLDER/$SUFFIX16S"
fi

#clean up symbolic links
echo "clean up symbolic links"
while read GENOME_ID INPUT_FASTA_FILE; do
	SYMBOLIC_LINK="$SYM_LINK_FOLDER/$INPUT_FASTA_FILE"
	rm "$SYMBOLIC_LINK"
done < "$FASTA_LIST_FILE"

#clean up folders
echo "clean up folders"
rm -rf "$TMP_DIR"
#rm -rf "$WORKING_DIR"
#rm -rf "$OUTPUT_DIR"

cd "$HOME_DIR"
exit;


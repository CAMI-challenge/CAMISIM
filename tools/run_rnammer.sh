EXECUTABLE="$1"
KINGDOM="$2"
MOLTYPE="$3"
OUT_GFF="$4"
OUT_FASTA="$5"
IN_FASTA="$6"
TEMP_DIR="$7"

cd "$TEMP_DIR"

"$EXECUTABLE" -S "$KINGDOM" -m "$MOLTYPE" -gff "$OUT_GFF" -f "$OUT_FASTA" - < "$IN_FASTA"
#"{0} -S {1} -m {2} -gff {3} -f {4} - < {5}".format(rnammer, kingdom, moltype, out_gff, out_fasta, in_fasta)


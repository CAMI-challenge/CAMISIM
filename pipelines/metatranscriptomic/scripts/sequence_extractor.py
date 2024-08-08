import gffutils
import pyfaidx

def extract_sequence(db, fasta, gene_identifier, child_feature_type):
    # Retrieve the gene feature from the database
    gene = db[gene_identifier]

    # Retrieve all child features (e.g., exons or CDS) for the gene
    children = list(db.children(gene, featuretype=child_feature_type, order_by='start'))

    child_sequences = []
    for child in children:
        # Extract the sequence for each child feature using gffutils
        # The sequence returned will be reverse-complemented for minus-strand features.
        child_seq = child.sequence(fasta, use_strand=True)

        if child_feature_type == 'CDS': # In this case consider the phase attribute
            # phase_attr = child.attributes.get('phase') if the 8th column is not named 'pahse', the return is 'None'
            phase_attr = child[7]
            phase = int(phase_attr) if phase_attr.isdigit() else 0

            # Adjust for phase: See gff3 documentation (http://gmod.org/wiki/GFF3) :
            # "For features of type "CDS", the phase indicates where the feature begins with reference to the reading frame. 
            # The phase is one of the integers 0, 1, or 2, indicating the number of bases that should be removed from the beginning of this feature to reach the first base of the next codon. 
            # In other words, a phase of "0" indicates that the next codon begins at the first base of the region described by the current line, a phase of "1" indicates that the next codon begins at the second base of this region, and a phase of "2" indicates that the codon begins at the third base of this region. 
            # (... )If there is no phase, put a "." (a period) in this field. 
            # For forward strand features, phase is counted from the start field. For reverse strand features, phase is counted from the end field."
            
            # Since the sequence is already reverse-complemented for minus-strand features, just cut remove from the beginning.
            child_seq = child_seq[phase:] if phase > 0 else child_seq
        
        child_sequences.append(child_seq)
    
    # Concatenate child sequences in the correct order for the strand
    if gene.strand == '+':
        combined_sequence = ''.join(child_sequences)
    elif gene.strand == '-':
        combined_sequence = ''.join(reversed(child_sequences))

    return gene.seqid, gene.start, combined_sequence

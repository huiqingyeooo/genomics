# This script counts the number of basepairs in a fasta file (excluding gaps)
from Bio import SeqIO

def count_basepairs(fasta_file):
    total_basepairs = 0
    # Iterate over all sequences in the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Remove gaps ('-' or 'N') and count base pairs
        sequence = record.seq
        cleaned_sequence = sequence.replace('-', '').replace('N', '')
        total_basepairs += len(cleaned_sequence)
    
    return total_basepairs

# Provide your FASTA file path here
fasta_file = "concatenated.fasta"
total_basepairs = count_basepairs(fasta_file)
print(f"Total number of base pairs (excluding gaps): {total_basepairs}")

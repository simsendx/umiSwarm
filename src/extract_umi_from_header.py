# Import libraries
import argparse
from Bio import SeqIO # Requires biopython

# Define UMI prcoessing function for fasta files
def process_umi(fasta, output):
    """Extract UMIs from fasta header and add to sequence."""
    # Initilialize variables
    result_file = open(output, "w")
    i=0 # Counter for no records
    j=0 # Counter for no. of records with UMI

    # Open fasta file
    for record in SeqIO.parse(fasta, "fasta"):
        i += 1
        header = record.id
        umi = header.split(':')[-1]
        sequence = record.seq

        record.seq = umi + sequence

        # Skip records that contain N bases in the UMI
        if "N" in umi:
            j += 1
            continue

        SeqIO.write(record, result_file, "fasta")
        
    print("Processed " + str(i) + " records. Found " + str(j) + "UMIs with Ns.")

#----------// Main //-------------

# Get command line arugments
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", type=str)
parser.add_argument("-o", "--output", type=str)

args = parser.parse_args()
fasta = args.fasta
output = args.output

# Run functions
process_umi(fasta, output)
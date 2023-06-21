# Import libraries
import argparse
from Bio import SeqIO # Requires biopython

def addUmiToHeader(infasta, outfasta, n_umi = 19):
    """ Add UMI from fasta sequence to header """

    # Initilialize variables
    i = 0 # Counter for no records

    outfile = open(outfasta, "w")

    # Open fasta file
    for record in SeqIO.parse(infasta, "fasta"):
        i += 1
        # Get header sequence
        header = record.id
        # Get read sequence
        sequence = record.seq

        # Extract UMI 
        umi = sequence[:n_umi]

        # Extract amplicon sequence without umi, primers and adapters
        new_sequence = sequence[n_umi:]

        # Generate new header sequence
        new_header = header + "umi=" + umi + ";"

        # Write new fasta record
        outfile.write(">" + str(new_header) + "\n")
        outfile.write(str(new_sequence) + "\n")

    # Output to console
    print("Processed " + str(i) + " records.")


#----------// Main //-------------

# Get command line arugments
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--fastain", type=str)
parser.add_argument("-o", "--fastaout", type=str)

args = parser.parse_args()
fasta = args.fastain
output = args.fastaout

# Run script
addUmiToHeader(fasta, output, n_umi = 19)
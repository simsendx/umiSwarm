
import argparse
import re
import pandas as pd
from Bio import SeqIO # Requires biopython

def fastaToOTU(fasta):
    """ Add UMI from fasta sequence to header """
    # Initilialize variables
    i = 0 # Counter for no records

    df = pd.DataFrame(columns = ['Read', 'Abundance'])

    # Open fasta file
    for record in SeqIO.parse(fasta, "fasta"):
        i += 1
        header = record.id
        sequence = record.seq

        # Extract UMI family size
        match = re.search(r'size=(\d+)', header)
        if match:
            size = match.group(1)

            new_record = pd.DataFrame([{'Read': str(sequence), 'Abundance': size}])
            
            # Add row to dataframe
            df = pd.concat([df, new_record], axis = 0, ignore_index=True)

    return(df)


#----------// Main //-------------
#
# Get command line arugments
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", type=str)
parser.add_argument("-o", "--otu", type=str)

args = parser.parse_args()
fasta = args.fasta
otu_out = args.otu

# Run functions
cons_table = fastaToOTU(fasta)

# Save consensus table as csv
cons_table.to_csv(otu_out, index = False)
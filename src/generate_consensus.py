## Generate consensus reads
#
# Merge UMI families and generaqte consensus reads.
#
# This scripts returns an abudance table with 3 columns:
#
# UMI sequence
# Consensus read
# Abudance/count of reads with that UMI
#
# Date: 2023-04-18
# Author: Stefan Filges

# Import libraries
import argparse
import re
import pandas as pd
from Bio import SeqIO # Requires biopython

def generate_abundance_table(fasta, n_umi = 19):
    """ Add UMI from fasta sequence to header """
    # Initilialize variables
    i = 0 # Counter for no records
    df = pd.DataFrame(columns = ['UMI', 'Read', 'Size'])

    # Open fasta file
    for record in SeqIO.parse(fasta, "fasta"):
        i += 1
        header = record.id
        sequence = record.seq
        umi = sequence[:n_umi]
        new_sequence = sequence[n_umi:]

        # Extract UMI family size
        match = re.search(r'size=(\d+)', header)
        if match:
            size = match.group(1)

        new_row = pd.DataFrame([{'UMI': str(umi), 'Read': str(new_sequence), 'Size': size}])

        # Add row to dataframe
        df = pd.concat([df, new_row], axis = 0, ignore_index=True)
    
    print("Processed " + str(i) + " records.")
    return df


# Define UMI prcoessing function for fasta files
def generate_consensus(fasta, n_umi = 19, threshold = 20):
    """Merge consensus reads."""

    # Get UMIs, read sequences and UMI counts
    df = generate_abundance_table(fasta, n_umi = n_umi)

    df_out = pd.DataFrame(columns = ['UMI', 'Read', 'Size'])

    umi_database = df['UMI'].unique()
    print(umi_database)

    for umi in umi_database:
        # Get UMI family = all reads with the same UMI
        family = df[df['UMI'] == umi]

        # Sort by size
        family = family.sort_values(by=['Size'], axis = 0, ascending = False)

        # Sum of UMI counts
        size = family['Size'].astype(int).sum(axis = 0)

        if(size <= threshold):
            next()
                
        # Get major read sequence
        centroid = family['Read'].iloc[0]
        
        # Define a new row for the output
        new_row = pd.DataFrame([{'UMI': str(umi), 'Read': str(centroid), 'Size': size}])
        
        # Append new row to output dataframe
        df_out = pd.concat([df_out, new_row], axis = 0, ignore_index=True)

    print(df_out.head())
    return df_out

#----------// Main //-------------

# Get command line arugments
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", type=str)
parser.add_argument("-o", "--output", type=str)
parser.add_argument("-u", "--umi_length", type=int)

args = parser.parse_args()
fasta = args.fasta
output = args.output
n_umi = args.umi_length

# Run functions
cons_table = generate_consensus(fasta, n_umi = n_umi)

# Save consensus table as csv
cons_table.to_csv(output)
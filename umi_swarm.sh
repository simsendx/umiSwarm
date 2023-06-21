#!/bin/bash

# A wrapper script for alignment free analysis of metabarcoding
# data containing UMIs.

#####################################################
#----------- Define command line options -----------#
#####################################################
display_help() {
    echo "$(basename "$0") [-h] [-i n] [options]" >&2
    echo
    echo "   -h, --help          Display this helpful screen."
    echo "   -i, --input-dir     Input directory for fastq files. Default is current working directory."
    echo "   -u --umi_length     UMI length, default is 19."
    echo "   -t --threads        Number of threads to use. Default is 12"
    echo "   -d --distance       Distance value to use for swarm. Default is 1."
    echo
    exit 1
}

#####################################################
#------------- Initilialize parameters -------------#
#####################################################

# Directories for python scripts
UMISWARM="/mnt/c/Users/Stefan/Documents/GitHub/umiSwarm/src/extract_umi_from_header.py"
CONS="/mnt/c/Users/Stefan/Documents/GitHub/umiSwarm/src/generate_consensus.py"
addUmiToHeader="/mnt/c/Users/Stefan/Documents/GitHub/umiSwarm/src/add_umi_to_header.py"

# Set defaults for general parameters
umi_length=19  # UMI-length to tranfer to header
umi=true       # Use UMI in swarm
ncore=12       # Number of threads to use
multiqc=false  # do not run multiqc unless it is installed
fastqc=false   # do not run fastqc unless it is installed
distance=1     # distance value to use with swarm

#####################################################
#----------------- Get user input ------------------#
#####################################################

# Use getops to get command line input supplied by the user
while getopts ':hfi:u:t:d:' option; do
  case "$option" in
    h | --help)
        display_help
        exit 0
        ;;
    i | --input-dir)
        FILES=$OPTARG
        ;;
    u | --umi_length)
        umi_length=$OPTARG
        ;;
    t | --threads)
        ncore=$OPTARG
        ;;
    d | --distance)
        distance=$OPTARG
        ;;
    :) printf "missing argument for -%i\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done
shift $((OPTIND - 1))

#####################################################
###-------------- Check dependencies -------------###
#####################################################

# Define shell output colors
YELLOW=$(tput setaf 3)
GREEN=$(tput setaf 2)
RED=$(tput setaf 1)
NC=$(tput sgr0)

# Check working directory
printf '\n%s %s %s\n' $GREEN "Checking working directory..." $NC
if [[ $FILES = "" ]]
  then
    printf '%s %s %s\n' $YELLOW "...specified input directory does not exist." $NC
    exit
else
  printf '%s %s %s %s\n' $YELLOW "...working directory is: " $NC $FILES 
  runDir=$FILES
fi

# Begin dependency check
printf '%s %s %s\n' $GREEN "Checking dependencies..." $NC

# Check if fastp is installed
if ! command -v fastp &> /dev/null
  then
    printf '%s %s\n' $YELLOW "...fastp could not be found."
    printf '%s %s\n' "Please install fastp: " "https://github.com/OpenGene/fastp. Exiting."
    exit
  else
    printf '%s %s %s\n' $YELLOW "...fastp is installed." $NC
  fi

# Check if cutadapt is installed
if ! command -v cutadapt &> /dev/null
  then
    printf '%s %s\n' $YELLOW "...cutadapt could not be found."
    printf '%s %s\n' "Please install cutadapt. Exiting."
    exit
  else
    printf '%s %s %s\n' $YELLOW "...cutadapt is installed." $NC
  fi

# Check if multiqc is installed
if ! command -v multiqc &> /dev/null
  then
    printf '%s %s %s\n' $YELLOW "...multiqc could not be found." $NC
    printf ' %s\n' "Please install multqic"
  else
    printf '%s %s %s\n' $YELLOW "...multiqc is installed." $NC
    multiqc=true
  fi

# Check if fastqc is installed
if ! command -v fastqc &> /dev/null
  then
    printf '%s %s %s\n' $YELLOW "...fastqc could not be found." $NC
    printf ' %s\n' "Please install fastqc"
  else
    printf '%s %s %s\n' $YELLOW "...fastqc is installed." $NC
    fastqc=true
  fi

# Check if vsearch is installed
if ! command -v vsearch &> /dev/null
  then
    printf '%s %s %s\n' $YELLOW "...vsearch could not be found. " $NC
    printf ' %s\n' "Please install vsearch"
    exit
  else
    printf '%s %s %s\n' $YELLOW "...vsearch is installed." $NC
  fi

# Check if swarm is installed
if ! command -v swarm &> /dev/null
  then
    printf '%s %s %s\n' $YELLOW "...swarm could not be found." $NC
    printf ' %s\n' "Please install swarm"
    exit
  else
    printf '%s %s %s\n' $YELLOW "...swarm is installed." $NC
  fi

# Check if pigz is installed
if ! command -v pigz &> /dev/null
  then
    printf '%s %s %s\n' $YELLOW "...pigz could not be found." $NC
    printf ' %s\n' "Please install pigz"
    exit
  else
    printf '%s %s %s\n' $YELLOW "...pigz is installed." $NC
  fi

# Report all dependencies installed
printf '%s %s %s\n' $GREEN "All dependencies installed." $NC

#####################################################
###------------ Processing fastq files -----------###
#####################################################

# Enter working directory
cd $FILES

# Process fastq files in input directory
for fastq in *.fastq.gz
do
  PREFIX=${fastq%%[_]*} # Get sample name
  if [[ $fastq =~ R1 ]] # If file name contains R1, proceed
  then
    # define R1
    fq1=$fastq 
    printf '%s \n' $fq1

    # define R2 by replacing R1 with R2
    fq2=${fastq//R1/R2} 
    printf '%s \n' $fq2
            
    # Perform adapter trimming, R1/R2 merging and UMI processing using fastp, quality filtering
    fastp --in1=$fq1 --in2=$fq2 \
      --merge \
      --unpaired1="${PREFIX}_unpaired.fastq.gz" \
      --unpaired2="${PREFIX}_unpaired.fastq.gz" \
      --merged_out="${PREFIX}_merged.fastq.gz" \
      --failed_out="${PREFIX}_failed.fastq.gz" \
      --detect_adapter_for_pe \
      --trim_poly_g \
      --disable_quality_filtering \
      --umi \
      --umi_loc=read1 \
      --umi_len=$umi_length \
      --thread=$ncore \
      --correction \
      --overlap_len_require=12  \
      --length_required=75 \
      --json="${PREFIX}.json" \
      --html="${PREFIX}.html" \
      --report_title=${PREFIX}

    # Trim primer sequences using cudadapt (forward primer contains SiMSen-Seq stem adapter)
    cutadapt -g ATGGGAAAGAGTGTCCGTCGGTAAAACTCGTGCCAGC \
      -a CAAACTGGGATTAGATACCCCACTATG \
      -o "${PREFIX}_merged_noPrimers.fastq.gz" \
      --discard-untrimmed \
      -j $ncore \
      --minimum-length 35 "${PREFIX}_merged.fastq.gz"

      # Convert to fasta
      # This does not work, because swarm doe not like the unknown base "N"
      #gunzip -c "${PREFIX}_merged_noPrimers.fastq.gz" | \
      #  awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' \
      #  > "${PREFIX}_merged_noPrimers.fasta"

      # Use vsearch to convert to fasta, removing "N" bases
      vsearch \
        --fastq_maxns 0 \
        --fastq_filter "${PREFIX}_merged_noPrimers.fastq.gz" \
        --fastaout "${PREFIX}_merged_noPrimers.fasta"

      # If UMI should be used:
      if $umi
        then
        # Add UMI back to main sequence using python script (requires biopython):
        ## conda install -c conda-forge biopython
        ## pip install biopython
        ## sudo apt-get install python-biopython-doc

        printf '%s \n' 'Running python to extract UMIs.'
        python3 "$UMISWARM" \
          --fasta "${PREFIX}_merged_noPrimers.fasta" \
          --output "${PREFIX}_merged_noPrimers_umi_in_sequence.fasta"

        # Dereplicate unique sequences using vsearch
        # This merges all unqiue sequences and adds the abundance to the header
        vsearch \
          --derep_fulllength "${PREFIX}_merged_noPrimers_umi_in_sequence.fasta" \
          --sizeout \
          --relabel_sha1 \
          --fasta_width 0 \
          --output "${PREFIX}_merged_noPrimers_dereplicated.fasta"

        # Compress input fasta
        pigz "${PREFIX}_merged_noPrimers_umi_in_sequence.fasta"
          
        # If umi should not be used:
        else
          # Dereplicate unique sequences using vsearch
          # This merges all unqiue sequences and add the abundance to the header
          vsearch \
            --derep_fulllength "${PREFIX}_merged_noPrimers.fasta" \
            --sizeout \
            --relabel_sha1 \
            --fasta_width 0 \
            --output "${PREFIX}_merged_noPrimers_dereplicated.fasta"
      fi

      # Compress input fasta file
      pigz "${PREFIX}_merged_noPrimers.fasta"

      # Run swarm
      swarm \
        -f \
        -t $ncore \
        -z \
        -d $distance \
        -s $PREFIX.stats \
        -w "${PREFIX}_swarm.fasta" \
        -o s$PREFIX.swarm "${PREFIX}_merged_noPrimers_dereplicated.fasta"

      # Add UMI back to header
      # UMI in sequence might interfere with chimera detection?
      python3 "$addUmiToHeader" \
        --fastain "${PREFIX}_swarm.fasta" \
        --fastaout "${PREFIX}_swarm_umi_in_header.fasta"

      # Compress output fasta files
      #pigz "${PREFIX}_merged_noPrimers_dereplicated.fasta"
      #pigz "${PREFIX}_swarm.fasta"

      # make dereplicate representatives, merge on only marker region.
      vsearch \
        --sizeout \
        --derep_fulllength "${PREFIX}_swarm_umi_in_header.fasta" \
        --output "${PREFIX}_swarm_umi_in_header_dereplicated.fasta" 

      # removes singletons and small clusters, and can be used directly in dada2 due to correct fastawidth
      vsearch \
        --fastx_filter "${PREFIX}_swarm_umi_in_header_dereplicated.fasta" \
        --fastaout "${PREFIX}_swarm_umi_in_header_dereplicated_filtered.fasta"  \
        --minsize 20 \
        --fasta_width 0 

      # Use chime2 to detect de novo chimeras
      vsearch \
        --uchime2_denovo "${PREFIX}_swarm_umi_in_header_dereplicated_filtered.fasta" \
        --chimeras "${PREFIX}_chim.fasta" \
        --nonchimeras "${PREFIX}_non_chim.fasta"

      # Merge reads by UMI and create abundance table
      #printf "%s \n" "Merging consensus families."

      #python3 "$CONS" \
      #  --fasta "${PREFIX}_non_chim.fasta" \
      #  --output "${PREFIX}_cons_otu.csv" \
      #  --umi_length 19

    fi
done

##############################################
###------ Merge consensus tables ----------###
##############################################

# TODO
# Write script to merge multiple consensus tables into a single matrix
# with sequences are rownames, sample names as column names and 
# UMI familiy sizes as values. This can then be used for downstream
# processing with phyloseq.

##############################################
###---------- Annotate sequences ----------###
##############################################

#download reference database
#wget https://raw.githubusercontent.com/EivindStensrud/ScandiFish/main/ScandiFish_12s_v1.2/ScandiFish_12s_v1.4_nf.fasta # Downloads database
#wget https://raw.githubusercontent.com/EivindStensrud/ScandiFish/main/ScandiFish_12s_v1.2/ScandiFish_12s_v1.4_nf.fasta.md5 # Downloads md5sum

#md5sum -c ScandiFish_12s_v1.4.fasta.md5 # Checks if database is correct.



# Run nucleotide BLAST
# BLASTN compares the sequences and keeps up to 100 reference sequences.
# Input is DADA2 output from function uniquesToFasta(seqtab.nochim).

#blastn \
#  -max_target_seqs 100 \
#  -evalue 0.01 \
#  -query /cluster/projects/nn9745k/02_results/40_simsenseq/umi_and_seq/rep-seqs.fna \
#  -out /cluster/projects/nn9745k/02_results/40_simsenseq/umi_and_seq/umi_and_seqs_output_blast_results \
#  -db /cluster/projects/nn9745k/03_databases/fish/ScandiFish_12s_v1.4/ScandiFish_12s_v1.4.fasta_db \
#  -outfmt 6 \
#  -num_threads $ncore 

##############################################
###------------ Phyloseq report -----------###
##############################################




##############################################
###------------- File cleanup -------------###
##############################################

# Check if fastqc folder exists in working directory
if [[ -d fastp_reports ]] 
then
    echo "Fastp report folder exists."
else
    echo "Fastp report folder does not exist. Creating a new one."
    mkdir fastp_reports
fi

# Move fastp reports to fastp report folder
mv *html *json fastp_reports

# Check if fastqc is installed before generating reports
if $fastqc
    then
    # Check if fastqc folder exists in working directory
    if [[ -d fastqc_reports ]] 
    then
        echo "Fastqc folder exists."
    else
        echo "Fastqc folder does not exist. Creating a new one."
        mkdir fastqc_reports
    fi

    # run fastqc and store reports in fastqc report folder
    fastqc -t 12 -o fastqc_reports *_merged_noPrimers.fastq.gz
fi

# Check if multiqc is installed before summarizing reports
if $multiqc
then
    # collate reports using multiqc
    multiqc fastqc_reports
fi

# Check if fastqc folder exists in working directory
if [[ -d sequence_archive ]] 
then
    echo "Folder for sequence archive exists."
else
    echo "Archive folder does not exist. Creating a new one."
    mkdir sequence_archive
fi

# Check if fastqc folder exists in working directory
if [[ -d reports ]] 
then
    echo "Folder for reports exists."
else
    echo "Reports folder does not exist. Creating a new one."
    mkdir reports
fi

# Move reports to single common folder
mv fastqc_reports fastp_reports multiqc* reports

# Move intermediate fastq and fasta to archive
mv *fasta.gz *unpaired* *failed* sequence_archive
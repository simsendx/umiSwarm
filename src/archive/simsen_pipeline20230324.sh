
######################################
### Fish eDNA simsen seq 20230313


RAWDIR=/cluster/projects/nn9745k/01_raw_data
RESDIR=/cluster/projects/nn9745k/02_results

# transfer files to scratch

mkdir $RESDIR/40_simsenseq
mkdir $RESDIR/40_simsenseq/figs


##############################################
### merge fwd and rev reads ###

cd $RAWDIR/40_simsenseq/Data


# RUN flash
# RUN filtering in dada2 first

module purge
module load Flash2/2.2.00-GCC-10.3.0

#optimalize the min-overlap?
for f in *.fastq.gz
    do
        FILE=${f#$DIRS}
        PREFIX=${FILE%%[_]*}
        if [[ $FILE =~ R1 ]] #If file name contains R1
        then
            flash2 $FILE ${FILE//R1/R2} --min-overlap=12 --compress --output-prefix=$PREFIX 2>&1 | tee flash.log
        fi
    done

# make fasta file, VSEARCH do not like fastq
module purge
module load VSEARCH/2.21.1-GCC-10.3.0
## with fastq_filter, size restriction may be applied.
for f in *.extendedFrags.fastq.gz
    do
        FILE=${f#$DIRS}
        PREFIX=${FILE%%[.]*}
        vsearch --fastq_maxns 0 --fastq_filter $FILE --fastaout $RESDIR/40_simsenseq/no_denoising/$PREFIX.filtered.fasta.gz
    done

cd $RESDIR/40_simsenseq/no_denoising

# dereplication
## Identical UMIs + sequences regions are merged. Get number of reads.

for f in *.filtered.fasta.gz
    do
        FILE=${f#$DIRS}
        PREFIX=${FILE%%[.]*}
        vsearch --sizeout --derep_fulllength $FILE --output $PREFIX.dereplicated.fasta.gz
    done



# RUN SWARM
## Running SWARM on UMI+marker region should take advantage of UMIs
## Most common/center cluster should be the correct one
## Cluster sequences on full length, detect sequencing errors using local threshold
## Important, only peaks with >100 amplicons are considered peaks!
## Consideres shallow sequencing? Differences, d = >1 can be used.
## -f --fastidious can increase the quality, but only with d=1, makes a "artificial" amplicon to connect networks. Increase quality.
## boundary, b, change what is regarded as a large cluster. Default >=3
## -w keeps only center cluster
module purge
module load swarm/3.0.0-GCC-9.3.0
for f in *.filtered.fasta.gz
    do
        FILE=${f#$DIRS}
        PREFIX=${FILE%%[.]*}
        swarm -z -s $PREFIX.stats -w $PREFIX.swarm.fasta -o s$PREFIX.swarm $FILE
    done


# RUN cutadapt

### cutadapt
module purge
module load cutadapt/1.18-foss-2018b-Python-3.6.6


for f in *.swarm.fasta
    do
        FILE=${f#$DIRS}
        PREFIX=${FILE%%[.]*}
        cutadapt -o $PREFIX.swarm.filtered.fasta --minimum-length 150 $FILE
    done

### removes UMIs and primers

for f in *.swarm.filtered.fasta
    do
        FILE=${f#$DIRS}
        PREFIX=${FILE%%[.]*}
        cutadapt -o $PREFIX.trimmed.fasta -g GTCGGTAAAACTCGTGCCAGC...CAAACTGGGATTAGATACCCCACTATG --discard-untrimmed $FILE
    done

# make dereplicate representatives, merge on only marker region.
##Sort by size for streamline for chimeara detection.
module purge
module load VSEARCH/2.21.1-GCC-10.3.0

for f in *.trimmed.fasta
    do
        FILE=${f#$DIRS}
        PREFIX=${FILE%%[.]*}
        vsearch --sizeout --derep_fulllength --sortbysize $FILE --output $PREFIX.dereplicated.trimmed.fasta
    done

#de novo chimearas detection 
## test with uchime2_denovo, uchime_denovo and maybe some more.

for f in *.dereplicated.trimmed.fasta
	do 
		FILE=${f#$DIRS}
		PREFIX=${FILE%%[.]*}
		vsearch --uchime2_denovo $FILE --chimeras $PREFIX.chim.fasta --nonchimeras $PREFIX.non_chim.fasta
	done

        


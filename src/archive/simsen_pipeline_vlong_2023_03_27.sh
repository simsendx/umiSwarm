####Script to translocate UMI's to FASTQ header
### https://unix.stackexchange.com/questions/510164/remove-and-add-sequence-information-at-specific-position-in-a-file

#currently only tested on one sample. 
#Gives a result, but not validated the result.
#Not streamlined for processing several samples at the same time yet.
							 
salloc --ntasks=1 --mem-per-cpu=4G --time=00:30:00 --qos=devel --account=nn9745k 

module purge
module load Flash2/2.2.00-GCC-10.3.0

RAWDIR=/cluster/projects/nn9745k/01_raw_data/40_simsenseq/test_Stensrud_SimSen
mkdir /cluster/projects/nn9745k/01_raw_data/40_simsenseq/test_Stensrud_SimSen/Merged
RESDIR=/cluster/projects/nn9745k/01_raw_data/40_simsenseq/test_Stensrud_SimSen/Merged
cd $RAWDIR

## Merges R1 and R2 reads
## mat alter min-overlap to larger numbers
for f in *.fastq.gz
    do
        FILE=${f#$DIRS}
        PREFIX=${FILE%%[_]*}
        if [[ $FILE =~ R1 ]] #If file name contains R1
        then
            flash2 $FILE ${FILE//R1/R2} --min-overlap=12 --compress --output-directory=$RESDIR --output-prefix=$PREFIX 2>&1 | tee flash.log

            
        fi
    done	
	

## Removal of UMI and translocate to header
##Need to change directory, and change input file to only take samples which want to be taken

#working
#MiFish-UF-v1, UMI region between adapter and primer sequence: NNNNAAANNANNAAANNNNATGGGAAAGAGTGTCC, length 35
#MiFish-UF-v2, UMI region between adapter and primer sequence: NNNNAAANNANNAAANNNNATGGGAAAGAGTGTCCNNNNNN, length 41

# FILE=${f#$DIRS} do not work properly
module purge

RAWDIR=/cluster/projects/nn9745k/01_raw_data/40_simsenseq/test_Stensrud_SimSen/Merged
RESDIR=/cluster/projects/nn9745k/01_raw_data/40_simsenseq/test_Stensrud_SimSen/UMI_header
cd $RAWDIR

for f in *extendedFrags.fastq.gz
    do
        FILE=${f#$DIRS} #fil
        if [[ $FILE =~ v1 ]]
        then
            zcat $FILE | awk '(FNR-1)%4 == 0 { illu=$1; adapt=$2; next }
                (FNR-2) %4 == 0 { UMI = substr($0,1,35)
                print name, chr, len, UMI ;
                seq = substr($0,36)
                print seq ; next }
                (FNR-3) %4 == 0 {print $0 ; next}
                (FNR-4) %4 == 0 {qseq = substr($0,36)
                print qseq ; next}
                    
                '| gzip -c >$RESDIR/UMI_$FILE
       
        fi
    done



for f in *extendedFrags.fastq.gz
    do
        FILE=${f#$DIRS} #fil
        if [[ $FILE =~ v2 ]]
        then
            zcat $FILE | awk '(FNR-1)%4 == 0 { illu=$1; adapt=$2; next }
                (FNR-2) %4 == 0 { UMI = substr($0,1,41)
                print illu, adapt, UMI ;
                seq = substr($0,42)
                print seq ; next }
                (FNR-3) %4 == 0 {print $0 ; next}
                (FNR-4) %4 == 0 {qseq = substr($0,42)
                print qseq ; next}
                    
                '| gzip -c >$RESDIR/UMI_$FILE
       
        fi
    done



module load cutadapt/1.18-foss-2018b-Python-3.6.6
# RUN cutadapt

##############################################
### run cutadapt to remove primer sequences
##############################################
RAWDIR=/cluster/projects/nn9745k/01_raw_data/40_simsenseq/test_Stensrud_SimSen/UMI_header
mkdir $RAWDIR/AdaptersRemoved
RESDIR=/cluster/projects/nn9745k/01_raw_data/40_simsenseq/test_Stensrud_SimSen/AdaptersRemoved

#cd $RAWDIR/40_simsenseq/Data
cd $RAWDIR

for f in *extendedFrags.fastq.gz
    do
        FILE=${f#$DIRS}
        PREFIX=${FILE%%[.]*}
        if [[ $FILE =~ UMI ]]
        then

            cutadapt -o $RESDIR/$PREFIX.trimmed.fastq.gz -g GTCGGTAAAACTCGTGCCAGC...CAAACTGGGATTAGATACCCCACTATG --discard-untrimmed --minimum-length 35 $FILE
        fi     
    done


RAWDIR=/cluster/projects/nn9745k/01_raw_data/40_simsenseq/test_Stensrud_SimSen/AdaptersRemoved
RESDIR=/cluster/projects/nn9745k/01_raw_data/40_simsenseq/test_Stensrud_SimSen/AdaptersRemoved
cd $RAWDIR

###
###Demultiplexes on UMI header. One file per UMI. Not sure what happens if UMI is not unique. May affect the result?
###Make Swarm for each file? Keep centroids? May result in several peaks per file? Keep centroids from Swarm

module purge
module load BBMap/38.98-GCC-11.2.0


for f in *.fastq.gz
	do
        FILE=${f#$DIRS}
        PREFIX=${FILE%%[_]*}
        if [[ $FILE =~ trimmed ]] #If file name contains R1
        then
            demuxbyname.sh in=$FILE out=demux%_$PREFIX.fastq.gz delimiter=space prefixmode=f
        fi
    done


#Dereplicate all files, no N's, report errors, sort, and output should be fasta.
module purge
module load VSEARCH/2.21.1-GCC-10.3.0

for f in *.fastq.gz
    do
        FILE=${f#$DIRS}
        PREFIX=${FILE%%[.]*}
        if [[ $FILE =~ demux ]]
        then
            vsearch --fastq_maxns 0 --fastq_filter $FILE --fastaout $PREFIX.dereplicated.fasta.gz
        fi
    done

#dereplicates on each UMI
#Then make Swarm on each file?

for f in *.dereplicated.fasta.gz
    do
        FILE=${f#$DIRS}
        PREFIX=${FILE%%[.]*}
        if [[ $FILE =~ demux ]]
        then
            vsearch --sizeout --derep_fulllength $FILE --output $PREFIX.dereplicated2.fasta.gz
        fi
    done

#For each fastafile, keep centroids. Should keep different sequences from same UMI.
module load swarm/3.0.0-GCC-9.3.0


#-w sort decreasing cluster abundance, easier for de novo chimera detection.
#-f fastidious, performes a second clustering, to reduce number of small clusters. Maybe be bad to have?
#-s, statistics file
#-z usearch abundance style.
#due to shallow sequencing per UMI, d=2 can be tested

for f in *.dereplicated2.fasta.gz
    do
        FILE=${f#$DIRS}
        PREFIX=${FILE%%[.]*}
        swarm -z -f -s Swarm_out/$PREFIX.fastidious.stats -w Swarm_out/$PREFIX.fastidious.swarm.fasta -o Swarm_out/$PREFIX.fastidious.swarm $FILE
    done

###
#Chimaera detection
###
###De Novo chimaera detection 
module purge
module load VSEARCH/2.21.1-GCC-10.3.0

#uchime2_denovo, parent is atleast 2 times more present than chimera, uchime3_denovo uses 16.
#do not detect chims atm, could be due to centroids filtering them out.
#Maybe need to collaps all fasta files.
#some fasta sequences have ~94% identity to herrings, something fishy here.


for f in *.fastidious.swarm.fasta
    do
        FILE=${f#$DIRS}
        PREFIX=${FILE%%[.]*}
        if [[ $FILE =~ demux ]]
        then
            vsearch --uchime2_denovo $FILE --chimeras $PREFIX.chim.fasta --nonchimeras $PREFIX.non_chim.fasta
        fi
    done

		



#Calib is a potential, but not sure if they are suitable, bad readme.




######################################
### Fish eDNA simsen seq 20230313
salloc --ntasks=1 --mem-per-cpu=4G --time=00:30:00 --qos=devel --account=nn9745k 

### cutadapt
module load cutadapt/1.18-foss-2018b-Python-3.6.6


## Amphibians

RAWDIR=/cluster/projects/nn9745k/01_raw_data
RESDIR=/cluster/projects/nn9745k/02_results

# transfer files to scratch

mkdir $RESDIR/40_simsenseq
mkdir $RESDIR/40_simsenseq/AdaptersRemoved
mkdir $RESDIR/40_simsenseq/figs


##############################################
### run cutadapt to remove primer sequences
##############################################

cd $RAWDIR/40_simsenseq/Data

# RUN cutadapt
# Standard eDNA approach

for f in *.fastq.gz
    do
        FILE=${f#$DIRS}
        if [[ $FILE =~ DNA*R1 ]] #If file name contains R1
        then
            cutadapt -g GTCGGTAAAACTCGTGCCAGC -G CATAGTGGGGTATCTAATCCCAGTTTG -o $RESDIR/40_simsenseq/AdaptersRemoved/$FILE -p $RESDIR/40_simsenseq/AdaptersRemoved/${FILE//R1/R2} --discard-untrimmed --minimum-length 35 $FILE ${FILE//R1/R2}

        fi
    done


#####
# To keep UMIs
# RUN filtering in dada2 first
# cutadapt removes only reverse primer.
## not working properly, continue w/ FLASH2


"""cd $RAWDIR/40_simsenseq/test_Stensrud_SimSen

for f in *.fastq.gz
    do
        FILE=${f#$DIRS}
        if [[ $FILE =~ R1 ]] #If file name contains R1
        then
            cutadapt -G CATAGTGGGGTATCTAATCCCAGTTTG -o $RESDIR/40_simsenseq/test/cutadapt_23_3_23/$FILE -p $RESDIR/40_simsenseq/test/cutadapt_23_3_23/${FILE//R1/R2} --discard-untrimmed --minimum-length 35 $FILE ${FILE//R1/R2}

        fi
    done"""
	
### Rearrange fastq file, so UMI is kept in header.
awk -v n=6 '(FNR-1) % 2 == 0 { name=$1; chr=$2; len=$3; next }
            (FNR-2) % 4 == 0 { seq=substr($0,1,n) }
                             { print name "." seq, chr, len
                               print substr($0,n+1) }'


	
cd $RAWDIR/40_simsenseq/test_Stensrud_SimSen

module purge
module load Flash2/2.2.00-GCC-10.3.0

### works, but needs to be done on raw reads. Otherwise UMIs is gone

for f in *.fastq.gz
    do
        FILE=${f#$DIRS}
        PREFIX=${FILE%%[_]*}
        if [[ $FILE =~ R1 ]] #If file name contains R1
        then
            flash2 $FILE ${FILE//R1/R2} --min-overlap=12 --compress --output-prefix=$PREFIX 2>&1 | tee flash.log

            # flash2 simsen-MiFish-v2-1ng-mock-6cyc-3_S12_L001_R1_001.fastq.gz simsen-MiFish-v2-1ng-mock-6cyc-3_S12_L001_R2_001.fastq.gz --min-overlap=12 --compress --output-prefix=$FILE 2>&1 | tee flash.log

        fi
    done




# RUN UCLUST
module purge
module load VSEARCH/2.21.1-GCC-10.3.0



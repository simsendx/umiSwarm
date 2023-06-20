
#UMI header, R1
for f in *.fastq.gz
	do
		FILE=${f#$DIRS}
		if [[ $FILE = *-v1-*R2_001.fastq.gz ]] #if file is merged, R1 read, v1
		then
			zcat $FILE | awk '(FNR-1)%4 == 0 { name=$1; chr=$2; len=$3; next }
						(FNR-2) %4 == 0 { head =substr($0,1,35)  
						 print name "." , chr, len, head ;
						 seq =substr($0,36)
						 print seq ; next}
						(FNR-3) %4 == 0 {print $0 ; next}
						(FNR-4) %4 == 0 {qseq = substr($0,36) 
						print qseq ; next}
			
						'| gzip -c >$RESDIR/UMI_$FILE
		elif [[ $FILE = *-v1-*R2_001.fastq.gz ]]
		then
			cp $FILE $RESDIR/UMI_$FILE
		
		fi
    done
	
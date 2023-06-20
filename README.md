# umiSwarm

Pipeline for eDNA metabarcoding analysis containing unique molecular identifiers.

## Basic usage

```
./umi_swarm -i <directory containing fastq files>
```

### Additional paramters

```
./umi_swarm -i <directory containing fastq files> -t nCPU -u 19 -d 1
```

- -t; --threads: No. of CPU cores to use.
- -u; --umi_length: Expected UMI length to trim.
- -d; --distance: Distance measure to use with swarm.

## Requirements

The pipeline depends on the following:

- Python3
- fastp
- cutadapt
- [vsearch](https://github.com/torognes/vsearch)
- [swarm](https://github.com/torognes/swarm) 
- BioPython

These can be installed using:

```
conda install -c bioconda fastp
conda install -c bioconda cutadapt
conda install -c conda-forge biopython
```

The following are optional to create QC reports. If they are not
installed, no reports will be generated.

- fastqc
- multiqc

## Installation


## Documentation
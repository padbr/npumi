# npumi
Derive amplicon sequence variants from nanopore amplicon sequencing done in a specific way

# Preamble
The scripts here are to document the derivation of Amplicon Sequence Variants (ASVs) from Nanopore amplicon data. The purpose of this repository is so that I can document how I have processed data in some peer reviewed publications. They are very specific to manner in which the amplicon libraries were prepared and are not intended to constitute a general pipeline to derive ASVs from Nanopore amplicon data.

# Requirements
### Hardware
These scripts can benefit from a machine with a high CPU core count if the work is split up appropriately (see below). For the most part, a modern workstation should have sufficient memory if not too many jobs are processed in parallel. One of the scripts, 03_npumi_demux.py, can run a LOT faster in a high memory machine, but it has a slower, lower-memory mode, if necessary. Memory usage and processing time will depend on the total amount of data and how much maximum data there is for any individual barcode. So very high volume datasets will have higher hardware requirements than a more modest dataset.

## Operating system
The scripts here have only been tested in ubuntu. Likely almost any current version of linux should suffice. The scripts here have not been tested in Apple nor Windows operating systems. The use of the subprocess module makes me think that these will not work outside of linux, though I am less certain of this regarding Apple machines.

## Dependencies
The scripts were written for python >=3.6. For the most part, standard python libraries are used, the versions of which should not be critical. Biopython is the only non-standard python library required. ('pip install --user biopython' if your distribution doesn't ship a biopython package). One of the scripts uses subprocess to call 'minimap2' and 'racon', both of which need to be installed appropriately.
        

# Installation
Perhaps the best way to use the scripts here is to simply download, extract, keep them in the same folder and add that to the system's PATH variable.
> export PATH="/path/to/these/scripts":$PATH

Alternatively, they may be copied to e.g. /usr/local/bin but that might be messy.


# Usage
The scripts are run in the order in which they are named, with the output of the first script being the input for the second and so on. For IO reasons, especially if working with HDDs, it is probably beneficial to gzip compress any input fastq files for these scripts.

### Basic workflow
The following is an overview of how it would look if the data was processed in a single thread.
Assume raw reads relative to the current directory are in '00_raw/reads.fq.gz'.
The following assumes that the barcodes have not been renamed in the 'primers.py' file.
> mkdir 01_framed 02_annotated 03_demux

> 01_frame_amplicons.py -i 00_raw/reads.fq.gz -l GTCTCGTGGGCTCGG -u 15 -b 24 -m 300 -M 1800 > 01_framed/reads.fq && gzip 01_framed/reads.fq

> 02_annotate_amplicons.py -i 01_framed/reads.fq.gz -o 02_annotated/reads.fq && gzip 02_annotated/reads.fq

> cd 02_annotated && 03_npumi_demux.py -e fq.gz && gzip BC*.fastq && mv BC*.fastq.gz ../03_demux

> cd ../03_demux && for file in BC*.fastq.gz; do 04_align_umid_clusters.py -i $file -o $(echo $file | sed s'/\.fastq\.gz/\.fasta/'); done

The derived ASV sequences are now in the resulting fasta files.
Consider deleting intermediate files.


### Process the data in chunks to speed it up
To take advantage of a server with a large number of CPU cores, it is probably best to split the data into smaller chunks to be processed in parallel. e.g.:
> zcat reads.fq.gz | split -l 4000000 -a 2 -d - split_

> for file in split_*; do mv $file ${file}".fq"; done

> pigz split*.fq # pigz is parallel implementation of gnuzip

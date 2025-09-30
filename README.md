<img src="docs/genasims_logo.jpg" alt="logo" title="GENAsims logo" height="400" align="right" />

# GENAsims

**GENAsims** generates realistic somatic mutations for different mutagens (e.g., APOBEC, UV) by sampling across genomic features such as replication timing bins, gene vs intergenic regions, and transcription/replication strands. 

## Installation

To install GENAsims from GitHub navigate to a directory where you want to install the program and clone the repository:
```
git clone https://github.com/KazanovLab/GENAsims
```

To compile GENAsims and install it system-wide navigate into the cloned directory and run:
```
cd GENAsims
make
sudo make install
```

## Quick start

This tool has two modes:

Index build — create a compact index of genome positions to enable fast sampling by multiple genomic features.

Mutation simulation — generate mutations according to user-defined degree of mutagenesis.

1) Build indices

Inputs (full/absolute paths):

-g — genome FASTA

-a — genome annotation (e.g., GFF3)

-s — system directory with information on replication timing regions & replication strands

-o — output directory

Example

```
genasims \
  -g   /humanGenome/hg19.fa \
  -a   /humanAnnotation/Homo_sapiens.GRCh37.87.chr.gff3 \
  -r   /RT/ESC_smooth_PC_corrected_average_10000_strand.txt \
  -o   /GENASIMS_indices/
```

2) Simulate mutations

Inputs:

-i — path to the index dir from step 1 (/GENASIMS_indices/)

-s — system directory with mutagen's distribution parameters

-n — total number of mutations to simulate

percentage per mutagen: currently APOBEC (-a) and UV (-u) are supported (percentages should sum to 100)

-o — output directory

Example

```
genasims \
  -i   /GENASIMS_indices/ \
  -s   /GENASIMS_system/ \
  -n 500000 -a 60 -u 40 \
  -o   /outd/
```

## Reporting Bugs and Feature Requests
Please use the [GitHub issue tracker](https://github.com/KazanovLab/GENAsims/issues) to report bugs or suggest features.

## Citing
*to be submitted*

## Funding
This study was supported by the Scientific and Technological Research Council of Turkey (TUBITAK) under Grant Number 123E476. The authors thank TUBITAK for their support. 


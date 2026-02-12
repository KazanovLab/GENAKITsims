<img src="docs/genasims_logo.jpg" alt="logo" title="GKsims logo" height="400" align="right" />

# GKsims

**GKsims** generates realistic somatic mutations for different mutagens (e.g., APOBEC, UV) by sampling across genomic features such as replication timing bins, gene vs intergenic regions, and transcription/replication strands. 

## Installation

To install GKsims from GitHub navigate to a directory where you want to install the program and clone the repository:
```
git clone https://github.com/KazanovLab/GKsims
```

To compile GKsims and install it system-wide navigate into the cloned directory and run:
```
cd GKsims
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

-r — system directory with information on replication timing regions & replication strands

-o — output directory

Example

```
gksims \
  -g   /humanGenome/hg19.fa \
  -a   /humanAnnotation/Homo_sapiens.GRCh37.87.chr.gff3 \
  -r   /RT/ESC_smooth_PC_corrected_average_10000_strand.txt \
  -o   /GKsims_indices/
```

2) Simulate mutations

Inputs:

-b — path to the index dir from step 1 (/GKsims_indices/)

-s — system directory with mutagen's distribution parameters

-n — total number of mutations to simulate

percentage per mutagen: currently APOBEC (-p) and UV (-u) are supported (percentages should sum to 100)

-o — output directory

Example

```
gksims \
  -b   /GKsims_indices/ \
  -s   /GKsims_system/ \
  -n 500000 -p 60 -u 40 \
  -o   /outd/
```

## Reporting Bugs and Feature Requests
Please use the [GitHub issue tracker](https://github.com/KazanovLab/GKsims/issues) to report bugs or suggest features.

## Citing
*to be submitted*


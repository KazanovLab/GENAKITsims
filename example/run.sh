../gksims -g /Users/mar/BIO/BIODATA/HUMAN/CH37/hg19.fa -a /Users/mar/BIO/BIODATA/HUMAN/CH37/Homo_sapiens.GRCh37.87.chr.gff3 -r ./RT/ESC_smooth_PC_corrected_average_10000_strand.txt -o ./GKsims_base

../gksims -b ./GKsims_base/hg19_Base.bin -s ../system -o ./results -n 10000 -p 100

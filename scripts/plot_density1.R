# install.packages(c("data.table","ggplot2"))
library(data.table)
library(ggplot2)

# --- 1. Read VCF ---
read_vcf_chr <- function(vcf_path, chr="1") {
  dt <- fread(vcf_path, skip="#CHROM", header=TRUE, sep="\t",
              select=c("#CHROM","POS"), col.names=c("CHROM","POS"))
  dt[CHROM %in% c(chr, paste0("chr", chr))]
}

# --- 2. Density calculation ---
density_bins <- function(positions, bin_size_bp = 1e6, chr_len = NULL) {
  stopifnot(length(positions) > 0)
  if (is.null(chr_len)) chr_len <- max(positions)
  
  breaks <- seq(0, chr_len, by = bin_size_bp)
  if (tail(breaks,1) < chr_len) breaks <- c(breaks, chr_len)
  
  h <- hist(positions, breaks = breaks, plot = FALSE)
  mids <- h$mids
  counts <- h$counts
  per_mb <- counts / (bin_size_bp / 1e6)   # per 1 MB
  
  data.table(pos = mids, density_per_mb = per_mb)
}

# --- 3. Multiple VCF ---
mut_density_bins_multi <- function(vcf_paths, chr="1", bin_size_bp=1e6) {
  rbindlist(lapply(vcf_paths, function(p){
    dt <- read_vcf_chr(p, chr=chr)
    dens <- density_bins(dt$POS, bin_size_bp=bin_size_bp)
    dens[, sample := basename(p)]
  }))
}

# --- 4. Usage example ---
vcfs <- c("/Users/mar/BIO/BIODATA/PCAWG/VCF/BLCA-US/cda1a403-16b6-487c-a82a-c377d1d0f89d.consensus.20160830.somatic.snv_mnv.vcf.gz",
          "/Users/mar/BIO/BIODATA/PCAWG/VCF/BLCA-US/448fe471-3f4e-4dc8-a4e0-6f147dc93abe.consensus.20160830.somatic.snv_mnv.vcf.gz",
          "/Users/mar/BIO/BIODATA/PCAWG/VCF/BLCA-US/301d6ce3-4099-4c1d-8e50-c04b7ce91450.consensus.20160830.somatic.snv_mnv.vcf.gz",
          "/Users/mar/BIO/BIODATA/PCAWG/VCF/BLCA-US/acc629cb-ad03-4cec-9b21-922e4932ef3e.consensus.20160830.somatic.snv_mnv.vcf.gz",
          "/Users/mar/BIO/BIODATA/PCAWG/VCF/BLCA-US/7d2a22eb-7344-4cba-ad7d-94c3f9ef3d7c.consensus.20160830.somatic.snv_mnv.vcf.gz",
          "/Users/mar/BIO/BIODATA/PCAWG/VCF/BLCA-US/2b142863-b963-4cc9-8f8f-c72503c93390.consensus.20160830.somatic.snv_mnv.vcf.gz",
          "/Users/mar/BIO/BIODATA/PCAWG/VCF/BLCA-US/b73523d7-f5a5-4140-8537-4df4d1ecf465.consensus.20160830.somatic.snv_mnv.vcf.gz",
          "/Users/mar/BIO/BIODATA/PCAWG/VCF/BLCA-US/8c619cbc-9e91-4716-9711-5236e55d8f46.consensus.20160830.somatic.snv_mnv.vcf.gz")  

# Any window size:
bin_size_bp <- 5e5 

dens_dt <- mut_density_bins_multi(vcfs, chr="1", bin_size_bp=bin_size_bp)

# --- 5. Visualization ---
ggplot(dens_dt, aes(pos/1e6, density_per_mb, color=sample)) +
  geom_line(linewidth=0.1) + geom_point(size=0.3) +
  labs(
    x = "Chromosome 1 position (Mb)",
    y = sprintf("Mutations per Mb (bin size = %.1f kb)", bin_size_bp/1e3),
    color = "VCF"
  ) +
  theme_minimal(base_size = 13)

ggsave("/Users/mar/BIO/PROJECTS/GKsims/plots/chr1.png",units="mm",width=1000,height=50,dpi=300,limitsize=F)

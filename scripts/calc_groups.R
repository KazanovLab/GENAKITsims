library(data.table)
library(ggplot2)

source("/Users/mar/BIO/PROJECTS/GKsims/R/parse_GKbase.R")

gene_cols <- c("Genes","InterGenes")
coding_cols_full <- c("Coding","NonCoding","UndefCoding")
coding_cols <- c("Coding","NonCoding")
RTstrand_cols <- c("Leading","Lagging","UndefStrand")
RT_cols <- c("rt1","rt2","rt3","rt4","rt5","rt6","rt7","UndefRT")

all0 <- CJ(
  gene_key     = gene_cols,
  coding_key   = coding_cols_full,
  RTstrand_key = RTstrand_cols,
  RT_key       = RT_cols
)

# Apply both constraints
all_combos <- all0[
  !(gene_key == "InterGenes" & coding_key != "UndefCoding") &
    !(gene_key == "Genes"      & coding_key == "UndefCoding") &
    !(RT_key == "UndefRT"      & RTstrand_key != "UndefStrand")
]

get_key <- function(dt, cols) {
  sub <- dt[, ..cols]
  mat <- as.matrix(sub == "+")
  idx <- max.col(mat, ties.method="first")
  idx[rowSums(mat) == 0] <- NA
  cols[idx]
}

dtall <- data.table()

for(f in list.files("/Users/mar/BIO/PROJECTS/GKsims/parameters/fact_trace")) {
  dt <- read.csv(paste0("/Users/mar/BIO/PROJECTS/GKsims/parameters/fact_trace/",f),sep='\t')
  dt <- data.table(dt)
  setnames(dt,"noCoding","NonCoding")
  setnames(dt,"interGenes","InterGenes")
  setnames(dt,"undefSTRAND","UndefStrand")
  setnames(dt,"undefRT","UndefRT")
  
  dt[, gene_key := get_key(.SD, gene_cols)]
  dt[, coding_key := get_key(.SD, coding_cols)]
  dt[, RTstrand_key := get_key(.SD, RTstrand_cols)]
  dt[, RT_key := get_key(.SD, RT_cols)]
  dt[is.na(coding_key), coding_key := "UndefCoding"]
  
  dtgr <- dt[,.("cnt"=.N),by=.(gene_key,coding_key,RTstrand_key,RT_key)]
  
  dtgroup <- merge(all_combos,dtgr,by=c("gene_key","coding_key","RTstrand_key","RT_key"),all.x=T)
  dtgroup[is.na(cnt), cnt := 0]
  dtgroup[, sample := f]
  dtgroup[, total_mut := sum(cnt)]
  dtgroup[, fraction := cnt/total_mut]
  
  dtgrp <- merge(dtgroup,gkbase_cnts,by=c("gene_key","coding_key","RTstrand_key","RT_key"),all.x=T)
  dtgrp[, density := cnt / target_cnt]
  dtgrp[, total_density := sum(density)]
  dtgrp[, fraction_density := density / total_density]
  
  if(!all.equal(sum(dtgrp$fraction_density), 1)){
   stop("Sum of fraction density is not 1")
  }
  
  dtall <- rbind(dtall, dtgrp)
}

dtall[, key := sprintf("%s|%s|%s|%s", gene_key, coding_key, RTstrand_key, RT_key)]

#dtall <- dtall[RTstrand_key != "UndefStrand" | RT_key == "UndefRT" ]

samples <- dtall[, .N,by=.(sample,total_mut)]
samples[, samples_mut := sum(total_mut)]
samples[, coef := total_mut / samples_mut]

dtall <- merge(dtall,samples[,.(sample,coef)],by="sample")

dtmean <- dtall[,.("weighted_fraction_density"=sum(coef*fraction_density)),by=key]

#dt <- dtall[sample %like% "^mutA_cda"]


dt1 <- dtall[,.(key,"fd"=fraction_density,sample)]
dt2 <- dtmean[,sample := "weighted_mean_real"]
setnames(dt2,"weighted_fraction_density","fd")

dt1[,trans_group := "pale"]
dt2[,trans_group := "vivid"]

dtplot1 <- rbind(dt1,dt2)

dtplot1[, sample_type := "real"]


# Simulated

dtall <- data.table()

for(f in list.files("/Users/mar/BIO/PROJECTS/GKsims/results/rndVCF")) {
  dt <- read.csv(paste0("/Users/mar/BIO/PROJECTS/GKsims/results/rndVCF/",f),sep='\t')
  dt <- data.table(dt)
  setnames(dt,"Genes","gene_key")
  setnames(dt,"Coding","coding_key")
  setnames(dt,"STRAND","RTstrand_key")
  setnames(dt,"RT","RT_key")
  
  dt[coding_key == "", coding_key:="UndefCoding"]
  dt[RTstrand_key == "leading", RTstrand_key := "Leading"]
  dt[RTstrand_key == "lagging", RTstrand_key := "Lagging"]
  dt[RTstrand_key == "", RTstrand_key := "UndefStrand"]
  dt[, RT_key := paste0("rt",RT_key)]
  dt[RT_key == "rt0", RT_key := "UndefRT"]
  
  dtgr <- dt[,.("cnt"=.N),by=.(gene_key,coding_key,RTstrand_key,RT_key)]
  
  dtgroup <- merge(all_combos,dtgr,by=c("gene_key","coding_key","RTstrand_key","RT_key"),all.x=T)
  dtgroup[is.na(cnt), cnt := 0]
  dtgroup[, sample := f]
  dtgroup[, total_mut := sum(cnt)]
  dtgroup[, fraction := cnt/total_mut]
  
  dtgrp <- merge(dtgroup,gkbase_cnts,by=c("gene_key","coding_key","RTstrand_key","RT_key"),all.x=T)
  dtgrp[, density := cnt / target_cnt]
  dtgrp[, total_density := sum(density)]
  dtgrp[, fraction_density := density / total_density]
  
  if(!all.equal(sum(dtgrp$fraction_density), 1)){
    stop("Sum of fraction density is not 1")
  }
  
  dtall <- rbind(dtall, dtgrp)
}

dtall[, key := sprintf("%s|%s|%s|%s", gene_key, coding_key, RTstrand_key, RT_key)]

#dtall <- dtall[RTstrand_key != "UndefStrand" | RT_key == "UndefRT" ]

samples <- dtall[, .N,by=.(sample,total_mut)]
samples[, samples_mut := sum(total_mut)]
samples[, coef := total_mut / samples_mut]

dtall <- merge(dtall,samples[,.(sample,coef)],by="sample")

dtmean <- dtall[,.("weighted_fraction_density"=sum(coef*fraction_density)),by=key]

dt1 <- dtall[,.(key,"fd"=fraction_density,sample)]
dt2 <- dtmean[,sample := "weighted_mean_simulated"]
setnames(dt2,"weighted_fraction_density","fd")

dt1[,trans_group := "pale"]
dt2[,trans_group := "vivid"]

dtplot2 <- rbind(dt1,dt2)

dtplot2[, sample_type := "simulated"]

dtplot <- rbind(dtplot1,dtplot2)

ggplot(dtplot, aes(x = key, y = fd, group = sample, color = sample, alpha = trans_group, linetype = sample_type)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.6) +
  scale_alpha_manual(values = c("vivid" = 1, "pale" = 0.1)) +
  scale_linetype_manual(values = c("real" = "solid", "simulated" = "dashed")) +
  labs(x = "gene|coding|RTstrand|RT", y = "Counts", color = "Sample") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.minor = element_blank()
  )


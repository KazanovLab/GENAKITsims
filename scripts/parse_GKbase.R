library(data.table)

# Read file - update the file path as needed
gkbase_count <- read.csv("/Users/mar/BIO/PROJECTS/GKsims/parameters/mutr_TC.txt", sep = "\t", header=F)
gkbase_count <- data.table(gkbase_count)

# Parse second column (starting from 4th character)
# 4th char: "g" = genes, "i" = intergenes
# 5th char: "c" = coding, "n" = noCoding, "-" = undefined
# 6th char: replication timing bin number
# 7th char: leading/lagging strands (<, >, -)

# Get the column name (assuming it's the second column)
col2_name <- names(gkbase_count)[2]

# Extract substring from 4th character onwards
gkbase_count[, pattern := substr(get(col2_name), 4, nchar(get(col2_name)))]

# Parse the components
gkbase_count[, gene_key := ifelse(substr(pattern, 1, 1) == "g", "Genes", 
                        ifelse(substr(pattern, 1, 1) == "i", "InterGenes", NA))]

gkbase_count[, coding_key := ifelse(substr(pattern, 2, 2) == "c", "Coding",
                          ifelse(substr(pattern, 2, 2) == "n", "NonCoding",
                                ifelse(substr(pattern, 2, 2) == "-", "UndefCoding", NA)))]

# Map RT bin (0-7) to RT_key labels
rt_labels <- c("UndefRT","rt1","rt2","rt3","rt4","rt5","rt6","rt7")
gkbase_count[, RT_bin := as.numeric(substr(pattern, 3, 3))]
gkbase_count[, RT_key := rt_labels[RT_bin + 1]]  # +1 because R indexing starts at 1

gkbase_count[, strand := substr(pattern, 4, 4)]  # <, >, or -

# Map strand to readable names
gkbase_count[, RTstrand_key := ifelse(strand == ">", "Leading",
                          ifelse(strand == "<", "Lagging",
                                ifelse(strand == "-", "UndefStrand", NA)))]

gkbase_cnts <- gkbase_count[,.("target_cnt"=sum(V3)),by=.(gene_key,coding_key,RT_key,RTstrand_key)]


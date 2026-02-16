# Kenneth B. Hoehn
# 2/10/26
# Format SNPs into concatenated sequences for BEAST

library(dplyr)
library(tidyr)

t = read.csv("083024_alltime_allpos_header_PAO1_lookup_SnpEff (1)-1.csv")
snps = read.csv("snps.csv")
times = c(128, 279, 613)

m = match(names(t),paste0("X",snps$ID))
t = rbind(t, snps[m,]$REF)
t[nrow(t),1] = "REF_REF_REF"

seqs = sapply(1:nrow(t), function(x)paste(t[x,2:ncol(t)], collapse=""))

# parse metadata from sequence ids
seqs = tibble(sequence_id = t[,1], sequence_alignment=seqs)
seqs$time = times[as.numeric(substr(seqs$sequence_id, 2, 2))]
seqs$location = sapply(strsplit(seqs$sequence_id, split="_"), function(x)x[2])

# remove reference strain
fseqs = filter(seqs, sequence_id != "REF_REF_REF")

st = paste0(">",fseqs$sequence_id, "-", fseqs$time, "-", substr(fseqs$location,2,2),
	"\n",fseqs$sequence_alignment)
writeLines(st, con="sample_full.fa")


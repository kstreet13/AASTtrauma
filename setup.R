require(readxl)

counts <- read_excel('data/Trauma neutrophil counts + metadata.xlsx')
meta <- t(counts[1:8,])
colnames(meta) <- meta[1,]; meta <- meta[-1,]
meta <- as.data.frame(meta)

counts <- as.matrix(counts[-c(1:8), ])
rownames(counts) <- counts[,1]
counts <- counts[,-1]
mode(counts) <- 'numeric'

all(colnames(counts) == rownames(meta))

require(SummarizedExperiment)
se <- SummarizedExperiment(assays = list(counts = counts), colData = meta)

rm(counts, meta)

# Filtering step for ALL analysis:  We only want to analyze genes with counts
# >=5 in >=20% of samples (15 samples).
# .2 * ncol(se)
filt <- which(rowSums(assay(se,'counts') >= 5) > 14.6)
se <- se[filt, ]
rm(filt)

se$Time[se$Patient.ID %in% c('HC37','HC38','CT1')] <- 'CTRL'
se$Time <- factor(se$Time, levels = c('CTRL','0h','8h','24h','72h'))

# se$lactate_group <- se$Phenotype
# se$lactate_group[se$lactate_group %in% c('Trauma Low Lactate','Trauma Intermediate Lactate')] <- 'Low/Int'
# se$lactate_group[se$lactate_group == 'Trauma High  Lactate'] <- 'High'

se$lactate_group <- c('Low/Int','High')[1+(as.numeric(se$`Admission lactate`) > 4)]
se$lactate_group[se$Phenotype == 'Healthy Control'] <- 'CTRL'
se$lactate_group <- factor(se$lactate_group, levels = c('CTRL','Low/Int','High'))

source('setup.R')

# It will be easier to explain to reviewers a comparison to healthy controls
# rather than comparing to 0h,n so I think this is what we should do.  If you do
# the differential expression and find no significant differences (as Jeremy
# did), then we will go for Plan B which is comparison to timepoint 0.  But for
# now let’s plan to define differentially-expressed genes as the union of genes
# significantly different from healthy controls at timepoints 0, 8, 24, and 72
# hrs.

require(DESeq2)
ds <- DESeqDataSetFromMatrix(countData = round(assay(se,'counts')),
                             colData = colData(se),
                             design= ~Time)
res <- DESeq(ds)
res1 <- results(res, alpha = .05, contrast = c('Time', '0h', 'CTRL'))
res2 <- results(res, alpha = .05, contrast = c('Time', '8h', 'CTRL'))
res3 <- results(res, alpha = .05, contrast = c('Time', '24h', 'CTRL'))
res4 <- results(res, alpha = .05, contrast = c('Time', '72h', 'CTRL'))

names(res1) <- paste0('0h_', names(res1))
names(res2) <- paste0('8h_', names(res2))
names(res3) <- paste0('24h_', names(res3))
names(res4) <- paste0('72h_', names(res4))

deseq <- cbind(res1,res2,res3,res4)

rm(res,res1,res2,res3,res4,ds)


# there are no DE genes
rownames(deseq)[which(deseq$`0h_padj` < .1)]
rownames(deseq)[which(deseq$`8h_padj` < .1)]
rownames(deseq)[which(deseq$`24h_padj` < .1)]
rownames(deseq)[which(deseq$`72h_padj` < .1)]



# then we will go for Plan B which is comparison to timepoint 0.

se$Time[se$Patient.ID %in% c('HC37','HC38','CT1')] <- 'CTRL'
se$Time <- factor(se$Time, levels = c('0h','8h','24h','72h','CTRL'))
ds <- DESeqDataSetFromMatrix(countData = round(assay(se,'counts')),
                             colData = colData(se),
                             design= ~Time)
res <- DESeq(ds)
res1 <- results(res, alpha = .05, contrast = c('Time', '8h', '0h'))
res2 <- results(res, alpha = .05, contrast = c('Time', '24h', '0h'))
res3 <- results(res, alpha = .05, contrast = c('Time', '72h', '0h'))

names(res1) <- paste0('8h_', names(res1))
names(res2) <- paste0('24h_', names(res2))
names(res3) <- paste0('72h_', names(res3))

deseq <- cbind(res1,res2,res3)

rm(res,res1,res2,res3,ds)


rownames(deseq)[which(deseq$`8h_padj` < .05)]
rownames(deseq)[which(deseq$`24h_padj` < .05)]
rownames(deseq)[which(deseq$`72h_padj` < .05)]

tab <- deseq[which(deseq$`72h_padj` < .05), c('72h_log2FoldChange','72h_pvalue','72h_padj')]
write.table(tab, file='~/Desktop/DEgenes_72h_vs_0h.csv', quote=FALSE)



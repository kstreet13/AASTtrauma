source('setup.R')
cpm <- 1000000 * t(t(assay(se,'counts')) / colSums(assay(se,'counts')))
assay(se,'cpm') <- cpm
rm(cpm)

# DE between lactate high/low groups at each timepoint
{
    require(DESeq2)
    
    ind1 <- which(se$Time == '0h')
    ds <- DESeqDataSetFromMatrix(countData = round(assay(se,'counts')[,ind1]),
                                 colData = colData(se)[ind1,],
                                 design= ~lactate_group)
    res <- DESeq(ds)
    res1 <- results(res, alpha = .05, contrast = c('lactate_group', 'High', 'Low/Int'))
    
    ind2 <- which(se$Time == '8h')
    ds <- DESeqDataSetFromMatrix(countData = round(assay(se,'counts')[,ind2]),
                                 colData = colData(se)[ind2,],
                                 design= ~lactate_group)
    res <- DESeq(ds)
    res2 <- results(res, alpha = .05, contrast = c('lactate_group', 'High', 'Low/Int'))

    ind3 <- which(se$Time == '24h')
    ds <- DESeqDataSetFromMatrix(countData = round(assay(se,'counts')[,ind3]),
                                 colData = colData(se)[ind3,],
                                 design= ~lactate_group)
    res <- DESeq(ds)
    res3 <- results(res, alpha = .05, contrast = c('lactate_group', 'High', 'Low/Int'))

    ind4 <- which(se$Time == '72h')
    ds <- DESeqDataSetFromMatrix(countData = round(assay(se,'counts')[,ind4]),
                                 colData = colData(se)[ind4,],
                                 design= ~lactate_group)
    res <- DESeq(ds)
    res4 <- results(res, alpha = .05, contrast = c('lactate_group', 'High', 'Low/Int'))
    
    names(res1) <- paste0(names(res1), '_0h')
    names(res2) <- paste0(names(res2), '_8h')
    names(res3) <- paste0(names(res3), '_24h')
    names(res4) <- paste0(names(res4), '_72h')
    
    deseq <- cbind(res1,res2,res3,res4)
    
    rm(res,res1,res2,res3,res4,ds)
} # deseq


# DE genes
g.idx <- which(abs(deseq$padj_0h) < .05 |
                   abs(deseq$padj_8h) < .05 |
                   abs(deseq$padj_24h) < .05 |
                   abs(deseq$padj_72h) < .05)

se2 <- se[,order(se$Time)]

heatdat <- log1p(assay(se2,'cpm')[g.idx, ])

# colored by timepoint
cc <- colorby(se2$Time, colors = c('gray',brewer.pal(9,'Blues')[c(2,4,6,8)]))
heatmap(heatdat, Rowv = TRUE, Colv = NA, ColSideColors = cc, labRow = NA, scale = 'row')

# colored by lactate high/low
cc <- colorby(se2$lactate_group, colors = c('gray',brewer.pal(9,'Paired')[5:6]))
heatmap(heatdat, Rowv = TRUE, Colv = NA, ColSideColors = cc, labRow = NA, scale = 'row')

# ordered by timepoint and lactate
ord <- order(se2$Time, se2$lactate_group)
cc <- colorby(se2$Time, colors = c('gray',brewer.pal(9,'Blues')[c(2,4,6,8)]))[ord]
cc <- colorby(se2$lactate_group, colors = c('gray',brewer.pal(9,'Paired')[5:6]))[ord]
heatmap(heatdat[,ord], Rowv = TRUE, Colv = NA, ColSideColors = cc, labRow = NA, scale = 'row')



# cluster genes?
genedat <- t(scale(t(heatdat)))

pca <- prcomp(genedat)

require(mclust)
mc <- Mclust(pca$x[,1:10])

table(mc$classification)

clusList <- lapply(1:5, function(i){
    rownames(deseq)[which(mc$classification==i)]
})
names(clusList) <- paste0('clus_',1:5)

max(lengths(clusList))
for(i in 1:5){
    clusList[[i]] <- c(clusList[[i]], rep('',64-length(clusList[[i]])))
}
clusList <- as.data.frame(clusList)

write.table(clusList, file='~/Desktop/fig4B_clusList.csv', quote=FALSE, sep=',')



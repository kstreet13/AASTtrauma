source('setup.R')
cpm <- 1000000 * t(t(assay(se,'counts')) / colSums(assay(se,'counts')))
assay(se,'cpm') <- cpm
rm(cpm)

# DE relative to healthy controls
{
    require(DESeq2)

    ds <- DESeqDataSetFromMatrix(countData = round(assay(se,'counts')),
                                 colData = colData(se),
                                 design= ~Time)
    res <- DESeq(ds)
    res1 <- results(res, alpha = .05, contrast = c('Time', '0h', 'CTRL'))
    res2 <- results(res, alpha = .05, contrast = c('Time', '8h', 'CTRL'))
    res3 <- results(res, alpha = .05, contrast = c('Time', '24h', 'CTRL'))
    res4 <- results(res, alpha = .05, contrast = c('Time', '72h', 'CTRL'))
    res5 <- results(res, alpha = .05, contrast = c('Time', '8h', '0h'))
    res6 <- results(res, alpha = .05, contrast = c('Time', '24h', '8h'))
    res7 <- results(res, alpha = .05, contrast = c('Time', '72h', '24h'))
    res8 <- results(res, alpha = .05, contrast = c('Time', '24h', '0h'))
    res9 <- results(res, alpha = .05, contrast = c('Time', '72h', '8h'))
    res10 <- results(res, alpha = .05, contrast = c('Time', '72h', '0h'))
    
        
    names(res1) <- paste0('0h_', names(res1))
    names(res2) <- paste0('8h_', names(res2))
    names(res3) <- paste0('24h_', names(res3))
    names(res4) <- paste0('72h_', names(res4))
    names(res5) <- paste0('8-0h_', names(res5))
    names(res6) <- paste0('24-8h_', names(res6))
    names(res7) <- paste0('72-24h_', names(res7))
    names(res8) <- paste0('24-0h_', names(res5))
    names(res9) <- paste0('72-8h_', names(res6))
    names(res10) <- paste0('72-0h_', names(res7))
    
    deseq <- cbind(res1,res2,res3,res4,res5,res6,res7,res8,res9,res10)
    
    rm(res,res1,res2,res3,res4,res5,res6,res7,res8,res9,res10,ds)
} # deseq

# DE genes
ind <- grep('padj',names(deseq))
g.lists <- lapply(ind, function(xi){ which(deseq[,xi] < .05) })
g.idx <- unique(do.call(c, g.lists))


se2 <- se[,order(se$Time)]

heatdat <- log1p(assay(se2,'cpm')[g.idx, ])

cc <- colorby(se2$Time, colors = c('gray',brewer.pal(9,'Blues')[c(2,4,6,8)]))

heatmap(heatdat, Rowv = TRUE, Colv = NA, ColSideColors = cc, labRow = NA, scale = 'row')

plot.new()
legend('left', legend = levels(se2$Time), fill = c('gray',brewer.pal(9,'Blues')[c(2,4,6,8)]), bty='n')

# cluster genes?
genedat <- t(scale(t(heatdat)))

pca <- prcomp(genedat)

require(mclust)
mc <- Mclust(pca$x[,1:10])

table(mc$classification)

clusList <- lapply(1:7, function(i){
    rownames(deseq)[which(mc$classification==i)]
})
names(clusList) <- paste0('clus_',1:7)

max(lengths(clusList))
for(i in 1:7){
    clusList[[i]] <- c(clusList[[i]], rep('',148-length(clusList[[i]])))
}
clusList <- as.data.frame(clusList)

write.table(clusList, file='~/Desktop/clusList.csv', quote=FALSE, sep=',')



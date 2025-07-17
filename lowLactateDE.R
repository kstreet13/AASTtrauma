source('setup.R')
cpm <- 1000000 * t(t(assay(se,'counts')) / colSums(assay(se,'counts')))
assay(se,'cpm') <- cpm
rm(cpm)

se <- se[,se$lactate_group %in% c('CTRL','Low/Int')]

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
    names(res8) <- paste0('24-0h_', names(res8))
    names(res9) <- paste0('72-8h_', names(res9))
    names(res10) <- paste0('72-0h_', names(res10))
    
    deseq <- cbind(res1,res2,res3,res4,res5,res6,res7,res8,res9,res10)
    
    rm(res,res1,res2,res3,res4,res5,res6,res7,res8,res9,res10,ds)
} # deseq

# DE genes
ind <- grep('padj',names(deseq))
g.lists <- lapply(ind, function(xi){ which(deseq[,xi] < .05) })
g.idx <- unique(do.call(c, g.lists))

write.csv(rownames(deseq)[g.idx], file='~/Desktop/fig2_LowLactate_genes.csv', quote=FALSE, row.names = FALSE, col.names = FALSE)

se2 <- se[,order(se$Time)]

heatdat <- log1p(assay(se2,'cpm')[g.idx, ])

cc <- colorby(se2$Time, colors = c('gray',brewer.pal(9,'Blues')[c(2,4,6,8)]))

heatmap(heatdat, Rowv = TRUE, Colv = NA, ColSideColors = cc, labRow = NA, scale = 'row')

plot.new()
legend('left', legend = levels(se2$Time), fill = c('gray',brewer.pal(9,'Blues')[c(2,4,6,8)]), bty='n')


require(pheatmap)
anncol <- data.frame(Timepoint = se2$Time)
rownames(anncol) <- colnames(heatdat)
ann_colors = list(
    Timepoint = c('gray',brewer.pal(9,'Blues')[c(2,4,6,8)])
)
names(ann_colors$Timepoint) <- levels(anncol$Timepoint)

pheatmap(heatdat, cluster_cols = FALSE, scale = 'row', show_rownames = FALSE,
         main='Low/Int Lactate', border_color = NA,
         annotation_col = anncol, annotation_colors = ann_colors)

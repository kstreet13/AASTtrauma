source('setup.R')
cpm <- 1000000 * t(t(assay(se,'counts')) / colSums(assay(se,'counts')))
assay(se,'cpm') <- cpm
rm(cpm)

# Pro-Invasive signature
sig <- read_excel('data/Pro-invasive signature.xlsx')
table(sig$Direction)
genes <- sig$Gene[sig$Gene %in% rownames(se)]

# Heatmap 1) Expression at each timepoint of the pro-invasive 51 gene signature
# in **high-lactate group**, log2FC relative to healthy controls.
{
    require(DESeq2)
    
    se2 <- se[, se$lactate_group == "High" | se$Time == "CTRL"]
    
    ds <- DESeqDataSetFromMatrix(countData = round(assay(se2,'counts')),
                                 colData = colData(se2),
                                 design= ~Time)
    res <- DESeq(ds)
    res1 <- results(res, alpha = .05, contrast = c('Time', '0h', 'CTRL'))
    res2 <- results(res, alpha = .05, contrast = c('Time', '8h', 'CTRL'))
    res3 <- results(res, alpha = .05, contrast = c('Time', '24h', 'CTRL'))
    res4 <- results(res, alpha = .05, contrast = c('Time', '72h', 'CTRL'))
    
    names(res1) <- paste0(names(res1), '_0h')
    names(res2) <- paste0(names(res2), '_8h')
    names(res3) <- paste0(names(res3), '_24h')
    names(res4) <- paste0(names(res4), '_72h')
    
    deseq <- cbind(res1,res2,res3,res4)
    
    rm(res,res1,res2,res3,res4,ds)
} # deseq

l2fcH <- cbind(deseq$log2FoldChange_0h,
              deseq$log2FoldChange_8h,
              deseq$log2FoldChange_24h,
              deseq$log2FoldChange_72h)
l2fcH <- l2fcH[which(rownames(deseq) %in% genes), ]
rownames(l2fcH) <- rownames(deseq)[which(rownames(deseq) %in% genes)]
colnames(l2fcH) <- paste0(c(0,8,24,72),'h')


cc <- colorby(factor(1:4), colors = c(brewer.pal(9,'Blues')[c(2,4,6,8)]))
heatmap(l2fcH, Rowv = TRUE, Colv = NA, ColSideColors = cc, labCol = NA,
        scale = 'row', main='High Lactate')

require(pheatmap)
mx <- max(abs(l2fcH))
pheatmap(l2fcH, cluster_cols = FALSE,
         color = colorRampPalette(c('darkred',2,'white',3,'darkgreen'))(100),
         breaks = seq(from=-mx, to=mx, length.out = 101),
         main='High Lactate')


write.table(l2fcH, file='~/Desktop/log2FC_high.csv', quote=FALSE, sep=',')

# Heatmap 2) Expression at each timepoint of the pro-invasive 51-gene signature
# in the low-lactate group.
{
    require(DESeq2)
    
    se2 <- se[, se$lactate_group == "Low/Int" | se$Time == "CTRL"]
    
    ds <- DESeqDataSetFromMatrix(countData = round(assay(se2,'counts')),
                                 colData = colData(se2),
                                 design= ~Time)
    res <- DESeq(ds)
    res1 <- results(res, alpha = .05, contrast = c('Time', '0h', 'CTRL'))
    res2 <- results(res, alpha = .05, contrast = c('Time', '8h', 'CTRL'))
    res3 <- results(res, alpha = .05, contrast = c('Time', '24h', 'CTRL'))
    res4 <- results(res, alpha = .05, contrast = c('Time', '72h', 'CTRL'))
    
    names(res1) <- paste0(names(res1), '_0h')
    names(res2) <- paste0(names(res2), '_8h')
    names(res3) <- paste0(names(res3), '_24h')
    names(res4) <- paste0(names(res4), '_72h')
    
    deseq <- cbind(res1,res2,res3,res4)
    
    rm(res,res1,res2,res3,res4,ds)
} # deseq

l2fcL <- cbind(deseq$log2FoldChange_0h,
              deseq$log2FoldChange_8h,
              deseq$log2FoldChange_24h,
              deseq$log2FoldChange_72h)
l2fcL <- l2fcL[which(rownames(deseq) %in% genes), ]
rownames(l2fcL) <- rownames(deseq)[which(rownames(deseq) %in% genes)]
colnames(l2fcL) <- paste0(c(0,8,24,72),'h')



cc <- colorby(factor(1:4), colors = c(brewer.pal(9,'Blues')[c(2,4,6,8)]))
heatmap(l2fcL, Rowv = TRUE, Colv = NA, ColSideColors = cc, labCol = NA,
        scale = 'none', main='Low/Int Lactate')

require(pheatmap)
mx <- max(abs(l2fcL))
pheatmap(l2fcL, cluster_cols = FALSE,
         color = colorRampPalette(c('darkred',2,'white',3,'darkgreen'))(100),
         breaks = seq(from=-mx, to=mx, length.out = 101),
         main='Low/Int Lactate')

write.table(l2fcL, file='~/Desktop/log2FC_low.csv', quote=FALSE, sep=',')


# order l2fcH and l2fcL
ord <- hclust(dist(cbind(l2fcL,l2fcH)))$order

require(pheatmap)
mx <- max(abs(cbind(l2fcL,l2fcH)))
pheatmap(l2fcL, cluster_cols = FALSE, cluster_rows = FALSE,
         color = colorRampPalette(c('darkred',2,'white',3,'darkgreen'))(100),
         breaks = seq(from=-mx, to=mx, length.out = 101),
         main='Low/Int Lactate')

pheatmap(l2fcH, cluster_cols = FALSE, cluster_rows = FALSE,
         color = colorRampPalette(c('darkred',2,'white',3,'darkgreen'))(100),
         breaks = seq(from=-mx, to=mx, length.out = 101),
         main='High Lactate')




# Score:  For each timepoint, let’s take the average of all the counts (KS: CPM)
# for all the genes in the high lactate group, and then all the counts in all
# the genes for the low lactate group.

genes2 <- genes
genes2 <- genes[genes != 'CD163']

gp <- paste(se$lactate_group, se$Time, sep='_')
sg <- colSums(assay(se,'counts')[genes2,])
sg <- log1p(1000000 * sg / colSums(assay(se,'counts')))

avg <- sapply(unique(gp), function(g){
    mean(sg[gp == g])
})

write.csv(avg, file='~/Desktop/supergene_AvglogCPM.csv', quote=FALSE)




# PCA using genes2 (signature minus CD163), averaged by timepoint/lactate
input <- assay(se,'cpm')[genes2, ]
input <- sapply(unique(gp), function(g){
    rowMeans(input[, gp == g])
})

pca <- BiocSingular::runPCA(t(input), rank=9)

#plot(pca$sdev^2)

plot(pca$x[,1:2], asp=1, col='white')
text(pca$x[,1:2], labels = rownames(pca$x))

# CT (text)
# 0h, 8h, 24h, 72h in dark/light red for lactate group
# % Var Exp on axes

varexp <- pca$sdev^2 / sum(pca$sdev^2)
tp <- rownames(pca$x)
tp <- gsub('^.*_','',tp)
tp[tp == 'CTRL'] <- 'CT'
hl <- rownames(pca$x)
hl <- gsub('_.*$','',hl)
hl[hl == 'CTRL'] <- 'CT'

cc <- hl
cc[hl == 'CT'] <- 'grey'
cc[hl == 'Low/Int'] <- brewer.pal(6,'Paired')[5]
cc[hl == 'High'] <- brewer.pal(6,'Paired')[6]

plot(pca$x[,1:2], asp=1, col='white',
     xlab = paste0('PC-1 (',format(100*varexp[1], digits=3),'%)'),
     ylab = paste0('PC-2 (',format(100*varexp[2], digits=3),'%)'))
text(pca$x[,1:2], labels = tp, col = cc)






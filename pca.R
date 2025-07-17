source('setup.R')
cpm <- 1000000 * t(t(assay(se,'counts')) / colSums(assay(se,'counts')))
assay(se,'cpm') <- cpm
rm(cpm)

filt <- which(rowMeans(assay(se,'cpm')) > 5)

pca <- BiocSingular::runPCA(t(log1p(assay(se,'cpm')[filt, ])), rank=73)

plot(pca$sdev^2)

cc <- colorby(se$Time, colors = c('gray',brewer.pal(9,'Blues')[c(2,4,6,8)]))
plot(pca$x, asp=1, col = cc,
     pch = c(1,16)[1+(se$lactate_group=='High')])

cc <- colorby(se$Time, colors = c('gray',brewer.pal(9,'Set1')[c(1,5,6,3)]))
plot(pca$x, asp=1, col = cc, cex=2,
     pch = c(1,16)[1+(se$lactate_group=='High')])

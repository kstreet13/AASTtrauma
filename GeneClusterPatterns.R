
# The DEG in high compared to low lactate trauma patients were clustered into 5
# clusters (.csv files attached for reference).  Are you able to plot these
# clusters, similar to the "line drawings" you did for the
# sepsis/vasoplegia paper, but distinguishing between high and low lactate gene
# expression levels over time perhaps by coloring expression red vs black?  I
# think this will help the audience understand what we mean by clustering these
# genes according to their dynamic behavior

clus1 <- read.csv('~/Downloads/high vs low lactate cluster 1.csv', header = FALSE)[,1]
clus2 <- read.csv('~/Downloads/high vs low lactate cluster 2.csv', header = FALSE)[,1]
clus3 <- read.csv('~/Downloads/high vs low lactate cluster 3.csv', header = FALSE)[,1]
clus4 <- read.csv('~/Downloads/high vs low lactate cluster 4.csv', header = FALSE)[,1]
clus5 <- read.csv('~/Downloads/high vs low lactate cluster 5csv.csv', header = FALSE)[,1]


source('setup.R')
se <- se[,se$lactate_group %in% c('Low/Int','High')]

cpm <- 1000000 * t(t(assay(se,'counts')) / colSums(assay(se,'counts')))




# messy line plot 
genelist <- unique(c(clus1,clus2,clus3,clus4,clus5))
linedat <- as.matrix(cpm[genelist, ])

# linedat <- linedat[,  c(paste0('SVS',1:4,'A'),paste0('SVS',1:4,'B'),paste0('SVS',1:4,'C'))]
zscores <- (linedat - rowMeans(linedat)) / rowSds(linedat)

meansLO <- cbind(
    rowMeans(zscores[,which(se$lactate_group=='Low/Int' & se$Time=='0h')]),
    rowMeans(zscores[,which(se$lactate_group=='Low/Int' & se$Time=='8h')]),
    rowMeans(zscores[,which(se$lactate_group=='Low/Int' & se$Time=='24h')]),
    rowMeans(zscores[,which(se$lactate_group=='Low/Int' & se$Time=='72h')])
)
meansHI <- cbind(
    rowMeans(zscores[,which(se$lactate_group=='High' & se$Time=='0h')]),
    rowMeans(zscores[,which(se$lactate_group=='High' & se$Time=='8h')]),
    rowMeans(zscores[,which(se$lactate_group=='High' & se$Time=='24h')]),
    rowMeans(zscores[,which(se$lactate_group=='High' & se$Time=='72h')])
)
lims <- range(cbind(meansLO,meansHI))
ind <- sample(nrow(linedat))

# separated, black = high, red = low/int
layout(matrix(1:6,nrow=3,byrow = TRUE))
# clus1
plot(c(1,4), lims, col='white', axes=FALSE,
     main = 'Cluster 1', xlab = 'Timepoint', ylab = 'Avg. Expression (z-score)')
box(); axis(2); axis(1, at=1:4, labels=c(0,8,24,72)); abline(h=0,lty=2)
for(g in clus1){
    lines(1:4, meansLO[g,], col=alpha(2, alpha=.3))
    lines(1:4, meansHI[g,], col=alpha(1, alpha=.3))
}
# clus2
plot(c(1,4), lims, col='white', axes=FALSE,
     main = 'Cluster 2', xlab = 'Timepoint', ylab = 'Avg. Expression (z-score)')
box(); axis(2); axis(1, at=1:4, labels=c(0,8,24,72)); abline(h=0,lty=2)
for(g in clus2){
    lines(1:4, meansLO[g,], col=alpha(2, alpha=.3))
    lines(1:4, meansHI[g,], col=alpha(1, alpha=.3))
}
# clus3
plot(c(1,4), lims, col='white', axes=FALSE,
     main = 'Cluster 3', xlab = 'Timepoint', ylab = 'Avg. Expression (z-score)')
box(); axis(2); axis(1, at=1:4, labels=c(0,8,24,72)); abline(h=0,lty=2)
for(g in clus3){
    lines(1:4, meansLO[g,], col=alpha(2, alpha=.3))
    lines(1:4, meansHI[g,], col=alpha(1, alpha=.3))
}
# clus4
plot(c(1,4), lims, col='white', axes=FALSE,
     main = 'Cluster 4', xlab = 'Timepoint', ylab = 'Avg. Expression (z-score)')
box(); axis(2); axis(1, at=1:4, labels=c(0,8,24,72)); abline(h=0,lty=2)
for(g in clus4){
    lines(1:4, meansLO[g,], col=alpha(2, alpha=.3))
    lines(1:4, meansHI[g,], col=alpha(1, alpha=.3))
}
# clus5
plot(c(1,4), lims, col='white', axes=FALSE,
     main = 'Cluster 5', xlab = 'Timepoint', ylab = 'Avg. Expression (z-score)')
box(); axis(2); axis(1, at=1:4, labels=c(0,8,24,72)); abline(h=0,lty=2)
for(g in clus5){
    lines(1:4, meansLO[g,], col=alpha(2, alpha=.3))
    lines(1:4, meansHI[g,], col=alpha(1, alpha=.3))
}

plot.new()
legend('left', legend=c('High Lactate','Low/Int Lactate'),
       lty=1, col = alpha(1:2, alpha=.35), bty = 'n')

layout(1)

# note: doesn't include controls

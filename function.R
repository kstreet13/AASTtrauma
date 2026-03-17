require(topGO)

de <- read.csv('data/genelists/metagen_two_datasets/meta_all_early.csv')


require(org.Hs.eg.db)
anno <- select(org.Hs.eg.db, keys=as.character(as.matrix(myGenes)), 
               columns="ENTREZID", keytype="SYMBOL")
names(anno) <- c('Gene','Entrez')

de <- merge(de, anno, all.x = TRUE)

sig <- de$Entrez[which(de$meta_fdr < .05)]

go <- goana(sig, universe = de$Entrez)
go <- go[order(go$P.DE), ]

topTerms <- function(res, topn = 10, ...){
    stopifnot(any(c('Pathway','Term') %in% names(res)))
    type <- c('Pathway','Term')[which.max(c('Pathway','Term') %in% names(res))]
    df <- res[topn:1, ]
    barplot(-log10(df$P.DE), horiz = TRUE, 
            xlab = '-log10 P-val', ylab = type, ...)
    t.x <- -log10(df$P.DE[nrow(df)])
    t.y <- seq(0.7, 0.7+1.2*(topn - 1), by = 1.2)
    text(t.x,t.y, df[[type]], pos = 2, cex=.8)
}

topTerms(go)

# all significant ones






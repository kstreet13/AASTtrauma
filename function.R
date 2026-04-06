require(topGO)

de <- read.csv('data/GSEA lists 3-17-26 2/core genes all timepoints/DEG core_genes_all_timepoints.csv')
univ <- read.csv('data/GSEA lists 3-17-26 2/core genes all timepoints/UNIVERSE core genes all timepoints.csv')

addEntrezIDs <- function(df){
    require(org.Hs.eg.db)
    anno <- select(org.Hs.eg.db, keys=df$Gene, 
                   columns="ENTREZID", keytype="SYMBOL")
    names(anno) <- c('Gene','Entrez')
    df <- merge(df, anno, all.x = TRUE)
    df
}

de <- addEntrezIDs(de)
univ <- addEntrezIDs(univ)

#sig <- de$Entrez[which(de$meta_fdr < .05)]

require(limma)
go <- goana(de$Entrez, universe = univ$Entrez)
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






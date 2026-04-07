require(topGO)

de <- read.csv('data/GSEA lists 3-17-26/core genes all timepoints/DEG core_genes_all_timepoints.csv')
univ <- read.csv('data/GSEA lists 3-17-26/core genes all timepoints/UNIVERSE core genes all timepoints.csv')

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
go$p.adj <- p.adjust(go$P.DE, method = 'fdr')

topTerms <- function(res){
    stopifnot(any(c('Pathway','Term') %in% names(res)))
    type <- c('Pathway','Term')[which.max(c('Pathway','Term') %in% names(res))]
    ind <- which(res$p.adj < .05)
    if(length(ind)==0){
        ind <- 1:10
    }
    return(res[ind, ])
}

topTerms(go)

topTermsPlot <- function(res, topn = 10, ...){
    stopifnot(any(c('Pathway','Term') %in% names(res)))
    type <- c('Pathway','Term')[which.max(c('Pathway','Term') %in% names(res))]
    df <- res[topn:1, ]
    barplot(-log10(df$P.DE), horiz = TRUE, 
            xlab = '-log10 P-val', ylab = type, ...)
    t.x <- -log10(df$P.DE[nrow(df)])
    t.y <- seq(0.7, 0.7+1.2*(topn - 1), by = 1.2)
    text(t.x,t.y, df[[type]], pos = 2, cex=.8)
}

topTermsPlot(go)

# all significant ones



# generic function for getting top GO terms from a given list of de genes and universe of possible genes
getTopTerms <- function(de, univ){
    de <- read.csv(de)
    univ <- read.csv(univ)
    
    de <- addEntrezIDs(de)
    univ <- addEntrezIDs(univ)
    
    #sig <- de$Entrez[which(de$meta_fdr < .05)]
    
    go <- goana(de$Entrez, universe = univ$Entrez)
    go <- go[order(go$P.DE), ]
    go$p.adj <- p.adjust(go$P.DE, method = 'fdr')
    
    return(topTerms(go))
}


# core genes all timepoints
res <- getTopTerms(de = 'data/GSEA lists 3-17-26/core genes all timepoints/DEG core_genes_all_timepoints.csv',
                    univ = 'data/GSEA lists 3-17-26/core genes all timepoints/UNIVERSE core genes all timepoints.csv')
write.csv(res, file = 'data/GSEA lists 3-17-26/core genes all timepoints/GOTERMS core_genes_all_timepoints.csv')

# overlap EMTAB5882 U GSE36809 U GLUEPMN
res <- getTopTerms(de = 'data/GSEA lists 3-17-26/overlap EMTAB5882 U GSE36809 U GLUEPMN/DEG EMTAB5882 U GSE36809 U GLUEPMN/DEG EMTAB5882 U GSE36809 U GLUEPMN_early.csv',
                   univ = 'data/GSEA lists 3-17-26/overlap EMTAB5882 U GSE36809 U GLUEPMN/UNIVERSE EMTAB5882 U GSE36809 U GLUEPMN.csv')
write.csv(res, file = 'data/GSEA lists 3-17-26/overlap EMTAB5882 U GSE36809 U GLUEPMN/GOTERMS EMTAB5882 U GSE36809 U GLUEPMN/GOTERMS EMTAB5882 U GSE36809 U GLUEPMN_early.csv')

res <- getTopTerms(de = 'data/GSEA lists 3-17-26/overlap EMTAB5882 U GSE36809 U GLUEPMN/DEG EMTAB5882 U GSE36809 U GLUEPMN/DEG EMTAB5882 U GSE36809 U GLUEPMN_middle.csv',
                   univ = 'data/GSEA lists 3-17-26/overlap EMTAB5882 U GSE36809 U GLUEPMN/UNIVERSE EMTAB5882 U GSE36809 U GLUEPMN.csv')
write.csv(res, file = 'data/GSEA lists 3-17-26/overlap EMTAB5882 U GSE36809 U GLUEPMN/GOTERMS EMTAB5882 U GSE36809 U GLUEPMN/GOTERMS EMTAB5882 U GSE36809 U GLUEPMN_middle.csv')

res <- getTopTerms(de = 'data/GSEA lists 3-17-26/overlap EMTAB5882 U GSE36809 U GLUEPMN/DEG EMTAB5882 U GSE36809 U GLUEPMN/DEG EMTAB5882 U GSE36809 U GLUEPMN_late.csv',
                   univ = 'data/GSEA lists 3-17-26/overlap EMTAB5882 U GSE36809 U GLUEPMN/UNIVERSE EMTAB5882 U GSE36809 U GLUEPMN.csv')
write.csv(res, file = 'data/GSEA lists 3-17-26/overlap EMTAB5882 U GSE36809 U GLUEPMN/GOTERMS EMTAB5882 U GSE36809 U GLUEPMN/GOTERMS EMTAB5882 U GSE36809 U GLUEPMN_late.csv')

# Metagen EMTAB5882 GSE36809 GLUEPMN AAST
res <- getTopTerms(de = 'data/GSEA lists 3-17-26/Metagen EMTAB5882 GSE36809 GLUEPMN AAST/DEG Metagen EMTAB5882 GSE36809 GLUEPMN AAST/DEG Metagen EMTAB5882 GSE36809 GLUEPMN AAST_early.csv',
                   univ = 'data/GSEA lists 3-17-26/Metagen EMTAB5882 GSE36809 GLUEPMN AAST/UNIVERSE Metagen EMTAB5882 GSE36809 GLUEPMN AAST.csv')
write.csv(res, file = 'data/GSEA lists 3-17-26/Metagen EMTAB5882 GSE36809 GLUEPMN AAST/GOTERMS Metagen EMTAB5882 GSE36809 GLUEPMN AAST/GOTERMS Metagen EMTAB5882 GSE36809 GLUEPMN AAST_early.csv')

res <- getTopTerms(de = 'data/GSEA lists 3-17-26/Metagen EMTAB5882 GSE36809 GLUEPMN AAST/DEG Metagen EMTAB5882 GSE36809 GLUEPMN AAST/DEG Metagen EMTAB5882 GSE36809 GLUEPMN AAST_middle.csv',
                   univ = 'data/GSEA lists 3-17-26/Metagen EMTAB5882 GSE36809 GLUEPMN AAST/UNIVERSE Metagen EMTAB5882 GSE36809 GLUEPMN AAST.csv')
write.csv(res, file = 'data/GSEA lists 3-17-26/Metagen EMTAB5882 GSE36809 GLUEPMN AAST/GOTERMS Metagen EMTAB5882 GSE36809 GLUEPMN AAST/GOTERMS Metagen EMTAB5882 GSE36809 GLUEPMN AAST_middle.csv')

res <- getTopTerms(de = 'data/GSEA lists 3-17-26/Metagen EMTAB5882 GSE36809 GLUEPMN AAST/DEG Metagen EMTAB5882 GSE36809 GLUEPMN AAST/DEG Metagen EMTAB5882 GSE36809 GLUEPMN AAST_late.csv',
                   univ = 'data/GSEA lists 3-17-26/Metagen EMTAB5882 GSE36809 GLUEPMN AAST/UNIVERSE Metagen EMTAB5882 GSE36809 GLUEPMN AAST.csv')
write.csv(res, file = 'data/GSEA lists 3-17-26/Metagen EMTAB5882 GSE36809 GLUEPMN AAST/GOTERMS Metagen EMTAB5882 GSE36809 GLUEPMN AAST/GOTERMS Metagen EMTAB5882 GSE36809 GLUEPMN AAST_late.csv')

# metagen EMTAB5882, GSE36809
res <- getTopTerms(de = 'data/GSEA lists 3-17-26/metagen EMTAB5882, GSE36809/DEG metagen EMTAB5882, GSE36809/DEG metagen EMTAB5882, GSE36809_early.csv',
                   univ = 'data/GSEA lists 3-17-26/metagen EMTAB5882, GSE36809/UNIVERSE metagen EMTAB5882, GSE36809.csv')
write.csv(res, file = 'data/GSEA lists 3-17-26/metagen EMTAB5882, GSE36809/GOTERMS metagen EMTAB5882, GSE36809/GOTERMS metagen EMTAB5882, GSE36809_early.csv')

res <- getTopTerms(de = 'data/GSEA lists 3-17-26/metagen EMTAB5882, GSE36809/DEG metagen EMTAB5882, GSE36809/DEG metagen EMTAB5882, GSE36809_middle.csv',
                   univ = 'data/GSEA lists 3-17-26/metagen EMTAB5882, GSE36809/UNIVERSE metagen EMTAB5882, GSE36809.csv')
write.csv(res, file = 'data/GSEA lists 3-17-26/metagen EMTAB5882, GSE36809/GOTERMS metagen EMTAB5882, GSE36809/GOTERMS metagen EMTAB5882, GSE36809_middle.csv')

res <- getTopTerms(de = 'data/GSEA lists 3-17-26/metagen EMTAB5882, GSE36809/DEG metagen EMTAB5882, GSE36809/DEG metagen EMTAB5882, GSE36809_late.csv',
                   univ = 'data/GSEA lists 3-17-26/metagen EMTAB5882, GSE36809/UNIVERSE metagen EMTAB5882, GSE36809.csv')
write.csv(res, file = 'data/GSEA lists 3-17-26/metagen EMTAB5882, GSE36809/GOTERMS metagen EMTAB5882, GSE36809/GOTERMS metagen EMTAB5882, GSE36809_late.csv')

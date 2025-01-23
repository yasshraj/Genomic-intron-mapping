#' Calculate Genome-Wide Splicing Scores
#'
#' This function computes splicing scores for all genes in a ballgown object, 
#' allowing filtering of low-expressed genes and introns based on user-defined criteria.
#'
#' @param bg A ballgown object.
#' @param gene.select A logical expression to filter genes based on gene-level coverage (default: keeps genes with coverage >= 1 in at least 95% of samples).
#' @param intron.select A logical expression to filter introns based on junction counts (default: keeps introns with counts >= 5 in at least 5% of samples).
#' @return A list with two elements:
#' - `score`: A matrix of splicing scores with introns as rows and samples as columns.
#' - `intron`: A GRanges object with intron structure.
#' @details The splicing score is calculated as junction count divided by gene-level per-base read coverage.

spliceGenome <- function(bg, gene.select = "rowQuantiles(x,probs = 0.05)>=1", 
    intron.select = "rowQuantiles(x,probs = 0.95)>=5") {
    options(stringsAsFactors = FALSE)
    cat("---Calculate gene-level read coverage:\n")
    t_cov <- texpr(bg, "cov")
    t_cov <- split.data.frame(t_cov, geneIDs(bg))
    g_cov <- lapply(t_cov, colSums)
    g_cov <- do.call("rbind", g_cov)
    if (!is.na(gene.select)) {
        g_cov <- subset(g_cov, subset = eval(parse(text = gene.select), 
            list(x = g_cov)))
    }
    cat("\t", nrow(g_cov), "genes selected.\n")
    
    cat("---Extract intron-level read ucount:\n")
    intron <- ballgown::structure(bg)$intron
    a <- vapply(intron$transcripts, function(x) {
        ((eval(parse(text = x))[1]))
    }, 1)
    a <- as.character(a)
    intron$gene_id <- geneIDs(bg)[a]
    intron <- intron[intron$gene_id %in% rownames(g_cov)]    
    ie <- iexpr(bg, "ucount")
    ie <- subset(ie, rownames(ie) %in% intron$id)
    if (!is.na(intron.select)) {
        ie <- subset(ie, eval(parse(text = intron.select), list(x = ie)))
    }
    intron <- intron[match(rownames(ie), intron$id)]
    g_cov <- g_cov[match(intron$gene_id, rownames(g_cov)), , drop = FALSE]
    cat("\t", nrow(g_cov), "introns in", length(unique(intron$gene_id)), 
        "genes selected.\n")
    score <- ifelse(g_cov == 0, NA, ie / g_cov)
    rownames(score) <- intron$id
    colnames(score) <- sampleNames(bg)
    list(score = score, intron = intron)
}

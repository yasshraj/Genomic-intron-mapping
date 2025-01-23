# Intron Analysis Utilities
#
# This file contains a collection of utility functions for analyzing intron data from a Ballgown object. 
# The functions are designed to create and filter intron data based on start positions, associate gene IDs with introns, and more.

# Create a DataFrame of Introns with the Same Start Position
#
# Groups introns by their start positions and assigns a cluster ID. Filters out unnecessary columns.
#
# @param bg Ballgown object
# @return Data table of introns with the same start position
# Example usage: same_start <- create_same_start_df(bg)
create_same_start_df <- function(bg) {
    expr_data  <- expr(bg)
    intron_data <- expr_data$intron

    data <- intron_data
    data <- as.data.table(data)

    setkey(data, start)
    data[, cluster_id := rleid(start)]
    
    exclude_cols <- grep("^rcount|^mrcount", names(data), value = TRUE)
    selected_data <- data[, !exclude_cols, with = FALSE]
    
    same_start <- selected_data[, if (.N > 1) .SD, by = cluster_id]
    return(same_start)
}

# Delete Data Based on Cluster ID or Intron ID
#
# Removes rows from the dataset based on a specified cluster ID or intron ID.
#
# @param data Data table containing intron data
# @param cluster_id Cluster ID to filter by (optional)
# @param i_id Intron ID to filter by (optional)
# @return Filtered data table
delete_data <- function(data, cluster_id = NULL, i_id = NULL) {
  if (!is.null(cluster_id) & !is.null(i_id)) {
    stop("Please provide either cluster_id or i_id, not both.")
  }
  
  if (!is.null(cluster_id)) {
    data <- data[data$cluster_id != cluster_id, ]
  }
  
  if (!is.null(i_id)) {
    data <- data[data$i_id != i_id, ]
  }
  
  return(data)
}

# Associate Gene IDs with Introns
#
# Matches intron IDs with their corresponding gene IDs and filters based on user-defined criteria.
#
# @param bg Ballgown object
# @param gene_select Logical expression to filter genes (optional)
# @param intron_select Logical expression to filter introns (optional)
# @return Data frame of gene IDs, intron IDs, and read counts
associate_gene_with_intron <- function(bg, gene_select = NA, intron_select = NA) {
  
  t_cov <- texpr(bg, "cov")
  t_cov <- split.data.frame(t_cov, geneIDs(bg))
  g_cov <- lapply(t_cov, colSums)
  g_cov <- do.call("rbind", g_cov)
  
  if (!is.na(gene_select)) {
    g_cov <- subset(g_cov, subset = eval(parse(text = gene_select), list(x = g_cov)))
  }
  
  intron <- ballgown::structure(bg)$intron
  
  intron$gene_id <- geneIDs(bg)[vapply(intron$transcripts, function(x) {
    as.numeric(eval(parse(text = x))[1])
  }, 1)]
  
  filtered_introns <- intron[intron$gene_id %in% rownames(g_cov), ]
  
  ie <- iexpr(bg, "ucount") 
  ie <- ie[rownames(ie) %in% filtered_introns$id, ]
  
  if (!is.na(intron_select)) {
    ie <- subset(ie, eval(parse(text = intron_select), list(x = ie)))
  }
  
  result <- data.frame(
    gene_id = filtered_introns$gene_id,
    i_id = filtered_introns$id,
    read_counts = rowSums(ie)
  )
  return(result)
}

# Filter Data by Intron ID and Add Gene ID
#
# Merges a filter list with gene IDs from a Ballgown object and orders the result by cluster ID if available.
#
# @param bg Ballgown object
# @param filter_list Data table with an 'i_id' column
# @return Data table with the added 'gene_id' column
# Example usage: filtered_same_start <- filter_by_iid(bg, same_start1)
filter_by_iid <- function(bg, filter_list) {
  filter_list <- as.data.table(filter_list)
  result <- associate_gene_with_intron(bg) 
  result <- as.data.table(result)

  if (!"i_id" %in% colnames(filter_list)) {
    stop("'i_id' column not found in the data.table")
  }
  
  if (!all(c("i_id", "gene_id") %in% colnames(result))) {
    stop("'result' must contain 'i_id' and 'gene_id' columns")
  }
  
  i_id_index <- which(colnames(filter_list) == "i_id")
  
  filter_list <- cbind(filter_list[, 1:i_id_index, with = FALSE], gene_id = NA_character_, filter_list[, (i_id_index + 1):ncol(filter_list), with = FALSE])
  
  setkey(filter_list, i_id)
  setkey(result, i_id)
  
  filter_list[result, gene_id := result$gene_id, on = .(i_id)]

  if ("cluster_id" %in% colnames(filter_list)) {
    setorder(filter_list, cluster_id)
  }
  
  return(filter_list)
}
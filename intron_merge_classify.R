# This file contains functions that performs several operations on intron data, including merging data frames, replacing ucount values, 
# identifying non-matching ucount rows, and comparing ucount values for different clusters.

# Merging Data and Replacing Ucount Values with Score
#
# This function merges two data frames, 'df1' and 'df2', based on the 'i_id' column from 'df2' 
# and the 'intron.id' column from 'df1'. It also replaces values in the 'ucount' columns of 'df2' 
# with the corresponding values from the 'score' columns of 'df1'.
#
# @param df1 Data frame containing 'intron.id' and 'score' columns
# @param df2 Data frame containing 'i_id' and 'ucount' columns
# @return Merged data frame with replaced 'ucount' valuesmerge_and_replace_names <- function(df1, df2) {
merge_and_replace_names <- function(df1, df2) {
    if (!'intron.id' %in% colnames(df1) || !'i_id' %in% colnames(df2)) {
        stop("df1 must contain 'intron.id' and df2 must contain 'i_id'")
    }

    merged_df <- merge(df2, df1, by.x = "i_id", by.y = "intron.id", all.x = TRUE)

    ucount_columns <- grep("^ucount\\.", colnames(merged_df), value = TRUE)
    score_columns <- grep("^score\\.", colnames(merged_df), value = TRUE)

    for (ucount_col in ucount_columns) {
        name_without_prefix <- sub("^ucount\\.", "", ucount_col)
        corresponding_score_col <- paste0("score.", name_without_prefix)

        if (corresponding_score_col %in% score_columns) {
            merged_df[[ucount_col]] <- ifelse(is.na(merged_df[[corresponding_score_col]]), merged_df[[ucount_col]], merged_df[[corresponding_score_col]])
        }
    }

    df1_cols <- colnames(df1)
    result_df <- subset(merged_df, select = setdiff(colnames(merged_df), df1_cols))

    if ("cluster_id" %in% colnames(result_df)) {
        setorder(result_df, cluster_id)
    }
    
    return(result_df)
}

# Finding Non-Matching Ucount Rows Between Two Data Frames
#
# This function compares 'ucount' columns between two data frames and returns rows where there are mismatches.
# It checks for differences between the 'ucount' columns in both data frames and returns the rows where mismatches occur.
#
# @param df1 First data frame containing 'ucount' columns
# @param df2 Second data frame containing 'ucount' columns
# @return Data frame with rows having non-matching 'ucount' values
non_matching_ucount_rows <- function(df1, df2) {
    df1 <- as.data.table(df1)
    df2 <- as.data.table(df2)
    
    ucount_cols1 <- grep("^ucount", names(df1), value = TRUE)
    ucount_cols2 <- grep("^ucount", names(df2), value = TRUE)
    
    common_ucount_cols <- intersect(ucount_cols1, ucount_cols2)
    
    non_matched_rows <- df1[, common_ucount_cols, with = FALSE] != df2[, common_ucount_cols, with = FALSE]
    non_matched_df <- df1[apply(non_matched_rows, 1, any), ]
    
    return(non_matched_df)
}

# 3. Function to remove introns that do not have a pair based on cluster_id
non_matched1 <- non_matched[, if (.N > 1) .SD, by = cluster_id]

# Comparing Ucount Values for Specific Cluster ID
#
# This function compares 'ucount' values between two rows within a specified cluster ID and classifies them 
# based on their relative values. It classifies the 'ucount' values into three categories:
# 'fivePrimeH' (first value is 10 times greater than the second), 'threePrimeH' (second value is 10 times greater), 
# and 'equalH' (neither value is significantly greater than the other).
#
# @param df Data frame containing 'ucount' columns and 'cluster_id'
# @param cluster_id_value Cluster ID for comparison
# @return List of classified results for 'ucount' comparisons
# example usage: compare_ucounts(non_matched1, 1)
compare_ucounts <- function(df, cluster_id_value) {
    filtered_df <- df[df$cluster_id == cluster_id_value, ]
    if (nrow(filtered_df) != 2) {
        stop("The filtered data does not contain exactly two rows.")
    }
    row1 <- filtered_df[1, ]
    row2 <- filtered_df[2, ]
    result <- list(
        fivePrimeH = list(),
        threePrimeH = list(),
        equalH = list()
    )
    ucount_columns <- grep("^ucount", colnames(df), value = TRUE)
    for (col in ucount_columns) {
        value1 <- row1[[col]]
        value2 <- row2[[col]]
        
        if (value1 >= 10 * value2) {
            result$fivePrimeH[[col]] <- c(value1, value2)
        } else if (value2 >= 10 * value1) {
            result$threePrimeH[[col]] <- c(value1, value2)
        } else {
            result$equalH[[col]] <- c(value1, value2)
        }
    }
    return(result)
}

# Comparing Ucount Values for All Cluster IDs
#
# This function compares 'ucount' values for all rows within each cluster and classifies them similarly to the 
# 'compare_ucounts' function, but this time for all unique cluster IDs in the data frame.
#
# @param df Data frame containing 'ucount' columns and 'cluster_id'
# @return List of classified results for all 'ucount' comparisons across cluster IDs
# example usage: compare_all_ucounts(non_matched1)
compare_all_ucounts <- function(df) {
    unique_cluster_ids <- unique(df$cluster_id)
    
    data_info <- list()
    for (cluster_id_value in unique_cluster_ids) {
        filtered_df <- df[df$cluster_id == cluster_id_value, ]
        if (nrow(filtered_df) != 2) {
            warning(paste("Cluster ID", cluster_id_value, "does not have exactly two rows. Skipping."))
            next
        }
        
        row1 <- filtered_df[1, ]
        row2 <- filtered_df[2, ]
        
        cluster_result <- list(
            fivePrimeH = list(),
            threePrimeH = list(),
            equalH = list()
        )
        ucount_columns <- grep("^ucount", colnames(df), value = TRUE)
        
        for (col in ucount_columns) {
            value1 <- row1[[col]]
            value2 <- row2[[col]]
            
            if (value1 >= 10 * value2) {
                cluster_result$fivePrimeH[[col]] <- c(value1, value2)
            } else if (value2 >= 10 * value1) {
                cluster_result$threePrimeH[[col]] <- c(value1, value2)
            } else {
                cluster_result$equalH[[col]] <- c(value1, value2)
            }
        }
        
        cluster_name <- paste0("yrow", cluster_id_value)
        data_info[[cluster_name]] <- cluster_result
    }
    
    return(data_info)
}

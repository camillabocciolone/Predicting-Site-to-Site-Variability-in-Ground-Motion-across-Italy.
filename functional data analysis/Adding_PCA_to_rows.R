# Load required library
library(tidyverse)

df_itaca <- read.csv("ITACAs2s_2_flatfile/ITACA_and_maps_V1.csv", sep=" ")

perform_pca <- function(df, columns_to_pca, n_components = NULL, scale. = TRUE, save_prefix = "pca_results") {
  
  # Extract data for PCA
  pca_data <- df %>% select(all_of(columns_to_pca))
  
  # Total number of rows before filtering
  total_rows <- nrow(df)
  
  # Filter out rows with NA or Inf values in PCA columns
  valid_rows <- pca_data %>%
    filter(if_all(everything(), ~ !is.na(.) & is.finite(.)))
  
  # Number of valid rows
  valid_row_count <- nrow(valid_rows)
  removed_row_count <- total_rows - valid_row_count
  
  # Stop if nothing remains
  if (valid_row_count == 0) {
    stop("No complete cases remain after removing NA/Inf values.")
  }
  
  # Perform PCA
  pca <- prcomp(valid_rows, center = TRUE, scale. = scale.)
  
  # Determine number of PCs to keep
  if (is.null(n_components)) {
    n_components <- ncol(pca$x)
  } else {
    n_components <- min(n_components, ncol(pca$x))
  }
  
  # Get PCA scores
  pca_scores <- as.data.frame(pca$x[, 1:n_components])
  colnames(pca_scores) <- paste0("PC", 1:n_components)
  
  # Apply the same filtering to the full dataframe
  df_clean <- df %>%
    filter(if_all(all_of(columns_to_pca), ~ !is.na(.) & is.finite(.)))
  
  # Add scores to the cleaned dataframe
  df_with_pcs <- bind_cols(df_clean, pca_scores)
  
  # Save loadings and explained variance
  loadings <- as.data.frame(pca$rotation[, 1:n_components])
  write.csv(loadings, paste0(save_prefix, "_loadings.csv"), row.names = TRUE)
  
  explained_variance <- summary(pca)$importance[2, 1:n_components]
  write.csv(data.frame(PC = paste0("PC", 1:n_components),
                       ExplainedVariance = explained_variance),
            paste0(save_prefix, "_explained_variance.csv"), row.names = FALSE)
  
  # Return both the data and some metadata
  return(list(
    df_with_pcs = df_with_pcs,
    total_rows = total_rows,
    rows_removed = removed_row_count,
    rows_used = valid_row_count
  ))
}



# We want to apply PCA on dS2S columns
dS2S_column_names = c("dS2S_SA_T0_010", "dS2S_SA_T0_025", "dS2S_SA_T0_040",
                      "dS2S_SA_T0_050", "dS2S_SA_T0_070", "dS2S_SA_T0_100",
                      "dS2S_SA_T0_150", "dS2S_SA_T0_200", "dS2S_SA_T0_250",
                      "dS2S_SA_T0_300", "dS2S_SA_T0_350", "dS2S_SA_T0_400",
                      "dS2S_SA_T0_450", "dS2S_SA_T0_500", "dS2S_SA_T0_600",
                      "dS2S_SA_T0_700", "dS2S_SA_T0_750", "dS2S_SA_T0_800",
                      "dS2S_SA_T0_900", "dS2S_SA_T1_000", "dS2S_SA_T1_200",
                      "dS2S_SA_T1_400", "dS2S_SA_T1_600", "dS2S_SA_T1_800",
                      "dS2S_SA_T2_000", "dS2S_SA_T2_500", "dS2S_SA_T3_000",
                      "dS2S_SA_T3_500", "dS2S_SA_T4_000", "dS2S_SA_T4_500",
                      "dS2S_SA_T5_000", "dS2S_SA_T6_000", "dS2S_SA_T7_000",
                      "dS2S_SA_T8_000", "dS2S_SA_T9_000", "dS2S_SA_T10_000")

result <- perform_pca(df_itaca, columns_to_pca = dS2S_column_names)

# Access results
df_pca <- result$df_with_pcs
cat("Total rows:", result$total_rows, "\n")
cat("Rows removed due to NA/Inf:", result$rows_removed, "\n")
cat("Rows used for PCA:", result$rows_used, "\n")


write.table(df_pca, "ITACAs2s_2_flatfile/ITACA_V3_maps_and_PCA.csv", col.names = TRUE, sep=";")

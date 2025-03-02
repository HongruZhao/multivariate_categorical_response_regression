  #------------------------------------------------
  # Minimal R script to create the final dataset
  #------------------------------------------------
  # setwd("/panfs/jay/groups/27/shenx/zhao1118/AD-real-data")

  # 1. Load the GEOquery package
  library("GEOquery")
  
  # 2. Download the GDS record and convert to ExpressionSet
  gds <- getGEO("GDS3268")
  eset <- GDS2eSet(gds, do.log2 = FALSE)
  
  # 3. Extract the expression matrix and phenotype data
  expr_matrix <- exprs(eset)
  pheno_data  <- pData(eset)
  
  # 4. Remove the first column in the phenotype data
  pheno_data_new <- pheno_data[, -1, drop = FALSE]
  
  # 5. Filter out rows in expr_matrix with too many NAs
  threshold_fraction <- 1
  n_cols            <- ncol(expr_matrix)
  non_na_counts     <- rowSums(!is.na(expr_matrix))
  keep_rows         <- non_na_counts >= threshold_fraction * n_cols
  expr_matrix_filtered <- expr_matrix[keep_rows, ]
  
  # 6. Align sample IDs and reorder if needed
  commonIDs <- intersect(rownames(pheno_data_new), colnames(expr_matrix_filtered))
  pheno_data_new2       <- pheno_data_new[commonIDs, , drop = FALSE]
  expr_matrix_filtered2 <- expr_matrix_filtered[, commonIDs, drop = FALSE]
  
  # 7. Transpose the filtered expression data and combine
  expr_matrix_filtered_t <- t(expr_matrix_filtered2)
  ultimate_dataset       <- cbind(pheno_data_new2, expr_matrix_filtered_t)
  
  # 8. Optionally remove the 4th column in ultimate_dataset (as in your code):
  dataset <- ultimate_dataset[, -4]
  
  # 'dataset' is now your final output (202 x (4 + number_of_genes - 1)).
  write.csv(dataset, "dataset.csv", row.names = FALSE)
  
  

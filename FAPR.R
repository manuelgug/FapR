##############################
##### PHASING WITH FAPR ######
##############################

library(future.apply)

FAPR <- function(MICROHAP_DATA, COI_DATA = NULL, verbose = TRUE) {
  
  unique_resmarkers <- unique(MICROHAP_DATA$resmarker)
  
  RESULTS_FINAL <- data.frame(
    SampleID = character(0), 
    matrix(character(0), nrow = 0, ncol = length(unique_resmarkers)),
    HAPLO_FREQ = numeric(0), 
    HAPLO_FREQ_RECALC = numeric(0)
  )
  
  # Rename the columns to include the resmarkers
  colnames(RESULTS_FINAL)[2:(1 + length(unique_resmarkers))] <- unique_resmarkers
  
  # Extract unique sample names
  unique_samples <- unique(MICROHAP_DATA$SampleID)
  
  # Set up parallel backend
  plan(multisession, workers = parallel::detectCores() - 2)  # Adjust the number of workers as needed
  
  # Use future_lapply to parallelize the loop
  results_list <- future_lapply(unique_samples, function(sample) {
    
    i_counter <- 0
    MOST_LIKELY_HAPLOS <- data.frame()
    MOST_LIKELY_HAPLOS_FREQS <- data.frame()
    
    unique_resmarkers <- unique(MICROHAP_DATA$resmarker)
    
    RESULTS <- data.frame(
      SampleID = character(0), 
      matrix(character(0), nrow = 0, ncol = length(unique_resmarkers)),
      HAPLO_FREQ = numeric(0), 
      HAPLO_FREQ_RECALC = numeric(0)
    )
    
    # Rename the columns to include the resmarkers
    colnames(RESULTS)[2:(1 + length(unique_resmarkers))] <- unique_resmarkers
    
    # 1) Select sample
    sID <- MICROHAP_DATA[MICROHAP_DATA$SampleID == sample,]
    sID$norm.reads.locus <- round(sID$norm.reads.locus, 4)
    
    # 2) Select sample's COI
    if (!is.null(COI_DATA)) {  
      COI <- COI_DATA[COI_DATA$SampleID == sample,]$COI  
    } else { 
      COI <- max(table(sID$resmarker)) 
    }
    
    # 3) Format data
    new_df <- data.frame(matrix(ncol = length(sID$resmarker), nrow = 1))
    colnames(new_df) <- sID$resmarker
    new_df[1, ] <- sID$Microhaplotype
    new_df <- rbind(new_df, sID$norm.reads.locus)
    
    unique_resmarkers <- unique(colnames(new_df))
    
    resulting_dataframes <- list()
    
    for (resmarker in unique_resmarkers) {
      columns <- which(names(new_df) == resmarker)
      df <- t(as.data.frame(new_df[, columns]))
      df <- as.data.frame(df)
      df$V2 <- as.numeric(df$V2)
      colnames(df) <- c(resmarker, "norm.reads.locus")
      rownames(df) <- NULL
      resulting_dataframes[[resmarker]] <- df
    }
    
    # Correct abundances
    resulting_dataframes <- lapply(resulting_dataframes, function(df) df[order(-df$norm.reads.locus), ])
    n_alleles_amps <- sapply(resulting_dataframes, nrow)
    
    adjust_by_allele_count <- function(df_list) {
      max_rows <- max(sapply(df_list, nrow))
      mean_matrix <- matrix(NA, nrow = max_rows, ncol = length(df_list))
      for (i in seq_along(df_list)) {
        norm_reads <- df_list[[i]]$norm.reads.locus
        mean_matrix[1:length(norm_reads), i] <- norm_reads
      }
      row_means <- rowMeans(mean_matrix, na.rm = TRUE)
      adjusted_dfs <- lapply(df_list, function(df) {
        df$norm.reads.locus <- row_means[1:nrow(df)]
        return(df)
      })
      return(adjusted_dfs)
    }
    
    allele_group_list <- split(resulting_dataframes, n_alleles_amps)
    adjusted_allele_groups <- lapply(allele_group_list, adjust_by_allele_count)
    resulting_dataframes <- do.call(c, adjusted_allele_groups)
    names(resulting_dataframes) <- sapply(resulting_dataframes, function(df) colnames(df)[[1]])
    
    alleles <- lapply(resulting_dataframes, function(df) df[[1]])
    freqs <- lapply(resulting_dataframes, function(df) df[[2]])
    
    comb_alleles_matrix <- expand.grid(alleles)
    comb_freqs_matrix <- expand.grid(freqs)
    
    if (nrow(comb_alleles_matrix) == 0) {
      return(NULL)
    }
    
    if (dim(comb_alleles_matrix)[1] != 1) {
      
      while (COI > i_counter) {
        i_counter <- i_counter + 1
        
        comb_freqs_matrix$probs <- Reduce(`*`, comb_freqs_matrix[unique_resmarkers])
        
        highest_prob <- as.numeric(which.max(comb_freqs_matrix$probs))
        
        most_likely_hap <- paste(as.matrix(comb_alleles_matrix[highest_prob, ]), collapse = "_")
        
        if (verbose) {
          print(paste0(sample, " #", i_counter, ": ", most_likely_hap, " is the most likely true haplotype."))  
        }
        
        MOST_LIKELY_HAPLOS <- rbind(MOST_LIKELY_HAPLOS, comb_alleles_matrix[highest_prob, ])
        temp <- comb_freqs_matrix[highest_prob, ]
        temp$HAPLO_FREQ <- min(comb_freqs_matrix[highest_prob, 1:length(unique_resmarkers)])
        temp$HAPLO_FREQ_RECALC <- NA
        MOST_LIKELY_HAPLOS_FREQS <- rbind(MOST_LIKELY_HAPLOS_FREQS, temp)
        
        min_allele_from_most_likely_hap <- min(comb_freqs_matrix[highest_prob, 1:length(unique_resmarkers)])
        
        row_to_match <- as.matrix(comb_alleles_matrix[highest_prob, ])
        mask <- sapply(colnames(comb_alleles_matrix), function(col_name) {
          comb_alleles_matrix[, col_name] == row_to_match[, col_name]
        })
        
        comb_freqs_matrix <- comb_freqs_matrix[, 1:length(unique_resmarkers)]
        comb_freqs_matrix[mask] <- comb_freqs_matrix[mask] - min_allele_from_most_likely_hap
        
        MOST_LIKELY_HAPLOS_FREQS$HAPLO_FREQ_RECALC <- MOST_LIKELY_HAPLOS_FREQS$HAPLO_FREQ / sum(MOST_LIKELY_HAPLOS_FREQS$HAPLO_FREQ)
        
        RESULTS <- cbind(SampleID = sample, MOST_LIKELY_HAPLOS, HAPLO_FREQ = MOST_LIKELY_HAPLOS_FREQS$HAPLO_FREQ, HAPLO_FREQ_RECALC = MOST_LIKELY_HAPLOS_FREQS$HAPLO_FREQ_RECALC)
      }  
      
    } else { 
      RESULTS <- cbind(SampleID = sample, comb_alleles_matrix, HAPLO_FREQ = 1, HAPLO_FREQ_RECALC = 1)
    }
    
    return(RESULTS)
  })
  
  # Combine the results from all parallel computations
  RESULTS_FINAL <- do.call(rbind, results_list)
  
  # Add haplotype name
  RESULTS_FINAL$haplotype <- apply(RESULTS_FINAL[, 2:(ncol(RESULTS_FINAL)-2)], 1, paste, collapse = "_")
  
  return(RESULTS_FINAL)
}

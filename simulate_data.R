

library(ggbeeswarm)
library(ggplot2)
library(dplyr)

############# GENERATE HAPLOTYPE PROPORTIONS

# Function to generate haplotype proportions
generate_numbers_summing_to_1 <- function(n) {
  # Generate n-1 random numbers between 0 and 1
  random_numbers <- runif(n - 1, 0, 1)
  
  # Sort the random numbers in ascending order
  sorted_numbers <- sort(random_numbers)
  
  # Calculate the differences between consecutive numbers
  differences <- c(sorted_numbers[1], diff(sorted_numbers), 1 - sorted_numbers[length(sorted_numbers)])
  
  return(differences)
}


# BENCHMARK

RESULTS_BENCH_ALL <- data.frame()
RESULTS_BENCH_ALL_FREQS <- list()

for (n in c(2:7)){
  
  n_haplos <- n
  
  print(paste0("Looping through ", n, " haplotypes"))
  
  for (i in 1:100) {
  
  proportions <- generate_numbers_summing_to_1(n_haplos)
  
  cumulative_sums <- cumsum(proportions)
  
  in_sample_freqs <- list()
  
  # Iterate over the indices of the proportions vector
  for (i in seq_along(proportions)) {
    # Calculate the cumulative sum up to index i
    cum_sum <- sum(proportions[1:i])
    
    # Check if the cumulative sum is equal to 1
    if (cum_sum == 1) {
      in_sample_freqs[[i]] <- 1
    } else {
      # Concatenate the cumulative sum up to index i with the remaining proportions
      in_sample_freqs[[i]] <- c(cum_sum, proportions[(i + 1):length(proportions)])
    }
  }
  
  names(in_sample_freqs) <- rev(seq(1:length(in_sample_freqs)))
  
  in_sample_freqs
  
  
  
  ############# GENERATE HAPLOTYPES (((NAMES DON'T MEAN ANYTHING!! IT'S JUST FOR FAPR COMPATBITLITY)))
  
  res_loci <- c("dhps_431", "dhps_437", "dhps_540", "dhps_581", "dhfr_51","dhfr_59", "dhfr_108")
  res_variants <- c("A", "B", "C", "D", "E", "F", "G")
  
  
  genos <- data.frame()
  
  for (locus in res_loci){
    
    for (res_variant in res_variants){
      
      # set variants
      nvar <- sample(1:n_haplos, 1)
      genotype <- res_variants[1:nvar]
      
      geno <- as.data.frame(cbind(locus, genotype))
      
    }
    
    #concat loci in the same df
    genos <- rbind(genos, geno)
    
  }
  
  #make the last locus always n_haplos number of variants 
  genos <- genos[genos$locus != "dhfr_108",]
  dhfr_108 <- as.data.frame(cbind(locus = "dhfr_108", genotype = res_variants[1:n_haplos]))
  genos<- rbind(genos, dhfr_108)
  
  
  ############# MERGE HAPLOTYPES WITH FREQUENCIES
  
  locus_counts <- table(genos$locus)
  
  # Create a new column filled with NA
  genos$freq <- NA
  
  # Iterate over unique loci
  for (locus in unique(genos$locus)) {
    # Get count for current locus
    count <- locus_counts[locus]
    # Check if count exists in in_sample_freqs
    if (count %in% seq_along(in_sample_freqs)) {
      # Index into in_sample_freqs and assign to corresponding rows
      genos$freq[genos$locus == locus] <- in_sample_freqs[[as.character(count)]]
    }
  }
  
  colnames(genos) <- c("resmarker", "AA", "norm.reads.locus")
  
  genos <- cbind(SampleID = "shit_sample", genos)
  
  resmarkers_table <- genos
  
  
  
  
  
  
  ################## MAIN LOOP ################## 
  
  unique_samples <- unique(resmarkers_table$SampleID)
  RESULTS_FINAL <- data.frame(SampleID = character(0), dhps_431 = character(0), dhps_437 = character(0), dhps_540 = character(0), dhps_581 = character(0), dhfr_51 = character(0), dhfr_59 = character(0), dhfr_108 = character(0), HAPLO_FREQ = numeric(0), HAPLO_FREQ_RECALC = numeric(0))
  
  # INIT LOOP HERE!
  for (sample in unique_samples){
    
    #sample <-"NQS1_-_QC2_10k_S249"
    
    i_counter <- 0
    MOST_LIKELY_HAPLOS <- data.frame()
    MOST_LIKELY_HAPLOS_FREQS <- data.frame()
    RESULTS <- data.frame(SampleID = character(0), dhps_431 = character(0), dhps_437 = character(0), dhps_540 = character(0), dhps_581 = character(0), dhfr_51 = character(0), dhfr_59 = character(0), dhfr_108 = character(0), HAPLO_FREQ = numeric(0), HAPLO_FREQ_RECALC = numeric(0), HAPLO_FREQ = numeric(0), HAPLO_FREQ_RECALC = numeric(0))
    
    # 1) select sample
    sID <- resmarkers_table[resmarkers_table$SampleID == sample,]
    
    # 2) select sample's COI
    #COI<- round(moire_output[moire_output$sample_id == sample,]["post_coi_med"]) #truncated post_coi_mean seems to work best for controls. however, needs more testing
    
    # 3) format data
    new_df <- data.frame(matrix(ncol = length(sID$resmarker), nrow=1))
    colnames(new_df) <- sID$resmarker
    new_df[1,] <-sID$AA
    new_df <- rbind(new_df, sID$norm.reads.locus )
    
    unique_resmarkers <- unique(colnames(new_df))
    
    resulting_dataframes <- list()
    # Loop through each unique resmarker (colname)
    for (resmarker in unique_resmarkers) {
      columns <- which(names(new_df) == resmarker)
      df <- t(as.data.frame(new_df[, columns]))
      df <- as.data.frame(df)
      df$V2 <- as.numeric(df$V2)
      colnames(df) <- c(resmarker, "norm.reads.locus")
      rownames(df) <- NULL
      resulting_dataframes[[resmarker]] <- df
    }
    
    alleles<-list(resulting_dataframes$dhps_431$dhps_431,
                  resulting_dataframes$dhps_437$dhps_437, 
                  resulting_dataframes$dhps_540$dhps_540, 
                  resulting_dataframes$dhps_581$dhps_581,
                  resulting_dataframes$dhfr_51$dhfr_51,
                  resulting_dataframes$dhfr_59$dhfr_59,
                  resulting_dataframes$dhfr_108$dhfr_108) #order is important
    
    freqs<-list(resulting_dataframes$dhps_431$norm.reads.locus,
                resulting_dataframes$dhps_437$norm.reads.locus, 
                resulting_dataframes$dhps_540$norm.reads.locus, 
                resulting_dataframes$dhps_581$norm.reads.locus,
                resulting_dataframes$dhfr_51$norm.reads.locus,
                resulting_dataframes$dhfr_59$norm.reads.locus,
                resulting_dataframes$dhfr_108$norm.reads.locus) #order is important
    
    comb_alleles <- expand.grid(alleles)
    comb_freqs <- expand.grid(freqs)
    
    # Check if comb_alleles is empty, and if so, skip to the next sample
    if (nrow(comb_alleles) == 0) {
      cat("Skipping sample", sample, "\n")
      next
    }
    
    comb_alleles_matrix <- as.data.frame(comb_alleles)
    colnames(comb_alleles_matrix) <- c("dhps_431", "dhps_437", "dhps_540", "dhps_581", "dhfr_51","dhfr_59", "dhfr_108")
    comb_freqs_matrix <- as.data.frame(comb_freqs)
    colnames(comb_freqs_matrix) <- c("dhps_431", "dhps_437", "dhps_540", "dhps_581", "dhfr_51","dhfr_59", "dhfr_108")
    
    # 4) phase
    if (dim(comb_alleles_matrix)[1] != 1){ #basically, don't process monoallelic samples 'cause they make the loop crash
      
      while (dim(MOST_LIKELY_HAPLOS_FREQS)[1] == 0 || 1-sum(RESULTS$HAPLO_FREQ) > 0.01) { ## PULIR CONDICIÓN? (previous condition: i_counter != COI && 1-sum(RESULTS$HAPLO_FREQ) > 0.0001)
        
        i_counter <- i_counter + 1
        
        # Calculate probs if all haplotypes were present
        comb_freqs_matrix$probs <- comb_freqs_matrix$dhps_431 * comb_freqs_matrix$dhps_437 * comb_freqs_matrix$dhps_540 * comb_freqs_matrix$dhps_581 * comb_freqs_matrix$dhfr_51 * comb_freqs_matrix$dhfr_59  * comb_freqs_matrix$dhfr_108
        
        #remove haplotypes with prob = 0
        #comb_freqs_matrix <- subset(comb_freqs_matrix, probs != 0)
        
        # Calculate SD and CV
        comb_freqs_matrix$freq_mean <- rowMeans(comb_freqs_matrix, na.rm = TRUE)
        comb_freqs_matrix$SD <- apply(comb_freqs_matrix[, 1:7], 1, sd)
        comb_freqs_matrix$CV <- (comb_freqs_matrix$SD / comb_freqs_matrix$freq_mean)
        
        ## Select the "BEST" haplo: highest prob and lowest CV
        lowest_CV <- which.min(comb_freqs_matrix$CV)
        highest_prob <- as.numeric(which.max(comb_freqs_matrix$probs))
        
        # #do CV and probs agree with each other?
        # if (lowest_CV == highest_prob) {
        #   most_likely_hap <- paste(as.matrix(comb_alleles_matrix[highest_prob, ]), collapse = "_")
        #   print(paste(sample, "#", i_counter, ":", most_likely_hap, "is the most likely true haplotype.", collapse = " "))
        # } else {
        #   most_likely_hap1 <- paste(as.matrix(comb_alleles_matrix[highest_prob, ]), collapse = "_")
        #   most_likely_hap2 <- paste(as.matrix(comb_alleles_matrix[lowest_CV, ]), collapse = "_")
        #   print(paste(sample, "#", i_counter, ": One of", most_likely_hap1, "and", most_likely_hap2, "is the most likely true haplotype. Visually examine the plot."))
        # }
        
        # Append most likely haplo
        MOST_LIKELY_HAPLOS <- rbind(MOST_LIKELY_HAPLOS, comb_alleles_matrix[highest_prob, ])
        temp <- comb_freqs_matrix[highest_prob, ]
        temp$HAPLO_FREQ <- min(comb_freqs_matrix[highest_prob, 1:7])
        temp$HAPLO_FREQ_RECALC <- NA
        MOST_LIKELY_HAPLOS_FREQS <- rbind(MOST_LIKELY_HAPLOS_FREQS, temp)
        
        # Select minimum allele freq from the most likely haplotype
        min_allele_from_most_lilely_hap <- min(comb_freqs_matrix[highest_prob, 1:7])
        
        # Boolean mask to detect alleles that are present on the most likely haplotype
        row_to_match <- as.matrix(comb_alleles_matrix[highest_prob, ])
        mask <- sapply(colnames(comb_alleles_matrix), function(col_name) {
          comb_alleles_matrix[, col_name] == row_to_match[, col_name]
        })
        
        # Subtract min_allele_from_most_likely_hap from the cells where mask is TRUE and ignore the specified column
        comb_freqs_matrix <- comb_freqs_matrix[, 1:7]
        comb_freqs_matrix[mask] <- comb_freqs_matrix[mask] - min_allele_from_most_lilely_hap
        
        #recalculate proportions of final haplos
        MOST_LIKELY_HAPLOS_FREQS$HAPLO_FREQ_RECALC <- MOST_LIKELY_HAPLOS_FREQS$HAPLO_FREQ / sum(MOST_LIKELY_HAPLOS_FREQS$HAPLO_FREQ)
        
        RESULTS <- cbind(SampleID = sample, MOST_LIKELY_HAPLOS, HAPLO_FREQ = MOST_LIKELY_HAPLOS_FREQS$HAPLO_FREQ, HAPLO_FREQ_RECALC = MOST_LIKELY_HAPLOS_FREQS$HAPLO_FREQ_RECALC)
      }  
      
    }else{ 
      
      #FORMAT AND ADD MONOALLELIC SAMPLES HERE
      RESULTS <- cbind(SampleID = sample, comb_alleles_matrix, HAPLO_FREQ = 1, HAPLO_FREQ_RECALC = 1)
    }
    
    #DONE
    RESULTS_FINAL <- rbind(RESULTS, RESULTS_FINAL)
  }
  
  
  RESULTS_FINAL$haplotype <- paste(RESULTS_FINAL$dhps_431, RESULTS_FINAL$dhps_437, RESULTS_FINAL$dhps_540, RESULTS_FINAL$dhps_581, RESULTS_FINAL$dhfr_51, RESULTS_FINAL$dhfr_59, RESULTS_FINAL$dhfr_108, sep = "_")
  RESULTS_FINAL_multiallelic <- RESULTS_FINAL[RESULTS_FINAL$HAPLO_FREQ_RECALC < 1, ]
  
  
  
  ############# EVALUATION
  
  # Extract observed and expected values
  observed <- RESULTS_FINAL$HAPLO_FREQ_RECALC  #USE RESULTS_FINAL_multiallelic IF NEEDED
  expected <- rev(sort(in_sample_freqs[[1]]))
  
  n_exp <- length(expected)
  n_obs <- length(observed)
  
  difference <- n_exp - n_obs
  
  # Ensure both vectors are of the same length by adding zeros
  max_len <- max(length(observed), length(expected))
  observed_padded <- c(observed, rep(0, max_len - length(observed)))
  expected_padded <- c(expected, rep(0, max_len - length(expected)))
  
  # Create dataframe for frequencies
  FREQS <- data.frame(observed = observed_padded, expected = expected_padded)
  
  # Append the dataframe to the list
  RESULTS_BENCH_ALL_FREQS <- c(RESULTS_BENCH_ALL_FREQS, list(FREQS))
  
  # Truncate both vectors to the minimum length (thus, calculate RMSE only for true positives)
  min_len <- min(length(observed), length(expected))
  observed_truncated <- observed[1:min_len]
  expected_truncated <- expected[1:min_len]
  
  # Calculate RMSE
  rmse <- sqrt(mean((observed_truncated - expected_truncated)^2))
  
  # Combine the RMSE, PERC, and freqs_list into a data frame
  RESULTS_BENCH <- data.frame(RMSE = rmse, PERC = difference, n_haplos = paste0(n_haplos, "_haplos"))
  RESULTS_BENCH_ALL <- rbind(RESULTS_BENCH_ALL, RESULTS_BENCH)
  }
}

#RESULTS_BENCH_ALL_FREQS
# Combine all dataframes in the list into a single dataframe
#RESULTS_BENCH_ALL_FREQS <- do.call(rbind, RESULTS_BENCH_ALL_FREQS)

#RESULTS_BENCH_ALL


stacked_input <- RESULTS_BENCH_ALL %>% 
  group_by(n_haplos, PERC) %>%
  summarize(counts = as.vector(table(PERC)))

# Define a lookup table
lookup_table <- data.frame(PERC = c("0", "1", "2", "3", "-1", "-2", "-3"),
                           acc = c("exact_haplos", "missing_1_haplo", "missing_2_haplos", "missing_3_haplos",
                                   "extra_1_haplo", "extra_2_haplos", "extra_3_haplos"))

# Merge the lookup table with stacked_input based on PERC
stacked_input <- merge(stacked_input, lookup_table, by = "PERC", all.x = TRUE)

stacked_input$acc <- factor(stacked_input$acc,
                            levels = c("extra_3_haplos", "extra_2_haplos", "extra_1_haplo",
                                       "missing_3_haplos", "exact_haplos", "missing_1_haplo", "missing_2_haplos"))

# Define custom colors
custom_colors <- c("exact_haplos" = "limegreen", "missing_1_haplo" = "orange", 
                   "missing_2_haplos" = "red2", "missing_3_haplos" = "black", 
                   "extra_1_haplo" = "pink", "extra_2_haplos" = "violet", 
                   "extra_3_haplos" = "purple")

# Create the stacked barplot
perc_plot <- ggplot(stacked_input, aes(x = n_haplos, y = counts, fill = acc)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "", y = "Accuracy") +
  scale_fill_manual(values = custom_colors) +  # Set custom colors
  theme_minimal()

perc_plot


max_rmse <- max(RESULTS_BENCH_ALL$RMSE)

rmse_plot <- ggplot(RESULTS_BENCH_ALL, aes(x = n_haplos, y = RMSE, color = n_haplos)) +
  geom_boxplot() +
  geom_quasirandom(alpha = 0.5)+
  labs(title = "", x = "", y = "RMSE for true positives") +
  theme_minimal() +
  ylim(0, max_rmse)+
  guides(color = "none")

rmse_plot

## AHORA METER RUIDO!!!! 


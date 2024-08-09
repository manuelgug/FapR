
library(gridExtra)
library(ggbeeswarm)
library(ggplot2)
library(dplyr)
library(reshape2)
library(zoo)
library(corrplot)


######----------------------------------------------------------------------------------------

set.seed(42069)

##### GENERATE HAPLOTYPE PROPORTIONS


# Function to generate haplotype proportions
generate_numbers_summing_to_1 <- function(n) {
  
  # Special case for n == 1
  if (n == 1) {
    return(1)
  }
  
  # Generate n-1 random numbers between 0.01 and 0.99
  random_numbers <- runif(n - 1, 0.01, 0.99) # 0.01 min because this is the most common MAF
  
  # Sort the random numbers in ascending order
  sorted_numbers <- sort(random_numbers)
  
  # Calculate the differences between consecutive numbers
  differences <- c(sorted_numbers[1], diff(sorted_numbers), 1 - sorted_numbers[length(sorted_numbers)])
  
  # Ensure the minimum value of 0.01
  differences <- pmax(differences, 0.01)
  
  # Adjust to ensure they sum to 1 after enforcing minimum
  total <- sum(differences)
  differences <- differences / total
  
  return(differences)
}


##### GENERATE DATA

max_haplos <- 5
individuals <- 1000

SIM_DATA <- list()

for (n in c(2:max_haplos)){
  
  n_haplos <- n
  
  print(paste0("Looping through ", n, " haplotypes"))
  
  for (i in 1:individuals) {

    proportions <- generate_numbers_summing_to_1(n_haplos)
    
    cumulative_sums <- cumsum(proportions)
    
    in_sample_freqs <- list()
    
    # Iterate over the indices of the proportions vector
    for (j in seq_along(proportions)) {
      # Calculate the cumulative sum up to index i
      cum_sum <- sum(proportions[1:j])
      
      # Check if the cumulative sum is equal to 1
      if (cum_sum == 1) {
        in_sample_freqs[[j]] <- 1
      } else {
        # Concatenate the cumulative sum up to index i with the remaining proportions
        in_sample_freqs[[j]] <- c(cum_sum, proportions[(j + 1):length(proportions)])
      }
    }
    
    names(in_sample_freqs) <- rev(seq(1:length(in_sample_freqs)))
    
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
    
    genos <- cbind(SampleID = paste0("shit_sample_", i,"_", n_haplos), genos)
    
    resmarkers_table <- genos
    
    # Append to SIM_DATA list
    SIM_DATA[[paste0("shit_sample_", i, "_", n_haplos)]] <- resmarkers_table
    
  }
}

#convert to a single df
SIM_DATA <- bind_rows(SIM_DATA)

#check data!

#min should not be below 0.01
summary(SIM_DATA$norm.reads.locus)

check_SIM_DATA <- SIM_DATA %>%
  group_by(SampleID, resmarker) %>%
  summarise(total_freq = sum(norm.reads.locus))

# all markers there? must be 7
length(unique(check_SIM_DATA$resmarker))

# all freqs = 1? must be 2100 == 300 * 7
sum(check_SIM_DATA$total_freq) == length(unique(check_SIM_DATA$SampleID)) * length(unique(check_SIM_DATA$resmarker))


#CREATE COI DATA
unique_sample_ids <- unique(SIM_DATA$SampleID)
coi_values <- sapply(unique_sample_ids, function(id) {
  substr(id, nchar(id), nchar(id))
})
COI_SIM_DATA <- data.frame(SampleID = unique_sample_ids, COI = as.numeric(coi_values), row.names = NULL)


#save data
saveRDS(SIM_DATA, paste0("SIM_DATA_clean_", "haps_", max_haplos,  "_ind_", individuals, ".RDS"))
saveRDS(COI_SIM_DATA, paste0("COI_SIM_DATA_clean", "haps_", max_haplos,  "_ind_", individuals, ".RDS"))


######----------------------------------------------------------------------------------------

### ADD NOISE TO SIMULATED DATA

# adding noise means that freqs for each resmarker vary a bit; for instance, it may not be 55-45 for all in a 2 haps sample, but rather dhfr_59: 57-43, dhfr_437: 50-50, dhps_540: 52-48
# This will increase the errors. How will FAPR perform?

SIM_DATA <- readRDS(paste0("SIM_DATA_clean_", "haps_", max_haplos,  "_ind_", individuals, ".RDS"))

# Function to introduce noise to norm.reads.locus values
introduce_noise <- function(data, max_change = 0, min_change = 0) {
  
  # Group data by SampleID and resmarker
  multi_aa_data <- data %>%
    group_by(SampleID, resmarker) %>%
    # Filter for rows with more than one AA
    filter(n_distinct(AA) > 1) %>%
    ungroup() %>%
    mutate(
      # Generate random factors to multiply the norm.reads.locus values
      random_factor = runif(n(), max(min_change, 1 - max_change), 1 + max_change),
      # Normalize the factors so that they sum up to 1
      random_factor = random_factor / sum(random_factor),
      # Update the norm.reads.locus values by multiplying them with the random factors
      norm.reads.locus = norm.reads.locus * random_factor
    ) %>%
    # Group by SampleID and resmarker again to calculate the sum of norm.reads.locus
    group_by(SampleID, resmarker) %>%
    mutate(sum_norm_reads = sum(norm.reads.locus)) %>%
    ungroup() %>%
    # Adjust the norm.reads.locus values to ensure they sum up to 1
    mutate(norm.reads.locus = norm.reads.locus / sum_norm_reads) %>%
    select(-sum_norm_reads)
  
  # Filter rows where there is only one AA and combine with the noisy data
  noisy_data <- bind_rows(
    data %>% group_by(SampleID, resmarker) %>% filter(n_distinct(AA) == 1),
    multi_aa_data
  )
  
  # Remove the random_factor column
  noisy_data <- as.data.frame(noisy_data[,-5])
  
  return(noisy_data)
}


# Define the range of max_change values
max_change_values <- seq(0, 0.4, by = 0.05)
min_change <- 0

# Loop through each max_change value
for (max_change in max_change_values) {
  
  # Apply the introduce_noise function to SIM_DATA
  SIM_DATA_noisy <- introduce_noise(SIM_DATA, max_change = max_change, min_change = min_change)
  
  # Get the corresponding indices of rows in noisy_SIM_DATA based on SampleID, resmarker, and AA
  indices <- match(apply(SIM_DATA[c("SampleID", "resmarker", "AA")], 1, paste, collapse = "_"), 
                   apply(SIM_DATA_noisy[c("SampleID", "resmarker", "AA")], 1, paste, collapse = "_"))
  
  # Order noisy_SIM_DATA according to the indices
  SIM_DATA_noisy <- SIM_DATA_noisy[indices, ]
  
  #CREATE COI DATA
  unique_sample_ids <- unique(SIM_DATA_noisy$SampleID)
  coi_values <- sapply(unique_sample_ids, function(id) {
    substr(id, nchar(id), nchar(id))
  })
  COI_SIM_DATA <- data.frame(SampleID = unique_sample_ids, COI = as.numeric(coi_values), row.names = NULL)
  
  # Save the noisy data frame to an RDS file
  filename <- paste0("SIM_DATA_", "haps_", max_haplos,  "_ind_", individuals, "_max_change_", max_change, ".RDS")
  saveRDS(SIM_DATA_noisy, filename)
  
  filename_coi <- paste0("SIM_DATA_COI_", "haps_", max_haplos,  "_ind_", individuals, "_max_change_", max_change, ".RDS")
  saveRDS(COI_SIM_DATA, filename_coi)
  
  # Optionally print a message to indicate progress
  print(paste("Saved file:", filename))
}


##############################----------------------------------------------------------------------------------------
##### PHASING WITH FAPR ######
##############################

FAPR <- function(SIM_DATA, COI_SIM_DATA,  verbose = TRUE) {
  
  # Initialize the final results data frame
  RESULTS_FINAL <- data.frame(SampleID = character(0), dhps_431 = character(0), dhps_437 = character(0), 
                              dhps_540 = character(0), dhps_581 = character(0), dhfr_51 = character(0), 
                              dhfr_59 = character(0), dhfr_108 = character(0), HAPLO_FREQ = numeric(0), 
                              HAPLO_FREQ_RECALC = numeric(0))
  
  # Extract unique sample names
  unique_samples <- unique(SIM_DATA$SampleID)
  
  # Loop through each sample
  for (sample in unique_samples) {
    
    i_counter <- 0
    MOST_LIKELY_HAPLOS <- data.frame()
    MOST_LIKELY_HAPLOS_FREQS <- data.frame()
    RESULTS <- data.frame(SampleID = character(0), dhps_431 = character(0), dhps_437 = character(0), 
                          dhps_540 = character(0), dhps_581 = character(0), dhfr_51 = character(0), 
                          dhfr_59 = character(0), dhfr_108 = character(0), HAPLO_FREQ = numeric(0), 
                          HAPLO_FREQ_RECALC = numeric(0))
    
    # 1) Select sample
    sID <- SIM_DATA[SIM_DATA$SampleID == sample,]
    sID$norm.reads.locus <- round(sID$norm.reads.locus, 4)
    
    # 2) Select sample's COI
    COI <- COI_SIM_DATA[COI_SIM_DATA$SampleID == sample,]$COI  
    
    # 3) Format data
    new_df <- data.frame(matrix(ncol = length(sID$resmarker), nrow = 1))
    colnames(new_df) <- sID$resmarker
    new_df[1, ] <- sID$AA
    new_df <- rbind(new_df, sID$norm.reads.locus)
    
    unique_resmarkers <- unique(colnames(new_df))
    
    # Initialize list to store data frames
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
    
    ####### Correct abundances #######

    # Function to order each data frame by norm.reads.locus in descending order
    order_by_norm_reads <- function(df) {
      df[order(-df$norm.reads.locus), ]
    }

    resulting_dataframes <- lapply(resulting_dataframes, order_by_norm_reads)
    n_alleles_amps <- sapply(resulting_dataframes, nrow)

    # Function to get the mean of norm.reads.locus for each row position within a group
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

    ########################################################3
    
    alleles <- list(resulting_dataframes$dhps_431$dhps_431,
                    resulting_dataframes$dhps_437$dhps_437, 
                    resulting_dataframes$dhps_540$dhps_540, 
                    resulting_dataframes$dhps_581$dhps_581,
                    resulting_dataframes$dhfr_51$dhfr_51,
                    resulting_dataframes$dhfr_59$dhfr_59,
                    resulting_dataframes$dhfr_108$dhfr_108) # Order is important
    
    freqs <- list(resulting_dataframes$dhps_431$norm.reads.locus,
                  resulting_dataframes$dhps_437$norm.reads.locus, 
                  resulting_dataframes$dhps_540$norm.reads.locus, 
                  resulting_dataframes$dhps_581$norm.reads.locus,
                  resulting_dataframes$dhfr_51$norm.reads.locus,
                  resulting_dataframes$dhfr_59$norm.reads.locus,
                  resulting_dataframes$dhfr_108$norm.reads.locus) # Order is important
    
    comb_alleles <- expand.grid(alleles)
    comb_freqs <- expand.grid(freqs)
    
    # Check if comb_alleles is empty, and if so, skip to the next sample
    if (nrow(comb_alleles) == 0) {
      cat("Skipping sample", sample, "\n")
      next
    }
    
    comb_alleles_matrix <- as.data.frame(comb_alleles)
    colnames(comb_alleles_matrix) <- c("dhps_431", "dhps_437", "dhps_540", "dhps_581", "dhfr_51", "dhfr_59", "dhfr_108")
    comb_freqs_matrix <- as.data.frame(comb_freqs)
    colnames(comb_freqs_matrix) <- c("dhps_431", "dhps_437", "dhps_540", "dhps_581", "dhfr_51", "dhfr_59", "dhfr_108")
    
    # 4) Phase
    if (dim(comb_alleles_matrix)[1] != 1) { # Skip monoallelic samples as they make the loop crash
      
      while (COI > i_counter) {
        i_counter <- i_counter + 1
        
        # Calculate probs if all haplotypes were present
        comb_freqs_matrix$probs <- comb_freqs_matrix$dhps_431 * comb_freqs_matrix$dhps_437 * 
          comb_freqs_matrix$dhps_540 * comb_freqs_matrix$dhps_581 * 
          comb_freqs_matrix$dhfr_51 * comb_freqs_matrix$dhfr_59 * 
          comb_freqs_matrix$dhfr_108
        highest_prob <- as.numeric(which.max(comb_freqs_matrix$probs))
        
        # Extract most probable haplo
        most_likely_hap <- paste(as.matrix(comb_alleles_matrix[highest_prob, ]), collapse = "_")
        
        if (verbose){
          print(paste0(sample, " #", i_counter, ": ", most_likely_hap, " is the most likely true haplotype."))  
        }
        
        # Append most likely haplo
        MOST_LIKELY_HAPLOS <- rbind(MOST_LIKELY_HAPLOS, comb_alleles_matrix[highest_prob, ])
        temp <- comb_freqs_matrix[highest_prob, ]
        temp$HAPLO_FREQ <- min(comb_freqs_matrix[highest_prob, 1:7])
        temp$HAPLO_FREQ_RECALC <- NA
        MOST_LIKELY_HAPLOS_FREQS <- rbind(MOST_LIKELY_HAPLOS_FREQS, temp)
        
        # Select minimum allele freq from the most likely haplotype
        min_allele_from_most_likely_hap <- min(comb_freqs_matrix[highest_prob, 1:7])
        
        # Boolean mask to detect alleles that are present on the most likely haplotype
        row_to_match <- as.matrix(comb_alleles_matrix[highest_prob, ])
        mask <- sapply(colnames(comb_alleles_matrix), function(col_name) {
          comb_alleles_matrix[, col_name] == row_to_match[, col_name]
        })
        
        # Subtract min_allele_from_most_likely_hap from the cells where mask is TRUE and ignore the specified column
        comb_freqs_matrix <- comb_freqs_matrix[, 1:7]
        comb_freqs_matrix[mask] <- comb_freqs_matrix[mask] - min_allele_from_most_likely_hap
        
        # Recalculate proportions of final haplos
        MOST_LIKELY_HAPLOS_FREQS$HAPLO_FREQ_RECALC <- MOST_LIKELY_HAPLOS_FREQS$HAPLO_FREQ / sum(MOST_LIKELY_HAPLOS_FREQS$HAPLO_FREQ)
        
        RESULTS <- cbind(SampleID = sample, MOST_LIKELY_HAPLOS, HAPLO_FREQ = MOST_LIKELY_HAPLOS_FREQS$HAPLO_FREQ, HAPLO_FREQ_RECALC = MOST_LIKELY_HAPLOS_FREQS$HAPLO_FREQ_RECALC)
      }  
      
    } else { 
      
      # Format and add monoallelic samples here
      RESULTS <- cbind(SampleID = sample, comb_alleles_matrix, HAPLO_FREQ = 1, HAPLO_FREQ_RECALC = 1)
    }
    
    # Append results to the final dataframe
    RESULTS_FINAL <- rbind(RESULTS, RESULTS_FINAL)
  }
  
  # Add haplotype name
  RESULTS_FINAL$haplotype <- paste(RESULTS_FINAL$dhps_431, RESULTS_FINAL$dhps_437, 
                                   RESULTS_FINAL$dhps_540, RESULTS_FINAL$dhps_581, 
                                   RESULTS_FINAL$dhfr_51, RESULTS_FINAL$dhfr_59, 
                                   RESULTS_FINAL$dhfr_108, sep = "_")
  
  return(RESULTS_FINAL)
}

# Get all SIM_DATA and COI_SIM_DATA files
sim_data_files <- list.files(pattern = "^SIM_DATA_haps_.*_ind_.*_max_change_.*\\.RDS")
SIM_DATA_list <- setNames(lapply(sim_data_files, readRDS), sim_data_files)
coi_data_files <- list.files(pattern = "^SIM_DATA_COI_haps_.*_ind_.*_max_change_.*\\.RDS")
COI_SIM_DATA_list <- setNames(lapply(coi_data_files, readRDS), coi_data_files)

print(paste("Loaded", length(SIM_DATA_list), "SIM_DATA files and", length(COI_SIM_DATA_list), "COI_SIM_DATA files into named lists."))


## RUN FAPR!!
for (data in 1:length(SIM_DATA_list)){
  
  data_name <- names(SIM_DATA_list)[data]
  print(paste0("------ Processing Dataset: ", data_name, " -------"))
  
  # 1) subset data
  SIM_DATA <- SIM_DATA_list[[data]]
  COI_SIM_DATA <- COI_SIM_DATA_list[[data]]
  
  # 2) run fapr
  RESULTS_FAPR <- FAPR(SIM_DATA, COI_SIM_DATA, verbose = FALSE)
  
  # 3) save results
  saveRDS(RESULTS_FAPR, paste0("FAPR_RESULTS_", data_name))

}


######----------------------------------------------------------------------------------------

############# EVALUATION OF noisy o clean data aquí

# 1) arrange data

# Get all SIM_DATA and COI_SIM_DATA files
sim_data_files <- list.files(pattern = "^SIM_DATA_haps_.*_ind_.*_max_change_.*\\.RDS")
SIM_DATA_list <- setNames(lapply(sim_data_files, readRDS), sim_data_files)
coi_data_files <- list.files(pattern = "^SIM_DATA_COI_haps_.*_ind_.*_max_change_.*\\.RDS")
COI_SIM_DATA_list <- setNames(lapply(coi_data_files, readRDS), coi_data_files)

# import fapr results
fapr_results_files <- list.files(pattern = "^FAPR_RESULTS_SIM_DATA_haps_.*_ind_.*_max_change_.*\\.RDS")
FAPR_RESULTS_list <- setNames(lapply(fapr_results_files, readRDS), fapr_results_files)

# Remove columns 2 to 8 from all data frames in FAPR_RESULTS_list
FAPR_RESULTS_list <- lapply(FAPR_RESULTS_list, function(df) {
  df[, -c(2:8)] # Keep all columns except 2 to 8
})

# Order rows by HAPLO_FREQ_RECALC for each unique SampleID in FAPR_RESULTS_list
FAPR_RESULTS_list <- lapply(FAPR_RESULTS_list, function(df) {
  df[order(df$SampleID, -df$HAPLO_FREQ_RECALC), ]
})

# extract clean data to compare everything else with
clean_data <- FAPR_RESULTS_list$FAPR_RESULTS_SIM_DATA_haps_5_ind_1000_max_change_0.RDS

# Remove clean data from list
FAPR_RESULTS_list <- FAPR_RESULTS_list[!names(FAPR_RESULTS_list) %in% "FAPR_RESULTS_SIM_DATA_haps_5_ind_1000_max_change_0.RDS"]


# 2) Precision against noise. How does precision behaves when introducing noise to the data for each COI?
comparison_results_list <- list()
unique_samples <- unique(clean_data$SampleID)

# Loop through each element in FAPR_RESULTS_list
for (i in seq_along(FAPR_RESULTS_list)) {
  
  element_name <- names(FAPR_RESULTS_list)[i]
  cat("Processing:", element_name, "\n")
  
  # Initialize an empty data frame to store the results for the current FAPR_RESULTS_list element
  comparison_results <- data.frame(SampleID = character(),
                                   unique_haplotypes_clean = integer(),
                                   unique_haplotypes_FAPR = integer(),
                                   TP = integer(),
                                   FP = integer(),
                                   #FN = integer(),
                                   Precision = numeric(),
                                   FDR = numeric(),
                                   stringsAsFactors = FALSE)
  
  # Get the current FAPR_RESULTS_list element
  FAPR_results <- FAPR_RESULTS_list[[i]]
  
  # Get unique SampleIDs from the current FAPR_RESULTS_list element
  unique_fapr_samples <- unique(FAPR_results$SampleID)
  
  # Find common SampleIDs between clean_data and the current FAPR_RESULTS_list element
  common_samples <- intersect(unique_samples, unique_fapr_samples)
  
  # Loop through each common sample and calculate the desired counts
  for (sample in common_samples) {
    
    # Subset the data for the current sample
    clean_haplotypes <- unique(clean_data$haplotype[clean_data$SampleID == sample])
    FAPR_haplotypes <- unique(FAPR_results$haplotype[FAPR_results$SampleID == sample])
    
    # Count unique haplotypes
    unique_haplotypes_clean <- length(clean_haplotypes)
    unique_haplotypes_FAPR <- length(FAPR_haplotypes)
    
    # Count the number of intersecting haplotypes (True Positives)
    TP <- length(intersect(clean_haplotypes, FAPR_haplotypes))
    
    # Calculate False Positives and False Negatives
    FP <- length(setdiff(FAPR_haplotypes, clean_haplotypes))
    #FN <- length(setdiff(clean_haplotypes, FAPR_haplotypes))
    
    # Calculate Precision and False Discovery Rate (FDR)
    Precision <- ifelse((TP + FP) > 0, TP / (TP + FP), NA)
    FDR <- ifelse((TP + FP) > 0, FP / (TP + FP), NA)
    
    # Store the results in the data frame
    comparison_results <- rbind(comparison_results, 
                                data.frame(SampleID = sample,
                                           unique_haplotypes_clean = unique_haplotypes_clean,
                                           unique_haplotypes_FAPR = unique_haplotypes_FAPR,
                                           TP = TP,
                                           FP = FP,
                                           #FN = FN,
                                           Precision = Precision,
                                           FDR = FDR,
                                           stringsAsFactors = FALSE))
  }

  comparison_results_list[[element_name]] <- comparison_results
}

comparison_results_list

# Add max_change column and bind rows
comparison_results_list_with_max_change <- lapply(names(comparison_results_list), function(name) {
  df <- comparison_results_list[[name]]
  df$max_change <- name
  return(df)
})

# Combine all data frames into one
combined_comparison_results <- bind_rows(comparison_results_list_with_max_change)
combined_comparison_results$max_change <- gsub(".*_max_change_", "", gsub(".RDS$", "", combined_comparison_results$max_change))

# Convert unique_haplotypes_clean to factor and sort it alphabetically
combined_comparison_results$unique_haplotypes_clean <- factor(combined_comparison_results$unique_haplotypes_clean, levels = sort(unique(combined_comparison_results$unique_haplotypes_clean)))

# Plot
eval_adj_relabun <- ggplot(combined_comparison_results, aes(x = unique_haplotypes_clean, y = Precision, color = max_change, fill = max_change)) +
  geom_boxplot(alpha = 0.50, outlier.shape = NA) + 
  geom_jitter(width = 0.25, size = 1, alpha = 0.15) + 
  facet_wrap(~ max_change, scales = "free_x", nrow = 2) + 
  labs(
    x = "Expected Haplotypes",
    y = "Precision",
    title = "Precision by Expected Haplotypes",
    subtitle = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom")

eval_adj_relabun

ggsave("eval_adj_relabun.png", eval_adj_relabun, dpi = 300, width = 14, height = 10, bg = "white")



# 3) check evenness from allele data to see the correlation between precision and evenness. is precision given by the evennes of the individual resmarker freqs?

# Function to calculate various evenness metrics
calculate_evenness_metrics <- function(freqs) {
  # Number of alleles (S)
  S <- length(freqs)
  
  # Shannon-Wiener Index (H')
  H <- -sum(freqs * log(freqs))
  
  # Pielou's Evenness Index (J')
  J <- H / log(S)
  
  # Simpson's Diversity Index (D)
  D <- sum(freqs^2)
  
  # Simpson's Evenness Index (E1/D)
  E1_D <- (1 / D) / S
  
  # Shannon's Equitability Index (EH)
  EH <- H / log(S)
  
  # Berger-Parker Index
  p_max <- max(freqs)
  Berger_Parker <- 1 / p_max
  
  # Smith & Wilson’s Evenness Index (Evar)
  Evar <- 1 - (2 / pi) * atan(var(freqs) / mean(freqs))
  
  # Return a list of evenness metrics
  return(list(J = J, E1_D = E1_D, EH = EH, Berger_Parker = Berger_Parker, Evar = Evar))
}

# Apply this to each element in SIM_DATA_list
evenness_results_list <- lapply(SIM_DATA_list, function(df) {
  df %>%
    group_by(SampleID, resmarker) %>%
    summarise(
      Pielou_J = calculate_evenness_metrics(norm.reads.locus)$J,
      Simpson_E1_D = calculate_evenness_metrics(norm.reads.locus)$E1_D,
      Shannon_EH = calculate_evenness_metrics(norm.reads.locus)$EH,
      Berger_Parker = calculate_evenness_metrics(norm.reads.locus)$Berger_Parker,
      Evar = calculate_evenness_metrics(norm.reads.locus)$Evar
    )
})

# Add max_change column and bind rows
evenness_results_list <- lapply(names(evenness_results_list), function(name) {
  df <- evenness_results_list[[name]]
  df$max_change <- name
  return(df)
})

# Combine all data frames into one
evenness_results_list <- bind_rows(evenness_results_list)
evenness_results_list$max_change <- gsub(".*_max_change_", "", gsub(".RDS$", "", evenness_results_list$max_change))

#calculate mean metrics for each sample (mean from non NA results for each resmarker)
evenness_means <- evenness_results_list %>%
  group_by(SampleID, max_change) %>%
  summarise(
    mean_Pielou_J = mean(Pielou_J[is.finite(Pielou_J)], na.rm = TRUE),
    mean_Simpson_E1_D = mean(Simpson_E1_D[is.finite(Simpson_E1_D)], na.rm = TRUE),
    mean_Shannon_EH = mean(Shannon_EH[is.finite(Shannon_EH)], na.rm = TRUE),
    mean_Berger_Parker = mean(Berger_Parker[is.finite(Berger_Parker)], na.rm = TRUE),
    mean_Evar = mean(Evar[is.finite(Evar)], na.rm = TRUE),
    max_change = first(max_change)
  )

sapply(evenness_means[,c(-1,-2)], function(x) sd(x, na.rm = TRUE))

# Calculate the correlation matrix
cor_matrix <- cor(evenness_means[,c(-1,-2)])

corrplot(cor_matrix, method = "circle", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         addCoef.col = "black", 
         order = "hclust") 

#remove pielou a(same as shannon EH apparently) and Evar
evenness_means <- evenness_means[!colnames(evenness_means) %in% c("mean_Pielou_J", "mean_Evar")]

cor_matrix <- cor(evenness_means[,c(-1,-2)])

corrplot(cor_matrix, method = "circle", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         addCoef.col = "black", 
         order = "hclust") 

# merge with precision

evenness_precision <- merge(evenness_means, combined_comparison_results[c("SampleID", "Precision", "unique_haplotypes_clean", "max_change")], by = c("SampleID", "max_change"))

#mean_Shannon_EH looks good
ggplot(evenness_precision, aes(x = mean_Shannon_EH, y = Precision, color = max_change))+
  #geom_boxplot(alpha = 0.50, outlier.shape = NA) + 
  geom_jitter(width = 0.1, size = 0.5, alpha = 0.15) + 
  facet_wrap(~ unique_haplotypes_clean, scales = "free_x", nrow = 2) + 
  geom_smooth(method = "lm", se = T)+
  theme_minimal()

##

















############# EVALUATION OF noisy o clean data aquí  (PROBLEAMS CON NOISY. CLEAN FUNCIONA. CHECAR) # 
# SIM_DATA <- readRDS(paste0("SIM_DATA_noisy_", "haps_", max_haplos,  "_ind_", individuals, ".RDS")) # 

#SIM_DATA <- as.data.frame(SIM_DATA)

#unique_samples <- unique(SIM_DATA$SampleID)

RESULTS_BENCH_ALL <- data.frame()
RESULTS_BENCH_ALL_FREQS <- list()

for (sample in unique_samples){
  
  print(sample)
  
  SIM_DATA_subset <- SIM_DATA[SIM_DATA$SampleID == sample,]
  SIM_DATA_subset <- SIM_DATA_subset[SIM_DATA_subset$norm.reads.locus != 1,] #remove monoallelic loci
  
  # Count the number of unique AA values for each resmarker
  resmarker_counts <- SIM_DATA_subset %>%
    group_by(resmarker) %>%
    summarize(unique_AA_count = n_distinct(AA))
  
  # Find the maximum unique AA count
  n_haplos <- max(resmarker_counts$unique_AA_count)
  
  # Filter for rows where the unique AA count is equal to the maximum
  SIM_DATA_subset <- SIM_DATA_subset %>%
    group_by(resmarker) %>%
    filter(n_distinct(AA) == n_haplos)
  
  
  RESULTS_FINAL_subet <- RESULTS_FINAL[RESULTS_FINAL$SampleID == sample,]
  
  # Extract observed and expected values
  observed <- RESULTS_FINAL_subet$HAPLO_FREQ_RECALC  
  #expected <- rev(sort(unique(SIM_DATA_subset$norm.reads.locus)))
  
  # Calculate the means in windows of size n_haplos and use as expected haplo freq
  by_ <- nrow(SIM_DATA_subset)/n_haplos
  
  if (by_ == 1){ # this means the sample is monoallelic for all loci except 1, so no need to average out anything
    
    expected <-  rev(sort(SIM_DATA_subset$norm.reads.locus))
    
  }else{ # if the sample is polyallelic for multiple loci, average freq is needed
    
    expected <- rollapply(rev(sort(SIM_DATA_subset$norm.reads.locus)), width = 1, by = by_, FUN = mean, align = "left", partial = TRUE)

  }

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
  
  # Combine the RMSE, diff_haplos, and freqs_list into a data frame
  RESULTS_BENCH <- data.frame(RMSE = rmse, diff_haplos = difference, n_haplos = paste0(n_haplos, "_haplos"))
  RESULTS_BENCH_ALL <- rbind(RESULTS_BENCH_ALL, RESULTS_BENCH)
}


#checkpoint
# saveRDS(RESULTS_BENCH_ALL_FREQS,  paste0("RESULTS_BENCH_ALL_FREQS_", "haps_", max_haplos,  "_ind_", individuals, ".RDS"))
# saveRDS(RESULTS_BENCH_ALL, paste0("RESULTS_BENCH_ALL_", "haps_", max_haplos,  "_ind_", individuals, ".RDS"))


######----------------------------------------------------------------------------------------


#### VISUALIZATION

# RESULTS_BENCH_ALL_FREQS <- readRDS(paste0("RESULTS_BENCH_ALL_FREQS_", "haps_", max_haplos,  "_ind_", individuals, ".RDS"))
# RESULTS_BENCH_ALL <- readRDS(paste0("RESULTS_BENCH_ALL_", "haps_", max_haplos,  "_ind_", individuals, ".RDS"))
# 

RESULTS_BENCH_ALL$Individual <- rownames(RESULTS_BENCH_ALL)


# DESCRIBE THE DATASET: EXPECTED_VALUES FOR EACH N HAPLOTYPES

chunk_size <- individuals # Define the chunk size (it's the same as individuals simulated)
num_chunks <- ceiling(length(RESULTS_BENCH_ALL_FREQS) / chunk_size)

# #shannon diversity as "evennes" index
# shannon_diversity <- function(x) {
#   -sum(x * log(x))
# }

# Loop through each chunk
for (i in 1:num_chunks) {
  start_index <- (i - 1) * chunk_size + 1
  end_index <- min(i * chunk_size, length(RESULTS_BENCH_ALL_FREQS))
  
  chunk_data <- RESULTS_BENCH_ALL_FREQS[start_index:end_index]
  
  names(chunk_data)<- c(start_index:end_index)
  
  # DESCRIBE THE DATASET: EXPECTED_VALUES FOR EACH N HAPLOTYPES
  expected_values <- lapply(chunk_data, function(x) x$expected[x$expected != 0])
  expected_values <- as.data.frame(do.call(rbind, expected_values))
  
  data_for_histogram_long <- reshape2::melt(expected_values)
  
  a <- ggplot(data_for_histogram_long, aes(x = value, color = variable)) +
    geom_density(alpha = 0.5, linewidth = 2) +
    labs(x = "Expected Frequency", y = "Density", title = "Plots of Expected Values from Simulated Data") +
    theme_minimal() +
    guides(color = FALSE)+
    theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 14),
          axis.text.y = element_text(size = 20),
          axis.title = element_text(size = 20),
          legend.text = element_text(size = 20))
  
  expected_values <- melt(t(expected_values))
  colnames(expected_values) <- c("Haplo", "Individual", "Freq")
  
  # Merge the two data frames by Individual column
  expected_values <- merge(RESULTS_BENCH_ALL, expected_values, by = c("Individual"), all.y = TRUE)
  
  # Identify levels of Individual for Haplo == "V1" and sort them
  levels_v1 <- unique(expected_values$Individual[expected_values$Haplo == "V1"])
  levels_v1 <- levels_v1[order(expected_values$Freq[expected_values$Haplo == "V1"], decreasing = FALSE)]
  
  # Reorder levels of Individual based on Freq for Haplo == "V1"
  expected_values$Individual <- factor(expected_values$Individual, levels = levels_v1)
  
  #shannon reorder
  # ordered_individual <- expected_values %>%
  #   group_by(Individual) %>%
  #   summarise(shannon = shannon_diversity(Freq)) %>%
  #   arrange(shannon) %>%
  #   pull(Individual)
  
  #min haplo freq reorder
  # Calculate the minimum Freq for each Individual
  min_freq <- expected_values %>%
    group_by(Individual) %>%
    summarise(min_freq = min(Freq))
  
  ordered_individual <- min_freq %>%
    arrange(min_freq) %>%
    pull(Individual)
  
  
  # Convert to factor with levels in the order they appear
  ordered_individual <- factor(ordered_individual, levels = unique(expected_values$Individual))
  
  expected_values$Individual <- factor(expected_values$Individual, levels = ordered_individual)
  
  # b <- ggplot(expected_values, aes(fill = as.factor(Haplo), y = Freq, x = Individual)) +
  #   geom_bar(position = "stack", stat = "identity") +
  #   labs(x = "Individual", y = "Expected Frequency", title = "") +
  #   #guides(fill = FALSE) +
  #   theme(axis.text.x = element_blank())
  
  b <- ggplot(expected_values, aes(fill = as.factor(-diff_haplos), y = Freq, x = Individual)) +
    geom_bar(position = "stack", stat = "identity", color = "black") +
    labs(x = "Individual", y = "Expected Frequency", title = "") +
    scale_fill_discrete(name = "Extra haplotypes")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          axis.text.y = element_text(size = 20),
          axis.title = element_text(size = 20),
          legend.text = element_text(size = 20))
  
  ggsave(paste0("haplos_", i+1, "_plots.png"), arrangeGrob(a, b, ncol = 1, nrow=2), dpi = 300, height= 30, width = 40)
}


# BENCHMARK PLOTS YOOO

create_benchmark_plots <- function(RESULTS_BENCH_ALL, bench_name = "") {
  # Calculate stacked input
  stacked_input <- RESULTS_BENCH_ALL %>% 
    group_by(n_haplos, diff_haplos) %>%
    summarize(counts = as.vector(table(diff_haplos)))
  
  # Define a lookup table
  lookup_table <- data.frame(diff_haplos = c("0", "1", "2", "3", "-1", "-2", "-3", "-4"),
                             acc = c("exact_haplos", "missing_1_haplo", "missing_2_haplos", "missing_3_haplos",
                                     "extra_1_haplo", "extra_2_haplos", "extra_3_haplos", "extra_4_haplos"))
  
  # Merge the lookup table with stacked_input based on diff_haplos
  stacked_input <- merge(stacked_input, lookup_table, by = "diff_haplos", all.x = TRUE)
  
  stacked_input$acc <- factor(stacked_input$acc,
                              levels = c("extra_4_haplos", "extra_3_haplos", "extra_2_haplos", "extra_1_haplo",
                                         "exact_haplos", "missing_1_haplo", "missing_2_haplos", "missing_3_haplos"))
  
  # Define custom colors
  custom_colors <- c("exact_haplos" = "limegreen", "missing_1_haplo" = "orange", 
                     "missing_2_haplos" = "red2", "missing_3_haplos" = "black", 
                     "extra_1_haplo" = "pink", "extra_2_haplos" = "violet", 
                     "extra_3_haplos" = "purple", "extra_4_haplos" = "blue")
  
  # Create the stacked barplot
  perc_plot <- ggplot(stacked_input, aes(x = n_haplos, y = counts, fill = acc)) +
    geom_bar(stat = "identity") +
    labs(title = "", x = "", y = "Samples") +
    scale_fill_manual(values = custom_colors) +  # Set custom colors
    theme_minimal()
  
  ggsave(paste0(bench_name,"__haplo_discovery_benchmark.png"), perc_plot, bg = "white", height = 6, width = 10)
  
  
  max_rmse <- max(RESULTS_BENCH_ALL$RMSE)
  
  rmse_plot <- ggplot(RESULTS_BENCH_ALL, aes(x = n_haplos, y = RMSE, color = n_haplos)) +
    geom_boxplot() +
    geom_quasirandom(alpha = 0.5)+
    labs(title = "", x = "", y = "RMSE for true positives") +
    theme_minimal() +
    ylim(0, max_rmse)+
    guides(color = "none")
  
  ggsave(paste0(bench_name,"__haplo_freqs_benchmark.png"), rmse_plot, bg = "white", height = 6, width = 10)
}


create_benchmark_plots(RESULTS_BENCH_ALL, "testing")


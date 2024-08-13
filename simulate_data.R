
### SIMULATE HAPLOTYPE DATA TO PHASE 

library(dplyr)

# SET PRIORS FROM EDA -----

max_haplos <- 5 
max_amps <- 4 

individuals <- 5000 # for each amount of haplos

# placeholder names
res_loci <- paste0("HAP", seq(1, max_amps))
res_variants <- c(LETTERS[1:max_haplos])

# Define noise
max_change_values <- seq(0, 0.4, by = 0.05)
min_change <- 0


# GENERATE HAPLOTYPE PROPORTIONS ------

set.seed(42069)

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


# GENERATE DATA -----

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
    genos <- genos[genos$locus != last(res_loci),]
    HAP7 <- as.data.frame(cbind(locus = last(res_loci), genotype = res_variants[1:n_haplos]))
    genos<- rbind(genos, HAP7)
    
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
    
    colnames(genos) <- c("resmarker", "Microhaplotype", "norm.reads.locus")
    
    genos <- cbind(SampleID = paste0("shit_sample_", i,"_", n_haplos), genos)
    
    resmarkers_table <- genos
    
    # Append to SIM_DATA list
    SIM_DATA[[paste0("shit_sample_", i, "_", n_haplos)]] <- resmarkers_table
    
  }
}

#convert to a single df
SIM_DATA <- bind_rows(SIM_DATA)

#check data!

#min should not be < ~0.01
summary(SIM_DATA$norm.reads.locus)

check_SIM_DATA <- SIM_DATA %>%
  group_by(SampleID, resmarker) %>%
  summarise(total_freq = sum(norm.reads.locus))

# all markers there? must be = to max_amps
if (length(unique(check_SIM_DATA$resmarker)) == max_amps){
  print("Simulated amps: good")
} else {
  print("well shit.")
}

# all freqs = 1?
if (sum(check_SIM_DATA$total_freq) == length(unique(check_SIM_DATA$SampleID)) * length(unique(check_SIM_DATA$resmarker))){
  print("Simulated freqs: good")
} else {
  print("well shit.")
}


# CREATE COI DATA: this is gonna be an optional input for FapR
unique_sample_ids <- unique(SIM_DATA$SampleID)
coi_values <- sapply(unique_sample_ids, function(id) {
  substr(id, nchar(id), nchar(id))
})
COI_SIM_DATA <- data.frame(SampleID = unique_sample_ids, COI = as.numeric(coi_values), row.names = NULL)


#save data
saveRDS(SIM_DATA, paste0("SIM_DATA_clean_", "haps_", max_haplos,  "_ind_", individuals, ".RDS"))
saveRDS(COI_SIM_DATA, paste0("COI_SIM_DATA_clean", "haps_", max_haplos,  "_ind_", individuals, ".RDS"))


# ADD NOISE TO SIMULATED DATA ------

# adding noise means that freqs for each resmarker vary a bit; for instance, it may not be 55-45 for all in a 2 haps sample, but rather dhfr_59: 57-43, dhfr_437: 50-50, dhps_540: 52-48
# This will increase the errors. How will FAPR perform?

introduce_noise <- function(data, max_change = 0, min_change = 0) { # Function to introduce noise to norm.reads.locus values
  
  # Group data by SampleID and resmarker
  multi_aa_data <- data %>%
    group_by(SampleID, resmarker) %>%
    filter(n_distinct(Microhaplotype) > 1) %>%  # Filter for rows with more than one Microhaplotype
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
  
  # Filter rows where there is only one Microhaplotype and combine with the noisy data
  noisy_data <- bind_rows(
    data %>% group_by(SampleID, resmarker) %>% filter(n_distinct(Microhaplotype) == 1),
    multi_aa_data
  )
  
  # Remove the random_factor column
  noisy_data <- as.data.frame(noisy_data[,-5])
  
  return(noisy_data)
}


# Introduce noise
for (max_change in max_change_values) {
  
  # Apply the introduce_noise function to SIM_DATA
  SIM_DATA_noisy <- introduce_noise(SIM_DATA, max_change = max_change, min_change = min_change)
  
  # Get the corresponding indices of rows in noisy_SIM_DATA based on SampleID, resmarker, and Microhaplotype
  indices <- match(apply(SIM_DATA[c("SampleID", "resmarker", "Microhaplotype")], 1, paste, collapse = "_"), 
                   apply(SIM_DATA_noisy[c("SampleID", "resmarker", "Microhaplotype")], 1, paste, collapse = "_"))
  
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


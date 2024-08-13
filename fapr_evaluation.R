
library(gridExtra)
library(ggbeeswarm)
library(ggplot2)
library(dplyr)
library(reshape2)
library(zoo)
library(corrplot)
library(broom)
library(boot)

source("calculate_diversity_metrics.R")
source("FAPR.R")

# IMPORT SIM DATA ------
sim_data_files <- list.files(pattern = "^SIM_DATA_haps_.*_ind_.*_max_change_.*\\.RDS")
SIM_DATA_list <- setNames(lapply(sim_data_files, readRDS), sim_data_files)
coi_data_files <- list.files(pattern = "^SIM_DATA_COI_haps_.*_ind_.*_max_change_.*\\.RDS")
COI_SIM_DATA_list <- setNames(lapply(coi_data_files, readRDS), coi_data_files)

print(paste("Loaded", length(SIM_DATA_list), "SIM_DATA files and", length(COI_SIM_DATA_list), "COI_SIM_DATA files into named lists."))


# RUN FAPR ON SIM DATA -----
for (data in 1:length(SIM_DATA_list)){
  
  data_name <- names(SIM_DATA_list)[data]
  print(paste0("------ Processing Dataset: ", data_name, " -------"))
  
  # 1) subset data
  SIM_DATA <- SIM_DATA_list[[data]]
  COI_SIM_DATA <- COI_SIM_DATA_list[[data]]
  
  # 2) run fapr
  RESULTS_FAPR <- FAPR(SIM_DATA, COI_SIM_DATA, verbose = F)
  
  # 3) save results
  saveRDS(RESULTS_FAPR, paste0("FAPR_RESULTS_", data_name))

}


# EVALUATION 1:  ARRANGE DATA -----------

# Get all SIM_DATA and COI_SIM_DATA files
sim_data_files <- list.files(pattern = "^SIM_DATA_haps_.*_ind_.*_max_change_.*\\.RDS")
SIM_DATA_list <- setNames(lapply(sim_data_files, readRDS), sim_data_files)
coi_data_files <- list.files(pattern = "^SIM_DATA_COI_haps_.*_ind_.*_max_change_.*\\.RDS")
COI_SIM_DATA_list <- setNames(lapply(coi_data_files, readRDS), coi_data_files)

# import fapr results
fapr_results_files <- list.files(pattern = "^FAPR_RESULTS_SIM_DATA_haps_.*_ind_.*_max_change_.*\\.RDS")
FAPR_RESULTS_list <- setNames(lapply(fapr_results_files, readRDS), fapr_results_files)

# Remove individual hap variants from all data frames in FAPR_RESULTS_list
FAPR_RESULTS_list <- lapply(FAPR_RESULTS_list, function(df) {
  df[, c(1, (ncol(df)-2):ncol(df))]
})

# Order rows by HAPLO_FREQ_RECALC for each unique SampleID in FAPR_RESULTS_list
FAPR_RESULTS_list <- lapply(FAPR_RESULTS_list, function(df) {
  df[order(df$SampleID, -df$HAPLO_FREQ_RECALC), ]
})

#extract clean data filename
clean_data_filename <- names(FAPR_RESULTS_list)[grepl("max_change_0.RDS", names(FAPR_RESULTS_list))]

# extract clean data to compare everything else with
clean_data <- FAPR_RESULTS_list[[clean_data_filename]]

# Remove clean data from list
FAPR_RESULTS_list <- FAPR_RESULTS_list[!names(FAPR_RESULTS_list) %in% clean_data_filename]


# EVALUATION 2:  CALCULATE PRECISION -----------

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
                                           Precision = Precision,
                                           FDR = FDR,
                                           stringsAsFactors = FALSE))
  }

  comparison_results_list[[element_name]] <- comparison_results
}


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
  geom_jitter(width = 0.25, size = 1, alpha = 0.02) + 
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


# same plot as above but simpler and w/ 95% CI

# Define a function to calculate the median of a sample
median_ci <- function(data, indices) {
  sample <- data[indices]
  return(median(sample))
}

# Function to calculate median and CI for a given subset of data
compute_median_ci <- function(data) {
  # Perform bootstrapping
  results <- boot(data$Precision, statistic = median_ci, R = 10000)
  # Extract the 95% confidence intervals
  ci <- boot.ci(results, type = "perc")
  # Return a tibble with median and CI
  tibble(
    median_precision = median(data$Precision, na.rm = TRUE),
    ci_lower = ci$percent[4],
    ci_upper = ci$percent[5]
  )
}

# Calculate median precision and 95% CI for each max_change of each SampleID
median_precision_ci <- combined_comparison_results %>%
  group_by(unique_haplotypes_clean, max_change) %>%
  do(compute_median_ci(.))

# Customize facet labels with "MOI = " and "noise = " prefixes
facet_labels_median <- labeller(
  unique_haplotypes_clean = function(x) paste("MOI = ", x)
)

# Create the plot
median_precision_plot <- ggplot(median_precision_ci, aes(x = max_change, y = median_precision, color = unique_haplotypes_clean)) +
  geom_point(size = 3) +  # Plot the median precision as points
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, size = 1) +  # Add error bars for CI
  facet_wrap(~ unique_haplotypes_clean, scales = "free_x", nrow = 2, labeller = facet_labels_median) +  # Facet by unique_haplotypes_clean
  theme_minimal(base_size = 12) +
  ylim(0,1)+
  labs(
    x = "Noise",
    y = "Median Precision",
    title = "Median Precision and 95% CI"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 10, face = "bold", hjust = 0),  
    strip.placement = "top",  
    panel.background = element_rect(fill = "white", color = "black"), 
    panel.grid.major = element_line(color = "grey90"), 
    panel.grid.minor = element_line(color = "grey90"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none"
  )

median_precision_plot

ggsave("median_CI_precision_plot.png", median_precision_plot, dpi = 300, width = 8, height = 6, bg = "white")



# EVALUATION 3:  ANALYSE EVENNESS -----------

# 3) check evenness from allele data to see the correlation between precision and evenness. is precision given by the evennes of the individual resmarker freqs?

# Apply this to each element in SIM_DATA_list
evenness_results_list <- lapply(SIM_DATA_list, function(df) {
  df %>%
    group_by(SampleID, resmarker) %>%
    summarise(
      Pielou_J = calculate_evenness_metrics(norm.reads.locus)$J,
      Simpson_E1_D = calculate_evenness_metrics(norm.reads.locus)$E1_D,
      Shannon_EH = calculate_evenness_metrics(norm.reads.locus)$EH,
      Berger_Parker = calculate_evenness_metrics(norm.reads.locus)$Berger_Parker,
      Evar = calculate_evenness_metrics(norm.reads.locus)$Evar,
      highest_freq = max(norm.reads.locus)
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
    mean_max_freq = mean(highest_freq[highest_freq != 1], na.rm = TRUE), # Exclude values of 1
    max_change = first(max_change)
  )

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


# EVALUATION 4:  EVENNESS VS PRECISION -----------

# merge with precision
evenness_precision <- merge(evenness_means, combined_comparison_results[c("SampleID", "unique_haplotypes_clean", "Precision", "max_change")], by = c("SampleID", "max_change"))


# Customize facet labels with "MOI = " and "noise = " prefixes
facet_labels <- labeller(
  unique_haplotypes_clean = function(x) paste("MOI = ", x),
  max_change = function(x) paste("noise = ", x)
)

# Create the plot
ev_pres_plot_shannon <- ggplot(evenness_precision, aes(x = mean_Shannon_EH, y = Precision, color = max_change)) +
  geom_jitter(width = 0.05, height = 0.05, size = 1, alpha = 0.05) + 
  geom_smooth(method = "lm", color = "grey60") +
  facet_wrap(~ unique_haplotypes_clean + max_change, scales = "free_x", nrow = 4, labeller = facet_labels) +
  xlim(0, NA) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 10, face = "bold", hjust = 0),  
    strip.placement = "top",  
    panel.background = element_rect(fill = "white", color = "black"), 
    panel.grid.major = element_line(color = "grey90"), 
    panel.grid.minor = element_line(color = "grey90"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none"
  ) +
  labs(
    x = "Mean Shannon's Equitability Index (EH)",
    y = "Precision",
    title = ""
  )

ev_pres_plot_shannon

ggsave("evenness_presision_shannon_EH_plot.png", ev_pres_plot_shannon, dpi = 300, width = 20, height = 14, bg = "white")


# Create the plot
ev_pres_plot_max_Freq <- ggplot(evenness_precision, aes(x = mean_max_freq, y = Precision, color = max_change)) +
  geom_jitter(width = 0.05, height = 0.05, size = 1, alpha = 0.05) + 
  geom_smooth(method = "lm", color = "grey60") +
  facet_wrap(~ unique_haplotypes_clean + max_change, scales = "free_x", nrow = 4, labeller = facet_labels) +
  xlim(0, NA) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 10, face = "bold", hjust = 0),  
    strip.placement = "top",  
    panel.background = element_rect(fill = "white", color = "black"), 
    panel.grid.major = element_line(color = "grey90"), 
    panel.grid.minor = element_line(color = "grey90"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none"
  ) +
  labs(
    x = "Mean Major Allele Frequency",
    y = "Precision",
    title = ""
  )

ev_pres_plot_max_Freq

ggsave("evenness_presicion_max_freq_plot.png", ev_pres_plot_max_Freq, dpi = 300, width = 20, height = 14, bg = "white")

##

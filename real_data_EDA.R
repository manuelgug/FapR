
# Exploratory Data Analysis (EDA)

# Objective: Know your priors from real mad4hatter ampseq data and then simulate tons of data based on them.

library(dplyr)
library(ggplot2)
library(progress)
library(fs)
library(stringr)
library(boot)
library(tidyr)

source("calculate_diversity_metrics.R")
source("FAPR.R")

## IMPORT AND CLEAN REAL DATA ---------------

directory_path <- "../FILTERED_DATA_FROM_CLUSTER/"
data_all <- data.frame()

pb <- progress_bar$new(
  format = "[:bar] :percent ETA: :eta",
  total = length(dir_ls(path = directory_path, regexp = "_FILTERED", ignore.case = TRUE))
)

for (folder_path in dir_ls(path = directory_path, regexp = "_FILTERED", ignore.case = TRUE)) {
  pb$tick() 
  
  folder_name <- path_file(folder_path)
  file_path <- file.path(folder_path, "resmarker_microhap_table_global_max_0_filtered.csv")
  
  cat("\n")
  print(folder_name)
  
  if (file.exists(file_path)) {
    data <- read.csv(file_path)
    data$run <- sub("^.*/([^/]+)/[^/]+$", "\\1", file_path)
    data_all <- rbind(data_all, data)
  }
}

data_all$run <- sub("_RESULTS_v0.1.8_FILTERED$", "", data_all$run)

# Create resmarker column
data_all$resmarker <- paste(data_all$Gene, data_all$MicrohapIndex, sep = "_")

# Keep only resmarkers of interest
data_all <- data_all[data_all$resmarker %in% c("dhps_431/436/437", "dhps_540/581", "dhfr_16/51/59", "dhfr_108/164"),]

# remove controls
controls <- c("dd2", "pm", "3d7", "hb3", "neg", "undetermined")

data_all <- data_all %>%
  filter(!str_detect(tolower(SampleID), paste(controls, collapse = "|")))

# keep samples with all haplos. if not, fapr collapses
total_haplos <- length(unique(data_all$resmarker))

data_all <- data_all %>% 
  group_by(SampleID) %>%
  filter(length(unique(resmarker)) == total_haplos) %>%
  ungroup()


# CLEAN DATA
# Count the number of unique SampleID for each unique run
run_sample_counts <- data_all %>%
  group_by(run) %>%
  summarise(unique_sample_count = n_distinct(SampleID)) %>%
  arrange(desc(unique_sample_count))

run_sample_counts %>%
  mutate(run_info = paste("Run:", run, "|| Samples:", unique_sample_count)) %>%
  select(run_info) %>%
  print(n=100)

print(paste("Total samples:", length(unique(data_all$SampleID))))



#### Q1) HOW MANY UNIQUE HAPLOTYPES ARE THERE PER AMPLICON OF INTEREST? -----

data_all %>%
  group_by(resmarker) %>%
  summarise(max(n.alleles))

# RESULT: from 3 to 5. 5 is the max. 5 would be a good upper max for sim data


#### Q2) IS THERE A CORRELATION BETWEEN EVENNESS AND MAX ALLELE COUNT? -----

# calculate evenness metrics
evenness_real_data <- data_all %>%
  group_by(SampleID, resmarker) %>%
  summarise(
    Pielou_J = calculate_evenness_metrics(norm.reads.locus)$J,
    Simpson_E1_D = calculate_evenness_metrics(norm.reads.locus)$E1_D,
    Shannon_EH = calculate_evenness_metrics(norm.reads.locus)$EH,
    Berger_Parker = calculate_evenness_metrics(norm.reads.locus)$Berger_Parker,
    Evar = calculate_evenness_metrics(norm.reads.locus)$Evar,
    highest_freq = max(norm.reads.locus),
    max_alleles = max(n.alleles)
  )

#evenness across haplos
evenness_real_data_means <- evenness_real_data %>%
  group_by(SampleID) %>%
  summarise(
    mean_Pielou_J = mean(Pielou_J[is.finite(Pielou_J)], na.rm = TRUE),
    mean_Simpson_E1_D = mean(Simpson_E1_D[is.finite(Simpson_E1_D)], na.rm = TRUE),
    mean_Shannon_EH = mean(Shannon_EH[is.finite(Shannon_EH)], na.rm = TRUE),
    mean_Berger_Parker = mean(Berger_Parker[is.finite(Berger_Parker)], na.rm = TRUE),
    mean_Evar = mean(Evar[is.finite(Evar)], na.rm = TRUE),
    mean_max_freq = mean(highest_freq[highest_freq != 1], na.rm = TRUE), # Exclude values of 1
    max_alleles = max(max_alleles)
  )

#remove monoallelic samples
evenness_real_data_means<- evenness_real_data_means[!evenness_real_data_means$max_alleles == 1,]

# viz evenness and max allele count
hist(evenness_real_data_means$mean_Shannon_EH)
hist(evenness_real_data_means$max_alleles, breaks = 5)

# Perform the correlation test
cor_test <- cor.test(evenness_real_data_means$max_alleles, evenness_real_data_means$mean_Shannon_EH, method = "spearman", exact = F)

# Extract the correlation coefficient and p-value
cor_coef <- round(cor_test$estimate, 2)
p_val <- formatC(cor_test$p.value, format = "e", digits = 2)

ggplot(evenness_real_data_means, aes(x = max_alleles, y = mean_Shannon_EH)) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.3, color = "gray40", size = 4) +
  theme_minimal() +
  geom_smooth(method = "lm", color = "red", fill = "pink2") +
  annotate("text", x = Inf, y = -Inf, label = paste("r =", cor_coef, "\n", "p =", p_val),
           hjust = 1.1, vjust = -1.1, size = 4, color = "black")

# RESULT: the more more alleles, the more even the sample (allele freqs for a given resmarker are more even), although > 3 alleles are relatively rare


#### Q3) HOW ARE SEQ ERRORS / AMP YIELD IN IN-SAMPLE FREQS BETWEEN SAME ALLELE AMPLICONS ACROSS SAMPLES?-------------

# calculate error found between alleles from amplicons with the same amount of alleles. can i make an educated guess on what to expect?

#subset samples with > 1 multiallelic loci
many_multiallelic_loci_samples <- data_all %>%
  group_by(SampleID) %>%
  filter(sum(n.alleles > 1) > 2) %>%
  arrange(SampleID, resmarker, desc(norm.reads.locus)) %>%  
  ungroup()

#how much does the amplicons of interest vary in multiallelic samples?
mean_alleles <- many_multiallelic_loci_samples %>%
  group_by(resmarker) %>%
  summarise(mean_alleles = mean(n.alleles),
            sd = sd(n.alleles))

ggplot(mean_alleles, aes(x = resmarker, y = mean_alleles, fill = resmarker)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean_alleles - sd, ymax = mean_alleles + sd), width = 0.2) + 
  theme_minimal() +
  labs(
    title = "",
    x = "Resmarker Haplo",
    y = "Mean Alleles in Multiallelic Samples"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  
    legend.position = "none" 
  )

#remove monoallelic loci from each sample
many_multiallelic_loci_samples_nomono <- many_multiallelic_loci_samples[many_multiallelic_loci_samples$norm.reads.locus != 1,]

# rank alleles from major to minor
many_multiallelic_loci_samples_nomono <- many_multiallelic_loci_samples_nomono %>%
  group_by(SampleID, resmarker) %>%
  mutate(rank = row_number()) %>%
  ungroup() 

#calculate error: deviation from the mean freq
many_multiallelic_loci_samples_nomono <- many_multiallelic_loci_samples_nomono %>%
  group_by(SampleID, n.alleles, rank) %>%
  mutate(dev_from_mean = abs(mean(norm.reads.locus)- norm.reads.locus))

#format for plot
error_stats <- many_multiallelic_loci_samples_nomono %>%
  group_by(SampleID) %>%
  summarise(resmarkers = paste(unique(resmarker), collapse = "___"),
            dev_from_mean = first(dev_from_mean),
            n.alleles = (first(n.alleles))) %>%
  ungroup()

#remove resmarkers with no comparisons or zero error
error_stats <- error_stats[grepl("___", error_stats$resmarkers),]
error_stats <- error_stats[error_stats$dev_from_mean != 0,]


ggplot(error_stats, aes(x = as.factor(n.alleles), y = dev_from_mean, color = resmarkers, fill = resmarkers))+
  geom_boxplot(outliers = F, alpha = 0.35)+
  geom_jitter(alpha = 0.5, size = 3, width = 0.1)+
  facet_wrap(~resmarkers)+
  theme_minimal()+
  labs(
    title = "Error by amps with the same number of alleles",
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  
    strip.text = element_text(size = 10, face = "bold", hjust = 0),  
    strip.placement = "top",  
    panel.background = element_rect(fill = "white", color = "black"), 
    panel.grid.major = element_line(color = "grey90"), 
    panel.grid.minor = element_line(color = "grey90"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none"
  )

#RESULT: biallelic samples vary mostly below 0.05, with a max at about 0.25 when 3 amps have 2 alleles. there was only one instance of a sample with 3 alleles for more than 1 amplicon and it varied 0.2
# i guess a reasonable range of noise to add to simulated data would be from 0.05 to 0.4 just to stretch it a little


## FAPR TEST --------

#run fapr
r_real_data <- FAPR(data_all, verbose = F)

#order gene columns by locus
exclude_cols <- c("SampleID", "HAPLO_FREQ", "HAPLO_FREQ_RECALC")

cols_to_reorder <- setdiff(colnames(r_real_data), exclude_cols)

ordered_colnames <- cols_to_reorder[order(sapply(cols_to_reorder, function(x) {
  num <- as.numeric(unlist(regmatches(x, gregexpr("\\d+", x))))
  if (length(num) > 0) min(num) else Inf
}))]

r_real_data <- r_real_data[, c(exclude_cols[1], ordered_colnames, exclude_cols[2:3])]

r_real_data$haplotype <- apply(r_real_data[, ordered_colnames], 1, function(row) {paste(row, collapse = "_")})


#plot phased haplotypes freq
r_real_data_multiallelic <- r_real_data[r_real_data$HAPLO_FREQ_RECALC < 1,] #keep multiallelic
r_real_data_multiallelic <- r_real_data_multiallelic[r_real_data_multiallelic$HAPLO_FREQ_RECALC != 0,] #remove any ceros (WHY ARE THERE CEROS AFTER INCLUDING ANC DATA?)

# Function to calculate the median
median_fn <- function(data, indices) {
  return(median(data[indices], na.rm = TRUE))
}

haplotype_counts <- r_real_data_multiallelic %>% 
  group_by(haplotype) %>%
  summarise(
    hap_count = n(),
    median_HAPLO_FREQ_RECALC = median(HAPLO_FREQ_RECALC, na.rm = TRUE),
    ci_HAPLO_FREQ_RECALC = ifelse(hap_count > 1, {
      # Bootstrap for CI
      boot_result <- boot(HAPLO_FREQ_RECALC, median_fn, R = 1000)
      boot_ci <- boot.ci(boot_result, type = "perc")$percent[4:5] # 95% CI
      paste0(boot_ci[1], "-", boot_ci[2])
    }, NA)
  ) %>%
  arrange(desc(hap_count))

haplotype_counts <- haplotype_counts %>%
  separate(ci_HAPLO_FREQ_RECALC, into = c("low_CI", "high_CI"), sep = "-", convert = TRUE)

haplotype_counts <- haplotype_counts %>%
  mutate(freq = hap_count/sum(hap_count))

# Perform the correlation test
cor_test <- cor.test(haplotype_counts$freq, haplotype_counts$median_HAPLO_FREQ_RECALC, method = "spearman", exact = F)

# Extract the correlation coefficient and p-value
cor_coef <- round(cor_test$estimate, 2)
p_val <- formatC(cor_test$p.value, format = "e", digits = 2)

ggplot(haplotype_counts, aes(x = freq, y = median_HAPLO_FREQ_RECALC)) +
  geom_jitter(width = 0.025, height = 0, alpha = 0.3, color = "gray40", size = 4) +
  theme_minimal() +
  geom_smooth(method = "lm", color = "red", fill = "pink2") +
  annotate("text", x = Inf, y = -Inf, label = paste("r =", cor_coef, "\n", "p =", p_val),
           hjust = 1.1, vjust = -1.1, size = 4, color = "black")

thresh <- 0
haplotype_df <- haplotype_counts[haplotype_counts$hap_count > thresh,] #filter > n for better vis

ggplot(haplotype_df, aes(x = reorder(haplotype, freq), y = freq)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = sprintf("%.2f", median_HAPLO_FREQ_RECALC)), 
            hjust = -0.1, size = 3.5, color = "black") + 
  coord_flip() +  
  theme_minimal() +
  labs(x = "Phased Haplotype", 
       y = paste0("Frequency"), 
       title = paste0("Order: ", paste(ordered_colnames, collapse = ", ")), 
       subtitle = paste0("Present in > ", thresh, " multiallelic samples"))

ggplot(haplotype_df, aes(x = reorder(haplotype, freq), y = median_HAPLO_FREQ_RECALC)) +
  geom_errorbar(aes(ymin = low_CI, ymax = high_CI, color = freq), width = 0.8) +
  geom_point(aes(color = freq), size = 4, shape = 16) +  
  scale_color_gradient(low = "steelblue", high = "red") +  
  scale_fill_gradient(low = "steelblue", high = "red") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    x = "Phased Haplotype",
    y = "Median Within-Sample Haplotype Frequency (95% CI)",
    title = paste0("Order: ", paste(ordered_colnames, collapse = ", ")),
    subtitle = paste0("Present in > ", thresh, " multiallelic samples")) +
  coord_flip()


#run moire with only  dhps-dhfr amplicons

library(moire)
library(optparse)

# Define the command-line arguments
option_list <- list(
  make_option(c("--allele_table", "-i"), type = "character", help = "FILTERED allele table from mad4hatter v0.1.8", default = NULL),
  make_option(c("--output_prefix", "-o"), type = "character", help = "Distinctive prefix for your output files", default = NULL)
)

# Parse the command-line arguments
opt <- parse_args(OptionParser(option_list = option_list))

alleles<-read.csv(opt$allele_table)
colnames(alleles)[1] <- "sample_id"

#subset amplicons of interest
amps_names <- c("Pf3D7_04_v3-748105-748359-1B", "Pf3D7_04_v3-748374-748611-1B") #dhfr
amps_names <- c(amps_names, "Pf3D7_08_v3-549583-549807-1B", "Pf3D7_08_v3-549960-550215-1B") #dhps

condition <- alleles$locus %in% amps_names
subsetted_df <- subset(alleles, condition)

#corroborate
unique(subsetted_df$locus) == amps_names

# set MOIRE parameters
dat_filter <- moire::load_long_form_data(subsetted_df)
burnin <- 1e4
num_samples <- 1e4
pt_chains <- seq(1, .5, length.out = 20)

#run moire
mcmc_results <- moire::run_mcmc(
  dat_filter, is_missing = dat_filter$is_missing,
  verbose = T, burnin = burnin, samples_per_chain = num_samples,
  pt_chains = pt_chains, pt_num_threads = length(pt_chains),
  thin = 10
)

#checkpoint
saveRDS(mcmc_results, paste0(opt$output_prefix, "_MOIRE-RESULTS_dhfr_dhps_only.RDS"))

#resume checkpoint
mcmc_results <- readRDS(paste0(opt$output_prefix, "_MOIRE-RESULTS_dhfr_dhps_only.RDS"))

eff_coi <- moire::summarize_effective_coi(mcmc_results)
naive_coi <- moire::summarize_coi(mcmc_results)
relatedness <- moire::summarize_relatedness(mcmc_results)

input_df <- merge(naive_coi, eff_coi, by="sample_id")
input_df <- merge(input_df, relatedness, by="sample_id")

write.csv(input_df, paste0(opt$output_prefix, "_moire_output_dhfr_dhps.csv"))

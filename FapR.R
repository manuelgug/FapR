library(dplyr)
library(ggplot2)
library(gridExtra)
library(optparse)

# Define the command-line arguments
option_list <- list(
  make_option(c("--resmarkers_table", "-i"), type = "character", help = "FILTERED resmarkers table from mad4hatter v0.1.8", default = "../results_v0.1.8_RESMARKERS_FIX/SMC2_NextSeq03_221122_RESULTS_v0.1.8_FILTERED/resmarker_table_global_max_0_filtered.csv"),
  make_option(c("--output_prefix", "-o"), type = "character", help = "Distinctive prefix for your output files", default = "test")
  #, make_option(c("--moire_output", "-m"), type = "character", help = "Moire output calculated from the dhfr-dhps amplicons", default = "HFS22_01_moire_output_dhfr_dhps.csv")
)

# Parse the command-line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Check if the required arguments are provided
if (is.null(opt$resmarkers_table) || is.null(opt$output_prefix)) {
  stop("All --resmarkers_table and --output_prefix arguments are required.")
}

################## IMPORT AND FORMAT DATA ################## 

resmarkers_table <- read.csv(opt$resmarkers_table)

#subset relevant markers
markers_to_phase <- c("dhfr_51", "dhfr_59", "dhfr_108", "dhps_431", "dhps_437", "dhps_540", "dhps_581")

#create resmarker column
resmarkers_table$resmarker <-  paste0(resmarkers_table$Gene, "_", resmarkers_table$CodonID)

resmarkers_table <- resmarkers_table %>%
  filter(grepl(paste(markers_to_phase, collapse = "|"), resmarker))

#for resmarkers with 2 amplicons (as dhps_581), keep the one with the highest amount of reads for each variant
resmarkers_table <- resmarkers_table %>%
  group_by(SampleID, resmarker, AA) %>%
  filter(Reads == max(Reads)) %>%
  ungroup()

resmarkers_table <- as.data.frame(resmarkers_table)

#if there is no norm.read.locus column, create it
resmarkers_table<- resmarkers_table %>%
  group_by(SampleID,Gene, CodonID, CodonStart) %>%
  mutate(norm.reads.locus = Reads/sum(Reads)) %>%
  mutate(n.alleles = n())

resmarkers_table <- resmarkers_table[,c("SampleID", "resmarker", "AA", "norm.reads.locus")]

#moire_output <- read.csv(opt$moire_output) #dhfr-dhps specific moire run

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
      
      #do CV and probs agree with each other?
      if (lowest_CV == highest_prob) {
        most_likely_hap <- paste(as.matrix(comb_alleles_matrix[highest_prob, ]), collapse = "_")
        print(paste(sample, "#", i_counter, ":", most_likely_hap, "is the most likely true haplotype.", collapse = " "))
      } else {
        most_likely_hap1 <- paste(as.matrix(comb_alleles_matrix[highest_prob, ]), collapse = "_")
        most_likely_hap2 <- paste(as.matrix(comb_alleles_matrix[lowest_CV, ]), collapse = "_")
        print(paste(sample, "#", i_counter, ": One of", most_likely_hap1, "and", most_likely_hap2, "is the most likely true haplotype. Visually examine the plot."))
      }
      
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

################# LIMIT OF DETETION FLAGGING OF HAPLOS ##############

#thresholds may change. current ones work perfectly for the DD2 gradient, but more sequencing may be needed. ALSO, no dhps gradients atm
LOD_dhfr_51 <- 0.0372670807453416 # according to DD2 gradient: lowest correct value
LOD_dhfr_59 <- 0.0372670807453416 # according to DD2 gradient: lowest correct value
LOD_dhfr_108 <- 0.299350895389018 # according to DD2 gradient: lowest correct value; 0.164061810215753 in HB3 gradient and 0.133363977705396 in QS1QC2 mix, but using the more astringent tho. VARÍA ENTRE RUNS!! PENSAR QUÑÉ HACER.
LOD_dhps_431 <- 0.00991796253214094 # according to HB3 gradient: lowest correct value
LOD_dhps_437 <- 0.00991796253214094 # according to HB3 gradient: lowest correct value
#LOD_dhps_540 <- ? #no gradient for this position
#LOD_dhps_581 <- ? #no gradient for this position

#for each sample, if loci is multiallelic, flag haplos freq below LOD for each loci as "dubious" for each allele.
flag_haplotypes <- function(df, locus, lod_threshold) {
  # Create a new column with "dubious" for haplotypes with HAPLO_FREQ_RECALC < LOD threshold, "correct" otherwise
  df[paste0("flag_", locus)] <- ifelse(df$HAPLO_FREQ_RECALC < lod_threshold, "dubious", "correct")
  return(df)
}

RESULTS_FINAL_FLAGGED <- RESULTS_FINAL %>%
  group_by(SampleID) %>%
  do(flag_haplotypes(., "dhps_431", LOD_dhps_431)) %>%
  do(flag_haplotypes(., "dhps_437", LOD_dhps_437)) %>%
  do(flag_haplotypes(., "dhfr_51", LOD_dhfr_51)) %>%
  do(flag_haplotypes(., "dhfr_59", LOD_dhfr_59)) %>%
  do(flag_haplotypes(., "dhfr_108", LOD_dhfr_108))

#remove haplotypes with "dubious" on all flags
all_dubious_rows <- rowSums(RESULTS_FINAL_FLAGGED[, grepl("^flag_", names(RESULTS_FINAL_FLAGGED))] == "dubious") == length(RESULTS_FINAL_FLAGGED[, grepl("^flag_", names(RESULTS_FINAL_FLAGGED))])
RESULTS_FINAL_FLAGGED <- RESULTS_FINAL_FLAGGED[!all_dubious_rows, ]

#flag as PASSED if the haplotype is found to be correct in a certain amount of the multiallelic samples (pop_freq_thresh)
RESULTS_FINAL_FLAGGEDL_multiallelic <- RESULTS_FINAL_FLAGGED[RESULTS_FINAL_FLAGGED$HAPLO_FREQ_RECALC < 1, ]

haplos <- paste(RESULTS_FINAL_FLAGGEDL_multiallelic$dhps_431, RESULTS_FINAL_FLAGGEDL_multiallelic$dhps_437, RESULTS_FINAL_FLAGGEDL_multiallelic$dhps_540, RESULTS_FINAL_FLAGGEDL_multiallelic$dhps_581, RESULTS_FINAL_FLAGGEDL_multiallelic$dhfr_51, RESULTS_FINAL_FLAGGEDL_multiallelic$dhfr_59, RESULTS_FINAL_FLAGGEDL_multiallelic$dhfr_108, sep = "_")
haplo_counts <- table(haplos)
haplo_counts <- as.data.frame(haplo_counts)
haplo_counts$haplos <- factor(haplo_counts$haplos, levels = haplo_counts$haplos[order(-haplo_counts$Freq)])
haplo_counts$proportion <- haplo_counts$Freq / sum(haplo_counts$Freq)

thresh_pop_freq = mean(haplo_counts$proportion) # DECIDE ON HOW TO CALCUALTE THIS THRESHOLD!!!! MEAN IS JUST FOR TESTING PURPOSES

#if haplo_counts$proportion => threshold_pop_freq, add RESULTS_FINAL_FLAGGED$flag_pop_freq = "PASSED" whenever haplo_counts$haplo matches RESULTS_FINAL_FLAGGED$haplotype, else add NA

RESULTS_FINAL_FLAGGED <- RESULTS_FINAL_FLAGGED %>%
  left_join(haplo_counts, by = c("haplotype" = "haplos")) %>%
  mutate(flag_pop_freq = ifelse(proportion >= thresh_pop_freq, "PASSED", NA)) %>%
  select(-proportion)%>%
  select(-Freq)

write.csv(RESULTS_FINAL_FLAGGED, paste0(opt$output_prefix, "_phased_haplos.csv"), row.names =FALSE)

#subset the dataframe to remove rows that have at least one "dubious" flag and no "PASSED" flag
dubious_rows <- rowSums(RESULTS_FINAL_FLAGGED[, grepl("^flag_d", names(RESULTS_FINAL_FLAGGED))] == "dubious") > 0 & is.na(RESULTS_FINAL_FLAGGED$flag_pop_freq)
RESULTS_FINAL_FLAGGED_correct_only <- RESULTS_FINAL_FLAGGED[!dubious_rows, ]

write.csv(RESULTS_FINAL_FLAGGED_correct_only, paste0(opt$output_prefix, "_phased_haplos_correct_only.csv"), row.names =FALSE)


################## CHECKS, VISUALIZATIONS, VALIDATION ################## 

generate_haplo_summary_plots <- function(RESULTS_FINAL, props_plot, props_plot_multi, hist_plot, profile_plot) {
  # Formatting
  haplos <- paste(RESULTS_FINAL$dhps_431, RESULTS_FINAL$dhps_437, RESULTS_FINAL$dhps_540, RESULTS_FINAL$dhps_581, RESULTS_FINAL$dhfr_51, RESULTS_FINAL$dhfr_59, RESULTS_FINAL$dhfr_108, sep = "_")
  haplo_counts <- table(haplos)
  haplo_counts <- as.data.frame(haplo_counts)
  haplo_counts$haplos <- factor(haplo_counts$haplos, levels = haplo_counts$haplos[order(-haplo_counts$Freq)])
  haplo_counts$proportion <- haplo_counts$Freq / sum(haplo_counts$Freq)
  
  # Barplot of counts
  p <- ggplot(haplo_counts, aes(x = haplos, y = Freq, fill = haplos)) +
    geom_bar(stat = "identity") +
    labs(title = "Haplotype Counts", x = NULL, y = "Samples") +
    theme(axis.text.x = element_text(angle = 65, hjust = 1)) +
    guides(fill = "none")
  
  # Pie chart of proportions of haplos
  q <- ggplot(haplo_counts, aes(x = "", y = proportion, fill = haplos)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    labs(title = "Haplotype Frequency", x = "") +
    scale_x_discrete(labels = NULL)+
    theme(axis.title.x = element_blank()) 
  
  grid_plot1 <- grid.arrange(p, q, ncol = 2)
  
  RESULTS_FINAL$haplotype <- paste(RESULTS_FINAL$dhps_431, RESULTS_FINAL$dhps_437, RESULTS_FINAL$dhps_540, RESULTS_FINAL$dhps_581, RESULTS_FINAL$dhfr_51, RESULTS_FINAL$dhfr_59, RESULTS_FINAL$dhfr_108, sep = "_")
  RESULTS_FINAL_multiallelic <- RESULTS_FINAL[RESULTS_FINAL$HAPLO_FREQ_RECALC < 1, ]
  
  # Formatting
  haplos <- paste(RESULTS_FINAL_multiallelic$dhps_431, RESULTS_FINAL_multiallelic$dhps_437, RESULTS_FINAL_multiallelic$dhps_540, RESULTS_FINAL_multiallelic$dhps_581, RESULTS_FINAL_multiallelic$dhfr_51, RESULTS_FINAL_multiallelic$dhfr_59, RESULTS_FINAL_multiallelic$dhfr_108, sep = "_")
  haplo_counts <- table(haplos)
  haplo_counts <- as.data.frame(haplo_counts)
  haplo_counts$haplos <- factor(haplo_counts$haplos, levels = haplo_counts$haplos[order(-haplo_counts$Freq)])
  haplo_counts$proportion <- haplo_counts$Freq / sum(haplo_counts$Freq)
  
  # Barplot of counts
  r <- ggplot(haplo_counts, aes(x = haplos, y = Freq, fill = haplos)) +
    geom_bar(stat = "identity") +
    labs(title = "Haplotype Counts", x = NULL, y = "Samples") +
    theme(axis.text.x = element_text(angle = 65, hjust = 1)) +
    guides(fill = "none")
  
  # Pie chart of proportions of haplos
  s <- ggplot(haplo_counts, aes(x = "", y = proportion, fill = haplos)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    labs(title = "Haplotype Frequencies", x = "") +
    scale_x_discrete(labels = NULL)+
    theme(axis.title.x = element_blank()) 
  
  grid_plot1_1 <- grid.arrange(r, s, ncol = 2)
  
  # Histograms of frequency of each haplo
  histogram_plots <- list()
  
  for (haplo in unique(haplos)) {
    freq_vector <- RESULTS_FINAL$HAPLO_FREQ_RECALC[haplos == haplo]
    
    plot <- ggplot(data = data.frame(Frequency = freq_vector)) +
      geom_histogram(aes(x = Frequency), bins = 10, fill = "cadetblue3", color = "cadetblue3") +
      labs(title = haplo, x = "Frequency", y = "Samples") + coord_cartesian(xlim = c(0, 1))  # Set x-axis limits
    
    histogram_plots[[haplo]] <- plot
  }
  
  grid_plot2 <- grid.arrange(grobs = histogram_plots, ncol = 4)
  
  # Haplo profile (barplot) multiallelic only
  a <- ggplot(RESULTS_FINAL_multiallelic, aes(x = SampleID, y = HAPLO_FREQ_RECALC, fill = haplotype)) +
    geom_bar(stat = "identity") +
    labs(title = "Haplotype Profiles", x = NULL, y = "Frequency") +
    theme(axis.text.x = element_text(angle = 65, hjust = 1) , legend.position = "top") +
    guides(fill = guide_legend(title = "Haplotype"))
  
  grid_plot3 <- a
  
  # Save plots to files
  ggsave(props_plot, grid_plot1, width = 16, height = 9)
  ggsave(props_plot_multi, grid_plot1_1, width = 16, height = 9)
  ggsave(hist_plot, grid_plot2, width = 14, height = 10)
  ggsave(profile_plot, grid_plot3, width = 12, height = 9)
}

# plot all haplotypes
generate_haplo_summary_plots(RESULTS_FINAL, 
                             paste0(opt$output_prefix, "_haplo_counts_proportions.png"), 
                             paste0(opt$output_prefix, "_haplo_counts_proportions_multiallelic.png"), 
                             paste0(opt$output_prefix, "_haplos_histograms.png"), 
                             paste0(opt$output_prefix, "_haplo_profile_multiallelic.png"))

# plot only correct haplotypes
generate_haplo_summary_plots(RESULTS_FINAL_FLAGGED_correct_only, 
                             paste0(opt$output_prefix, "_haplo_counts_proportions_correct_only.png"), 
                             paste0(opt$output_prefix, "_haplo_counts_proportions_correct_only_multiallelic.png"), 
                             paste0(opt$output_prefix, "_haplos_histograms_correct_only.png"), 
                             paste0(opt$output_prefix, "_haplo_profile_multiallelic_correct_only.png"))


RESULTS_FINAL_FLAGGED_correct_only$Haplotype <- do.call(paste0, RESULTS_FINAL_FLAGGED_correct_only[, 2:8])
RESULTS_FINAL_FLAGGED_correct_only <- subset(RESULTS_FINAL_FLAGGED_correct_only, HAPLO_FREQ_RECALC != 1)

RESULTS_FINAL_FLAGGED_correct_only <- RESULTS_FINAL_FLAGGED_correct_only %>%
  arrange(SampleID, Haplotype)

coocurring_haplos <- RESULTS_FINAL_FLAGGED_correct_only %>%
  group_by(SampleID) %>%
  summarize(Coocurring_Haplotypes = paste(Haplotype, collapse = "_"))

coocurring_haplos_counts<-as.data.frame(colSums(table(coocurring_haplos)))

df<-cbind(rownames(coocurring_haplos_counts), coocurring_haplos_counts)

colnames(df)<-c("Haplotypes","Count")

df <- df %>%
  filter(grepl("_", Haplotypes))

df$Haplotypes <- reorder(df$Haplotypes, -df$Count)

cooc <- ggplot(df, aes(x = Haplotypes, y = Count)) +
  geom_bar(stat = "identity") +
  labs(title = "Co-ocurring haplotypes", x = "Haplotype combos", y = "Samples") +
  theme(axis.text.x = element_text(angle = 65, hjust = 1))

ggsave(paste0(opt$output_prefix, "_haplo_combos_multiallelic_correct_only.png"), cooc, width = 12, height = 9)

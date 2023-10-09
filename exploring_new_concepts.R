library(dplyr)
library(ggplot2)
library(gridExtra)


################## IMPORT AND FORMAT DATA ################## 

resmarkers_table <- read.csv("HSF22_01_resmarker_table_global_max_0_filtered_resmarkers_FIX_has_DD2.csv")

#subset relevant markers
markers_to_phase <- c("dhfr_51", "dhfr_59", "dhfr_108", "dhps_437", "dhps_540")

resmarkers_table <- resmarkers_table %>%
  filter(grepl(paste(markers_to_phase, collapse = "|"), resmarker))

resmarkers_table <- resmarkers_table[,c("SampleID", "resmarker", "AA", "norm.reads.locus")]


################## MAIN LOOP ################## 

unique_samples <- unique(resmarkers_table$SampleID)
RESULTS_FINAL <- data.frame(SampleID = character(0), dhps_437 = character(0), dhps_540 = character(0), dhfr_51 = character(0), dhfr_59 = character(0), dhfr_108 = character(0), HAPLO_FREQ = numeric(0), HAPLO_FREQ_RECALC = numeric(0))

# INIT LOOP HERE!

for (sample in unique_samples){
  
  eCOI_counter <- 0
  MOST_LIKELY_HAPLOS <- data.frame()
  MOST_LIKELY_HAPLOS_FREQS <- data.frame()
  RESULTS <- data.frame(SampleID = character(0), dhps_437 = character(0), dhps_540 = character(0), dhfr_51 = character(0), dhfr_59 = character(0), dhfr_108 = character(0), HAPLO_FREQ = numeric(0), HAPLO_FREQ_RECALC = numeric(0))
  
  # 1) select sample
  sID <- resmarkers_table[resmarkers_table$SampleID == sample,]
  
  # 2) select sample's eCOI
  #algo como: eCOI <- coi[coi$SampleID == sample,]
  eCOI<- 2 ### THIS IS A PLACE HOLDER FOR ECOI, WHICH SHOULD BE AUTOMATICALLY SET WHEN SELECTING THE SAMPLE (MAYBE ROUNDED ALSO? SHOULD BE INTEGER!)
  
  # 3) format data
  new_df <- data.frame(matrix(ncol = length(sID$resmarker), nrow=1))
  colnames(new_df) <- sID$resmarker
  new_df[1,] <-sID$AA
  new_df <- rbind(new_df, sID$norm.reads.locus )
  
  unique_resmarkers <- unique(colnames(new_df))
  
  # Initialize a list to store the resulting dataframes
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
  
  alleles<-list(resulting_dataframes$dhps_437$dhps_437, 
                resulting_dataframes$dhps_540$dhps_540, 
                resulting_dataframes$dhfr_51$dhfr_51,
                resulting_dataframes$dhfr_59$dhfr_59,
                resulting_dataframes$dhfr_108$dhfr_108) #order is important
  
  freqs<-list(resulting_dataframes$dhps_437$norm.reads.locus, 
              resulting_dataframes$dhps_540$norm.reads.locus, 
              resulting_dataframes$dhfr_51$norm.reads.locus,
              resulting_dataframes$dhfr_59$norm.reads.locus,
              resulting_dataframes$dhfr_108$norm.reads.locus) #order is important
  
  comb_alleles <- expand.grid(alleles)
  comb_freqs <- expand.grid(freqs)
  
  comb_alleles_matrix <- as.data.frame(comb_alleles)
  colnames(comb_alleles_matrix) <- c("dhps_437", "dhps_540", "dhfr_51","dhfr_59", "dhfr_108")
  comb_freqs_matrix <- as.data.frame(comb_freqs)
  colnames(comb_freqs_matrix) <- c("dhps_437", "dhps_540", "dhfr_51","dhfr_59", "dhfr_108")
  
  # 4) phase
  if (dim(comb_alleles_matrix)[1] != 1){ #basically, don't process monoallelic samples 'cause they make the loop crash
    
    while (dim(MOST_LIKELY_HAPLOS_FREQS)[1] == 0 ||  eCOI_counter != eCOI) {
      
      eCOI_counter <- eCOI_counter + 1
      
      # Calculate probs if all haplotypes were present
      comb_freqs_matrix$probs <- comb_freqs_matrix$dhps_437 * comb_freqs_matrix$dhps_540 * comb_freqs_matrix$dhfr_51 * comb_freqs_matrix$dhfr_59  * comb_freqs_matrix$dhfr_108
      
      # Calculate SD and CV
      comb_freqs_matrix$freq_mean <- rowMeans(comb_freqs_matrix, na.rm = TRUE)
      comb_freqs_matrix$SD <- apply(comb_freqs_matrix[, 1:5], 1, sd)
      comb_freqs_matrix$CV <- (comb_freqs_matrix$SD / comb_freqs_matrix$freq_mean)
      
      ## Select the "BEST" haplo: highest prob and lowest CV
      lowest_CV <- which.min(comb_freqs_matrix$CV)
      highest_prob <- as.numeric(which.max(comb_freqs_matrix$probs))
      
      #do CV and probs agree with each other?
      if (lowest_CV == highest_prob) {
        most_likely_hap <- paste(as.matrix(comb_alleles_matrix[highest_prob, ]), collapse = "_")
        print(paste(sample, "#", eCOI_counter, ":", most_likely_hap, "is the most likely true haplotype.", collapse = " "))
      } else {
        most_likely_hap1 <- paste(as.matrix(comb_alleles_matrix[highest_prob, ]), collapse = "_")
        most_likely_hap2 <- paste(as.matrix(comb_alleles_matrix[lowest_CV, ]), collapse = "_")
        print(paste(sample, "#", eCOI_counter, ": One of", most_likely_hap1, "and", most_likely_hap2, "is the most likely true haplotype. Visually examine the plot."))
      }
      
      # Append most likely haplo
      MOST_LIKELY_HAPLOS <- rbind(MOST_LIKELY_HAPLOS, comb_alleles_matrix[highest_prob, ])
      temp <- comb_freqs_matrix[highest_prob, ]
      temp$HAPLO_FREQ <- min(comb_freqs_matrix[highest_prob, 1:5])
      MOST_LIKELY_HAPLOS_FREQS <- rbind(MOST_LIKELY_HAPLOS_FREQS, temp)
      
      # Select minimum allele freq from the most likely haplotype
      min_allele_from_most_lilely_hap <- min(comb_freqs_matrix[highest_prob, 1:5])
      
      # Boolean mask to detect alleles that are present on the most likely haplotype
      row_to_match <- as.matrix(comb_alleles_matrix[highest_prob, ])
      mask <- sapply(colnames(comb_alleles_matrix), function(col_name) {
        comb_alleles_matrix[, col_name] == row_to_match[, col_name]
      })
      
      # Subtract min_allele_from_most_likely_hap from the cells where mask is TRUE and ignore the specified column
      comb_freqs_matrix <- comb_freqs_matrix[, 1:5]
      comb_freqs_matrix[mask] <- comb_freqs_matrix[mask] - min_allele_from_most_lilely_hap
    }
    
    #recalculate proportions of final haplos
    MOST_LIKELY_HAPLOS_FREQS$HAPLO_FREQ_RECALC <- MOST_LIKELY_HAPLOS_FREQS$HAPLO_FREQ / sum(MOST_LIKELY_HAPLOS_FREQS$HAPLO_FREQ)
    
    RESULTS <- cbind(SampleID = sample, MOST_LIKELY_HAPLOS, HAPLO_FREQ = MOST_LIKELY_HAPLOS_FREQS$HAPLO_FREQ, HAPLO_FREQ_RECALC = MOST_LIKELY_HAPLOS_FREQS$HAPLO_FREQ_RECALC)
    
  }else{ 
    
    #FORMAT AND ADD MONOALLELIC SAMPLES HERE
    RESULTS <- cbind(SampleID = sample, comb_alleles_matrix, HAPLO_FREQ = 1, HAPLO_FREQ_RECALC = 1)
    
  }
  
  #DONE
  RESULTS_FINAL <- rbind(RESULTS, RESULTS_FINAL)
  
}

write.csv(RESULTS_FINAL, "phased_haplos.csv", row.names =FALSE)


################## CHECKS, VISUALIZATIONS, VALIDATION ################## 


#multi <- RESULTS_FINAL$HAPLO_FREQ_RECALC < 1
#hist(RESULTS_FINAL[multi, c("HAPLO_FREQ")])
#hist(RESULTS_FINAL[multi, c("HAPLO_FREQ_RECALC")])

# formatting
haplos<-paste(RESULTS_FINAL$dhps_437, RESULTS_FINAL$dhps_540, RESULTS_FINAL$dhfr_51, RESULTS_FINAL$dhfr_59, RESULTS_FINAL$dhfr_108, sep = "_")

haplo_counts <- table(haplos)
haplo_counts <- as.data.frame(haplo_counts)
haplo_counts$haplos <- factor(haplo_counts$haplos, levels = haplo_counts$haplos[order(-haplo_counts$Freq)])

#barplot of counts
p <-ggplot(haplo_counts, aes(x = haplos, y = Freq)) +
  geom_bar(stat = "identity", fill = "cadetblue3") +
  labs(title = "Haplo Counts", x = "Haplo", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("counts.png", p, width = 12, height = 9)  # Adjust width and height as needed


# Histograms of frequency of each haplo
histogram_plots <- list()

for (haplo in unique(haplos)) {
  freq_vector <- RESULTS_FINAL$HAPLO_FREQ_RECALC[haplos == haplo]
  
  plot <- ggplot(data = data.frame(Frequency = freq_vector)) +
    geom_histogram(aes(x = Frequency), bins = 10, fill = "cadetblue3", color ="white") +
    labs(title = haplo, x = "Frequency", y = "Count")
  
  histogram_plots[[haplo]] <- plot
}

grid_plot <- grid.arrange(grobs = histogram_plots, ncol = 4)  # Change ncol as needed
ggsave("histograms_grid.png", grid_plot, width = 12, height = 9)  # Adjust width and height as needed





## with eCOI, determine how many MAX haplos are expected
## find the "BEST" complementary haplo(s) that add up to 100 as close as possible (without considering eCOI)
## find the "BEST" complementary haplo(s) that add up to 100 as close as possible considering eCOI
## same? PERFECT; not the same? WAT DO?

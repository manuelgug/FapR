
###ESTE MÉTODO ITERATIVO IGNORA COI. USANDO COI SERÍA UNA ESTRATEGIA DISTINTA, SUPONGO... PROBAR ASÍ CON TODAS LAS ITERACIONES POSIBLES PRIMERO?

library(dplyr)

################## IMPORT AND FORMAT DATA ################## 

resmarkers_table <- read.csv("HSF22_01_resmarker_table_global_max_0_filtered_resmarkers_FIX_has_DD2.csv")

#subset relevant snps
markers_to_phase <- c("dhfr_51", "dhfr_59", "dhfr_108", "dhps_437", "dhps_540")

resmarkers_table <- resmarkers_table %>%
  filter(grepl(paste(markers_to_phase, collapse = "|"), resmarker))

resmarkers_table <- resmarkers_table[,c("SampleID", "resmarker", "AA", "norm.reads.locus")]


#########################################
###test reformatting with one sample !!!!!!
##########################################

test<-resmarkers_table[resmarkers_table$SampleID =="N1933388_5_S115",] #select 1 sampleID

new_df <- data.frame(matrix(ncol = length(test$resmarker), nrow=1))
colnames(new_df) <- test$resmarker
new_df[1,] <-test$AA
new_df <- rbind(new_df, test$norm.reads.locus )

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
     resulting_dataframes$dhfr_108$dhfr_108)

freqs<-list(resulting_dataframes$dhps_437$norm.reads.locus, 
              resulting_dataframes$dhps_540$norm.reads.locus, 
              resulting_dataframes$dhfr_51$norm.reads.locus,
              resulting_dataframes$dhfr_59$norm.reads.locus,
              resulting_dataframes$dhfr_108$norm.reads.locus)

comb_alleles <- expand.grid(alleles)
comb_freqs <- expand.grid(freqs)

comb_alleles_matrix <- as.data.frame(comb_alleles)
colnames(comb_alleles_matrix) <- c("dhps_437", "dhps_540", "dhfr_51","dhfr_59", "dhfr_108")
comb_freqs_matrix <- as.data.frame(comb_freqs)
colnames(comb_freqs_matrix) <- c("dhps_437", "dhps_540", "dhfr_51","dhfr_59", "dhfr_108")


################## MAIN LOOP ################## 

MOST_LIKELY_HAPLOS <- data.frame()
MOST_LIKELY_HAPLOS_FREQS <- data.frame()


while (dim(MOST_LIKELY_HAPLOS_FREQS)[1] == 0 || sd(MOST_LIKELY_HAPLOS_FREQS[nrow(MOST_LIKELY_HAPLOS_FREQS), 1:5]) > 0.00000000001) {
  # Calculate mean freq for all possible haplotypes
  comb_freqs_matrix$freq_mean <- rowMeans(comb_freqs_matrix, na.rm = TRUE)
  
  # Calculate probs if all haplotypes were present
  comb_freqs_matrix$probs <- comb_freqs_matrix$dhps_437 * comb_freqs_matrix$dhps_540 * comb_freqs_matrix$dhfr_51 * comb_freqs_matrix$dhfr_59  * comb_freqs_matrix$dhfr_108
  
  # Calculate Standard Deviation (SD)
  comb_freqs_matrix$SD <- apply(comb_freqs_matrix[, 1:5], 1, sd)
  
  # Calculate Coefficient of Variation (CV)
  comb_freqs_matrix$CV <- (comb_freqs_matrix$SD / comb_freqs_matrix$freq_mean)
  
  ## Select the "BEST" haplo: highest prob and lowest CV
  lowest_CV <- which.min(comb_freqs_matrix$CV)
  highest_prob <- as.numeric(which.max(comb_freqs_matrix$probs))
  
  if (lowest_CV == highest_prob) {
    most_likely_hap <- paste(as.matrix(comb_alleles_matrix[highest_prob, ]), collapse = "_")
    print(paste(most_likely_hap, "is the most likely true haplotype.", collapse = " "))
  } else {
    most_likely_hap1 <- paste(as.matrix(comb_alleles_matrix[highest_prob, ]), collapse = "_")
    most_likely_hap2 <- paste(as.matrix(comb_alleles_matrix[lowest_CV, ]), collapse = "_")
    print(paste("One of", most_likely_hap1, "and", most_likely_hap2, "is the most likely true haplotype. Visually examine the plot."))
    
    # Plot CV vs. probs
    #plot(comb_freqs_matrix$CV, comb_freqs_matrix$probs)
  }
  
  # Append next best haplo
  MOST_LIKELY_HAPLOS <- rbind(MOST_LIKELY_HAPLOS, comb_alleles_matrix[highest_prob, ])
  temp <- comb_freqs_matrix[highest_prob, ]
  temp$HAPLO_FREQ <- min(comb_freqs_matrix[highest_prob, 1:5])
  MOST_LIKELY_HAPLOS_FREQS <- rbind(MOST_LIKELY_HAPLOS_FREQS, temp)
  
  # Select minimum allele freq of the most likely haplotype
  min_allele_from_most_lilely_hap <- min(comb_freqs_matrix[highest_prob, 1:5])
  
  # Boolean mask to detect alleles that are present on the most likely haplotype
  row_to_match <- as.matrix(comb_alleles_matrix[highest_prob, ])
  mask <- sapply(colnames(comb_alleles_matrix), function(col_name) {
    comb_alleles_matrix[, col_name] == row_to_match[, col_name]
  })
  
  # Subtract min_allele_from_most_lilely_hap from the cells where mask is TRUE and ignore the specified column
  comb_freqs_matrix <- comb_freqs_matrix[, 1:5]
  comb_freqs_matrix[mask] <- comb_freqs_matrix[mask] - min_allele_from_most_lilely_hap
}

#########################################################################################

#MOST_LIKELY_HAPLOS
#MOST_LIKELY_HAPLOS_FREQS

RESULT <- cbind(MOST_LIKELY_HAPLOS, MOST_LIKELY_HAPLOS_FREQS$HAPLO_FREQ)
colnames(RESULT)[5] <- "HAPLO_FREQ"










## with eCOI, determine how many MAX haplos are expected
## find the "BEST" complementary haplo(s) that add up to 100 as close as possible (without considering eCOI)
## find the "BEST" complementary haplo(s) that add up to 100 as close as possible considering eCOI
## same? PERFECT; not the same? WAT DO?

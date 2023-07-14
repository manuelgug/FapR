
library(ggplot2)
library(ggbeeswarm)
library(dplyr)
library(gridExtra)
library(purrr)


args = commandArgs(trailingOnly=T)
allele_data_filtered=args[1]
outfile=args[2]

############################### DATA PREPARATION ###########################################

filtered_allele.data<-read.csv(allele_data_filtered, sep ="\t")
#filtered_allele.data<-read.csv("allele_data_filtered_TESS2_run1.txt", sep ="\t") # has the DD2 mixes

#subset amplicons of interest
amps_names <- c("Pf3D7_04_v3-748105-748359-1B", "Pf3D7_04_v3-748374-748611-1B") #dhfr
amps_names <- c(amps_names, "Pf3D7_08_v3-549583-549807-1B", "Pf3D7_08_v3-549960-550215-1B") #dhps

amplicons_of_interest <- filtered_allele.data[filtered_allele.data$locus %in% amps_names, ]

#change names
new_names <- c(
  "Pf3D7_04_v3-748105-748359-1B" = "dhfr16_dhfr51_dhfr59",
  "Pf3D7_04_v3-748374-748611-1B" = "dhfr108_dhfr164",
  "Pf3D7_08_v3-549583-549807-1B" = "dhps436_dhps437_dhps431",
  "Pf3D7_08_v3-549960-550215-1B" = "dhps540_dhps581"
)

for (i in seq_along(new_names)) {
  pattern <- paste0("^", names(new_names)[i])
  amplicons_of_interest$locus <- sub(pattern, new_names[i], amplicons_of_interest$locus)
  amplicons_of_interest$allele <- sub(pattern, new_names[i], amplicons_of_interest$allele)
}

# Subset monoclonal samples (for the 4 given amplicons)
monoclonal_samples <- amplicons_of_interest %>%
  group_by(sampleID) %>%
  filter(all(n.alleles == 1))

# Subset polyclonal samples (for the 4 given amplicons)
polyclonal_samples <- anti_join(amplicons_of_interest, monoclonal_samples, by = c("locus", "sampleID"))


############################# EXPLORATORY DATA ANALYSIS ###########################################

ratio <- data.frame(
  SampleType = c("Monoclonal", "Polyclonal"),
  Count = c(length(unique(monoclonal_samples$sampleID)), length(unique(polyclonal_samples$sampleID)))
)

ratio$Percentage <- round(ratio$Count / sum(ratio$Count) * 100, 1)

pie_chart <- ggplot(ratio, aes(x = "", y = Count, fill = SampleType)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  labs(fill = "Sample Type", x = NULL, y = NULL, title = "Monoclonal to Polyclonal Ratio") +
  theme_void() +
  geom_text(aes(label = paste0(Count, " (", Percentage, "%)")), position = position_stack(vjust = 0.5))

poly_var<-polyclonal_samples %>%
  group_by(locus) %>%
  summarize(variance = var(norm.reads.locus))

variance<-ggplot(poly_var, aes(x = locus, y = variance, fill=locus)) +
  geom_bar(stat = "identity") +
  xlab("Amplicon") +
  ylab("Variance") +
  ggtitle("Variance in allele frequency for each amplicon") +
  theme_minimal()+
  guides(fill = "none")

p1<-ggplot(polyclonal_samples, aes(x = locus, y = n.alleles)) +
  ggbeeswarm::geom_quasirandom(aes(color = locus), dodge.width = 0.5, size = 2, show.legend = FALSE) +
  xlab("Amplicon") +
  ylab("Number of Alleles")+ 
  labs(title = "Number of alleles for each amplicon")

p2<-ggplot(polyclonal_samples, aes(x = n.alleles, fill = locus)) +
  geom_density(alpha = 0.5) +
  xlab("Number of Alleles") +
  ylab("Density") +
  scale_fill_discrete(name = "Amplicon")+ 
  labs(title = "Distribution of number of alleles")

p3 <- ggplot(polyclonal_samples, aes(x = locus, y = norm.reads.locus)) +
  ggbeeswarm::geom_quasirandom(aes(color = locus), dodge.width = 0.5, size = polyclonal_samples$n.alleles, show.legend = FALSE) +
  xlab("Amplicon") +
  ylab("Allele Frequency") +
  labs(title = "Allele frequencies for each amplicon") +
  scale_size(name = "Number of Alleles")

p4<-ggplot(polyclonal_samples, aes(x = norm.reads.locus, fill = locus)) +
  geom_density(alpha = 0.5) +
  xlab("Allele Frequency") +
  ylab("Density") +
  scale_fill_discrete(name = "Amplicon")+ 
  labs(title = "Distribution of allele frequencies")

grid <- gridExtra::arrangeGrob(pie_chart, variance, p1, p3, p2, p4, nrow = 3, ncol = 2)

ggsave(paste0(outfile,"_EDA.pdf"), grid, width = 18, height = 14, dpi = 300, device="pdf") 




###################################### PHASING ###############################################


#extract already phased samples (meaning, all amplicons minus 1 only have 1 allele [#if a sample has a value of 1 in all minus 1 unique locus, it's already phased, just split it by the number of alleles the remaining locus has.])
autophased_samples <- polyclonal_samples %>%
  group_by(sampleID) %>%
  filter(sum(n.alleles == 1) == length(unique(locus)) - 1) %>% #for these 4 loci, it means a sample has at least 3 monoallelic loci and 1 polyallelic
  ungroup()

#leave unphased samples
polyclonal_samples <- polyclonal_samples %>%
  anti_join(autophased_samples, by = c("sampleID", "locus"))

####FINAL FORMAT TO monoclonal_samples and autophased_samples
monoclonal_samples$haplotypes<-"HAPLOTYPE1" 

autophased_samples <- autophased_samples %>%
  arrange(sampleID, locus, n.alleles, desc(norm.reads.locus)) %>%
  group_by(locus) %>%
  arrange(desc(n.alleles))

autophased_samples <- autophased_samples %>% 
  group_by(sampleID) %>%
  mutate(haplotypes = ifelse(n.alleles == 1, paste0("HAPLOTYPE", 1:max(n.alleles), collapse = "_"), paste0("HAPLOTYPE", 1:max(n.alleles[n.alleles > 1])))) %>%
  arrange(sampleID)

PARTIAL_OUTPUT_1 <- rbind(monoclonal_samples, autophased_samples)

#################################################################################

#pick sample 
#sample_names[1] de ../SMC_retest_new_snps_new_output/mad4hatter-main/SMC2_NextSeq04_131222_results/allele_data.txt para probar con 1, 2, 3 y 4 alelos: GOOD
#sample_names[2] de ../SMC_retest_new_snps_new_output/mad4hatter-main/SMC2_NextSeq04_131222_results/allele_data.txt para probar con 1, 1, 2, 5 alelos GOOD
#sample_names[7] de ../SMC_retest_new_snps_new_output/mad4hatter-main/SMC2_NextSeq04_131222_results/allele_data.txt para probar con 1, 1, 3, 3 alelos: GOOD
#sample_names[9] de ../SMC_retest_new_snps_new_output/mad4hatter-main/SMC2_NextSeq04_131222_results/allele_data.txt para probar con 1, 2, 2, 3 alelos: GOOD

sample_names<-unique(polyclonal_samples$sampleID) #out of loop
PARTIAL_OUTPUT_2 <- data.frame()

for (sample in sample_names){
  
  muestra<-polyclonal_samples[grep(sample, polyclonal_samples$sampleID), ] 
  
  #sort
  sorted_polyclonal <- muestra %>%
    arrange(sampleID, locus, n.alleles, desc(norm.reads.locus)) %>%
    group_by(locus) %>%
    arrange(n.alleles)
  
  #split mono and polyallelic loci
  monoallelic_loci <- sorted_polyclonal %>%
    filter(n.alleles == 1)
  
  sorted_polyclonal <- sorted_polyclonal %>%
    filter(n.alleles > 1)
  
  #set max haplos expected
  n_haplos <- max(unique(sorted_polyclonal$n.alleles))
  
  #loci have the same number of alleles are phased just by sorting. give FINAL FORMAT as the action in the if statement and proceed with next sample
  if (length(unique(sorted_polyclonal$n.alleles)) == 1){ ##### ALL MULTIALLELIC LOCI HAVE THE SAME # OF ALLELES
    
    #GIVE FINAL FORMAT
    sorted_polyclonal<-rbind(monoallelic_loci, sorted_polyclonal)
    
    sorted_polyclonal <- sorted_polyclonal %>%
      arrange(sampleID, locus, n.alleles, desc(norm.reads.locus)) %>%
      group_by(locus) %>%
      arrange(desc(n.alleles))
    
    sorted_polyclonal <- sorted_polyclonal %>% 
      group_by(sampleID) %>%
      mutate(haplotypes = ifelse(n.alleles == 1, paste0("HAPLOTYPE", 1:max(n.alleles), collapse = "_"), paste0("HAPLOTYPE", 1:max(n.alleles[n.alleles > 1])))) %>%
      arrange(sampleID)
    
    PHASED_SAMPLE <- sorted_polyclonal
    
  } else { ###### AT LEAST 1 PAIR OF ALLELES HAS THE SAME # OF ALLELES (PERFORM ALLELE AVERAGING)
    
    ##loop to average frequencies of amplicons with the same number of alleles 
    unique_alleles <- unique(sorted_polyclonal$n.alleles)
    averaged_frequencies <- data.frame()
    
    for (allele in unique_alleles) {
      subset_rows <- sorted_polyclonal %>%
        filter(n.alleles == {{allele}})
      
      seq_unique_alleles <- seq(1, allele)
      
      if (length(unique(subset_rows$locus)) == 1) {
        
        break  # End the loop if there are no frequencies to average because number of alleles is not repeated across amplicons
        
      } else {
        
        for (num in seq_unique_alleles) {
          rows <- subset_rows %>%
            group_by(locus) %>%
            slice({{num}})
          
          #calculate average frequency and reads for each partial haplotype
          averaged_values <- mean(rows$norm.reads.locus)
          average_reads<-mean(rows$reads)
          
          #give proper format
          sampleID<-rows$sampleID[1]
          concatenated_locus <- paste(rows$locus, collapse = "@")
          concatenated_asv <- paste(rows$asv, collapse = "@")
          concatenated_alleles<- paste(rows$allele, collapse = "@")
          
          averaged_values <- cbind(sampleID=sampleID, locus=concatenated_locus, asv = concatenated_asv, reads = average_reads, allele=concatenated_alleles, norm.reads.allele=1, norm.reads.locus=averaged_values, n.alleles = allele)
          averaged_frequencies <- rbind(averaged_frequencies, averaged_values) 
          
        }
        #replace alleles with averaged partial haplotypes when needed
        sorted_polyclonal <- sorted_polyclonal[!(sorted_polyclonal$n.alleles %in% unique(averaged_frequencies$n.alleles)), ]
        sorted_polyclonal$norm.reads.allele <- 1
        sorted_polyclonal <- rbind(averaged_frequencies, as.data.frame(sorted_polyclonal))
      }
    }
      
  #REMAINING SAMPLES HAVE DIFFERENT # OF ALLELES ON EACH LOCI
  
  # Get adjacent pairs of loci and/or partial haplotypes
  adjacent_pairs <- lapply(1:(length(unique_alleles)-1), function(i) c(unique(sorted_polyclonal$locus)[i], unique(sorted_polyclonal$locus)[i+1]))
  adjacent_pairs <- discard(adjacent_pairs, ~any(is.na(.)))
  freqs <- data.frame(alleles = character(), sum_freqs = numeric())
  errors_list <- list() 
  
  # calcualte all sums of frequencies from the locus of the pair that has the most alleles
  for (pair in adjacent_pairs) {
    #cat("Locus Pair:", pair[1], "-", pair[2], "\n")
    
    combo_limit = abs(as.numeric(unique(sorted_polyclonal[sorted_polyclonal$locus==pair[1],][length(colnames(sorted_polyclonal))]))
                      -as.numeric(unique(sorted_polyclonal[sorted_polyclonal$locus==pair[2],][length(colnames(sorted_polyclonal))]))) +1 # combo_limit = max number of alleles that can be summed
    
    subset_rows_pair1 <- sorted_polyclonal[sorted_polyclonal$locus == pair[1], ]
    subset_rows_pair2 <- sorted_polyclonal[sorted_polyclonal$locus == pair[2], ]
    
    for (i in 1:combo_limit) {
      sums <- combn(as.numeric(subset_rows_pair2$norm.reads.locus), i, FUN = sum)
      
      # Get corresponding alleles
      allele_strings <- combn(subset_rows_pair2$allele, i, FUN = function(x) paste(x, collapse = "@"))
      
      # Create data frame with alleles and sum_freqs
      iteration_results <- data.frame(allele = allele_strings, sum_freqs = sums)
      freqs <- rbind(freqs, iteration_results)
    }
    
    #calculate errors
    errors<-abs(outer(as.numeric(subset_rows_pair1$norm.reads.locus), freqs$sum_freqs, "-"))
    colnames(errors)<- freqs$allele #pair[2] #highest amount of alleles
    rownames(errors)<-subset_rows_pair1$allele #pair[1], lowest amount of alleles
    
    # Prevent alleles from the most polyallellic amplicon from repeating on multiple haplotypes 
    for (i in 1:nrow(errors)) {
      
      # Fill columns matching row names with NAs: alleels from the same locus shouldn't be compared with each other
      matching_cols <- grepl(paste0("^", rownames(errors)[i], "(?:[^[:alnum:]]|$)"), colnames(errors))
      errors[i, matching_cols] <- NA
      errors[, apply(errors, 2, function(x) any(is.na(x)))] <- NA
      
      lowest_value <- min(errors[i, ], na.rm = TRUE)  # Exclude NA as the lowest value
      lowest_col <- colnames(errors)[which(errors[i, ] == lowest_value)][1] #### added [1] because happened to me that there was an instance of 2 columns with the exact same error. in the case i checked it got resolved by selecting the first. don't know how common this happens: muestra __N2011605_7_S83__ de ./SMC_retest_new_snps_new_output/mad4hatter-main/SMC2_NextSeq04_131222_results/allele_data.txt
      
      # Split lowest_col by '@' to handle multiple sub-strings
      lowest_cols <- unlist(strsplit(lowest_col, "@"))
      
      for (sub_col in lowest_cols) {
        pattern <- paste0("^", sub_col, "(?:[^[:alnum:]]|$)")
        
        # Grep the pattern in the rest of the columns and fill matched columns with NAs
        for (col in colnames(errors)[-which(colnames(errors) %in% lowest_col)]) {
          if (grepl(pattern, col)) {
            errors[, col] <- NA
            errors[, apply(errors, 2, function(x) any(is.na(x)))] <- NA
          }
        }
      }
    }
    errors_list[[paste(pair[1], "-", pair[2], sep = "")]] <- errors
  }
  
  
  ###linking partial haplotypes when needed (use sample that has multiple alleles in multiple loci to test)
  
  #extract partial haplotypes (those with the least error)
  partial_haplos <- list()
  
  for (i in 1:length(errors_list)) {
    df <- errors_list[[i]]
    haplos <- c()
    
    for (j in 1:nrow(df)) {
      rowname <- rownames(df)[j]
      lowest_value <- min(df[j, ], na.rm = TRUE)
      lowest_col <- colnames(df)[df[j, ] == lowest_value]
      lowest_col <- lowest_col[!is.na(lowest_col)]
      
      # Exclude the lowest column from further iterations
      df[, lowest_col] <- NA
      
      # Handle rows and columns with "@" symbol
      if (grepl("@", lowest_col)) {
        col_parts <- unlist(strsplit(lowest_col, "@"))
        for (part in col_parts) {
          result <- paste(rowname, part, sep = "@")
          #print(result)
          haplos <- c(haplos, result)
        }
      } else {
        result <- paste(rowname, lowest_col, sep = "@")
        #print(result)
        haplos <- c(haplos, result)
      }
    }
    
    partial_haplos[[i]] <- haplos
  }
  
  ###si sólo hay 1 dataframe en la lista partial_haplos (no hace falta linking, sólo acomodar en dataframe con @ como separador): IF
  if (length(partial_haplos) ==1 ){
    full_haplos<- data.frame(do.call(rbind, strsplit(partial_haplos[[1]], "@")))
  }else{
    ####si hay más de 1 dataframe en la lista partial_haplos: ELSE
    #formatting
    partial_haplos<-lapply(partial_haplos, function(x) {
      data.frame(strings = unlist(strsplit(x, "@")))
    })
    
    partial_haplos<-lapply(partial_haplos, function(df) {
      df$strings2 <- NA 
      df$strings2[c(FALSE, TRUE)] <- df$strings[c(FALSE, TRUE)] 
      df$strings[c(FALSE, TRUE)] <- NA  
      df <- fill(df, strings, .direction = "down") 
      df <- fill(df, strings2, .direction = "up")  
      df <- df[complete.cases(df), ]
      df <- df[seq(1, nrow(df), by = 2), ]
      return(df)
    })
    
    ##LINK PARTIAL HAPLOTYPES TOGETHER
    full_haplos <- data.frame()
    
    # Iterate through each pair of adjacent dataframes
    for (i in 1:(length(partial_haplos) - 1)) {
      
      if (i==0){
        i<-i+1
      }
      
      df1 <- partial_haplos[[i]]
      
      if (i + 1 <= length(partial_haplos)) {
        df2 <- partial_haplos[[i + 1]]
        
        for (j in 1:nrow(df1)) {
          element <- df1[j, ncol(df1)]
          
          if (nrow(df2) > 0) {
            matching_rows <- df2[grep(element, df2[, 1]), ]
            
            if (nrow(matching_rows) > 0) {
              df1_merged <- df1[j, -ncol(df1)]
              combined_rows <- cbind(df1_merged, matching_rows)
              full_haplos <- rbind(full_haplos, combined_rows)
            }
          }
        }
        
      } else {
        full_haplos<- as.data.frame(partial_haplos)
      }
    }
  }
  
  #label full haplotypes
  full_haplos$haplos <- paste0("HAPLOTYPE", 1:nrow(full_haplos))     
  
  ##if there are @ in sorted_polyclonal$locus... reformat with original alleles, removing the average alleles.
  if (any(grepl("@", sorted_polyclonal$locus))){
    #sort
    sorted_polyclonal <- muestra %>%
      arrange(sampleID, locus, n.alleles, desc(norm.reads.locus)) %>%
      group_by(locus) %>%
      arrange(n.alleles)
    
    #split mono and polyallelic loci
    monoallelic_loci <- sorted_polyclonal %>%
      filter(n.alleles == 1)
    
    sorted_polyclonal <- sorted_polyclonal %>%
      filter(n.alleles > 1)
  }
  
  #final format
  sorted_polyclonal$haplotypes<-NA
  
  for (i in 1:nrow(full_haplos)) {
    row_data <- full_haplos[i, ]
    match_indices <- sapply(row_data, function(element) grep(element, sorted_polyclonal$allele))
    match_indices <- unique(unlist(match_indices))
    if (length(match_indices) > 0) {
      for (index in match_indices) {
        if (is.na(sorted_polyclonal$haplotypes[index])) {
          sorted_polyclonal$haplotypes[index] <- row_data$haplos
        } else {
          existing_haplos <- sorted_polyclonal$haplotypes[index]
          new_haplos <- paste0(existing_haplos, "_", row_data$haplos)
          sorted_polyclonal$haplotypes[index] <- new_haplos
        }
      }
    }
  }
  
  monoallelic_loci$haplotypes<-paste0("HAPLOTYPE", 1:n_haplos, collapse = "_")
  
  #concatenate complete sample
  PHASED_SAMPLE<-rbind(monoallelic_loci, sorted_polyclonal)
  
  }
  PARTIAL_OUTPUT_2 <- rbind(PHASED_SAMPLE, PARTIAL_OUTPUT_2)
  
  print(paste(sample, ": PHASED!")) 
}

OUTPUT_FINAL <- rbind(PARTIAL_OUTPUT_1, PARTIAL_OUTPUT_2)

write.table(OUTPUT_FINAL, file=paste0(outfile,".csv"), quote=F, sep="\t", col.names=T, row.names=F)



###ESTE MÉTODO ITERATIVO IGNORA COI. USANDO COI SERÍA UNA ESTRATEGIA DISTINTA, SUPONGO... PROBAR ASÍ CON TODAS LAS ITERACIONES POSIBLES PRIMERO?

#4 amplicons, all biallelic

rm(a1, a2, a3, a4)

a1 <- data.frame(dhps_437 = c("A", "G"), 
                 freq = c(0.75, 0.25))
a2 <- data.frame(dhps_540 = c("K"), 
                 freq = c(1))
a3 <-data.frame(dhfr_51_59 = c("IR", "NC"), 
                 freq = c(0.84, 0.16))
a4 <- data.frame(dhfr_108 = c("N"), 
                 freq = c(1))


#create all possible combinations of alleles and freqs
alleles <- list(a1$dhps_437, a2$dhps_540, a3$dhfr_51_59, a4$dhfr_108)
freqs <- list(a1$freq, a2$freq, a3$freq, a4$freq)

comb_alleles <- expand.grid(alleles)
comb_freqs <- expand.grid(freqs)

comb_alleles_matrix <- as.data.frame(comb_alleles)
colnames(comb_alleles_matrix) <- c("dhps_437", "dhps_540", "dhfr_51_59", "dhfr_108")
comb_freqs_matrix <- as.data.frame(comb_freqs)
colnames(comb_freqs_matrix) <- c("dhps_437", "dhps_540", "dhfr_51_59", "dhfr_108")

comb_freqs_matrix$freq_mean <- rowMeans(comb_freqs_matrix, na.rm = TRUE) #calculate mean freq for all possible haplotypes
comb_freqs_matrix$probs<-comb_freqs_matrix$dhps_437 * comb_freqs_matrix$dhps_540 * comb_freqs_matrix$dhfr_51_59 *comb_freqs_matrix$dhfr_108 #calculate probs if all haploptypes were present
comb_freqs_matrix$SD <- apply(comb_freqs_matrix[, 1:4], 1, sd) # Calculate Standard Deviation (SD)
comb_freqs_matrix$CV <- (comb_freqs_matrix$SD / comb_freqs_matrix$freq_mean) # Calculate Coefficient of Variation (CV)

#plot shit
plot(comb_freqs_matrix$CV, comb_freqs_matrix$probs)

## select the "BEST" haplo: highest prob and lowest CV
lowest_CV <- which.min(comb_freqs_matrix$CV)
highest_prob <- as.numeric(which.max(comb_freqs_matrix$probs))

if (lowest_CV == highest_prob){
  
  most_likely_hap<-paste(as.matrix(comb_alleles_matrix[highest_prob, ]), collapse="_")
  print(paste( most_likely_hap, "is the most likely true haplotype.", collapse=" "))
  
}else{
  
  most_likely_hap1<-paste(as.matrix(comb_alleles_matrix[highest_prob, ]), collapse="_")
  most_likely_hap2<-paste(as.matrix(comb_alleles_matrix[lowest_CV, ]), collapse="_")
  print(paste("One of", most_likely_hap1, "and", most_likely_hap2, "is the most likely true haplotype. Visually examine the plot."))
  
  plot(comb_freqs_matrix$CV, comb_freqs_matrix$probs)
  
  ##BREAK??
}

#store most likely haplo1 in new dfs:

MOST_LIKELY_HAPLOS <- comb_alleles_matrix[highest_prob, ]
MOST_LIKELY_HAPLOS_FREQS <- comb_freqs_matrix[highest_prob, ]
MOST_LIKELY_HAPLOS_FREQS$HAPLO_FREQ<-min(comb_freqs_matrix[highest_prob,1:4])



#################################################################3

##asumiendo sólo un haplo claramente most likely... seeguir:

#select minimum allele freq of most likely haplotype
min_allele_from_most_lilely_hap <- min(comb_freqs_matrix[highest_prob,1:4])

#boolean mask to detect alleles that are present on the most likely haplotype
row_to_match <- as.matrix(comb_alleles_matrix[highest_prob, ])

  mask <- sapply(colnames(comb_alleles_matrix), function(col_name) {
  comb_alleles_matrix[, col_name] == row_to_match[,col_name]
})


modified_comb_freqs_matrix <- comb_freqs_matrix[,1:4]
  
# Subtract min_allele_from_most_lilely_hap from the cells where mask is TRUE and ignore the specified column
modified_comb_freqs_matrix[mask] <-  modified_comb_freqs_matrix[mask] - min_allele_from_most_lilely_hap



#NEXT ITERATION(?)
modified_comb_freqs_matrix$freq_mean <- rowMeans(modified_comb_freqs_matrix, na.rm = TRUE) #calculate mean freq for all possible haplotypes
modified_comb_freqs_matrix$probs<-modified_comb_freqs_matrix$dhps_437 * modified_comb_freqs_matrix$dhps_540 * modified_comb_freqs_matrix$dhfr_51_59 *modified_comb_freqs_matrix$dhfr_108 #calculate probs if all haploptypes were present
modified_comb_freqs_matrix$SD <- apply(modified_comb_freqs_matrix[, 1:4], 1, sd) # Calculate Standard Deviation (SD)
modified_comb_freqs_matrix$CV <- (modified_comb_freqs_matrix$SD / modified_comb_freqs_matrix$freq_mean) # Calculate Coefficient of Variation (CV)

#plot shit
plot(modified_comb_freqs_matrix$CV, modified_comb_freqs_matrix$probs)

## select the "BEST" haplo: highest prob and lowest CV
lowest_CV <- which.min(modified_comb_freqs_matrix$CV)
highest_prob <- as.numeric(which.max(modified_comb_freqs_matrix$probs))

if (lowest_CV == highest_prob){
  
  most_likely_hap<-paste(as.matrix(comb_alleles_matrix[highest_prob, ]), collapse="_")
  print(paste( most_likely_hap, "is the most likely true haplotype.", collapse=" "))
  
}else{
  
  most_likely_hap1<-paste(as.matrix(comb_alleles_matrix[highest_prob, ]), collapse="_")
  most_likely_hap2<-paste(as.matrix(comb_alleles_matrix[lowest_CV, ]), collapse="_")
  print(paste("One of", most_likely_hap1, "and", most_likely_hap2, "is the most likely true haplotype. Visually examine the plot."))
  
  plot(comb_freqs_matrix$CV, comb_freqs_matrix$probs)
  
  ##BREAK??
}

#append next best haplo:
MOST_LIKELY_HAPLOS <- rbind(MOST_LIKELY_HAPLOS,comb_alleles_matrix[highest_prob, ])
temp <-modified_comb_freqs_matrix[highest_prob, ]
temp$HAPLO_FREQ <- min(modified_comb_freqs_matrix[highest_prob,1:4])
MOST_LIKELY_HAPLOS_FREQS <- rbind(MOST_LIKELY_HAPLOS_FREQS, temp )


#################################################################3

##asumiendo sólo un haplo claramente most likely... seeguir:
##asumiendo sólo un haplo claramente most likely... seeguir:

#select minimum allele freq of most likely haplotype
min_allele_from_most_lilely_hap <- min(modified_comb_freqs_matrix[highest_prob,1:4])

#boolean mask to detect alleles that are present on the most likely haplotype
row_to_match <- as.matrix(comb_alleles_matrix[highest_prob, ])

mask <- sapply(colnames(comb_alleles_matrix), function(col_name) {
  comb_alleles_matrix[, col_name] == row_to_match[,col_name]
})


modified_comb_freqs_matrix <- modified_comb_freqs_matrix[,1:4]

# Subtract min_allele_from_most_lilely_hap from the cells where mask is TRUE and ignore the specified column
modified_comb_freqs_matrix[mask] <-  modified_comb_freqs_matrix[mask] - min_allele_from_most_lilely_hap



#NEXT ITERATION(?)
modified_comb_freqs_matrix$freq_mean <- rowMeans(modified_comb_freqs_matrix, na.rm = TRUE) #calculate mean freq for all possible haplotypes
modified_comb_freqs_matrix$probs<-modified_comb_freqs_matrix$dhps_437 * modified_comb_freqs_matrix$dhps_540 * modified_comb_freqs_matrix$dhfr_51_59 *modified_comb_freqs_matrix$dhfr_108 #calculate probs if all haploptypes were present
modified_comb_freqs_matrix$SD <- apply(modified_comb_freqs_matrix[, 1:4], 1, sd) # Calculate Standard Deviation (SD)
modified_comb_freqs_matrix$CV <- (modified_comb_freqs_matrix$SD / modified_comb_freqs_matrix$freq_mean) # Calculate Coefficient of Variation (CV)

#plot shit
plot(modified_comb_freqs_matrix$CV, modified_comb_freqs_matrix$probs)

## select the "BEST" haplo: highest prob and lowest CV
lowest_CV <- which.min(modified_comb_freqs_matrix$CV)
highest_prob <- as.numeric(which.max(modified_comb_freqs_matrix$probs))

if (lowest_CV == highest_prob){
  
  most_likely_hap<-paste(as.matrix(comb_alleles_matrix[highest_prob, ]), collapse="_")
  print(paste( most_likely_hap, "is the most likely true haplotype.", collapse=" "))
  
}else{
  
  most_likely_hap1<-paste(as.matrix(comb_alleles_matrix[highest_prob, ]), collapse="_")
  most_likely_hap2<-paste(as.matrix(comb_alleles_matrix[lowest_CV, ]), collapse="_")
  print(paste("One of", most_likely_hap1, "and", most_likely_hap2, "is the most likely true haplotype. Visually examine the plot."))
  
  plot(comb_freqs_matrix$CV, comb_freqs_matrix$probs)
  
  ##BREAK??
}

#append next best haplo:
MOST_LIKELY_HAPLOS <- rbind(MOST_LIKELY_HAPLOS,comb_alleles_matrix[highest_prob, ])
temp <-modified_comb_freqs_matrix[highest_prob, ]
temp$HAPLO_FREQ <- min(modified_comb_freqs_matrix[highest_prob,1:4])
MOST_LIKELY_HAPLOS_FREQS <- rbind(MOST_LIKELY_HAPLOS_FREQS, temp )

### DONE, CLEAN! poner chequeo: si current SD = 0, END.
















## with eCOI, determine how many MAX haplos are expected
## find the "BEST" complementary haplo(s) that add up to 100 as close as possible (without considering eCOI)
## find the "BEST" complementary haplo(s) that add up to 100 as close as possible considering eCOI
## same? PERFECT; not the same? WAT DO?

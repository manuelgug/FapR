
library(dplyr)
library(tidyr)


filename <- "../pf-haploatlas-PF3D7_0417200_population_summary.csv" #downloaded from https://pf-haploatlas.streamlit.app/
real_haplo_counts <- read.csv(filename)
gene_names <- read.csv("../pf_gene_names.csv") # downloaded from https://plasmodb.org/plasmo/app/step/434333253/download

## format table
row.names(real_haplo_counts)<-real_haplo_counts$ns_changes
real_haplo_counts <- real_haplo_counts[,4:14]
real_haplo_counts <- real_haplo_counts[rowSums(real_haplo_counts) > 0,] #remove haplos with zero counts


## complete the haplos using the mutations in the amplicons of interest only
# separate mutations
mutations <- strsplit(rownames(real_haplo_counts), "/")
unique_mutations <- unique(do.call(c, mutations))

# Function to order mutations ascendingly by locus 
order_mutations <- function(mutations) {
  numbers <- as.numeric(gsub("[A-Z]", "", mutations))
  mutations[order(numbers)]
}

# order mutations
mutations <- lapply(mutations, order_mutations)
unique_mutations <- order_mutations(unique_mutations)

# Function to replace the last letter with the first letter (wt locus)
replace_last_with_first <- function(mutation) {

  first_letter <- substr(mutation, 1, 1)
  middle_part <- substr(mutation, 2, nchar(mutation) - 1)
  last_letter <- substr(mutation, nchar(mutation), nchar(mutation))
  
  wt_locus <- paste0(first_letter, middle_part, first_letter)
  
  return(wt_locus)
}

# set wt loci
wt_loci <- sapply(unique_mutations, replace_last_with_first)


#complete the haplos with wt and ns mutations
mutations_all <- lapply(mutations, function(mut_vec) {
  # Replace each mutation in the vector with the corresponding wt_loci value
  wt_loci_mod <- wt_loci
  for (mut in mut_vec) {
    if (mut %in% names(wt_loci_mod)) {
      # Replace the last letter in wt_loci with the last letter of the mutation
      substr(wt_loci_mod[mut], nchar(wt_loci_mod[mut]), nchar(wt_loci_mod[mut])) <- substr(mut, nchar(mut), nchar(mut))
    }
  }
  return(wt_loci_mod)
})

#create haplo df
df_haplos <- as.data.frame(do.call(rbind,mutations_all))

#calculate relative frequencies by column: these are gonna be the WEIGHTS
real_haplo_counts <- real_haplo_counts %>%
  mutate(across(everything(), ~ . / sum(.)))

#merge tables
merged_haplos <- cbind(df_haplos, real_haplo_counts)


#Step 1: Identify all unique mutations
all_wt_ns_unique_mutations <- unique(do.call(c, mutations_all))

# Step 2: Create an empty result dataframe to store aggregated counts
columns_of_interest <- colnames(real_haplo_counts)
WEIGHTS <- data.frame(mutation = all_wt_ns_unique_mutations, matrix(0, nrow = length(all_wt_ns_unique_mutations), ncol = length(columns_of_interest)))
colnames(WEIGHTS)[2:ncol(WEIGHTS)] <- columns_of_interest

# Step 3: Loop through each unique mutation and sum the counts across all relevant rows
for (mut in all_wt_ns_unique_mutations) {
  # Create a boolean mask where the mutation occurs
  mask <- apply(df_haplos, 1, function(x) any(x == mut))
  
  # Sum the counts for this mutation across the relevant rows and store in the result dataframe
  WEIGHTS[WEIGHTS$mutation == mut, 2:ncol(WEIGHTS)] <- round(colSums(merged_haplos[mask, columns_of_interest]), 4)
}


#add gene_ID column
gene_ID <- sub(".*-(PF3D7_\\d+)_population_summary\\.csv$", "\\1",  basename(filename))

gene_names$Gene.Name.or.Symbol <- tolower(gene_names$Gene.Name.or.Symbol)

gene <- gene_names[gene_names$Gene.ID %in% gene_ID, ]$Gene.Name.or.Symbol
gene <- sub("-.*", "", gene)

WEIGHTS$locus <- paste0(gene, "_",  gsub("([A-Z])(\\d+)([A-Z])", "\\2", WEIGHTS$mutation)) # create locus column (gene + locus)
WEIGHTS$mutation <- substring(WEIGHTS$mutation, nchar(WEIGHTS$mutation), nchar(WEIGHTS$mutation)) # Keep only the actual mutation (regardless of wt or mut)

WEIGHTS <- WEIGHTS %>%
  mutate(locus_number = as.numeric(sub(".*_(\\d+)$", "\\1", WEIGHTS$locus))) %>%
  select(locus, everything()) %>%
  arrange(locus_number) %>%
  select(everything(), -locus_number)

output <- paste0("WEIGHTS_", gene, "_", basename(filename))

write.csv(WEIGHTS, output, row.names = F)

#
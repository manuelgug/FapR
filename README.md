# FapR: Frequency-based amplicon phasing in R

This script performs the phasing of `dfhr` and `dhps` amplicons from polyclonal infections of *Plasmodium falciparum*, utilizing only their allele frequencies.

## Dependencies

The following R packages must be installed: `ggplot2`, `dplyr`, `gridExtra`, and `optparse`.

## Usage

```shell
Rscript FapR.R -i [resmarker_table_global_max_0_filtered.csv] -o [output_prefix]
```

### Arguments

- `allele_data`: Path to the input allele data file.
- `outfile`: Prefix appended to the output files.

## Script Overview

The script can be divided into three main sections: data preparation, exploratory data analysis, and phasing.


### Phasing

#### Assumptions
1. All alleles in the input file are true alleles. Appropriate filtering is strongly suggested.
2. There are no copy number variants (CNV) on the amplicons.

#### Conditions


### Algorithm

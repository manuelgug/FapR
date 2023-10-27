# FapR: Frequency-based amplicon phasing in R

This script performs the phasing of resistance markers from polyclonal infections of *Plasmodium falciparum* utilizing only allele frequencies. As of now, 7 resistance markers across 4 [mad4hatter amplicons](https://www.protocols.io/view/mad4hatter-14egn779mv5d/v3) from the `dhfr` and `dhps` genes are supported.

## Dependencies

The following R packages must be installed: `ggplot2`, `dplyr`, `gridExtra`, and `optparse`.

## Usage

```shell
Rscript FapR.R -i [resmarker_table_global_max_0_filtered.csv] -o [output_prefix]
```

## Script Overview

#### Assumptions
1. All alleles in the input file are true alleles. Appropriate filtering is strongly suggested.
2. There are no copy number variants (CNV) on the amplicons.

### Phasing
FapR uses an iterative approach on which haplotypes are accepted based on 

1. **Probability of occurring in a sample**: haplotypes built from highly abundant resmarkers are more liekly to be true.
2. **Variance on the resmarker frequencies**: haplotypes built from similarly abundant resmarkers are more liekly to be true.

### Flagging
Phased haplotypes are flagged based on: 

1. **Frequency in the sequencing run** (assuming it is from a given population): allows to catch haplotypes that are frequent in the run, but have a low abundance in particular samples. This flag takes precedent over the following.
2. **Limit of detection of each amplicon** (experimentally tested): allows to catch haplotypes that are rare in the run, but moderate to highly abundant in particular samples. This also allows to build partial haplotypes.



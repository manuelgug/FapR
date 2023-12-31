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
FapR uses an iterative approach in which haplotypes are accepted based on: 

1. **Probability of occurring in a sample**: haplotypes built from highly abundant resmarkers are more likely to be true.
2. **Variance on the resmarker frequencies**: haplotypes built from similarly abundant resmarkers are more likely to be true.

![algo](https://github.com/manuelgug/FapR/blob/main/img/fapr_algo.png)

*Figure 1. FapR's phasing algorithm.*

![example](https://github.com/manuelgug/FapR/blob/main/img/fapr_example.png)

*Figure 2. Example of the phasing process. The best haplotype on each iteration (highest probability and lowest coefficient of variation) is highlighted in green and its assigned frequency in pink. Haplotypes with zero probability are highlighted in orange. This particular sample resulted in 3 haplotypes that add up to 99%.*

### Flagging
Phased haplotypes are flagged based on: 

1. **Frequency in the sequencing run** (assuming it is from a given population)
    + A single threshold derived from population frequency
    + Allows to catch haplotypes that are frequent in the run, but have low abundance in particular samples
    + This flag takes precedence over the following
      
2. **Limit of detection of each amplicon** (experimentally tested)
    + A threshold for each amplicon
    + Allows to catch haplotypes that are rare in the run, but moderate to highly abundant in particular samples
    + Allows to build partial haplotypes

![flagging](https://github.com/manuelgug/FapR/blob/main/img/fapr_flagging.png)

*Figure 3. FapR's flagging algorithm. Currently, population frequency threshold is the mean haplotype frequency.*

![flagging_example](https://github.com/manuelgug/FapR/blob/main/img/fapr_flagging_example.png)

*Figure 4. Example of the flagging results. Checkmarks are accepted haplotypes: green = correctly phased and also frequent in the population/run; blue = correctly phased but rare in the population/run; purple = inconclusive phasing but frequent in the population/run. Orange crosses are haplotypes with inconclusive phasing and also rare (or absent) in the population/run, thus being inconclusive haplotypes.*


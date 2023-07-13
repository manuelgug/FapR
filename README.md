# FapR: Frequency-based amplicon phasing in R

This script performs the phasing of `dfhr` and `dhpr` amplicons from polyclonal infections of *Plasmodium falciparum*, utilizing only their allele frequencies.

## Dependencies

The following R packages must be installed: `ggplot2`, `ggbeeswarm`, `dplyr`, `gridExtra`, and `purrr`.

## Usage

```shell
Rscript SCRIPT.R allele_data outfile
```

### Arguments

- `allele_data_filtered`: Path to the input allele data file.
- `outfile`: Prefix appended to the output files.

## Script Overview

The script can be divided into three main sections: data preparation, exploratory data analysis, and phasing.

### Data Preparation

The script first subsets the amplicons of interest (`dfhr` and `dhpr`) based on the specified amplicon names and updates their names accordingly. If you need to modify the amplicon names or add more amplicons, please modify these variables accordingly.

Next, the script divides the data into monoclonal and polyclonal samples based on the number of alleles. Monoclonal samples have only one allele for all specified amplicons, while polyclonal samples have multiple alleles for at least one amplicon.

### Exploratory Data Analysis

The script performs a brief exploratory data analysis on the polyclonal samples. It calculates and visualizes various statistics, such as the ratio of monoclonal to polyclonal samples, variance in allele frequency for each amplicon, the number of alleles for each amplicon, the distribution of the number of alleles, and allele frequencies for each amplicon.

FIGURE

### Phasing

#### Assumptions
1. All alleles in the input file are true alleles. (Filtering should be performed)
2. There are no copy number variants (CNV) on the amplicons.

#### Conditions
1. The amplicon with the most alleles determines the expected number of haplotypes.
2. Each allele from the amplicon with the most alleles should be used only once to build the haplotypes.
3. Alleles from the rest of the amplicons should be used at least once to build the haplotypes.

### Algorithm

#### Workflow

##### Samples that have a single polyallelic amplicon
No ambiguity exists in these samples. All of the monoallelic amplicons are paired with each of the alleles from the amplicon with the most alleles.

FIGURE

##### Samples that have the same number of alleles on all of their multiallelic loci
Alleles are sorted by frequency and grouped accordingly.

FIGURE

##### Samples with amplicons that have a different number of alleles

Amplicons are sorted in descending order by the number of alleles. Phasing is performed on adjacent pairs of amplicons, where each adjacent pair consists of a MajAMP (amplicon with the least number of alleles) and a MinAMP (amplicon with the most number of alleles).

FIGURE1

0. Frequencies from alleles belonging to amplicons with the same number of alleles are averaged and treated as a single allele from this point onward.

FIGURE2

1. Calculate the maximum number of alleles from MinAMP that can be paired with a single allele of MajAMP:
```shell
combo_limit = (#alleles MinAMP - #alleles MajAMP) + 1
```
2. Sum frequencies from all possible combinations of groups of `combo_limit` alleles from MinAMP.
 
```shell
combo_limit = 2

2A = 0.69
2B = 0.22
2C = 0.09
2A + 2B = 0.91
2A + 2C = 0.78
2B + 2C = 0.31
```

3. Subtract the resulting frequencies from the previous step, as well as individual alleles of MinAMP, from the most frequent allele of MajAMP to calculate the errors. The allele or combination of alleles resulting in the least error are then phased to the most frequent allele of MajAMP.
```shell
abs(allele 1 MajAMP - [combo1 MinAMP, combo 2 MinAMP, allele 1 MinAMP], ...)
```
```shell
abs(1A - [2A, 2B, 2C, 2A + 2B, 2A + 2C, 2B + 2C])

error calculations:
abs(0.725 - 0.69) = 0.035 <---- BEST
abs(0.725 - 0.22) = 0.505
abs(0.725 - 0.09) = 0.635
abs(0.725 - 0.91) = 0.255
abs(0.725 - 0.78) = 0.055 
abs(0.725 - 0.31) = 0.415
```
```shell
Partial haplotype 1 = 1A_2A
```

4. Remove all options, whether individual or combinations, that already contain phased alleles and repeat step 3 for each remaining allele of MajAMP.

```shell
abs(1A - [2B,  2C, 2B + 2C])

error calculations:
abs(0.275 - 0.22) = 0.055
abs(0.275 - 0.09) = 0.185
abs(0.275 - 0.31) = 0.035 <---- BEST
```
```shell
Partial haplotype 2 = 1B_2B
Partial haplotype 3 = 1B_2C
```
5. Repeat steps 1-5 for the next adjacent pair.

```shell
Partial haplotype 4 = 2A_3A
Partial haplotype 5 = 2A_3C
Partial haplotype 6 = 2B_3B
Partial haplotype 7 = 2C_3D
```
6. Link partial haplotypes by the common alleles to build the final haplotypes.

 ```shell
Final haplotype 1 = 1A_2A_3A
Final haplotype 2 = 1A_2A_3C
Final haplotype 3 = 1B_2B_3B
Final haplotype 4 = 1B_2C_3D
```
- The amplicon with the most alleles determines the expected number of haplotypes ✔
- Each allele from the amplicon with the most alleles should be used only once to build the haplotypes ✔
- Alleles from the rest of the amplicons should be used at least once to build the haplotypes ✔

7. Annotate the original allele data with a new column indicating the corresponding haplotype each allele belongs to.

